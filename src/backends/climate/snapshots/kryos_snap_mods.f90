module kryos_snap_mods
    ! Climate transforms (yelmox climate_mods): static changes applied to the
    ! adapter's climate_out after the snapshots and the terms, in a fixed order:
    !   taper -- the atmosphere towards the reference with a static weight
    !            field (regional taper)
    !   bias  -- bias correction of the atmosphere (static, or sliced at the
    !            forcing time when its groups have a time axis)
    ! Each is an optional group named by the instance keys `taper` and `bias`
    ! ("" = off; a named group is active); unset keys from the templates taper,
    ! bias of input/snapshots_defaults.nml. Their fields are varslice groups,
    ! read from the par file, else the first library file that has them.
    !
    ! Contract: each transform states what it touches in `now` (`ref` is never
    ! touched); climate_mods_init checks every active transform against what
    ! the backend supplies (climate_caps_class). Both change the atmosphere
    ! after the backend, so both need an atmosphere that drives the surface
    ! (caps_need_atm, the predicate of the uniform terms too).
    !   taper: now = ref + w*(now - ref) for tas, pr, ta_ann, pr_ann and
    !         tsl_ann (those the climate has), w in [0, 1] a static field on
    !         the climate grid, the reference at the current surface
    !         (taper_apply); with outside = "mean", ref + w*(now - ref) +
    !         (1 - w)*the mean anomaly (temperature; precipitation: ratio of the
    !         totals) over the cells with w = 1, per month.
    !   bias: tas, ta_ann, tsl_ann + dtas; pr * fpr, pr_ann by the ratio of the
    !         annual totals.

    use precision,        only: wp, dp
    use fesm_log,         only: LOG_INFO, log_msg, fesm_error, to_str
    use nml,              only: nml_read, nml_validate
    use ncio,             only: nc_write
    use coordinates,      only: grid_class
    use varslice,         only: varslice_class, varslice_init_nml, varslice_update, varslice_end
    use kryos_snap_out,   only: climate_out_class
    use kryos_snap_terms, only: climate_caps_class, caps_need_atm
    use kryos_snap_nml,   only: SNAP_DEFAULTS, FILE_LEN, group_file

    implicit none

    private

    ! The atmosphere towards the reference with a weight field w in [0, 1]
    ! (one 2D field), read from a varslice group on the climate grid (online
    ! remap allowed): now = ref + w*(now - ref); w = 1 keeps the climate, w = 0
    ! is the reference (outside = "ref") or the reference plus the mean anomaly
    ! of the cells with w = 1 (outside = "mean").
    type, public :: taper_class
        logical               :: active  = .false.
        character(len=64)     :: group   = ""
        character(len=256)    :: weight  = ""      ! the varslice group of w
        character(len=8)      :: outside = "ref"   ! ref | mean: the climate where w = 0
        real(wp), allocatable :: w(:,:)            ! [1] (nx,ny)
        logical,  allocatable :: full(:,:)         ! w = 1: the cells of the mean anomaly
    end type taper_class

    ! An additive temperature field and a precipitation factor field, monthly
    ! (12) or annual (1 field, all months), read from varslice groups on the
    ! climate grid (online remap allowed), each with a weight:
    !   tas = tas + f_tas*dtas,   pr = pr*(1 + f_pr*(fpr - 1))
    ! A group with a time axis is kept open and sliced at the forcing time at
    ! every update (bias_update).
    type, public :: bias_class
        logical               :: active = .false.
        character(len=64)     :: group  = ""
        real(wp)              :: f_tas = 1.0_wp, f_pr = 1.0_wp
        logical               :: has_dtas = .false., has_fpr = .false.
        logical               :: tdep = .false.   ! a group with a time axis
        type(varslice_class)  :: vs_dtas, vs_fpr
        real(wp), allocatable :: dtas(:,:,:)     ! [K]  (nx,ny,12), 0 without a dtas group
        real(wp), allocatable :: fpr(:,:,:)      ! [1]  (nx,ny,12), 1 without a fpr group
    end type bias_class

    type, public :: climate_mods_class
        type(climate_caps_class) :: caps
        type(taper_class)        :: taper
        type(bias_class)         :: bias
    end type climate_mods_class

    public :: climate_mods_init
    public :: climate_mods_apply
    public :: climate_mods_write_init, climate_mods_write_step

contains

    subroutine climate_mods_init(mods, par_file, library, taper_group, bias_group, caps, domain, grid, time)
        ! Read the transform groups named (an empty name: inactive) and check
        ! each active one's contract against what the backend supplies (caps).
        ! grid: the climate grid. time: the forcing time, for the bias groups
        ! with a time axis. Call after the backend init.
        type(climate_mods_class), intent(inout) :: mods
        character(len=*),         intent(in)    :: par_file, library(:), taper_group, bias_group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid
        real(dp),                 intent(in)    :: time

        mods%caps = caps

        call taper_init(mods%taper, par_file, library, taper_group, caps, domain, grid)
        call bias_init(mods%bias, par_file, library, bias_group, caps, domain, grid, time)

    end subroutine climate_mods_init

    subroutine climate_mods_apply(mods, out, time)
        ! The active transforms on `out`, in the fixed order; the bias fields
        ! sliced at `time` (the forcing time) first when they have a time axis.
        type(climate_mods_class), intent(inout) :: mods
        type(climate_out_class),  intent(inout) :: out
        real(dp),                 intent(in)    :: time

        if (mods%taper%active) call taper_apply(mods%taper, out)

        if (mods%bias%active) then
            if (mods%bias%tdep) call bias_update(mods%bias, time)
            call bias_apply(mods%bias, out)
        end if

    end subroutine climate_mods_apply

    subroutine climate_mods_write_init(mods, filename)
        ! The static fields of the active transforms, into the climate's 2D
        ! output file (xc, yc defined): the taper weight; the bias correction,
        ! weighted, as annual means (of the start, for a time-dependent bias).
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename

        if (mods%taper%active) &
            call nc_write(filename, "taper_w", mods%taper%w, dim1="xc", dim2="yc", units="1", &
                          long_name="Weight of the climate towards the reference (taper)")

        if (mods%bias%active) then
            associate (b => mods%bias)
            call nc_write(filename, "bias_dtas_ann", b%f_tas*sum(b%dtas, dim=3)/12.0_wp, &
                          dim1="xc", dim2="yc", units="K", &
                          long_name="Bias correction of the air temperature (ann)")
            call nc_write(filename, "bias_fpr_ann", &
                          1.0_wp + b%f_pr*(sum(b%fpr, dim=3)/12.0_wp - 1.0_wp), &
                          dim1="xc", dim2="yc", units="1", &
                          long_name="Bias correction factor of the precipitation (ann)")
            end associate
        end if

    end subroutine climate_mods_write_init

    subroutine climate_mods_write_step(mods, filename, ncid, n)
        ! Record n of the active transforms' diagnostics, into an open 2D
        ! output file: the weighted annual temperature correction of a
        ! time-dependent bias, bias_dtas_ann_now (the taper and a static bias
        ! are written once by climate_mods_write_init).
        type(climate_mods_class), intent(in) :: mods
        character(len=*),         intent(in) :: filename
        integer,                  intent(in) :: ncid, n

        if (mods%bias%active .and. mods%bias%tdep) then
            associate (b => mods%bias)
            call nc_write(filename, "bias_dtas_ann_now", b%f_tas*sum(b%dtas, dim=3)/12.0_wp, &
                          dim1="xc", dim2="yc", dim3="time", start=[1,1,n], &
                          count=[size(b%dtas,1),size(b%dtas,2),1], ncid=ncid, units="K", &
                          long_name="Bias correction of the air temperature at the forcing time (ann)", &
                          prec_abs=1e-3_wp)
            end associate
        end if

    end subroutine climate_mods_write_step

    ! ===== taper ===========================================================

    subroutine taper_init(tap, par_file, library, group, caps, domain, grid)
        ! Read the taper group (a named group is active), check its contract
        ! (caps_need_atm) and load the weight from its varslice group (static,
        ! one field) on the climate grid; a weight outside [0, 1] (or missing)
        ! is an error.
        type(taper_class),        intent(inout) :: tap
        character(len=*),         intent(in)    :: par_file, library(:), group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid

        type(varslice_class) :: vs
        character(len=FILE_LEN) :: gfile
        character(len=*), parameter :: proc = "kryos_snap_mods::taper_init"

        tap%group  = group
        tap%active = len_trim(group) > 0
        if (.not. tap%active) return

        gfile = group_file(par_file, library, trim(group), "the taper transform")
        call nml_validate(trim(gfile), SNAP_DEFAULTS, trim(group), defaults_group="taper")
        call caps_need_atm(caps, "["//trim(group)//"]")

        call nml_read(trim(gfile), trim(group), "weight",  tap%weight,  defaults_file=SNAP_DEFAULTS, defaults_group="taper")
        call nml_read(trim(gfile), trim(group), "outside", tap%outside, defaults_file=SNAP_DEFAULTS, defaults_group="taper")
        select case (trim(tap%outside))
        case ("ref", "mean")
        case default
            call fesm_error(proc, "["//trim(group)//"] outside must be ref or mean; got '"//trim(tap%outside)//"'")
        end select

        call varslice_init_nml(vs, trim(group_file(par_file, library, trim(tap%weight), "["//trim(group)//"] weight")), &
                               trim(tap%weight), domain=domain, grid_name=trim(grid%name), grid=grid)
        call varslice_update(vs)
        if (size(vs%var, 3) /= 1 .or. size(vs%var, 4) /= 1) then
            call fesm_error(proc, "["//trim(group)//"] weight "//trim(tap%weight)//" must be one 2D field", &
                            "got "//to_str(size(vs%var, 3))//" x "//to_str(size(vs%var, 4))//" slices")
        end if
        tap%w = vs%var(:,:,1,1)
        call varslice_end(vs)

        if (any(tap%w < 0.0_wp .or. tap%w > 1.0_wp)) then
            call fesm_error(proc, "["//trim(group)//"] weight "//trim(tap%weight)//" must lie in [0, 1] "// &
                            "(missing values: set fill in the group)", &
                            "range: "//to_str(real(minval(tap%w), dp))//" "//to_str(real(maxval(tap%w), dp)))
        end if

        tap%full = tap%w >= 1.0_wp
        if (trim(tap%outside) == "mean" .and. .not. any(tap%full)) then
            call fesm_error(proc, "["//trim(group)//"] outside = mean needs cells with w = 1 (the mean anomaly)", &
                            "weight "//trim(tap%weight)//" has none")
        end if

        call log_msg(LOG_INFO, "kryos_snap_mods", trim(group)//" active, weight = "//trim(tap%weight)// &
                     ", outside = "//trim(tap%outside)//", cells with w = 1: "//to_str(count(tap%full)))

    end subroutine taper_init

    subroutine taper_apply(tap, out)
        ! The atmosphere of `out` towards its reference: now = ref + w*(now - ref),
        ! for each field the climate has (now and ref), with the reference at
        ! the surface of `now` (now%z_srf): a reference given at its own surface
        ! (ref%z_srf: the elevation of the reference snapshot) is moved there
        ! with the monthly lapse rate and beta_p, as the climate is moved
        ! (tas - lapse*dz, pr*exp(-beta_p*lapse*dz)), so that w = 0 gives the
        ! reference climate at the model surface. The sea-level temperature
        ! tsl_ann needs no move.
        ! outside = "mean": where w < 1 the anomaly tends to the mean anomaly
        ! over the cells with w = 1 instead of zero, per month: the temperature
        ! difference dT(m) and the precipitation ratio of the totals f(m);
        ! ta_ann and tsl_ann with the mean of dT, pr_ann with the ratio of the
        ! annual totals.
        type(taper_class),       intent(in)    :: tap
        type(climate_out_class), intent(inout) :: out

        real(wp), allocatable :: dz(:,:), tas_r(:,:,:), pr_r(:,:,:), pr_ann_r(:,:)
        real(wp), allocatable :: dT(:), f(:)
        real(wp) :: dT_ann, f_ann
        integer  :: m, nx, ny, nm

        associate (now => out%now, ref => out%ref, w => tap%w, full => tap%full)

        nx = size(w, 1)
        ny = size(w, 2)
        allocate(dz(nx,ny))
        dz = 0.0_wp
        if (allocated(ref%z_srf)) dz = now%z_srf - ref%z_srf

        dT_ann = 0.0_wp
        if (allocated(now%tas) .and. allocated(ref%tas)) then
            nm = size(now%tas, 3)
            allocate(tas_r(nx,ny,nm), dT(nm))
            dT = 0.0_wp
            do m = 1, nm
                tas_r(:,:,m) = ref%tas(:,:,m) - out%lapse(m)*dz
                if (trim(tap%outside) == "mean") &
                    dT(m) = sum(now%tas(:,:,m) - tas_r(:,:,m), mask=full) / real(count(full), wp)
                now%tas(:,:,m) = tas_r(:,:,m) + w*(now%tas(:,:,m) - tas_r(:,:,m)) + (1.0_wp - w)*dT(m)
            end do
            dT_ann = sum(dT) / real(nm, wp)
            if (allocated(now%ta_ann)) &
                now%ta_ann = sum(tas_r, dim=3)/real(nm, wp) + &
                             w*(now%ta_ann - sum(tas_r, dim=3)/real(nm, wp)) + (1.0_wp - w)*dT_ann
        end if
        if (allocated(now%pr) .and. allocated(ref%pr)) then
            nm = size(now%pr, 3)
            allocate(pr_r(nx,ny,nm), f(nm))
            f = 1.0_wp
            do m = 1, nm
                pr_r(:,:,m) = ref%pr(:,:,m) * exp(-now%beta_p*out%lapse(m)*dz)
            end do
            if (trim(tap%outside) == "mean") then
                do m = 1, nm
                    if (sum(pr_r(:,:,m), mask=full) > 0.0_wp) &
                        f(m) = sum(now%pr(:,:,m), mask=full) / sum(pr_r(:,:,m), mask=full)
                end do
            end if
            f_ann = 1.0_wp
            if (trim(tap%outside) == "mean" .and. sum(sum(pr_r, dim=3), mask=full) > 0.0_wp) &
                f_ann = sum(sum(now%pr, dim=3), mask=full) / sum(sum(pr_r, dim=3), mask=full)
            do m = 1, nm
                now%pr(:,:,m) = w*now%pr(:,:,m) + (1.0_wp - w)*pr_r(:,:,m)*f(m)
            end do
            if (allocated(now%pr_ann)) then
                pr_ann_r   = sum(pr_r, dim=3)/real(nm, wp)*365.0_wp
                now%pr_ann = w*now%pr_ann + (1.0_wp - w)*pr_ann_r*f_ann
            end if
        end if
        if (allocated(now%tsl_ann) .and. allocated(ref%tsl_ann)) &
            now%tsl_ann = ref%tsl_ann + w*(now%tsl_ann - ref%tsl_ann) + (1.0_wp - w)*dT_ann

        end associate

    end subroutine taper_apply

    ! ===== bias ============================================================

    subroutine bias_init(bias, par_file, library, group, caps, domain, grid, time)
        ! Read the bias group (a named group is active), check its contract
        ! (caps_need_atm) and open its fields from the varslice groups it names
        ! (1 or 12 slices) on the climate grid, sliced at `time` when they
        ! have a time axis.
        type(bias_class),         intent(inout) :: bias
        character(len=*),         intent(in)    :: par_file, library(:), group, domain
        type(climate_caps_class), intent(in)    :: caps
        type(grid_class),         intent(in)    :: grid
        real(dp),                 intent(in)    :: time

        character(len=256) :: grp_dtas, grp_fpr
        character(len=FILE_LEN) :: gfile
        integer :: nx, ny

        bias%group  = group
        bias%active = len_trim(group) > 0
        if (.not. bias%active) return

        gfile = group_file(par_file, library, trim(group), "the bias transform")
        call nml_validate(trim(gfile), SNAP_DEFAULTS, trim(group), defaults_group="bias")
        call caps_need_atm(caps, "["//trim(group)//"]")

        call nml_read(trim(gfile), trim(group), "dtas",  grp_dtas, init=.true., &
                      defaults_file=SNAP_DEFAULTS, defaults_group="bias")
        call nml_read(trim(gfile), trim(group), "fpr",   grp_fpr,  init=.true., &
                      defaults_file=SNAP_DEFAULTS, defaults_group="bias")
        call nml_read(trim(gfile), trim(group), "f_tas", bias%f_tas, defaults_file=SNAP_DEFAULTS, defaults_group="bias")
        call nml_read(trim(gfile), trim(group), "f_pr",  bias%f_pr,  defaults_file=SNAP_DEFAULTS, defaults_group="bias")

        nx = grid%G%nx
        ny = grid%G%ny
        allocate(bias%dtas(nx,ny,12), bias%fpr(nx,ny,12))
        bias%dtas = 0.0_wp
        bias%fpr  = 1.0_wp
        bias%has_dtas = len_trim(grp_dtas) > 0
        bias%has_fpr  = len_trim(grp_fpr)  > 0
        if (bias%has_dtas) call open_field(bias%vs_dtas, trim(grp_dtas))
        if (bias%has_fpr)  call open_field(bias%vs_fpr,  trim(grp_fpr))
        bias%tdep = (bias%has_dtas .and. bias%vs_dtas%par%with_time) .or. &
                    (bias%has_fpr  .and. bias%vs_fpr%par%with_time)
        call bias_update(bias, time)

        call log_msg(LOG_INFO, "kryos_snap_mods", trim(group)//" active, dtas = "//trim(grp_dtas)// &
                     ", fpr = "//trim(grp_fpr)//trim(merge(", sliced at the forcing time", "                            ", bias%tdep)), &
                     "range dtas [K]: "//to_str(real(minval(bias%dtas), dp))//" "//to_str(real(maxval(bias%dtas), dp))// &
                     ", fpr: "//to_str(real(minval(bias%fpr), dp))//" "//to_str(real(maxval(bias%fpr), dp)))

    contains

        subroutine open_field(vs, vgroup)
            type(varslice_class), intent(inout) :: vs
            character(len=*),     intent(in)    :: vgroup
            call varslice_init_nml(vs, trim(group_file(par_file, library, vgroup, "["//trim(group)//"]")), &
                                   vgroup, domain=domain, grid_name=trim(grid%name), grid=grid)
            if (.not. vs%par%with_time) call varslice_update(vs)
        end subroutine open_field

    end subroutine bias_init

    subroutine bias_update(bias, time)
        ! The bias fields at `time`: the groups with a time axis sliced there
        ! (linear between the bracketing records, the first or last outside),
        ! then one field for all months, or 12 monthly fields.
        type(bias_class), intent(inout) :: bias
        real(dp),         intent(in)    :: time

        if (bias%has_dtas) then
            if (bias%vs_dtas%par%with_time) call varslice_update(bias%vs_dtas, [time], method="extrap")
            call fill_monthly(bias%dtas, bias%vs_dtas, "dtas")
        end if
        if (bias%has_fpr) then
            if (bias%vs_fpr%par%with_time) call varslice_update(bias%vs_fpr, [time], method="extrap")
            call fill_monthly(bias%fpr, bias%vs_fpr, "fpr")
        end if

    contains

        subroutine fill_monthly(var, vs, what)
            real(wp),             intent(out) :: var(:,:,:)
            type(varslice_class), intent(in)  :: vs
            character(len=*),     intent(in)  :: what
            integer :: m, nm
            nm = size(vs%var, 3)
            if (nm /= 1 .and. nm /= 12) then
                call fesm_error("kryos_snap_mods::bias_update", "["//trim(bias%group)//"] "//what// &
                                " must have 1 or 12 slices", "got "//to_str(nm))
            end if
            do m = 1, 12
                var(:,:,m) = vs%var(:,:,min(m,nm),1)
            end do
        end subroutine fill_monthly

    end subroutine bias_update

    subroutine bias_apply(bias, out)
        ! Correct the monthly air temperature and precipitation of `out`, and their
        ! annual means (the temperature additively, the precipitation by the
        ! ratio of the annual totals). The reference climate is unchanged.
        type(bias_class),        intent(in)    :: bias
        type(climate_out_class), intent(inout) :: out

        real(wp), allocatable :: pr_sum0(:,:)
        real(wp) :: dtas_ann(size(bias%dtas,1),size(bias%dtas,2))

        dtas_ann = sum(bias%dtas, dim=3) / 12.0_wp

        out%now%tas    = out%now%tas + bias%f_tas*bias%dtas
        out%now%ta_ann = out%now%ta_ann + bias%f_tas*dtas_ann
        if (allocated(out%now%tsl_ann)) out%now%tsl_ann = out%now%tsl_ann + bias%f_tas*dtas_ann

        pr_sum0     = sum(out%now%pr, dim=3)
        out%now%pr  = out%now%pr * (1.0_wp + bias%f_pr*(bias%fpr - 1.0_wp))
        if (allocated(out%now%pr_ann)) then
            where (pr_sum0 > 0.0_wp) out%now%pr_ann = out%now%pr_ann * sum(out%now%pr, dim=3) / pr_sum0
        end if

    end subroutine bias_apply

end module kryos_snap_mods
