module kryos_climate_snapshots
    ! climate backend snapshots (adapters-design section 4, contracts-design
    ! 2.3 and 7): yelmox_climate.f90 (snapshots branch only), snapshots.f90,
    ! climate_terms.f90, climate_mods.f90 and climate_out.f90 moved into kryos
    ! as the private modules kryos_snap_* (src/backends/climate/snapshots/).
    ! The backend computes in sp (fesm-utils wp); the exchange is dp: imports
    ! are narrowed once, exports widened from the sp values (P1), in the
    ! operand order of yelmox (yelmox_climate_update, kryos_coupling step_climate).
    !
    ! Forcing scalars: the terms name coefficients; the adapter keeps them in a
    ! private table (kryos_snap_coefs), refreshed at each update from the
    ! imports coef_<name> (the forcing's value, lag 0). A derived coefficient
    ! (from = air_anomaly) is produced here (export <name>, read by the core's
    ! climate -> forcing row) and used inside the same update, as yelmox
    ! coefs_set + terms_eval(derived=.true.).
    !
    ! The backend is stateless (the coefficients restart with the forcing, the
    ! exports with the core): no restart file. After a restore PASS_RESTART
    ! only re-derives the internal climate (for output) and logs the
    ! difference to the restored exports.

    use precision,            only: sp, dp, wp
    use fesm_log,             only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,          only: fields_get, fields_put
    use fesm_units,           only: units_conv
    use nml,                  only: nml_read, nml_has_param
    use ncio,                 only: nc_create, nc_write_dim, nc_write, nc_open, nc_close, nc_time_index
    use regions,              only: regions_select
    use kryos_component,      only: init_ctx_t, step_time_t, constraints_t, axis_t, &
                                    component_write_output, sec_year_check, NAME_LEN, PATH_LEN, &
                                    PASS_BOUNDARY, PASS_RESTART
    use kryos_climate,        only: climate_t
    use kryos_snap_out,       only: climate_out_class, clim_state_class, summer_mean
    use kryos_snap_coefs,     only: snap_coef_table_t, derived_spec, coefs_derived, coefs_is_derived, &
                                    coefs_set, coefs_value
    use kryos_snap_terms,     only: climate_caps_class, terms_class, terms_init, terms_check, terms_eval, &
                                    terms_save_ocean, terms_apply_atm, terms_apply_ocn, terms_write_step
    use kryos_snap_mods,      only: climate_mods_class, climate_mods_init, climate_mods_apply, &
                                    climate_mods_write_init, climate_mods_write_step
    use kryos_snap_snapshots, only: snapshots_class, climate_state_class, snapshots_init, snapshots_update, &
                                    snapshots_end, snapshots_write_init, snapshots_write_step, &
                                    snapshots_lapse_monthly

    implicit none

    private
    public :: snapshots_climate_t, snapshots_climate_constraints

    integer, parameter :: TERMS_MAX = 32

    ! Copy of the fields a state supplies (allocated) into the climate record.
    interface copy_allocated
        module procedure copy_allocated_1D, copy_allocated_2D, copy_allocated_3D
    end interface copy_allocated

    ! A derived coefficient this instance computes: the mean sea-level air
    ! anomaly of its climate update over a region (factor * mean(tsl_ann - ref)).
    type :: derived_air_t
        character(len=64)    :: name = ""
        real(wp)             :: factor = 1.0_wp
        logical, allocatable :: mask(:,:)
        integer              :: n = 0
    end type derived_air_t

    type, extends(climate_t) :: snapshots_climate_t
        ! Instance keys (defaults as the yelmox [forcing] group)
        character(len=NAME_LEN) :: group     = "snapshots"   ! the snapshots group
        character(len=NAME_LEN) :: reference = "pd"          ! the reference snapshot [snapshot_<name>]
        character(len=NAME_LEN) :: taper = "", bias = ""     ! transform groups, "" = off
        ! The backend
        logical                 :: south = .false.           ! southern hemisphere (seasons, lapse rates)
        type(snapshots_class)   :: snap
        type(terms_class)       :: terms
        type(climate_mods_class):: mods
        type(climate_out_class) :: out                        ! now/ref after each update (output, exports)
        type(snap_coef_table_t) :: coefs                      ! the forcing scalars the terms name
        type(derived_air_t), allocatable :: dair(:)
        real(dp)                :: pr_scale = 0.0_dp          ! mm d-1 -> kg m-2 yr-1: sec_year/86400
        real(dp), allocatable   :: xc(:), yc(:)               ! output axes [km]
    contains
        procedure :: init          => snapshots_init_
        procedure :: init_state    => snapshots_init_state_
        procedure :: update        => snapshots_update_
        procedure :: restart_write => snapshots_restart_write_
        procedure :: restart_read  => snapshots_restart_read_
        procedure :: write_output  => snapshots_write_output_
        procedure :: end           => snapshots_end_
    end type snapshots_climate_t

contains

    function snapshots_climate_constraints() result(con)
        ! Regional, any grid, several instances; stateless, so no restart file.
        type(constraints_t) :: con
        con = constraints_t(restart=.false.)
    end function snapshots_climate_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine snapshots_init_(self, ctx)
        ! climate_init of yelmox_climate (snapshots branch): the terms, the
        ! snapshots, the contracts of the terms, the derived air-anomaly
        ! coefficients, the transforms.
        class(snapshots_climate_t), intent(inout) :: self
        type(init_ctx_t),           intent(in)    :: ctx

        type(climate_caps_class) :: caps
        character(len=NAME_LEN)  :: names(TERMS_MAX)
        real(dp) :: offset

        call group_read(self, ctx, names)
        self%prec = sp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny

        ! Hemisphere: south when the climate grid lies mostly south (kryos.f90:365).
        self%south = (sum(ctx%grid%lat) / size(ctx%grid%lat) < 0.0_wp)

        call terms_init(self%terms, trim(ctx%par_file), ctx%library, names, self%coefs, trim(self%name))
        call snapshots_init(self%snap, trim(ctx%par_file), ctx%library, trim(ctx%domain), ctx%grid, ctx%time_rel, &
                            self%south, self%terms, trim(self%group), trim(self%reference))
        ! The atmosphere is required (yelmox_climate.f90:168-173).
        if (.not. self%snap%par%with_atm) then
            call fesm_error("kryos_climate_snapshots::init", "climate = snapshots needs an atmosphere", &
                            "instance "//trim(self%name)//": the reference snapshot `"//trim(self%reference)// &
                            "` binds no tas/pr")
        end if
        caps%name           = "snapshots"
        caps%atm_drives_srf = .true.
        caps%has_tsl_ann    = .true.
        caps%has_snapshots  = .true.
        caps%has_ocn_prof   = self%snap%par%with_ocn
        self%caps%has_ocn   = self%snap%par%with_ocn

        call terms_check(self%terms, caps)
        call dair_init(self, ctx)
        call climate_mods_init(self%mods, trim(ctx%par_file), ctx%library, trim(self%taper), trim(self%bias), &
                               caps, trim(ctx%domain), ctx%grid, ctx%time_rel)

        call coef_declarations(self)
        if (self%snap%par%with_ocn) call depth_axis(self)

        ! pr [mm d-1] -> [kg m-2 yr-1]: sec_year/86400 via fesm_units with the year of the
        ! run calendar (ctx%sec_year), the same call as the chion adapter (P3).
        call sec_year_check(ctx, self%name)
        call units_conv("mm d-1", "kg m-2 yr-1", self%pr_scale, offset, sec_year=ctx%sec_year)
        call log_msg(LOG_INFO, "kryos_climate_snapshots", "pr [mm d-1] -> [kg m-2 yr-1]: x "//to_str(self%pr_scale)// &
                     " (sec_year = "//to_str(ctx%sec_year)//" s / 86400 s)", "instance "//trim(self%name))

        self%xc = ctx%grid%G%x
        self%yc = ctx%grid%G%y
    end subroutine snapshots_init_

    subroutine snapshots_init_state_(self, st, pass)
        ! PASS_BOUNDARY: the climate at the start (yelmox step_climate with
        ! init=.true., kryos_startup.f90:101, 162): the same body as update
        ! (the backend is stateless; init and force change nothing for
        ! snapshots). PASS_RESTART: the exports are restored by the core; only
        ! the internal climate is re-derived (output, term values).
        class(snapshots_climate_t), intent(inout) :: self
        type(step_time_t),          intent(in)    :: st
        integer,                    intent(in)    :: pass

        select case (pass)
        case (PASS_BOUNDARY)
            call climate_step(self, st%t_eval)
            call exports_all(self, st%t1)
        case (PASS_RESTART)
            call climate_step(self, st%t_eval)
            call restore_check(self)
        end select
    end subroutine snapshots_init_state_

    subroutine snapshots_update_(self, st)
        ! climate_update (yelmox_climate.f90:262-345) at the forcing time
        ! st%t_eval = time_rel (clock = "forcing").
        class(snapshots_climate_t), intent(inout) :: self
        type(step_time_t),          intent(in)    :: st

        call climate_step(self, st%t_eval)
        call exports_all(self, st%t1)
    end subroutine snapshots_update_

    subroutine snapshots_restart_write_(self, dir, t)
        ! Stateless: no file (constraints restart = .false., the wrapper does not call this).
        class(snapshots_climate_t), intent(inout) :: self
        character(len=*),           intent(in)    :: dir
        real(dp),                   intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine snapshots_restart_write_

    subroutine snapshots_restart_read_(self, dir, t)
        class(snapshots_climate_t), intent(inout) :: self
        character(len=*),           intent(in)    :: dir
        real(dp),                   intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine snapshots_restart_read_

    subroutine snapshots_write_output_(self, file, t, stream, init)
        ! Stream 2D: <instance>.nc (yelmox snap.nc, climate_write_init and
        ! climate_write_2D). No 1D stream; any other: the default writer.
        class(snapshots_climate_t), intent(inout) :: self
        character(len=*),           intent(in)    :: file, stream
        real(dp),                   intent(in)    :: t
        logical,                    intent(in)    :: init

        integer :: ncid, n

        select case (trim(stream))
        case ("2D")
            if (init) call write_init(self, file, t)
            call nc_open(file, ncid, writable=.true.)
            n = nc_time_index(file, "time", t, ncid)
            call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)
            call write_2D(self, file, ncid, n)
            call nc_close(ncid)
        case ("1D")
            continue
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine snapshots_write_output_

    subroutine snapshots_end_(self)
        class(snapshots_climate_t), intent(inout) :: self
        call snapshots_end(self%snap)
    end subroutine snapshots_end_

    ! =========================================================================
    ! init helpers
    ! =========================================================================

    subroutine group_read(self, ctx, names)
        ! Adapter keys of the instance group; the core has validated them.
        class(snapshots_climate_t), intent(inout) :: self
        type(init_ctx_t),           intent(in)    :: ctx
        character(len=NAME_LEN),    intent(out)   :: names(TERMS_MAX)

        names = ""
        if (nml_has_param(ctx%par_file, ctx%instance, "snapshots")) &
            call nml_read(ctx%par_file, ctx%instance, "snapshots", self%group)
        if (nml_has_param(ctx%par_file, ctx%instance, "reference")) &
            call nml_read(ctx%par_file, ctx%instance, "reference", self%reference)
        if (nml_has_param(ctx%par_file, ctx%instance, "terms")) &
            call nml_read(ctx%par_file, ctx%instance, "terms", names, init=.true.)
        if (nml_has_param(ctx%par_file, ctx%instance, "taper")) &
            call nml_read(ctx%par_file, ctx%instance, "taper", self%taper)
        if (nml_has_param(ctx%par_file, ctx%instance, "bias")) &
            call nml_read(ctx%par_file, ctx%instance, "bias", self%bias)
    end subroutine group_read

    subroutine dair_init(self, ctx)
        ! The derived air-anomaly coefficients this instance computes (one per
        ! name used by a term), their region resolved on the climate grid
        ! (yelmox_climate.f90:211-244).
        class(snapshots_climate_t), intent(inout) :: self
        type(init_ctx_t),           intent(in)    :: ctx

        type(derived_spec) :: drv
        type(derived_air_t), allocatable :: d(:)
        integer :: k, n

        allocate(d(self%terms%n))
        n = 0
        do k = 1, self%terms%n
            associate (t => self%terms%term(k))
            if (len_trim(t%coef) == 0) cycle
            if (.not. coefs_is_derived(self%coefs, t%coef)) cycle
            drv = coefs_derived(self%coefs, t%coef)
            if (trim(drv%from) /= "air_anomaly") cycle
            if (any(d(1:n)%name == t%coef)) cycle
            n = n + 1
            d(n)%name   = t%coef
            d(n)%factor = drv%factor
            d(n)%mask   = regions_select(ctx%reg, trim(drv%region))
            d(n)%n      = count(d(n)%mask)
            if (d(n)%n == 0) then
                call fesm_error("kryos_climate_snapshots::init", "[coef_"//trim(t%coef)//"] region = '"// &
                                trim(drv%region)//"' selects no cell of the climate grid", &
                                "instance "//trim(self%name))
            end if
            call log_msg(LOG_INFO, "kryos_climate_snapshots", "coefficient "//trim(t%coef)//" = "// &
                         to_str(real(drv%factor, dp))//" x the mean sea-level air anomaly over '"//trim(drv%region)// &
                         "' ("//to_str(d(n)%n)//" of "//to_str(size(d(n)%mask))//" cells), set by this climate", &
                         "instance "//trim(self%name))
            end associate
        end do
        self%dair = d(1:n)
    end subroutine dair_init

    subroutine coef_declarations(self)
        ! The forcing scalars for the component layer (imports coef_<name>,
        ! exports of the derived ones).
        class(snapshots_climate_t), intent(inout) :: self
        integer :: k, n

        n = size(self%coefs%coef)
        allocate(self%coef_names(n), self%coef_units(n))
        do k = 1, n
            self%coef_names(k) = self%coefs%coef(k)%name
            self%coef_units(k) = self%coefs%coef(k)%units
        end do
        allocate(self%derived_names(size(self%dair)), self%derived_units(size(self%dair)))
        do k = 1, size(self%dair)
            self%derived_names(k) = self%dair(k)%name
            self%derived_units(k) = self%coefs%coef(index_of(self, self%dair(k)%name))%units
        end do
    end subroutine coef_declarations

    integer function index_of(self, name) result(k)
        class(snapshots_climate_t), intent(in) :: self
        character(len=*),           intent(in) :: name
        do k = 1, size(self%coefs%coef)
            if (trim(self%coefs%coef(k)%name) == trim(name)) return
        end do
        k = 0
    end function index_of

    subroutine depth_axis(self)
        ! The producer's ocean axis `depth` [m, positive down] (contracts-design
        ! 1, A8): the snapshots' model levels, widened.
        class(snapshots_climate_t), intent(inout) :: self
        type(axis_t) :: ax

        ax%name  = "depth"
        ax%units = "m"
        allocate(ax%v(size(self%snap%ref%depth)))
        ax%v = real(self%snap%ref%depth, dp)
        if (.not. allocated(self%axes)) allocate(self%axes(0))
        self%axes = [self%axes, ax]
    end subroutine depth_axis

    ! =========================================================================
    ! The climate update (yelmox_climate.f90:262-345, snapshots branch)
    ! =========================================================================

    subroutine climate_step(self, time)
        ! terms_eval -> snapshots_update -> out -> atmosphere terms -> derived
        ! coefficients -> ocean terms -> transforms. `time` = forcing time.
        class(snapshots_climate_t), intent(inout) :: self
        real(dp),                   intent(in)    :: time

        real(wp), allocatable :: z_srf(:,:)
        integer :: k

        allocate(z_srf(self%nx, self%ny))
        call fields_get(self%imports, self%hi_z_srf, z_srf)
        call coefs_from_imports(self)

        ! The value of every term's coefficient for this update.
        call terms_eval(self%terms, self%coefs)

        call snapshots_update(self%snap, z_srf=z_srf, time=time, terms=self%terms)
        call state_copy(self%out%now, self%snap%now)
        call state_copy(self%out%ref, self%snap%ref)

        ! The elevation dependence of the atmosphere (the wiring moves it to another surface).
        self%out%now%z_srf = z_srf
        if (.not. allocated(self%out%now%beta_p)) allocate(self%out%now%beta_p(size(z_srf,1), size(z_srf,2)))
        self%out%lapse      = snapshots_lapse_monthly(self%snap)
        self%out%now%beta_p = self%snap%now%beta_p

        ! The terms: the snapshots' ocean kept as the base of the ocean terms; the
        ! uniform terms on the atmosphere; the derived air anomaly of this update
        ! and the terms it drives re-evaluated; the ocean terms; the transforms.
        call terms_save_ocean(self%terms, self%out)
        call terms_apply_atm(self%terms, self%out)
        do k = 1, size(self%dair)
            call coefs_set(self%coefs, self%dair(k)%name, self%dair(k)%factor * &
                 sum(self%out%now%tsl_ann - self%out%ref%tsl_ann, mask=self%dair(k)%mask) / real(self%dair(k)%n, wp))
        end do
        call terms_eval(self%terms, self%coefs, derived=.true.)
        call terms_apply_ocn(self%terms, self%out)
        call climate_mods_apply(self%mods, self%out, time)
    end subroutine climate_step

    subroutine coefs_from_imports(self)
        ! The forcing scalars of this update: the imports coef_<name>, narrowed
        ! (P1; coef_class%f is sp). The derived ones carry the value last set
        ! (the forcing's lag-1 copy of this climate's export).
        class(snapshots_climate_t), intent(inout) :: self
        integer :: k
        do k = 1, size(self%coefs%coef)
            call fields_get(self%imports, self%hi_coef(k), self%coefs%coef(k)%f)
        end do
    end subroutine coefs_from_imports

    subroutine state_copy(dst, src)
        ! The fields of a snapshots state that the climate carries, those it built
        ! (yelmox_climate snapshots_state_copy).
        type(clim_state_class),    intent(inout) :: dst
        type(climate_state_class), intent(in)    :: src

        call copy_allocated(dst%z_srf,   src%z_srf)
        call copy_allocated(dst%tas,     src%tas)
        call copy_allocated(dst%pr,      src%pr)
        call copy_allocated(dst%tsl_ann, src%tsl_ann)
        call copy_allocated(dst%ta_ann,  src%ta_ann)
        call copy_allocated(dst%pr_ann,  src%pr_ann)
        call copy_allocated(dst%to_ann,  src%to_ann)
        call copy_allocated(dst%so_ann,  src%so_ann)
        call copy_allocated(dst%depth,   src%depth)
    end subroutine state_copy

    subroutine copy_allocated_1D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:)
        real(wp), allocatable, intent(in)    :: src(:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_1D

    subroutine copy_allocated_2D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:,:)
        real(wp), allocatable, intent(in)    :: src(:,:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_2D

    subroutine copy_allocated_3D(dst, src)
        real(wp), allocatable, intent(inout) :: dst(:,:,:)
        real(wp), allocatable, intent(in)    :: src(:,:,:)
        if (allocated(src)) dst = src
    end subroutine copy_allocated_3D

    ! =========================================================================
    ! Exports (adapters-design 4.4, contracts-design 2.3): INSTANT at t1
    ! =========================================================================

    subroutine exports_all(self, t1)
        class(snapshots_climate_t), intent(inout) :: self
        real(dp),                   intent(in)    :: t1

        real(dp), allocatable :: lapse3(:,:,:)
        integer :: k, m

        associate (now => self%out%now, ref => self%out%ref)
            call fields_put(self%exports, self%h_tas, now%tas, t1=t1)
            ! pr [mm d-1] -> [kg m-2 yr-1]: one dp product; the surface divides by the same factor (P3).
            call fields_put(self%exports, self%h_pr, real(now%pr, dp) * self%pr_scale, t1=t1)
            call fields_put(self%exports, self%h_z_srf_clim, now%z_srf, t1=t1)

            allocate(lapse3(self%nx, self%ny, 12))
            do m = 1, 12
                lapse3(:,:,m) = real(self%out%lapse(m), dp)
            end do
            call fields_put(self%exports, self%h_lapse, lapse3, t1=t1)
            call fields_put(self%exports, self%h_beta_p, now%beta_p, t1=t1)

            if (self%caps%has_ocn) then
                call fields_put(self%exports, self%h_to_ann, now%to_ann, t1=t1)
                call fields_put(self%exports, self%h_so_ann, now%so_ann, t1=t1)
                ! now - ref in sp on the climate grid, before any remap (coupling:612)
                call fields_put(self%exports, self%h_dto_ann, now%to_ann - ref%to_ann, t1=t1)
            end if
        end associate

        do k = 1, size(self%dair)
            call fields_put(self%exports, self%h_derived(k), coefs_value(self%coefs, self%dair(k)%name), t1=t1)
        end do
    end subroutine exports_all

    subroutine restore_check(self)
        ! After a restore: the re-derived atmosphere against the restored exports
        ! (equal when the restart time lies on the climate cadence).
        class(snapshots_climate_t), intent(inout) :: self

        real(dp), allocatable :: v(:,:,:)
        real(dp) :: dtas, dpr

        allocate(v(self%nx, self%ny, 12))
        call fields_get(self%exports, self%h_tas, v)
        dtas = maxval(abs(v - real(self%out%now%tas, dp)))
        call fields_get(self%exports, self%h_pr, v)
        dpr = maxval(abs(v - real(self%out%now%pr, dp) * self%pr_scale))
        call log_msg(LOG_INFO, "kryos_climate_snapshots", "restore: re-derived climate against the restored exports: "// &
                     "max |tas| = "//to_str(dtas)//" K, max |pr| = "//to_str(dpr)//" kg m-2 yr-1", &
                     "instance "//trim(self%name))
    end subroutine restore_check

    ! =========================================================================
    ! Output (yelmox_climate climate_write_init, climate_write_2D)
    ! =========================================================================

    subroutine write_init(self, filename, time)
        ! Create the 2D file: xc, yc, time (kryos_output io_dims_init), the
        ! ocean depth axis, the month axis, the transforms' static fields.
        class(snapshots_climate_t), intent(in) :: self
        character(len=*),           intent(in) :: filename
        real(dp),                   intent(in) :: time

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=self%xc, units="km")
        call nc_write_dim(filename, "yc", x=self%yc, units="km")
        call nc_write_dim(filename, "time", x=time, dx=1.0_dp, nx=1, units="year", unlimited=.true.)

        if (self%snap%par%with_ocn) call nc_write_dim(filename, "depth", x=self%snap%ref%depth, units="m")
        call snapshots_write_init(self%snap, filename)
        call climate_mods_write_init(self%mods, filename)
    end subroutine write_init

    subroutine write_2D(self, filename, ncid, n)
        ! One record, as applied (after the transforms): the annual and summer
        ! near-surface air temperature, the annual precipitation, the annual
        ! ocean profiles, the monthly atmosphere, the snapshots' elevation, the
        ! transforms' diagnostics and the value of every term's coefficient.
        class(snapshots_climate_t), intent(in) :: self
        character(len=*),           intent(in) :: filename
        integer,                    intent(in) :: ncid, n

        associate (out => self%out)

        call write2D("t2m_ann", out%now%ta_ann, "K", "Near-surface air temperature (ann)", &
                     prec_abs=1e-3_wp)
        if (allocated(out%now%tas)) &
            call write2D("t2m_sum", summer_mean(out%now%tas, self%south), "K", &
                         "Near-surface air temperature (sum)", prec_abs=1e-3_wp)
        if (allocated(out%now%pr_ann)) &
            call write2D("pr_ann", out%now%pr_ann, "mm/a", "Precipitation (ann)", &
                         prec_rel=1e-4_wp)
        if (allocated(out%now%to_ann)) &
            call write3D("to_ann", out%now%to_ann, "depth", "K", "Ocean temperature (ann)", &
                         prec_abs=1e-3_wp)
        if (allocated(out%now%so_ann)) &
            call write3D("so_ann", out%now%so_ann, "depth", "PSU", "Ocean salinity (ann)", &
                         prec_rel=1e-4_wp)

        ! The monthly atmosphere, then the snapshots' own state.
        call write3D("tas",     out%now%tas,     "month", "K",    "Near-surface air temperature", &
                     prec_abs=1e-3_wp)
        call write3D("pr",      out%now%pr,      "month", "mm/d", "Precipitation", &
                     prec_rel=1e-4_wp)
        call write2D("tsl_ann", out%now%tsl_ann, "K", "Sea-level air temperature (ann)", &
                     prec_abs=1e-3_wp)
        call snapshots_write_step(self%snap, filename, ncid, n)

        end associate

        call climate_mods_write_step(self%mods, filename, ncid, n)
        call terms_write_step(self%terms, filename, ncid, n)

    contains

        ! prec_abs/prec_rel: precision of the written values (ncio bit
        ! rounding), absent for lossless.
        subroutine write2D(vnm, var, units, long_name, prec_abs, prec_rel)
            character(len=*), intent(in) :: vnm, units, long_name
            real(wp),         intent(in) :: var(:,:)
            real(wp),         intent(in), optional :: prec_abs, prec_rel
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1),size(var,2),1], ncid=ncid, &
                          units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end subroutine write2D

        subroutine write3D(vnm, var, dim3, units, long_name, prec_abs, prec_rel)
            character(len=*), intent(in) :: vnm, dim3, units, long_name
            real(wp),         intent(in) :: var(:,:,:)
            real(wp),         intent(in), optional :: prec_abs, prec_rel
            call nc_write(filename, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                          start=[1,1,1,n], count=[size(var,1),size(var,2),size(var,3),1], &
                          ncid=ncid, units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end subroutine write3D

    end subroutine write_2D

end module kryos_climate_snapshots
