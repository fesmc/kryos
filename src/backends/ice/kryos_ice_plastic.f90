module kryos_ice_plastic
    ! ice backend `plastic` (docs/dev/plastic-backend.md): a perfectly-plastic ice sheet, the library
    ! PlasticIceSheet (module plastic_ice), as Gowan et al. (2021) built PaleoMIST, run online: at every
    ! update the ice mask of the margin chronology at t1 (kryos_plastic_margin), the basal shear stress
    ! at t1 (kryos_plastic_taub) and the imported bed and sea level give the thickness (plastic_solve,
    ! warm-started from the last solution), or, with drive = "smb", the imported smb does
    ! (plastic_solve_smb). Built with plastic=1; otherwise kryos_ice_plastic_stub.f90 provides this module.
    ! Everything is dp, as the exchange.
    !
    ! The ice is imposed, not integrated: dHidt = (H - H_prev)/dt and mb_net = rho_i dHidt, the mass
    ! change the margin and the bed impose, so the sea-level budget closes. smb = bmb = cmb = 0, or with
    ! mass_balance = "implied" the mass balance the geometry needs (plastic_smb_flux_divergence: dHidt
    ! + div q of the steady solver's face fluxes; the outflow to land ablated in the margin cells, the
    ! outflow to the ocean calved), with the surface velocity of the same closure (ux_s, uy_s). The ice
    ! is grounded by construction (f_grnd = 1 on it). Optional imports: the imports the backend does not read (T_srf, bmb_shlf, and smb with drive = "taub") may stay unwired; so may
    ! z_bed and z_sl (a domain without geo): until they are published the backend takes the reference
    ! bed and z_sl = 0.

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, LOG_WARN, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use ncio,            only: nc_read, nc_write, nc_size, nc_close
    use nml,             only: nml_read, nml_has_param, nml_has_group, nml_group_keys, nml_key_len, nml_replace
    use phys_constants,  only: phys_const_isset
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, component_write_output, component_optional, &
                               NAME_LEN, PATH_LEN, PASS_INIT, PASS_GEOMETRY, PASS_SPINUP, PASS_RESTART
    use kryos_outfile,   only: outfile_dims_2d, outfile_open, outfile_var2d
    use kryos_ice,       only: ice_t, ice_opt_t
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, margin_footprint, margin_interp_code
    use kryos_plastic_taub,   only: taub_t, taub_method_code, TAUB_M_CONST, TAUB_M_FILE, TAUB_M_DOMAINS, taub_set_const, &
                                    taub_read_field, taub_read_domains, taub_read_scale, taub_at, &
                                    taub_read_domain_factor, taub_set_smooth
    use kryos_plastic_friction, only: friction_t, friction_law, friction_map
    use plastic_ice,     only: plastic_par_t, plastic_t, plastic_glen_t, plastic_sliding_t, plastic_init, &
                               plastic_solve, plastic_solve_smb, plastic_end, plastic_smb_flux_divergence, &
                               plastic_stress_residual, plastic_margin_frac_from_sdist, PLASTIC_OK, &
                               PLASTIC_NOT_CONVERGED

    implicit none

    private
    public :: plastic_ice_t, plastic_ice_constraints

    character(len=*), parameter :: TABLE = "ice-plastic"      ! input/kryos-variables-ice-plastic.md
    real(dp), parameter :: KG_TO_GT = 1.0e-12_dp

    type, extends(ice_t) :: plastic_ice_t
        type(plastic_t)         :: ps
        type(plastic_par_t)     :: par
        type(plastic_glen_t)    :: rheo
        type(plastic_sliding_t) :: slide
        type(margin_t)          :: mg
        type(taub_t)            :: tb
        character(len=8)  :: drive = "taub"              ! "taub": plastic_solve | "smb": plastic_solve_smb
        logical           :: closure = .false.           ! rheology and sliding set (drive = "smb", mass_balance = "implied")
        logical           :: implied = .false.           ! mass_balance = "implied": smb, cmb, ux_s, uy_s from the closure
        real(dp)          :: mb_smooth = 0.0_dp          ! [m] mb_smooth_km: Gaussian sd of the conservative smoothing (0: none)
        type(friction_t)  :: fr                          ! slide_law "yield", "coulomb_n": the coefficient per cell
        real(dp) :: dx = 0.0_dp, dy = 0.0_dp             ! [m]
        real(dp), allocatable :: z_bed_ref(:,:), H_ice_pd(:,:), footprint(:,:)   ! static
        real(dp), allocatable :: z_bed(:,:), z_sl(:,:)   ! [m] the bed and sea level of the last solve
        real(dp), allocatable :: tau_b(:,:), sd(:,:)     ! [Pa], [m] its basal shear stress and signed distance
        logical,  allocatable :: mask(:,:)               ! its ice mask
        real(dp), allocatable :: mfrac(:,:)              ! margin_position = "sdist": the margin's distance from the
                                                         ! ice-free cell centres [fraction of the spacing]; else
                                                         ! unallocated (absent in the library calls: par%margin)
        real(dp), allocatable :: H_prev(:,:)             ! [m] the thickness at the previous update (dHidt)
        real(dp), allocatable :: divq(:,:), q_land(:,:), q_marine(:,:)   ! [m ice yr-1] (implied): flux divergence,
                                                         ! the outflow through the margin faces to land and to the ocean
        real(dp), allocatable :: us_x(:,:), us_y(:,:)    ! [m yr-1] (implied) surface velocity of the closure
        logical :: bed_ref = .false., sl_zero = .false.  ! logged once: the bed or sea level was not imported
    contains
        procedure :: init          => plastic_init_
        procedure :: init_state    => plastic_init_state_
        procedure :: update        => plastic_update_
        procedure :: flotation     => plastic_flotation_
        procedure :: restart_write => plastic_restart_write_
        procedure :: restart_read  => plastic_restart_read_
        procedure :: write_output  => plastic_write_output_
        procedure :: end           => plastic_end_
    end type plastic_ice_t

contains

    function plastic_ice_constraints() result(con)
        ! Any grid, several instances, regional, with a restart file.
        type(constraints_t) :: con
        con = constraints_t(build_flag="plastic=1")
    end function plastic_ice_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine plastic_init_(self, ctx)
        class(plastic_ice_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        integer :: stat

        self%prec = dp
        self%caps%has_flotation = .true.
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%dx = ctx%grid%G%dx*ctx%grid%cs%xy_conv
        self%dy = ctx%grid%G%dy*ctx%grid%cs%xy_conv

        call group_read(self, ctx)
        self%opt = ice_opt_t(z_bed_ref=.true., H_ice_ref=.true., footprint=.true., H_ice_pd=.true., H_grnd=.true., &
                             uxy_s=self%implied, mb_tot=self%implied)
        call component_optional(self, "z_bed")
        call component_optional(self, "z_sl")
        call component_optional(self, "T_srf")
        call component_optional(self, "bmb_shlf")
        if (trim(self%drive) /= "smb") call component_optional(self, "smb")

        if (.not. (allocated(ctx%ref%z_bed) .and. allocated(ctx%ref%H_ice))) then
            call fesm_error("kryos_ice_plastic::init", "reference topography incomplete", "instance "// &
                            trim(self%name)//": z_bed and H_ice are needed on grid `"//trim(self%grid)//"`")
        end if
        self%z_bed_ref = ctx%ref%z_bed
        self%H_ice_pd  = ctx%ref%H_ice
        call constants_check(self, ctx)

        self%footprint = margin_footprint(self%mg)
        allocate(self%z_bed(self%nx, self%ny), self%z_sl(self%nx, self%ny), self%tau_b(self%nx, self%ny), &
                 self%sd(self%nx, self%ny), self%mask(self%nx, self%ny), self%H_prev(self%nx, self%ny))
        self%z_bed  = self%z_bed_ref
        self%z_sl   = 0.0_dp
        self%H_prev = 0.0_dp
        if (self%implied) then
            allocate(self%divq(self%nx, self%ny), self%q_land(self%nx, self%ny), self%q_marine(self%nx, self%ny), &
                     self%us_x(self%nx, self%ny), self%us_y(self%nx, self%ny))
            self%divq = 0.0_dp;  self%q_land = 0.0_dp;  self%q_marine = 0.0_dp;  self%us_x = 0.0_dp;  self%us_y = 0.0_dp
        end if

        call plastic_init(self%ps, self%nx, self%ny, self%dx, self%dy, self%par, stat=stat)
        call status_check(self, stat, "plastic_init")
        call log_msg(LOG_INFO, "kryos_ice_plastic", "drive "//trim(self%drive)//", mode "//trim(self%par%mode)// &
                     ", rho_i "//to_str(self%par%rho_i)//", rho_w "//to_str(self%par%rho_w)//", footprint "// &
                     to_str(count(self%footprint > 0.5_dp))//" of "//to_str(self%nx*self%ny)//" cells", &
                     "instance "//trim(self%name))
    end subroutine plastic_init_

    subroutine plastic_init_state_(self, st, pass)
        ! PASS_INIT: the first solve at time_init on the reference bed (z_sl = 0; the imports are not
        ! published yet) and every INSTANT export, the static ones included. PASS_GEOMETRY (cold start,
        ! after geo; twice with a host-resolved geo): the solve on the bed and sea level geo gave, and
        ! the geometry. PASS_SPINUP: the MEAN exports of the virtual previous step (no change: zero).
        ! PASS_RESTART (after restart_read, which restores H, w, H_prev and z_srf, and the core's restore
        ! of the imports and exports): the inputs of the restored time (bed, sea level, mask, tau_b) and
        ! the diagnostics (ux_s, uy_s); no state export is put. Other passes: nothing.
        class(plastic_ice_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        integer,              intent(in)    :: pass

        select case (pass)
        case (PASS_INIT)
            call solve(self, st%t1, warm=.false.)
            self%H_prev = self%ps%H
            call exports_static(self, st%t1)
            call exports_geometry(self, st%t1)
            call exports_mass(self, st%t1)
        case (PASS_GEOMETRY)
            call boundary(self)
            call solve(self, st%t1, warm=.true.)
            self%H_prev = self%ps%H
            call exports_geometry(self, st%t1)
            call exports_mass(self, st%t1)
        case (PASS_SPINUP)
            call exports_rates(self, st)
        case (PASS_RESTART)
            call boundary(self)
            call mask_at(self, st%t1)
            if (trim(self%drive) == "taub") then
                call forcing_taub(self, st%t1)
            else
                self%tau_b = driving_stress(self)
            end if
            call diagnostics(self)
            call exports_velocity(self, st%t1)
        end select
    end subroutine plastic_init_state_

    subroutine plastic_update_(self, st)
        ! The thickness at t1 on the imported bed and sea level; the rates over [t0, t1].
        class(plastic_ice_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        call boundary(self)
        call solve(self, st%t1, warm=.true.)
        call exports_geometry(self, st%t1)
        call exports_mass(self, st%t1)
        call exports_rates(self, st)
        self%H_prev = self%ps%H
    end subroutine plastic_update_

    subroutine plastic_flotation_(self, H, z_bed, z_sl, f_grnd)
        ! The backend's rule: grounded where rho_i H >= rho_w (z_sl - z_bed), with its own densities.
        class(plastic_ice_t), intent(in)  :: self
        real(dp),             intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),             intent(out) :: f_grnd(:,:)
        f_grnd = merge(1.0_dp, 0.0_dp, self%par%rho_i*H >= self%par%rho_w*(z_sl - z_bed))
    end subroutine plastic_flotation_

    subroutine plastic_restart_write_(self, dir, t)
        ! <instance>_restart.nc: the thickness, the warm start (w = H**2, the mask it belongs to, whether it
        ! is usable: the library's warm-start rule), the previous thickness (dHidt of the next step) and
        ! z_srf, lossless (no precision given).
        class(plastic_ice_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t

        character(len=:), allocatable :: file

        file = restart_file(self, dir)
        call outfile_dims_2d(file, self%xkm, self%ykm, t)
        call nc_write(file, "H",      self%ps%H,     dim1="xc", dim2="yc", units="m",  long_name="Ice thickness")
        call nc_write(file, "w",      self%ps%w,     dim1="xc", dim2="yc", units="m2", long_name="H**2, warm start")
        call nc_write(file, "H_prev", self%H_prev,   dim1="xc", dim2="yc", units="m",  &
                      long_name="Ice thickness at the previous update")
        call nc_write(file, "z_srf",  self%ps%z_srf, dim1="xc", dim2="yc", units="m",  long_name="Surface elevation")
        call nc_write(file, "mask_prev", merge(1, 0, self%ps%wk%mask_prev), dim1="xc", dim2="yc", units="1", &
                      long_name="Ice mask of the stored w (warm start)")
        call nc_write(file, "has_warm", [merge(1, 0, self%ps%has_warm)], dim1="time", units="1", &
                      long_name="w is a usable warm start")
    end subroutine plastic_restart_write_

    subroutine plastic_restart_read_(self, dir, t)
        class(plastic_ice_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t

        character(len=*), parameter :: proc = "kryos_ice_plastic::restart_read"
        character(len=:), allocatable :: file
        real(dp) :: tf(1)
        integer, allocatable :: m(:,:)
        integer  :: hw(1)
        logical  :: exists

        file = restart_file(self, dir)
        inquire(file=file, exist=exists)
        if (.not. exists) call fesm_error(proc, "restart file missing", "instance "//trim(self%name)//": "//file)
        if (nc_size(file, "xc") /= self%nx .or. nc_size(file, "yc") /= self%ny) then
            call fesm_error(proc, "the restart file is not on the instance grid", "instance "//trim(self%name)// &
                            ": "//file)
        end if
        call nc_read(file, "time", tf)
        if (abs(tf(1) - t) > 1.0e-6_dp*max(1.0_dp, abs(t))) then
            call fesm_error(proc, "the restart file is of another time", "instance "//trim(self%name)//": "//file// &
                            " at t = "//to_str(tf(1))//" yr, the restore at "//to_str(t)//" yr")
        end if
        call nc_read(file, "H",      self%ps%H)
        call nc_read(file, "w",      self%ps%w)
        call nc_read(file, "H_prev", self%H_prev)
        call nc_read(file, "z_srf",  self%ps%z_srf)
        allocate(m(self%nx, self%ny))
        call nc_read(file, "mask_prev", m)
        self%ps%wk%mask_prev = m == 1
        call nc_read(file, "has_warm", hw)
        self%ps%has_warm = hw(1) == 1
    end subroutine plastic_restart_read_

    subroutine plastic_write_output_(self, file, t, stream, init)
        ! Stream `plastic` (out_plastic, <instance>_plastic.nc): the state of the last solve and its
        ! inputs (table ice-plastic), its convergence and stress residual (0-D); any other: the default
        ! writer.
        class(plastic_ice_t), intent(inout) :: self
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init

        integer  :: ncid, n
        real(dp) :: rmed, rmax

        select case (trim(stream))
        case ("plastic")
            if (init) call outfile_dims_2d(file, self%xkm, self%ykm, t)
            call outfile_open(file, t, ncid, n)
            call outfile_var2d(file, "H_ice", self%ps%H,  n, ncid, "m",  "Ice thickness", table=TABLE)
            call outfile_var2d(file, "z_srf", surface(self), n, ncid, "m", "Surface elevation", table=TABLE)
            call outfile_var2d(file, "z_bed", self%z_bed, n, ncid, "m",  "Bedrock elevation", table=TABLE)
            call outfile_var2d(file, "z_sl",  self%z_sl,  n, ncid, "m",  "Sea level", table=TABLE)
            call outfile_var2d(file, "mask",  merge(1.0_dp, 0.0_dp, self%mask), n, ncid, "1", "Ice mask", table=TABLE)
            call outfile_var2d(file, "sdist", self%sd,    n, ncid, "m",  "Signed distance to the margin", table=TABLE)
            call outfile_var2d(file, "tau_b", self%tau_b, n, ncid, "Pa", "Basal shear stress", table=TABLE)
            if (self%tb%method == TAUB_M_DOMAINS) &
                call outfile_var2d(file, "taub_domain", real(self%tb%dom, dp), n, ncid, "1", "Basal shear stress domain", &
                                   table=TABLE)
            if (self%implied) then
                call outfile_var2d(file, "div_q", self%divq, n, ncid, "m yr-1", &
                                   "Divergence of the ice flux (ice equivalent)", table=TABLE)
                call outfile_var2d(file, "q_land", self%q_land, n, ncid, "m yr-1", &
                                   "Outflow through the margin faces to land (ice equivalent)", table=TABLE)
                call outfile_var2d(file, "q_marine", self%q_marine, n, ncid, "m yr-1", &
                                   "Outflow through the margin faces to the ocean (ice equivalent)", table=TABLE)
                call nc_write(file, "abl_margin", -self%par%rho_i*sum(self%q_land)*self%dx*self%dy*KG_TO_GT, &
                              dim1="time", start=[n], count=[1], ncid=ncid, units="Gt yr-1", &
                              long_name="Ablation of the outflow to land at the margin (last solve)")
                call nc_write(file, "discharge", self%par%rho_i*sum(self%q_marine)*self%dx*self%dy*KG_TO_GT, &
                              dim1="time", start=[n], count=[1], ncid=ncid, units="Gt yr-1", &
                              long_name="Outflow to the ocean at the margin (last solve)")
            end if
            call nc_write(file, "converged", merge(1, 0, self%ps%converged), dim1="time", start=[n], count=[1], &
                          ncid=ncid, units="1", long_name="The last solve converged")
            call nc_write(file, "n_outer_used", self%ps%n_outer_used, dim1="time", start=[n], count=[1], ncid=ncid, &
                          units="1", long_name="Outer iterations of the last solve")
            call stress_residual(self, rmed, rmax)
            call nc_write(file, "stress_res_median", rmed, dim1="time", start=[n], count=[1], ncid=ncid, units="1", &
                          long_name="Median over the ice of |rho_i g H |grad z_srf| - tau_b| / tau_b (last solve)")
            call nc_write(file, "stress_res_max", rmax, dim1="time", start=[n], count=[1], ncid=ncid, units="1", &
                          long_name="Maximum over the ice of |rho_i g H |grad z_srf| - tau_b| / tau_b (last solve)")
            call nc_close(ncid)
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine plastic_write_output_

    subroutine plastic_end_(self)
        class(plastic_ice_t), intent(inout) :: self
        call plastic_end(self%ps)
    end subroutine plastic_end_

    ! =========================================================================
    ! The solve
    ! =========================================================================

    subroutine boundary(self)
        ! The bed and sea level of the imports once published; the reference bed and z_sl = 0 until then
        ! (a domain without geo: always), logged once.
        class(plastic_ice_t), intent(inout) :: self

        if (published(self, self%hi_z_bed)) then
            call fields_get(self%imports, self%hi_z_bed, self%z_bed)
        else
            self%z_bed = self%z_bed_ref
            if (.not. self%bed_ref) call log_msg(LOG_INFO, "kryos_ice_plastic", "z_bed is not imported: the "// &
                                                 "reference bed", "instance "//trim(self%name))
            self%bed_ref = .true.
        end if
        if (published(self, self%hi_z_sl)) then
            call fields_get(self%imports, self%hi_z_sl, self%z_sl)
        else
            self%z_sl = 0.0_dp
            if (.not. self%sl_zero) call log_msg(LOG_INFO, "kryos_ice_plastic", "z_sl is not imported: 0 m", &
                                                 "instance "//trim(self%name))
            self%sl_zero = .true.
        end if
    end subroutine boundary

    subroutine solve(self, t, warm)
        ! The plastic thickness at model time t on self%z_bed, self%z_sl: the mask and tau_b at t, then
        ! plastic_solve (or plastic_solve_smb with the imported smb). Non-convergence is a warning.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: warm

        real(dp), allocatable :: smb(:,:)
        integer :: stat

        call mask_at(self, t)
        select case (trim(self%drive))
        case ("taub")
            call forcing_taub(self, t)
            call plastic_solve(self%ps, self%z_bed, self%z_sl, self%tau_b, self%mask, warm=warm, stat=stat, &
                               margin_frac=self%mfrac)
        case ("smb")
            if (.not. published(self, self%hi_smb)) then
                call fesm_error("kryos_ice_plastic::solve", "drive = 'smb' needs the smb import", "instance "// &
                                trim(self%name)//": no smb was published (a surface instance in the domain)")
            end if
            allocate(smb(self%nx, self%ny))
            call fields_get(self%imports, self%hi_smb, smb)
            smb = smb/self%par%rho_i                    ! kg m-2 yr-1 -> m ice yr-1
            call plastic_solve_smb(self%ps, self%z_bed, self%z_sl, smb, self%mask, self%rheo, self%slide, warm=warm, &
                                   stat=stat, margin_frac=self%mfrac)
            self%tau_b = driving_stress(self)
        end select
        if (stat == PLASTIC_NOT_CONVERGED) then
            call log_msg(LOG_WARN, "kryos_ice_plastic", "the solve at t = "//to_str(t)//" yr did not converge ("// &
                         to_str(self%ps%n_outer_used)//" outer iterations, last change "// &
                         to_str(self%ps%outer_change)//" m): "//trim(self%ps%msg), "instance "//trim(self%name))
        else
            call status_check(self, stat, "the solve at t = "//to_str(t)//" yr")
        end if
        call diagnostics(self)
    end subroutine solve

    subroutine mask_at(self, t)
        ! The ice mask and signed distance at t (margin chronology) and, with margin_position = "sdist",
        ! the margin's position on the ice-free cells next to the ice: frac = -sd / sqrt(dx dy) (library
        ! helper). The backend's sd at a slice's time is the distance to the staircase of its cell faces
        ! (centre-to-centre distance to the nearest opposite cell less half a cell), so frac = 0.5 there
        ! ("face"); between the slices the interpolated sd moves the margin within the cells.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t
        call margin_at(self%mg, t, self%mask, self%sd)
        if (allocated(self%mfrac)) call plastic_margin_frac_from_sdist(self%sd, self%mask, self%dx, self%dy, self%mfrac)
    end subroutine mask_at

    subroutine forcing_taub(self, t)
        ! tau_b at t (drive = "taub"; with "smb" the solve derives it): positive on every ice cell.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t

        integer :: n

        call taub_at(self%tb, t, self%tau_b)
        n = count(self%mask .and. .not. self%tau_b > 0.0_dp)
        if (n > 0) then
            call fesm_error("kryos_ice_plastic::solve", to_str(n)//" ice cells have no positive tau_b at t = "// &
                            to_str(t)//" yr", "instance "//trim(self%name)//": cells outside every tau_b domain "// &
                            "take taub_default (0 = none)")
        end if
    end subroutine forcing_taub

    subroutine diagnostics(self)
        ! mass_balance = "implied": on the current geometry, the divergence of the steady solver's face fluxes,
        ! the outflow through the margin faces (to land, to the ocean) and the surface velocity of the closure.
        ! drive = "taub": the driving stress of the fluxes at most tau_b (the bed's yield stress; the margin faces'
        ! step to the margin surface would otherwise drive the outflow far above it).
        class(plastic_ice_t), intent(inout) :: self
        integer :: stat
        if (.not. self%implied) return
        if (len_trim(self%fr%kind) > 0) call friction_map_set(self)
        if (trim(self%drive) == "taub") then
            call plastic_smb_flux_divergence(self%ps, self%z_bed, self%z_sl, self%mask, self%rheo, self%slide, self%divq, &
                                             margin_frac=self%mfrac, stat=stat, q_land=self%q_land, &
                                             q_marine=self%q_marine, us_x=self%us_x, us_y=self%us_y, tau_cap=self%tau_b)
        else
            call plastic_smb_flux_divergence(self%ps, self%z_bed, self%z_sl, self%mask, self%rheo, self%slide, self%divq, &
                                             margin_frac=self%mfrac, stat=stat, q_land=self%q_land, &
                                             q_marine=self%q_marine, us_x=self%us_x, us_y=self%us_y)
        end if
        call status_check(self, stat, "plastic_smb_flux_divergence")
    end subroutine diagnostics

    subroutine friction_map_set(self)
        ! The sliding coefficient per cell of slide_law "yield" or "coulomb_n" (kryos_plastic_friction) on the current
        ! geometry and tau_b.
        class(plastic_ice_t), intent(inout) :: self
        call friction_map(self%fr, self%slide%m, self%tau_b, self%ps%H, self%z_bed, self%z_sl, self%mask, &
                          self%par%rho_i, self%par%rho_w, self%par%g, self%slide%fric_map)
    end subroutine friction_map_set

    subroutine status_check(self, stat, what)
        ! A library call that did not return PLASTIC_OK (bad input: nothing computed) stops; a closure's
        ! PLASTIC_NOT_CONVERGED (none expected) is a warning.
        class(plastic_ice_t), intent(in) :: self
        integer,              intent(in) :: stat
        character(len=*),     intent(in) :: what
        if (stat == PLASTIC_OK) return
        if (stat == PLASTIC_NOT_CONVERGED) then
            call log_msg(LOG_WARN, "kryos_ice_plastic", what//": not converged: "//trim(self%ps%msg), &
                         "instance "//trim(self%name))
            return
        end if
        call fesm_error("kryos_ice_plastic", what//" failed (stat "//to_str(stat)//"): "//trim(self%ps%msg), &
                        "instance "//trim(self%name))
    end subroutine status_check

    subroutine stress_residual(self, rmed, rmax)
        ! The relative stress residual |rho_i g H |grad z_srf| - tau_b| / tau_b of the last solve
        ! (plastic_stress_residual, cell-centred: 1 at the first ice cell by construction), its median and
        ! maximum over the ice cells; computed when the stream is written, not every step.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(out)   :: rmed, rmax

        real(dp), allocatable :: r(:,:)
        integer :: stat

        allocate(r(self%nx, self%ny))
        call plastic_stress_residual(self%ps, self%z_bed, self%z_sl, self%tau_b, self%mask, rmax, r=r, stat=stat, &
                                     margin_frac=self%mfrac)
        call status_check(self, stat, "plastic_stress_residual")
        rmed = median(pack(r, self%mask .and. self%ps%H > 0.0_dp .and. self%tau_b > 0.0_dp))
    end subroutine stress_residual

    function median(v) result(m)
        ! The median of v (0 for an empty v): quickselect on a copy.
        real(dp), intent(in) :: v(:)
        real(dp) :: m

        real(dp), allocatable :: a(:)
        real(dp) :: p, t
        integer  :: n, k, lo, hi, i, j

        n = size(v)
        m = 0.0_dp
        if (n == 0) return
        a = v
        k = (n + 1)/2
        lo = 1;  hi = n
        do while (lo < hi)
            p = a((lo + hi)/2)
            i = lo;  j = hi
            do while (i <= j)
                do while (a(i) < p);  i = i + 1;  end do
                do while (a(j) > p);  j = j - 1;  end do
                if (i <= j) then
                    t = a(i);  a(i) = a(j);  a(j) = t
                    i = i + 1;  j = j - 1
                end if
            end do
            if (k <= j) then
                hi = j
            else if (k >= i) then
                lo = i
            else
                exit
            end if
        end do
        m = a(k)
        if (mod(n, 2) == 0) m = 0.5_dp*(m + minval(a(k+1:n)))
    end function median

    function driving_stress(self) result(tau)
        ! [Pa] rho_i g H |grad z_srf| (centred differences, one-sided at the edges): tau_b of a plastic
        ! ice sheet at yield, for the diagnostics of drive = "smb".
        class(plastic_ice_t), intent(in) :: self
        real(dp) :: tau(self%nx, self%ny)

        real(dp), allocatable :: zs(:,:)
        real(dp) :: gx, gy
        integer  :: i, j, im, ip, jm, jp

        zs = surface(self)
        do j = 1, self%ny
            jm = max(j - 1, 1);  jp = min(j + 1, self%ny)
            do i = 1, self%nx
                im = max(i - 1, 1);  ip = min(i + 1, self%nx)
                gx = (zs(ip,j) - zs(im,j)) / (real(ip - im, dp)*self%dx)
                gy = (zs(i,jp) - zs(i,jm)) / (real(jp - jm, dp)*self%dy)
                tau(i,j) = self%par%rho_i*self%par%g*self%ps%H(i,j)*sqrt(gx*gx + gy*gy)
            end do
        end do
    end function driving_stress

    function surface(self) result(zs)
        ! z_srf: the library's on the ice, the bed or the sea surface elsewhere.
        class(plastic_ice_t), intent(in) :: self
        real(dp) :: zs(self%nx, self%ny)
        where (self%ps%H > 0.0_dp)
            zs = self%ps%z_srf
        elsewhere
            zs = max(self%z_bed, self%z_sl)
        end where
    end function surface

    ! =========================================================================
    ! Exports (contracts-design 2.1)
    ! =========================================================================

    subroutine exports_static(self, t1)
        ! Put once (PASS_INIT): the reference bed, the footprint (the union of the margin slices) and
        ! the present-day ice (the reference's).
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        call fields_put(self%exports, self%h_z_bed_ref, self%z_bed_ref, t1=t1)
        call fields_put(self%exports, self%h_footprint, self%footprint, t1=t1)
        call fields_put(self%exports, self%h_H_ice_pd,  self%H_ice_pd,  t1=t1)
    end subroutine exports_static

    subroutine exports_geometry(self, t1)
        ! INSTANT: the state, the bed and sea level used, f_ice, f_grnd (1 on the ice, the flotation
        ! rule elsewhere), H_grnd, and the reference ice H_ice_ref (with the geometry, as yelmo: a
        ! consumer reading it each step finds it current).
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1

        real(dp), allocatable :: f_grnd(:,:)

        associate (H => self%ps%H)
            allocate(f_grnd(self%nx, self%ny))
            call plastic_flotation_(self, H, self%z_bed, self%z_sl, f_grnd)
            where (H > 0.0_dp) f_grnd = 1.0_dp
            call fields_put(self%exports, self%h_H_ice,  H,             t1=t1)
            call fields_put(self%exports, self%h_z_srf,  surface(self), t1=t1)
            call fields_put(self%exports, self%h_z_bed,  self%z_bed,    t1=t1)
            call fields_put(self%exports, self%h_z_sl,   self%z_sl,     t1=t1)
            call fields_put(self%exports, self%h_f_ice,  merge(1.0_dp, 0.0_dp, H > 0.0_dp), t1=t1)
            call fields_put(self%exports, self%h_f_grnd, f_grnd,        t1=t1)
            call fields_put(self%exports, self%h_H_grnd, &
                            H - (self%par%rho_w/self%par%rho_i)*(self%z_sl - self%z_bed), t1=t1)
            call fields_put(self%exports, self%h_H_ice_ref, self%H_ice_pd, t1=t1)
        end associate
        call exports_velocity(self, t1)
    end subroutine exports_geometry

    subroutine exports_velocity(self, t1)
        ! INSTANT diagnostics (restart = .false.): the surface velocity of the closure at the cell centres.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        if (.not. self%implied) return
        call fields_put(self%exports, self%h_ux_s, self%us_x, t1=t1)
        call fields_put(self%exports, self%h_uy_s, self%us_y, t1=t1)
    end subroutine exports_velocity

    subroutine exports_mass(self, t1)
        ! INSTANT, M_ice [Gt] = rho_i sum(H) dx dy * 1e-12.
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        call fields_put(self%exports, self%h_M_ice, self%par%rho_i*sum(self%ps%H)*self%dx*self%dy*KG_TO_GT, t1=t1)
    end subroutine exports_mass

    subroutine exports_rates(self, st)
        ! MEAN over [t0, t1]: dHidt = (H - H_prev)/(t1 - t0), mb_net = rho_i dHidt [kg m-2 yr-1], bmb = 0.
        ! smb = cmb = 0 (the ice is imposed), or with mass_balance = "implied" the mass balance the geometry
        ! needs [kg m-2 yr-1]: on the ice smb = rho_i (dHidt + div q - q_land) (the outflow to land ablated in
        ! the margin cell), cmb = -rho_i q_marine; on the cells the ice left over the step their loss rho_i
        ! dHidt, as cmb where the bed lies below sea level, else as smb. Over the grid smb + cmb = mb_net
        ! (the margin outflow sums to div q); smb_tot, cmb_tot [Gt yr-1] the sums. div q and the outflow are
        ! those of the geometry at t1. With mb_smooth_km > 0, div q - q_land is smoothed conservatively over the
        ! ice first (smooth_conservative: totals, and so the closure, unchanged).
        class(plastic_ice_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        real(dp), allocatable :: dHdt(:,:), smb(:,:), cmb(:,:), x(:,:)

        dHdt = (self%ps%H - self%H_prev)/(st%t1 - st%t0)
        allocate(smb(self%nx, self%ny), cmb(self%nx, self%ny))
        smb = 0.0_dp;  cmb = 0.0_dp
        if (self%implied) then
            x = self%divq - self%q_land
            if (self%mb_smooth > 0.0_dp) call smooth_conservative(x, self%mask, self%dx, self%dy, self%mb_smooth)
            where (self%mask)
                smb = self%par%rho_i*(dHdt + x)
                cmb = -self%par%rho_i*self%q_marine
            elsewhere (self%z_sl > self%z_bed)
                cmb = self%par%rho_i*dHdt
            elsewhere
                smb = self%par%rho_i*dHdt
            end where
            call fields_put(self%exports, self%h_smb_tot, sum(smb)*self%dx*self%dy*KG_TO_GT, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_cmb_tot, sum(cmb)*self%dx*self%dy*KG_TO_GT, t1=st%t1, t0=st%t0)
        end if
        call fields_put(self%exports, self%h_dHidt,  dHdt,                 t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_smb,    smb,                  t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_bmb,    0.0_dp*dHdt,          t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_cmb,    cmb,                  t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_mb_net, self%par%rho_i*dHdt,  t1=st%t1, t0=st%t0)
    end subroutine exports_rates

    subroutine smooth_conservative(f, mask, dx, dy, sd)
        ! Gaussian smoothing of f over the cells of mask that keeps its sum: each cell spreads its value over
        ! the mask cells within 3 sd, weights exp(-r^2 / (2 sd^2)) normalised by their sum (scatter form), so
        ! sum(f) over the mask is unchanged; f off the mask is left as it is.
        real(dp), intent(inout) :: f(:,:)
        logical,  intent(in)    :: mask(:,:)
        real(dp), intent(in)    :: dx, dy, sd

        real(dp), allocatable :: g(:,:), wsum(:,:)
        integer  :: nx, ny, i, j, p, q, ri, rj
        real(dp) :: w

        nx = size(f, 1);  ny = size(f, 2)
        ri = ceiling(3.0_dp*sd/dx);  rj = ceiling(3.0_dp*sd/dy)
        allocate(g(nx, ny), wsum(nx, ny))
        wsum = 0.0_dp
        !$omp parallel do private(i, j, p, q, w) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (.not. mask(i, j)) cycle
                do q = max(j - rj, 1), min(j + rj, ny)
                    do p = max(i - ri, 1), min(i + ri, nx)
                        if (mask(p, q)) wsum(i, j) = wsum(i, j) + kern(p - i, q - j)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        g = 0.0_dp
        !$omp parallel do private(i, j, p, q, w) schedule(static)
        do q = 1, ny
            do p = 1, nx
                if (.not. mask(p, q)) cycle
                do j = max(q - rj, 1), min(q + rj, ny)
                    do i = max(p - ri, 1), min(p + ri, nx)
                        if (mask(i, j) .and. wsum(i, j) > 0.0_dp) g(p, q) = g(p, q) + f(i, j)*kern(p - i, q - j)/wsum(i, j)
                    end do
                end do
            end do
        end do
        !$omp end parallel do
        where (mask) f = g
    contains
        pure real(dp) function kern(a, b)
            integer, intent(in) :: a, b
            kern = exp(-((a*dx)**2 + (b*dy)**2)/(2.0_dp*sd*sd))
        end function kern
    end subroutine smooth_conservative

    ! =========================================================================
    ! init helpers
    ! =========================================================================

    subroutine group_read(self, ctx)
        ! The keys of the instance group (template [ice_plastic] of input/kryos_defaults.nml; the core
        ! has validated them; a key left out has the template's value, which the defaults here equal),
        ! then the margin, tau_b, the solver parameters, rheology and sliding.
        class(plastic_ice_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        character(len=*), parameter :: proc = "kryos_ice_plastic::init"
        character(len=PATH_LEN) :: margin_file, taub_file, taub_domain_file, taub_scale_file, taub_domain_factor_file
        character(len=NAME_LEN) :: margin_var, margin_interp, margin_extrapolate, mode, taub_method, taub_var
        character(len=NAME_LEN) :: taub_domain_var, taub_id_var, taub_value_var, taub_time_var, taub_override_grp
        character(len=NAME_LEN) :: slide_law, margin_position, mass_balance
        real(dp) :: mb_smooth_km, taub_smooth_km
        character(len=:), allocatable :: who
        real(dp) :: margin_H_min, taub_const, taub_default, rheo_A, rheo_n, slide_beta, slide_m, slide_C, slide_u0
        real(dp) :: slide_umax
        type(friction_t) :: fr
        integer,  allocatable :: codes(:)
        real(dp), allocatable :: vals(:)
        logical  :: hold
        integer  :: interp

        who = "instance "//trim(self%name)
        margin_file = "";  margin_var = "f_ice";  margin_H_min = 0.5_dp;  margin_interp = "sdist"
        margin_extrapolate = "error";  margin_position = "sdist"
        taub_method = "const";  taub_const = 0.0_dp;  taub_file = "";  taub_var = "taub"
        taub_domain_file = "";  taub_domain_var = "taub_domain";  taub_id_var = "domain_id";  taub_value_var = "taub"
        taub_time_var = "taub_t";  taub_override_grp = "";  taub_default = 0.0_dp;  taub_scale_file = ""
        taub_domain_factor_file = "";  taub_smooth_km = 50.0_dp
        mode = "surface"
        rheo_A = 0.0_dp;  rheo_n = 3.0_dp
        slide_law = "none";  slide_beta = 0.0_dp;  slide_m = 1.0_dp;  slide_C = 0.0_dp;  slide_u0 = 0.0_dp
        slide_umax = 0.0_dp
        mass_balance = "none";  mb_smooth_km = 0.0_dp

        call rd_s("drive", self%drive)
        ! solver (plastic_par_t; densities: see constants_check)
        self%par%rho_i = 910.0_dp;  self%par%rho_w = 1028.0_dp;  self%par%g = 9.81_dp
        call rd_s("mode", mode);                     self%par%mode = trim(mode)
        call rd_r("relax", self%par%relax)
        call rd_r("warm_frac_max", self%par%warm_frac_max)
        call rd_i("n_outer", self%par%n_outer)
        call rd_i("max_sweeps", self%par%max_sweeps)
        call rd_r("tol", self%par%tol)
        call rd_r("outer_tol", self%par%outer_tol)
        call rd_r("rho_i", self%par%rho_i)
        call rd_r("rho_w", self%par%rho_w)
        call rd_r("g", self%par%g)
        call rd_r("H_min", self%par%H_min)
        ! margin
        call rd_s("margin_file", margin_file)
        call rd_s("margin_var", margin_var)
        call rd_r("margin_H_min", margin_H_min)
        call rd_s("margin_interp", margin_interp)
        call rd_s("margin_extrapolate", margin_extrapolate)
        call rd_s("margin_position", margin_position)
        ! tau_b
        call rd_s("taub_method", taub_method)
        call rd_r("taub_const", taub_const)
        call rd_s("taub_file", taub_file)
        call rd_s("taub_var", taub_var)
        call rd_s("taub_domain_file", taub_domain_file)
        call rd_s("taub_domain_var", taub_domain_var)
        call rd_s("taub_id_var", taub_id_var)
        call rd_s("taub_value_var", taub_value_var)
        call rd_s("taub_time_var", taub_time_var)
        call rd_s("taub_override", taub_override_grp)
        call rd_r("taub_default", taub_default)
        call rd_s("taub_domain_factor_file", taub_domain_factor_file)
        call rd_s("taub_scale_file", taub_scale_file)
        call rd_r("taub_smooth_km", taub_smooth_km)
        ! rheology and sliding (drive = "smb"; mass_balance = "implied")
        call rd_r("rheo_A", rheo_A)
        call rd_r("rheo_n", rheo_n)
        call rd_s("slide_law", slide_law)
        call rd_r("slide_beta", slide_beta)
        call rd_r("slide_uy", fr%uy)
        call rd_r("slide_uy_b", fr%uy_b)
        call rd_r("slide_umax", slide_umax)
        call rd_r("slide_Hc", fr%Hc)
        call rd_r("slide_lambda_min", fr%lam_min)
        call rd_r("slide_zb0", fr%zb0)
        call rd_r("slide_zb1", fr%zb1)
        call rd_r("slide_p", fr%p)
        call rd_r("slide_m", slide_m)
        call rd_r("slide_C", slide_C)
        call rd_r("slide_u0", slide_u0)
        call rd_s("mass_balance", mass_balance)
        call rd_r("mb_smooth_km", mb_smooth_km)

        select case (trim(self%drive))
        case ("taub", "smb")
        case default
            call fesm_error(proc, "drive = '"//trim(self%drive)//"' is not recognised", who//"; use taub or smb")
        end select

        ! margin chronology
        interp = margin_interp_code(margin_interp)
        if (interp == 0) call fesm_error(proc, "margin_interp = '"//trim(margin_interp)//"' is not recognised", &
                                         who//"; use sdist or nearest")
        select case (trim(margin_extrapolate))
        case ("error");  hold = .false.
        case ("hold");   hold = .true.
        case default
            call fesm_error(proc, "margin_extrapolate = '"//trim(margin_extrapolate)//"' is not recognised", &
                            who//"; use error or hold")
        end select
        select case (trim(margin_position))
        case ("sdist")
            allocate(self%mfrac(self%nx, self%ny))
            self%mfrac = 0.0_dp
        case ("face", "centre")
            self%par%margin = trim(margin_position)
        case default
            call fesm_error(proc, "margin_position = '"//trim(margin_position)//"' is not recognised", &
                            who//"; use sdist, face or centre")
        end select
        if (len_trim(margin_file) == 0) call fesm_error(proc, "margin_file is not set", who)
        call substitute(margin_file, ctx)
        call margin_init(self%mg, trim(margin_file), trim(margin_var), margin_H_min, self%nx, self%ny, &
                         self%dx, self%dy, interp, hold)

        ! tau_b
        if (trim(self%drive) == "taub") then
            select case (taub_method_code(taub_method))
            case (TAUB_M_CONST)
                call taub_set_const(self%tb, taub_const)
            case (TAUB_M_FILE)
                if (len_trim(taub_file) == 0) call fesm_error(proc, "taub_method = 'file' needs taub_file", who)
                call substitute(taub_file, ctx)
                call taub_read_field(self%tb, trim(taub_file), trim(taub_var), self%nx, self%ny)
            case (TAUB_M_DOMAINS)
                if (len_trim(taub_domain_file) == 0) &
                    call fesm_error(proc, "taub_method = 'domains' needs taub_domain_file", who)
                call substitute(taub_domain_file, ctx)
                if (len_trim(taub_override_grp) > 0) then
                    call override_read(self, ctx, trim(taub_override_grp), codes, vals)
                    call taub_read_domains(self%tb, trim(taub_domain_file), trim(taub_domain_var), trim(taub_id_var), &
                                           trim(taub_value_var), trim(taub_time_var), taub_default, hold, self%nx, &
                                           self%ny, codes=codes, values=vals, table="&"//trim(taub_override_grp))
                else
                    call taub_read_domains(self%tb, trim(taub_domain_file), trim(taub_domain_var), trim(taub_id_var), &
                                           trim(taub_value_var), trim(taub_time_var), taub_default, hold, self%nx, &
                                           self%ny)
                end if
            case default
                call fesm_error(proc, "taub_method = '"//trim(taub_method)//"' is not recognised", &
                                who//"; use const, file or domains")
            end select
            if (len_trim(taub_override_grp) > 0 .and. taub_method_code(taub_method) /= TAUB_M_DOMAINS) &
                call fesm_error(proc, "taub_override needs taub_method = 'domains'", who)
            if (len_trim(taub_domain_factor_file) > 0) then
                if (taub_method_code(taub_method) /= TAUB_M_DOMAINS) &
                    call fesm_error(proc, "taub_domain_factor_file needs taub_method = 'domains'", who)
                call substitute(taub_domain_factor_file, ctx)
                call taub_read_domain_factor(self%tb, trim(taub_domain_factor_file))
                call log_msg(LOG_INFO, "kryos_ice_plastic", "taub_domain_factor_file: "// &
                             to_str(count(self%tb%dfac /= 1.0_dp))//" of "//to_str(size(self%tb%dfac))// &
                             " domains scaled, factor "//to_str(minval(self%tb%dfac))//" to "// &
                             to_str(maxval(self%tb%dfac)), who//": "//trim(taub_domain_factor_file))
            end if
            if (len_trim(taub_scale_file) > 0) then
                if (taub_method_code(taub_method) /= TAUB_M_DOMAINS) &
                    call fesm_error(proc, "taub_scale_file needs taub_method = 'domains'", who)
                call substitute(taub_scale_file, ctx)
                call taub_read_scale(self%tb, trim(taub_scale_file))
                call log_msg(LOG_INFO, "kryos_ice_plastic", "taub_scale_file: "//to_str(size(self%tb%slf, 1))// &
                             " groups ("//to_str(count(self%tb%sgrp > 0))//" of "//to_str(size(self%tb%sgrp))// &
                             " domains scaled), "//to_str(size(self%tb%stime))//" nodes, factor "// &
                             to_str(exp(minval(self%tb%slf)))//" to "//to_str(exp(maxval(self%tb%slf))), &
                             who//": "//trim(taub_scale_file))
            end if
            call taub_set_smooth(self%tb, 1000.0_dp*taub_smooth_km, self%dx, self%dy)
        end if

        ! rheology and sliding: drive = "smb" and mass_balance = "implied" need both
        self%closure = rheo_A > 0.0_dp .and. trim(slide_law) /= "none"
        select case (trim(mass_balance))
        case ("none");     self%implied = .false.
        case ("implied");  self%implied = .true.
        case default
            call fesm_error(proc, "mass_balance = '"//trim(mass_balance)//"' is not recognised", who//"; use none or implied")
        end select
        if (mb_smooth_km < 0.0_dp) call fesm_error(proc, "mb_smooth_km must be >= 0", who)
        self%mb_smooth = 1000.0_dp*mb_smooth_km
        if ((trim(self%drive) == "smb" .or. self%implied) .and. .not. self%closure) then
            call fesm_error(proc, "drive = 'smb' and mass_balance = 'implied' need the rheology and sliding", who// &
                            ": set rheo_A > 0 and slide_law (linear, weertman, coulomb, yield or coulomb_n)")
        end if
        if (self%closure) call rheology_set(self, rheo_A, rheo_n, slide_law, slide_beta, slide_m, slide_C, slide_u0, &
                                            slide_umax, fr)
        if (len_trim(self%fr%kind) > 0 .and. trim(self%drive) == "smb") &
            call fesm_error(proc, "slide_law = '"//trim(self%fr%kind)//"' needs drive = 'taub'", who// &
                            ": its coefficient follows the imposed tau_b")

    contains

        subroutine rd_s(key, v)
            character(len=*), intent(in)    :: key
            character(len=*), intent(inout) :: v
            if (nml_has_param(ctx%par_file, ctx%instance, key)) call nml_read(ctx%par_file, ctx%instance, key, v)
        end subroutine rd_s

        subroutine rd_r(key, v)
            character(len=*), intent(in)    :: key
            real(dp),         intent(inout) :: v
            if (nml_has_param(ctx%par_file, ctx%instance, key)) call nml_read(ctx%par_file, ctx%instance, key, v)
        end subroutine rd_r

        subroutine rd_i(key, v)
            character(len=*), intent(in)    :: key
            integer,          intent(inout) :: v
            if (nml_has_param(ctx%par_file, ctx%instance, key)) call nml_read(ctx%par_file, ctx%instance, key, v)
        end subroutine rd_i

    end subroutine group_read

    subroutine rheology_set(self, A, n, law, beta, m, C, u0, umax, fr)
        ! The Glen rheology and the sliding law of the flux closure, from the keys
        ! rheo_* and slide_*. The one place that names their components. "yield", "coulomb_n": the library's
        ! Weertman, regularised Coulomb law with a coefficient per cell (friction_map, every step).
        class(plastic_ice_t), intent(inout) :: self
        real(dp),             intent(in)    :: A, n, beta, m, C, u0, umax
        character(len=*),     intent(in)    :: law
        type(friction_t),     intent(in)    :: fr

        character(len=*), parameter :: proc = "kryos_ice_plastic::init"

        self%rheo%A = A
        self%rheo%n = n
        select case (trim(law))
        case ("linear", "weertman")
            if (.not. beta > 0.0_dp) call fesm_error(proc, "slide_law = '"//trim(law)//"' needs slide_beta > 0", &
                                                    "instance "//trim(self%name))
        case ("coulomb")
            if (.not. (C > 0.0_dp .and. u0 > 0.0_dp)) call fesm_error(proc, "slide_law = 'coulomb' needs slide_C > 0 "// &
                                                                     "and slide_u0 > 0", "instance "//trim(self%name))
        case ("yield")
            if (.not. (fr%uy > 0.0_dp .and. m > 0.0_dp)) call fesm_error(proc, "slide_law = 'yield' needs slide_uy > 0 "// &
                                                                        "and slide_m > 0", "instance "//trim(self%name))
        case ("coulomb_n")
            if (.not. (fr%Hc > 0.0_dp .and. u0 > 0.0_dp .and. m > 0.0_dp .and. umax > 0.0_dp)) &
                call fesm_error(proc, "slide_law = 'coulomb_n' needs slide_Hc, slide_u0, slide_m and slide_umax > 0", &
                                "instance "//trim(self%name))
            if (.not. (fr%lam_min > 0.0_dp .and. fr%lam_min <= 1.0_dp)) &
                call fesm_error(proc, "slide_lambda_min must be in (0, 1]", "instance "//trim(self%name))
            if (fr%lam_min < 1.0_dp .and. .not. fr%zb1 > fr%zb0) &
                call fesm_error(proc, "slide_zb1 must be above slide_zb0", "instance "//trim(self%name))
            if (fr%p < 0.0_dp) call fesm_error(proc, "slide_p must be >= 0", "instance "//trim(self%name))
        case default
            call fesm_error(proc, "slide_law = '"//trim(law)//"' is not recognised", "instance "//trim(self%name)// &
                            "; use none, linear, weertman, coulomb, yield or coulomb_n")
        end select
        if (umax < 0.0_dp) call fesm_error(proc, "slide_umax must be >= 0", "instance "//trim(self%name))
        self%fr = fr
        self%fr%kind = trim(law)
        self%slide%law = friction_law(self%fr)
        if (len_trim(self%slide%law) > 0) then
            allocate(self%slide%fric_map(self%nx, self%ny))
            self%slide%fric_map = huge(1.0_dp)
        else
            self%fr%kind = ""
            self%slide%law = trim(law)
        end if
        self%slide%beta  = beta
        self%slide%m     = m
        self%slide%C     = C
        self%slide%u_0   = u0
        self%slide%u_max = umax
    end subroutine rheology_set

    subroutine override_read(self, ctx, grp, codes, vals)
        ! The group named by taub_override: codes (domain ids) and taub [Pa], one value per code (as the
        ! till_cf group of the yelmo backend).
        class(plastic_ice_t),  intent(in)  :: self
        type(init_ctx_t),      intent(in)  :: ctx
        character(len=*),      intent(in)  :: grp
        integer,  allocatable, intent(out) :: codes(:)
        real(dp), allocatable, intent(out) :: vals(:)

        character(len=*), parameter :: proc = "kryos_ice_plastic::init"
        integer, parameter :: n_max = 256               ! entries counted up to n_max + 1
        character(len=8), parameter :: req(2) = [character(len=8) :: "codes", "taub"]
        character(len=64) :: sv(n_max+1)
        character(len=:), allocatable :: who
        character(len=nml_key_len), allocatable :: keys(:)
        integer :: k, n, m

        who = "instance "//trim(self%name)//": taub_override = '"//grp//"'"
        if (.not. nml_has_group(ctx%par_file, grp)) &
            call fesm_error(proc, "group &"//grp//" of taub_override is missing from the par file", who)
        keys = nml_group_keys(ctx%par_file, grp)
        do k = 1, size(keys)
            if (.not. any(req == keys(k))) &
                call fesm_error(proc, "&"//grp//": unknown key '"//trim(keys(k))//"'", who//"; keys: codes, taub")
        end do
        do k = 1, size(req)
            if (.not. any(keys == req(k))) call fesm_error(proc, "&"//grp//": key '"//trim(req(k))//"' is missing", who)
        end do

        sv = ""
        call nml_read(ctx%par_file, grp, "codes", sv)
        n = count(len_trim(sv) > 0)
        sv = ""
        call nml_read(ctx%par_file, grp, "taub", sv)
        m = count(len_trim(sv) > 0)
        if (n /= m .or. n == 0 .or. n > n_max) then
            call fesm_error(proc, "&"//grp//": codes and taub must have the same number of entries (1 to "// &
                            to_str(n_max)//")", who//"; got "//to_str(n)//" codes, "//to_str(m)//" taub")
        end if
        allocate(codes(n), vals(n))
        call nml_read(ctx%par_file, grp, "codes", codes)
        call nml_read(ctx%par_file, grp, "taub",  vals)
        do k = 1, n
            call log_msg(LOG_INFO, "kryos_ice_plastic", "taub_override: domain "//to_str(codes(k))//": "// &
                         to_str(vals(k))//" Pa", who)
        end do
    end subroutine override_read

    subroutine constants_check(self, ctx)
        ! The densities and g of the group against the run's constants: the globe's ice load (VILMA2)
        ! and the hub use the run's rho_ice, so a difference is logged as a warning.
        class(plastic_ice_t), intent(in) :: self
        type(init_ctx_t),     intent(in) :: ctx

        call differs("rho_i", self%par%rho_i, real(ctx%cnst%rho_ice, dp), "rho_ice")
        call differs("rho_w", self%par%rho_w, real(ctx%cnst%rho_sw, dp), "rho_sw")
        call differs("g",     self%par%g,     real(ctx%cnst%g, dp), "g")

    contains

        subroutine differs(key, v, c, cname)
            character(len=*), intent(in) :: key, cname
            real(dp),         intent(in) :: v, c
            if (.not. phys_const_isset(c)) return
            if (abs(v - c) > 1.0e-9_dp*abs(c)) then
                call log_msg(LOG_WARN, "kryos_ice_plastic", key//" = "//to_str(v)//" differs from the run's "// &
                             cname//" = "//to_str(c), "instance "//trim(self%name)//": the ice mass the globe "// &
                             "gathers uses the run's constants")
            end if
        end subroutine differs

    end subroutine constants_check

    subroutine substitute(path, ctx)
        ! {domain} and {grid_name} of a file name (as recon_path).
        character(len=*), intent(inout) :: path
        type(init_ctx_t), intent(in)    :: ctx
        call nml_replace(path, "{domain}",    trim(ctx%domain))
        call nml_replace(path, "{grid_name}", trim(ctx%grid_name))
    end subroutine substitute

    logical function published(self, h)
        ! Has the import with handle h been published?
        class(plastic_ice_t), intent(in) :: self
        integer,              intent(in) :: h
        published = h > 0
        if (published) published = self%imports%f(h)%published
    end function published

    function restart_file(self, dir) result(file)
        class(plastic_ice_t), intent(in) :: self
        character(len=*),     intent(in) :: dir
        character(len=:), allocatable    :: file
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
    end function restart_file

end module kryos_ice_plastic
