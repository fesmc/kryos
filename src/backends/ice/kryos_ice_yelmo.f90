module kryos_ice_yelmo
    ! ice backend Yelmo (adapters-design section 2 and its decision of
    ! 2026-10-24), built with yelmo=1; otherwise kryos_ice_yelmo_stub.f90
    ! provides this module. Yelmo computes in sp and the exchange is dp: imports
    ! are narrowed once, exports widened from Yelmo's sp values (P1-P5), in
    ! the operand order of yelmox (kryos_coupling.f90, kryos_startup.f90).

    use, intrinsic :: iso_fortran_env, only: int64
    use precision,       only: sp, dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_has_param
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, grid_ref_t, &
                               component_write_output, NAME_LEN, PATH_LEN, &
                               PASS_INIT, PASS_GEOMETRY, PASS_STATE, PASS_SPINUP, PASS_RESTART
    use kryos_ice,       only: ice_t, ice_opt_t
    use yelmo_defs,      only: MASK_ICE_NONE
    use yelmo,           only: yelmo_class, ytopo_input_class, yelmo_init, yelmo_init_grid, yelmo_init_state, &
                               yelmo_update, yelmo_update_equil, yelmo_end, yelmo_print_bound, &
                               yelmo_restart_write, yelmo_restart_init, yelmo_regions_update, &
                               yelmo_regions_write, yelmo_write_init, yelmo_write_step
    use yelmo_io,         only: yelmo_restart_read_topo_bnd
    use yelmo_topography, only: calc_ytopo_diagnostic
    use topography,       only: calc_H_grnd, calc_f_grnd_subgrid_linear, calc_f_grnd_subgrid_area

    implicit none

    private
    public :: yelmo_ice_t, yelmo_ice_constraints

    ! Variable sets of the streams yelmo and yelmo_sm (yelmox kryos_output.f90)
    character(len=56), parameter :: YELMO_VARS_2D(24) = [ character(len=56) :: &
        "H_ice","z_srf","z_bed","mask_bed","uxy_b","uxy_s","uxy_bar", &
        "ux_bar","uy_bar","cb_ref","N_eff","beta","taub","taud","visc_bar", &
        "T_prime_b","hyd_W_til","mb_net","smb","bmb","cmb","z_sl","bmb_shlf", &
        "H_ice_ref" ]
    character(len=56), parameter :: YELMO_VARS_2D_SM(10) = [ character(len=56) :: &
        "H_ice","z_srf","z_bed","mask_bed","uxy_s","smb","bmb","z_sl", &
        "hyd_W_til","T_prime_b" ]

    type, extends(ice_t) :: yelmo_ice_t
        type(yelmo_class) :: ylmo
        character(len=NAME_LEN) :: group = "yelmo"        ! Yelmo's top group in the par file
        real(sp) :: init_equil_time = 0.0_sp              ! [yr] cold start: yelmo_update_equil, 0 = none
        real(sp) :: init_time_thrm  = 0.0_sp              ! [yr] then equilibrate with topography fixed, 0 = none
        real(dp) :: t_Q_geo = -huge(1.0_dp), t_H_sed = -huge(1.0_dp)   ! publication times last copied
        integer(int64) :: tick_prepared = -1                  ! clock tick (st%tick) of the last prepare
    contains
        procedure :: init          => yelmo_init_
        procedure :: init_state    => yelmo_init_state_
        procedure :: prepare       => yelmo_prepare_
        procedure :: update        => yelmo_update_
        procedure :: flotation     => yelmo_flotation_
        procedure :: restart_write => yelmo_restart_write_
        procedure :: restart_read  => yelmo_restart_read_
        procedure :: restart_files => yelmo_restart_files_
        procedure :: write_output  => yelmo_write_output_
        procedure :: end           => yelmo_end_
    end type yelmo_ice_t

contains

    function yelmo_ice_constraints() result(con)
        ! Any grid, several instances, regional, with a restart file.
        type(constraints_t) :: con
        con = constraints_t(build_flag="yelmo=1")
    end function yelmo_ice_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine yelmo_init_(self, ctx)
        ! kryos_startup.f90 / kryos.f90:303-313 (yelmo_init with the prepared
        ! reference) and domain_regions_init (region file names).
        class(yelmo_ice_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx

        type(ytopo_input_class) :: topo

        call group_read(self, ctx)
        self%prec = sp
        self%caps%has_flotation = .true.
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = ice_opt_t(T_shlf=.true., tf_shlf=.false., Q_geo=.true., H_sed=.true., Qd=.false., &
                              uxy_s=.false., dz_bed_dt_corr=.false., z_bed_ref=.true., H_ice_ref=.true., &
                              footprint=.true., H_ice_pd=.true.)

        call yelmo_init_grid(self%ylmo%grd, ctx%grid)
        call reference_narrow(self, ctx%ref, topo)
        call yelmo_init(self%ylmo, filename=trim(ctx%par_file), grid_def="none", time=ctx%time, &
                        domain=trim(ctx%domain), grid_name=trim(ctx%grid_name), group=trim(self%group), &
                        reg=ctx%reg, topo_init=topo, topo_pd=topo, cnst=ctx%cnst, prepared=.true.)

        call region_files(self)
        call reference_check(self, topo)
    end subroutine yelmo_init_

    subroutine yelmo_init_state_(self, st, pass)
        ! Passes of core-design 3.2; others: nothing. PASS_RESTART (after
        ! restart_read and the core's restore of the imports and exports):
        ! Yelmo's diagnostics on the restored state; no export is put, the
        ! core restored them as the bundle holds them (adapters-design A5).
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        integer,            intent(in)    :: pass

        select case (pass)
        case (PASS_INIT)                                ! the initialised state; MEAN exports not put
            call exports_geometry(self, st%t1)
            call exports_reference(self, st%t1)
            call exports_static(self, st%t1)
            call exports_mass(self, st%t1)
        case (PASS_GEOMETRY)                            ! domain_init_isostasy tail (kryos_startup.f90:246-249)
            call bnd_geometry(self)
            call topo_diagnostic(self)
            call exports_geometry(self, st%t1)
        case (PASS_STATE)                               ! kryos_startup.f90:178-179 (bnd: prepare)
            call yelmo_print_bound(self%ylmo%bnd)
            call yelmo_init_state(self%ylmo, time=st%t1, thrm_method="robin-cold")
        case (PASS_SPINUP)                              ! domain_init_ice 252-287, the experiment hook before it
            if (self%init_equil_time > 0.0_sp) &
                call yelmo_update_equil(self%ylmo, st%t1, time_tot=self%init_equil_time, dt=1.0_sp, &
                                        topo_fixed=.false.)
            if (self%init_time_thrm > 0.0_sp) &
                call yelmo_update_equil(self%ylmo, st%t1, time_tot=self%init_time_thrm, dt=real(st%dtt, sp), &
                                        topo_fixed=.true.)
            call exports_all(self, st)
        case (PASS_RESTART)                             ! domain_restart_read tail (kryos_startup.f90:521-527)
            call bnd_geometry(self)                     ! the restored imports; exports: restored by the core
            call topo_diagnostic(self)
            call yelmo_regions_update(self%ylmo)
        end select
    end subroutine yelmo_init_state_

    subroutine yelmo_prepare_(self, st)
        ! The imports of the step into Yelmo's boundary (couple_to_yelmo), where
        ! an experiment hook may change them before update or, at the cold
        ! start, PASS_STATE. The tick is recorded for the check in update.
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st

        call bnd_from_imports(self)
        self%tick_prepared = st%tick
    end subroutine yelmo_prepare_

    subroutine yelmo_update_(self, st)
        ! step_icesheet (kryos_coupling.f90:386-396) on the boundary prepared for the step.
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st

        if (st%tick /= self%tick_prepared) then
            call fesm_error("kryos_ice_yelmo::update", "the boundary was not prepared for this step", &
                            "instance "//trim(self%name)//": clock tick "//to_str(st%tick)//", last prepared at tick "// &
                            to_str(self%tick_prepared)//"; the run sequence must call domain_prepare before the update")
        end if
        call yelmo_update(self%ylmo, st%t1)
        call exports_geometry(self, st%t1)
        call exports_mass(self, st%t1)
        call exports_rates(self, st)
    end subroutine yelmo_update_

    subroutine yelmo_flotation_(self, H, z_bed, z_sl, f_grnd)
        ! Yelmo's own grounded fraction of ice H over z_bed and z_sl, as
        ! calc_ytopo_diagnostic forms f_grnd: calc_H_grnd (use_f_ice = .false.:
        ! the thickness at face value, so no f_ice) with the run's densities, then
        ! the subgrid rule of ytopo.gl_sep on the instance grid. Narrowed once,
        ! computed in sp, widened: on Yelmo's own geometry it is tpo%now%f_grnd
        ! bit for bit (pinning points enter f_grnd_pin, not f_grnd).
        class(yelmo_ice_t), intent(in)  :: self
        real(dp),           intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),           intent(out) :: f_grnd(:,:)

        real(sp), allocatable :: H_s(:,:), f_ice(:,:), H_grnd(:,:), f_g(:,:), f_x(:,:), f_y(:,:)

        allocate(H_s(self%nx, self%ny), f_ice(self%nx, self%ny), H_grnd(self%nx, self%ny), f_g(self%nx, self%ny))
        H_s   = real(H, sp)
        f_ice = 1.0_sp                                  ! not read: use_f_ice = .false.
        call calc_H_grnd(H_grnd, H_s, f_ice, real(z_bed, sp), real(z_sl, sp), &
                         self%ylmo%bnd%c%rho_ice, self%ylmo%bnd%c%rho_sw, use_f_ice=.false.)
        select case (self%ylmo%tpo%par%gl_sep)
        case (1)                                        ! binary f_grnd
            allocate(f_x(self%nx, self%ny), f_y(self%nx, self%ny))
            call calc_f_grnd_subgrid_linear(f_g, f_x, f_y, H_grnd, self%ylmo%tpo%par%boundaries)
        case (3)                                        ! grounded area of the interpolated H_grnd
            call calc_f_grnd_subgrid_area(f_g, H_grnd=H_grnd, boundaries=self%ylmo%tpo%par%boundaries)
        case default
            call fesm_error("kryos_ice_yelmo::flotation", "ytopo.gl_sep = "//to_str(self%ylmo%tpo%par%gl_sep)// &
                            " is not 1 or 3", "instance "//trim(self%name))
        end select
        f_grnd = real(f_g, dp)
    end subroutine yelmo_flotation_

    subroutine yelmo_restart_write_(self, dir, t)
        ! kryos_startup.f90:427
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        call yelmo_restart_write(self%ylmo, restart_file(self, dir), t)
    end subroutine yelmo_restart_write_

    function yelmo_restart_files_(self) result(files)
        ! The files yelmo_restart_write creates: <instance>_restart.nc and, for
        ! each Lagrangian tracer backend that is on and has state (ytrc_restart_write),
        ! the sidecar <instance>_restart_elsa.nc / <instance>_restart_tracer.nc.
        class(yelmo_ice_t), intent(in) :: self
        character(len=PATH_LEN), allocatable :: files(:)

        associate (trc => self%ylmo%trc)
            files = [character(len=PATH_LEN) :: trim(self%name)//"_restart.nc"]
            if (trc%par%use_elsa) then
                if (allocated(trc%elsa%now%d_iso)) &
                    files = [character(len=PATH_LEN) :: files, trim(self%name)//"_restart_elsa.nc"]
            end if
            if (trc%par%use_tracer) then
                if (allocated(trc%trc%now%x)) &
                    files = [character(len=PATH_LEN) :: files, trim(self%name)//"_restart_tracer.nc"]
            end if
        end associate
    end function yelmo_restart_files_

    subroutine yelmo_restart_read_(self, dir, t)
        ! kryos_startup.f90:493-497
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t

        character(len=:), allocatable :: file

        file = restart_file(self, dir)
        call yelmo_restart_read_topo_bnd(self%ylmo%tpo, self%ylmo%bnd, self%ylmo%time, &
                                         self%ylmo%par%restart_interpolated, self%ylmo%grd, &
                                         self%ylmo%par%domain, self%ylmo%par%grid_name, file, t)
        call yelmo_restart_init(self%ylmo, file, t)
        self%ylmo%par%use_restart = .true.
    end subroutine yelmo_restart_read_

    subroutine yelmo_write_output_(self, file, t, stream, init)
        ! Yelmo's own streams (kryos_output.f90:55-190); any other: the default writer.
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init

        select case (trim(stream))
        case ("yelmo")
            if (init) call yelmo_write_init(self%ylmo, file, time_init=t, units="years")
            call yelmo_write_step(self%ylmo, file, t, nms=YELMO_VARS_2D, compare_pd=.false.)
        case ("yelmo_sm")
            if (init) call yelmo_write_init(self%ylmo, file, time_init=t, units="years")
            call yelmo_write_step(self%ylmo, file, t, nms=YELMO_VARS_2D_SM, compare_pd=.false.)
        case ("yelmo_ts")
            if (trim(file) /= trim(self%ylmo%reg%fnm)) then
                call fesm_error("kryos_ice_yelmo::write_output", "stream `yelmo_ts` writes Yelmo's own region files", &
                                "instance "//trim(self%name)//": the core's file `"//trim(file)// &
                                "` is not Yelmo's global region file `"//trim(self%ylmo%reg%fnm)//"`")
            end if
            call yelmo_regions_write(self%ylmo, t, init=init, units="years")
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine yelmo_write_output_

    subroutine yelmo_end_(self)
        class(yelmo_ice_t), intent(inout) :: self
        call yelmo_end(self%ylmo, time=self%ylmo%tpo%par%time)
    end subroutine yelmo_end_

    ! =========================================================================
    ! init helpers
    ! =========================================================================

    subroutine group_read(self, ctx)
        ! Adapter keys of the instance group; the core has validated them.
        class(yelmo_ice_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx

        if (nml_has_param(ctx%par_file, ctx%instance, "yelmo")) &
            call nml_read(ctx%par_file, ctx%instance, "yelmo", self%group)
        if (nml_has_param(ctx%par_file, ctx%instance, "init_equil_time")) &
            call nml_read(ctx%par_file, ctx%instance, "init_equil_time", self%init_equil_time)
        if (nml_has_param(ctx%par_file, ctx%instance, "init_time_thrm")) &
            call nml_read(ctx%par_file, ctx%instance, "init_time_thrm", self%init_time_thrm)
    end subroutine group_read

    subroutine reference_narrow(self, ref, topo)
        ! The core's reference on this grid as Yelmo's topography input (sp).
        class(yelmo_ice_t),      intent(in)  :: self
        type(grid_ref_t),        intent(in)  :: ref
        type(ytopo_input_class), intent(out) :: topo

        if (.not. (allocated(ref%z_bed) .and. allocated(ref%H_ice) .and. &
                   allocated(ref%z_srf) .and. allocated(ref%z_bed_sd))) then
            call fesm_error("kryos_ice_yelmo::init", "reference topography incomplete", &
                            "instance "//trim(self%name)//": z_bed, H_ice, z_srf and z_bed_sd are needed on grid `"// &
                            trim(self%grid)//"`")
        end if
        topo%z_bed    = real(ref%z_bed,    sp)
        topo%H_ice    = real(ref%H_ice,    sp)
        topo%z_srf    = real(ref%z_srf,    sp)
        topo%z_bed_sd = real(ref%z_bed_sd, sp)
    end subroutine reference_narrow

    subroutine region_files(self)
        ! Yelmo writes its region files itself: <instance>_yelmo_ts[_<region>].nc
        ! in the run directory (domain_regions_init, kryos_regions.f90:20-46).
        class(yelmo_ice_t), intent(inout) :: self
        integer :: i

        self%ylmo%outfldr  = ""
        self%ylmo%reg%fnm  = trim(self%name)//"_yelmo_ts.nc"
        do i = 1, self%ylmo%par%n_reg
            self%ylmo%regs(i)%fnm = trim(self%name)//"_yelmo_ts_"//trim(self%ylmo%regs(i)%name)//".nc"
        end do
    end subroutine region_files

    subroutine reference_check(self, topo)
        ! Yelmo takes the prepared reference as given (adapters-design 2.5): its
        ! bed reference is the core's bit for bit, its ice reference differs only
        ! where mask_ice removes ice.
        class(yelmo_ice_t),      intent(in) :: self
        type(ytopo_input_class), intent(in) :: topo

        integer  :: n
        real(dp) :: dmax

        n = count(self%ylmo%bnd%z_bed_ref /= topo%z_bed)
        if (n > 0) then
            dmax = maxval(abs(real(self%ylmo%bnd%z_bed_ref, dp) - real(topo%z_bed, dp)))
            call fesm_error("kryos_ice_yelmo::init", "z_bed_ref differs from the core's reference", &
                            "instance "//trim(self%name)//": "//to_str(n)//" cells differ, max |diff| = "// &
                            to_str(dmax)//" m")
        end if

        n = count(self%ylmo%bnd%H_ice_ref /= topo%H_ice)
        dmax = 0.0_dp
        if (n > 0) dmax = maxval(abs(real(self%ylmo%bnd%H_ice_ref, dp) - real(topo%H_ice, dp)))
        call log_msg(LOG_INFO, "kryos_ice_yelmo", "H_ice_ref against the core's reference: "// &
                     to_str(n)//" cells differ, max |diff| = "//to_str(dmax)//" m", &
                     "instance "//trim(self%name)//" (mask_ice removes ice outside the dynamic region)")
    end subroutine reference_check

    function restart_file(self, dir) result(file)
        class(yelmo_ice_t), intent(in) :: self
        character(len=*),   intent(in) :: dir
        character(len=:), allocatable  :: file
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
    end function restart_file

    ! =========================================================================
    ! Imports into bnd (adapters-design 2.1; yelmox couple_*_to_yelmo)
    ! =========================================================================

    subroutine bnd_geometry(self)
        ! Bed and sea level (couple_isostasy_to_yelmo): narrowed.
        class(yelmo_ice_t), intent(inout) :: self
        call fields_get(self%imports, self%hi_z_bed, self%ylmo%bnd%z_bed)
        call fields_get(self%imports, self%hi_z_sl,  self%ylmo%bnd%z_sl)
    end subroutine bnd_geometry

    subroutine bnd_from_imports(self)
        ! Every wired import, converted once. An optional import that is not
        ! wired keeps Yelmo's default; the static ones (Q_geo, H_sed) are
        ! copied when their publication time has changed.
        class(yelmo_ice_t), intent(inout) :: self

        real(dp), allocatable :: v(:,:)

        allocate(v(self%nx, self%ny))
        associate (bnd => self%ylmo%bnd)
            call bnd_geometry(self)

            call fields_get(self%imports, self%hi_smb, bnd%smb)                  ! couple_surface_to_yelmo:321
            bnd%smb = bnd%smb * bnd%c%conv_we_ie * 1e-3
            call fields_get(self%imports, self%hi_T_srf, bnd%T_srf)

            call fields_get(self%imports, self%hi_bmb_shlf, v)                   ! inverse of the shelf's * rho_ice
            bnd%bmb_shlf = real(v / real(bnd%c%rho_ice, dp), sp)
            if (wired(self, self%hi_T_shlf)) call fields_get(self%imports, self%hi_T_shlf, bnd%T_shlf)

            if (changed(self, self%hi_Q_geo, self%t_Q_geo)) then                 ! W m-2 -> mW m-2
                call fields_get(self%imports, self%hi_Q_geo, v)
                bnd%Q_geo = real(v * 1000.0_dp, sp)
            end if
            if (changed(self, self%hi_H_sed, self%t_H_sed)) call fields_get(self%imports, self%hi_H_sed, bnd%H_sed)
        end associate
    end subroutine bnd_from_imports

    logical function wired(self, h)
        ! Has the optional import with handle h been published?
        class(yelmo_ice_t), intent(in) :: self
        integer,            intent(in) :: h
        wired = h > 0
        if (wired) wired = self%imports%f(h)%published
    end function wired

    logical function changed(self, h, t_seen)
        ! Wired, and published at a time other than t_seen (updated to it).
        class(yelmo_ice_t), intent(in)    :: self
        integer,            intent(in)    :: h
        real(dp),           intent(inout) :: t_seen
        changed = wired(self, h)
        if (changed) changed = self%imports%f(h)%t1 /= t_seen
        if (changed) t_seen = self%imports%f(h)%t1
    end function changed

    subroutine topo_diagnostic(self)
        ! Reconcile f_ice, f_grnd, H_grnd, z_srf with the bed and sea level of bnd.
        class(yelmo_ice_t), intent(inout) :: self
        call calc_ytopo_diagnostic(self%ylmo%tpo, self%ylmo%dyn, self%ylmo%mat, self%ylmo%thrm, self%ylmo%bnd)
    end subroutine topo_diagnostic

    ! =========================================================================
    ! Exports (adapters-design 2.1): sp values widened, no dp arithmetic first,
    ! except the mass fluxes (* rho_ice, Yelmo's own density)
    ! =========================================================================

    subroutine exports_all(self, st)
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        call exports_geometry(self, st%t1)
        call exports_reference(self, st%t1)
        call exports_mass(self, st%t1)
        call exports_rates(self, st)
    end subroutine exports_all

    subroutine exports_geometry(self, t1)
        ! INSTANT: the state and the bed and sea level used.
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t1

        associate (tpo => self%ylmo%tpo%now, bnd => self%ylmo%bnd)
            call fields_put(self%exports, self%h_H_ice,  tpo%H_ice,  t1=t1)
            call fields_put(self%exports, self%h_z_srf,  tpo%z_srf,  t1=t1)
            call fields_put(self%exports, self%h_z_bed,  bnd%z_bed,  t1=t1)
            call fields_put(self%exports, self%h_z_sl,   bnd%z_sl,   t1=t1)
            call fields_put(self%exports, self%h_f_ice,  tpo%f_ice,  t1=t1)
            call fields_put(self%exports, self%h_f_grnd, tpo%f_grnd, t1=t1)
        end associate
    end subroutine exports_geometry

    subroutine exports_reference(self, t1)
        ! INSTANT, static: the reference Yelmo holds (processed by mask_ice).
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t1

        call fields_put(self%exports, self%h_z_bed_ref, self%ylmo%bnd%z_bed_ref, t1=t1)
        call fields_put(self%exports, self%h_H_ice_ref, self%ylmo%bnd%H_ice_ref, t1=t1)
    end subroutine exports_reference

    subroutine exports_static(self, t1)
        ! INSTANT, put once (PASS_INIT, after yelmo_init: mask_ice and the present
        ! day are final, with or without a restart): footprint = 1 where Yelmo
        ! holds ice (mask_ice /= MASK_ICE_NONE), and H_ice_pd = dta%pd%H_ice, the
        ! present-day ice, not bnd%H_ice_ref, which recon and ice_bnd overwrite.
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t1

        call fields_put(self%exports, self%h_footprint, &
                        merge(1.0_sp, 0.0_sp, self%ylmo%bnd%mask_ice /= MASK_ICE_NONE), t1=t1)
        call fields_put(self%exports, self%h_H_ice_pd, self%ylmo%dta%pd%H_ice, t1=t1)
    end subroutine exports_static

    subroutine exports_mass(self, t1)
        ! INSTANT, M_ice [Gt] = V_ice * rho_ice * 1e-3 in sp (yelmox.f90:185).
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t1

        real(sp) :: m_ice

        m_ice = self%ylmo%reg%V_ice * self%ylmo%bnd%c%rho_ice * 1e-3_sp
        call fields_put(self%exports, self%h_M_ice, m_ice, t1=t1)
    end subroutine exports_mass

    subroutine exports_rates(self, st)
        ! MEAN over [t0, t1], one put each: Yelmo's rates of the call [m ie yr-1];
        ! the mass fluxes in kg m-2 yr-1.
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st

        real(dp) :: rho

        rho = real(self%ylmo%bnd%c%rho_ice, dp)
        associate (tpo => self%ylmo%tpo%now)
            call fields_put(self%exports, self%h_dHidt,  tpo%dHidt, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_smb,    real(tpo%smb,    dp) * rho, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_bmb,    real(tpo%bmb,    dp) * rho, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_cmb,    real(tpo%cmb,    dp) * rho, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_mb_net, real(tpo%mb_net, dp) * rho, t1=st%t1, t0=st%t0)
        end associate
    end subroutine exports_rates

end module kryos_ice_yelmo
