module kryos_shelf_marine_shelf
    ! shelf backend marine_shelf, including bmb_method = "pico" (adapters-design
    ! section 6). No external library and no stub: the code is compiled here.
    !
    ! The yelmox marine_shelf and pico code lives in src/backends/shelf/marine_shelf/,
    ! copied from yelmox 74f3a2b (branch topo-prep), unchanged except the module
    ! names:
    !   libs/marine_shelf.f90     -> kryos_mshlf.f90              (module kryos_mshlf)
    !   libs/pico/pico.f90        -> kryos_mshlf_pico.f90         (module kryos_mshlf_pico)
    !   libs/pico/pico_geometry.f90 -> kryos_mshlf_pico_geometry.f90
    !   libs/pico/pico_physics.f90  -> kryos_mshlf_pico_physics.f90
    !
    ! marine_shelf computes in sp, the exchange is dp: imports are narrowed once,
    ! exports widened from the sp values (P1); bmb_shlf [m ice yr-1, melt < 0]
    ! leaves as kg m-2 yr-1 = bmb * rho_ice, the density marine_shelf holds from
    ! the run's constants (Yelmo divides by its own copy of the same number).
    ! Phase 1 takes ocean profiles on the depth axis: the two-call update.

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use ncio,            only: nc_close
    use nml,             only: nml_read, nml_has_param
    use kryos_outfile,   only: outfile_dims_2d, outfile_open, outfile_var2d
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, component_write_output, &
                               axis_find, NAME_LEN, PASS_BOUNDARY, PASS_RESTART
    use kryos_shelf,     only: shelf_t, shelf_opt_t
    use kryos_mshlf,     only: marshelf_class, marshelf_init, marshelf_update_shelf, marshelf_update, &
                               marshelf_end, marshelf_restart_write, marshelf_restart_read

    implicit none

    private
    public :: marshelf_shelf_t, marshelf_shelf_constraints

    ! Precisions of the 2D stream (input/kryos-variables-<TABLE>.md)
    character(len=*), parameter :: TABLE = "shelf-marine_shelf"

    type, extends(shelf_t) :: marshelf_shelf_t
        type(marshelf_class) :: mshlf                ! the backend object (with pico)
        character(len=NAME_LEN) :: group = "marine_shelf"   ! its group in the par file
        real(sp) :: dx = 0.0_sp                      ! [m] grid spacing (marshelf_update, calc_shelf_slope)
        real(sp), allocatable :: depth(:)            ! [m] depth axis of the ocean profiles
        real(sp), allocatable :: H_ice(:,:), z_bed(:,:), f_grnd(:,:), z_sl(:,:), z_srf(:,:)   ! narrowed imports
        real(sp), allocatable :: to(:,:,:), so(:,:,:), dto(:,:,:)
    contains
        procedure :: init          => marshelf_init_
        procedure :: init_state    => marshelf_init_state_
        procedure :: update        => marshelf_update_
        procedure :: restart_write => marshelf_restart_write_
        procedure :: restart_read  => marshelf_restart_read_
        procedure :: write_output  => marshelf_write_output_
        procedure :: end           => marshelf_end_
    end type marshelf_shelf_t

contains

    function marshelf_shelf_constraints() result(con)
        ! Regional, any grid, several instances, with a restart file.
        type(constraints_t) :: con
        con = constraints_t()
    end function marshelf_shelf_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine marshelf_init_(self, ctx)
        ! marshelf_init (kryos.f90:425-427); the regions are the instance grid's.
        class(marshelf_shelf_t), intent(inout) :: self
        type(init_ctx_t),        intent(in)    :: ctx

        integer :: k

        call group_read(self, ctx)
        self%prec = sp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = shelf_opt_t(ocn_prof=.true., ocn_shelf=.false.)

        ! Phase 1: ocean profiles; the producer's depth axis comes in ctx%axes.
        k = 0
        if (allocated(ctx%axes)) k = axis_find(ctx%axes, "depth")
        if (k == 0) then
            call fesm_error("kryos_shelf_marine_shelf::init", "shelf `"//trim(self%name)// &
                            "` needs ocean profiles on a depth axis", &
                            "the climate of domain `"//trim(self%domain)//"` declares none (has_ocn); "// &
                            "shelf-base ocean values (has_ocn_shelf) arrive in phase 2")
        end if
        self%depth = real(ctx%axes(k)%v, sp)

        ! Fields come from the bundle: marshelf_init would read a restart file itself.
        call restart_key_check(self, ctx)

        call marshelf_init(self%mshlf, trim(ctx%par_file), trim(self%group), self%nx, self%ny, &
                           trim(ctx%domain), trim(ctx%grid_name), ctx%reg, self%corr%basins, self%corr%bmb_corr, &
                           xc=real(ctx%grid%G%x, sp), yc=real(ctx%grid%G%y, sp), cnst=ctx%cnst)

        self%dx = real(ctx%grid%G%dx * ctx%grid%cs%xy_conv, sp)

        allocate(self%H_ice(self%nx, self%ny), self%z_bed(self%nx, self%ny), self%f_grnd(self%nx, self%ny), &
                 self%z_sl(self%nx, self%ny),  self%z_srf(self%nx, self%ny))
        allocate(self%to(self%nx, self%ny, size(self%depth)), self%so(self%nx, self%ny, size(self%depth)), &
                 self%dto(self%nx, self%ny, size(self%depth)))
    end subroutine marshelf_init_

    subroutine marshelf_init_state_(self, st, pass)
        ! The shelf on the initial boundary fields (kryos_startup.f90:166-169):
        ! the same two calls as an update. PASS_RESTART: the diagnostic tf_shlf
        ! again, from the restored now (the state exports come with the
        ! bundle). Other passes: nothing.
        class(marshelf_shelf_t), intent(inout) :: self
        type(step_time_t),       intent(in)    :: st
        integer,                 intent(in)    :: pass

        select case (pass)
        case (PASS_BOUNDARY)
            call calc_shelf(self)
            call exports_put(self, st)
        case (PASS_RESTART)
            call exports_diagnostics(self, st)
        end select
    end subroutine marshelf_init_state_

    subroutine marshelf_update_(self, st)
        ! step_shelf (kryos_coupling.f90:610-619): the two-call update on the
        ! hub geometry and the climate's ocean profiles of this step.
        class(marshelf_shelf_t), intent(inout) :: self
        type(step_time_t),       intent(in)    :: st

        call calc_shelf(self)
        call exports_put(self, st)
    end subroutine marshelf_update_

    subroutine marshelf_restart_write_(self, dir, t)
        ! kryos_startup.f90:428; the file holds now. tf_corr and tf_corr_basin are the shelf
        ! component's (its <instance>_corr_restart.nc).
        class(marshelf_shelf_t), intent(inout) :: self
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t
        call marshelf_restart_write(self%mshlf, restart_file(self, dir), t)
    end subroutine marshelf_restart_write_

    subroutine marshelf_restart_read_(self, dir, t)
        ! kryos_startup.f90:514
        class(marshelf_shelf_t), intent(inout) :: self
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t
        associate (unused => t); end associate
        call marshelf_restart_read(self%mshlf, restart_file(self, dir))
    end subroutine marshelf_restart_read_

    subroutine marshelf_write_output_(self, file, t, stream, init)
        ! Stream 2D: bmb_shlf [m/yr], T_shlf, tf_shlf of mshlf%now (mshlf_write_step,
        ! kryos_output.f90:279-294), same units; precisions from the table. Others: default.
        class(marshelf_shelf_t), intent(inout) :: self
        character(len=*),        intent(in)    :: file, stream
        real(dp),                intent(in)    :: t
        logical,                 intent(in)    :: init

        integer :: ncid, n

        if (trim(stream) /= "2D") then
            call component_write_output(self, file, t, stream, init)
            return
        end if

        if (init) call outfile_dims_2d(file, self%xkm, self%ykm, t)
        call outfile_open(file, t, ncid, n)
        associate (now => self%mshlf%now)
            call outfile_var2d(file, "bmb_shlf", now%bmb_shlf, n, ncid, "m/yr", "Shelf basal mass balance", table=TABLE)
            call outfile_var2d(file, "T_shlf",   now%T_shlf,   n, ncid, "K",    "Shelf temperature",        table=TABLE)
            call outfile_var2d(file, "tf_shlf",  now%tf_shlf,  n, ncid, "K",    "Thermal forcing",          table=TABLE)
        end associate
        call nc_close(ncid)
    end subroutine marshelf_write_output_

    subroutine marshelf_end_(self)
        class(marshelf_shelf_t), intent(inout) :: self
        call marshelf_end(self%mshlf)
    end subroutine marshelf_end_

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    subroutine group_read(self, ctx)
        ! Adapter keys of the instance group; the core has validated them.
        class(marshelf_shelf_t), intent(inout) :: self
        type(init_ctx_t),        intent(in)    :: ctx

        if (nml_has_param(ctx%par_file, ctx%instance, "marine_shelf")) &
            call nml_read(ctx%par_file, ctx%instance, "marine_shelf", self%group)
    end subroutine group_read

    subroutine restart_key_check(self, ctx)
        ! [marine_shelf] restart must be "none" (restarts come from the bundle).
        class(marshelf_shelf_t), intent(in) :: self
        type(init_ctx_t),        intent(in) :: ctx

        character(len=512) :: rst

        rst = "none"
        call nml_read(ctx%par_file, trim(self%group), "restart", rst)
        select case (trim(rst))
        case ("none", "None", "")
        case default
            call fesm_error("kryos_shelf_marine_shelf::init", "["//trim(self%group)//"] restart must be 'None'", &
                            "instance "//trim(self%name)//": restarts come from the bundle, got '"//trim(rst)//"'")
        end select
    end subroutine restart_key_check

    function restart_file(self, dir) result(file)
        class(marshelf_shelf_t), intent(in) :: self
        character(len=*),        intent(in) :: dir
        character(len=:), allocatable       :: file
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
    end function restart_file

    subroutine calc_shelf(self)
        ! Imports narrowed once (P1), then the two calls of
        ! kryos_coupling.f90:616-619: ocean profiles to the shelf base, then the melt.
        class(marshelf_shelf_t), intent(inout) :: self

        call fields_get(self%imports, self%hi_H_ice,  self%H_ice)
        call fields_get(self%imports, self%hi_z_bed,  self%z_bed)
        call fields_get(self%imports, self%hi_f_grnd, self%f_grnd)
        call fields_get(self%imports, self%hi_z_sl,   self%z_sl)
        call fields_get(self%imports, self%hi_z_srf,  self%z_srf)
        call fields_get(self%imports, self%hi_to_ann,  self%to)
        call fields_get(self%imports, self%hi_so_ann,  self%so)
        call fields_get(self%imports, self%hi_dto_ann, self%dto)

        call marshelf_update_shelf(self%mshlf, self%H_ice, self%z_bed, self%f_grnd, self%z_sl, &
                                   self%dx, self%depth, self%to, self%so, dto_ann=self%dto)
        call marshelf_update(self%mshlf, self%H_ice, self%z_bed, self%f_grnd, self%z_sl, &
                             dx=self%dx, tf_corr=self%corr%tf_corr, tf_corr_basin=self%corr%tf_corr_basin, &
                             z_srf=self%z_srf)
    end subroutine calc_shelf

    subroutine exports_put(self, st)
        ! bmb_shlf: MEAN over [t0, t1], one put, kg m-2 yr-1 = bmb [m ice yr-1] * rho_ice
        ! (exact product of two sp numbers in dp). T_shlf: INSTANT at t1; then
        ! the diagnostics.
        class(marshelf_shelf_t), intent(inout) :: self
        type(step_time_t),       intent(in)    :: st

        associate (now => self%mshlf%now)
            call fields_put(self%exports, self%h_bmb_shlf, &
                            real(now%bmb_shlf, dp) * real(self%mshlf%par%rho_ice, dp), t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_T_shlf,  now%T_shlf,  t1=st%t1)
        end associate
        call exports_diagnostics(self, st)
    end subroutine exports_put

    subroutine exports_diagnostics(self, st)
        ! The diagnostic (restart = .false.) tf_shlf: INSTANT at t1. Also on
        ! PASS_RESTART, from now as marshelf_restart_read restored it.
        class(marshelf_shelf_t), intent(inout) :: self
        type(step_time_t),       intent(in)    :: st
        call fields_put(self%exports, self%h_tf_shlf, self%mshlf%now%tf_shlf, t1=st%t1)
    end subroutine exports_diagnostics

end module kryos_shelf_marine_shelf
