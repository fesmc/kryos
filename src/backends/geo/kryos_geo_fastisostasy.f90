module kryos_geo_fastisostasy
    ! geo backend FastIsostasy (adapters-design section 3), built with
    ! fastisostasy=1, else the stub: fastiso_geo_t (regional, isos_class) and
    ! fastiso_bsl_t (global part, bsl_class); the factory picks by scope.
    ! FastIsostasy and barysealevel compute in sp and the exchange is dp: imports
    ! are narrowed once, exports are the library's sp values widened (no dp
    ! arithmetic first), in the operand order of yelmox (kryos_coupling.f90,
    ! kryos_startup.f90). Decision of 2026-10-24: the regional instance is
    ! initialised from the core's reference on its grid (ctx%ref), not from the
    ! ice exports.

    use precision,       only: sp, dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_has_param
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_open, nc_close, nc_time_index
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, component_write_output, NAME_LEN, &
                               PASS_INIT, PASS_GEOMETRY, PASS_RESTART
    use kryos_geo,       only: geo_t, geo_opt_t
    use fastisostasy,    only: isos_class, isos_init, isos_init_ref, isos_init_state, isos_update, isos_end, &
                               isos_restart_write, bsl_class, bsl_init, bsl_update, bsl_restart_write, &
                               bsl_restart_read, bsl_write_init, bsl_write_step
    use barysealevel,    only: bsl_end

    implicit none

    private
    public :: fastiso_geo_t, fastiso_bsl_t, fastiso_geo_constraints, fastiso_bsl_constraints

    character(len=*), parameter :: BSL_GROUP = "barysealevel"   ! fixed in bsl_init (barysealevel.f90:59-61)

    type, extends(geo_t) :: fastiso_geo_t
        type(isos_class) :: isos                       ! the domain's FastIsostasy state
        type(bsl_class)  :: bsl_loc                    ! copy of the global bsl state (method, bsl_now, A_ocean_now)
        character(len=NAME_LEN) :: group = "isos"      ! FastIsostasy's group in the par file
        real(sp) :: dx = 0.0_sp                        ! [m] cell size handed to isos_init
        real(sp), allocatable :: z_bed_ref_g(:,:), H_ice_ref_g(:,:)   ! ctx%ref on the geo grid, narrowed
        real(dp), allocatable :: xc(:), yc(:)          ! [km] axes of the output files
    contains
        procedure :: init          => fastiso_geo_init_
        procedure :: init_state    => fastiso_geo_init_state_
        procedure :: update        => fastiso_geo_update_
        procedure :: restart_write => fastiso_geo_restart_write_
        procedure :: restart_read  => fastiso_geo_restart_read_
        procedure :: write_output  => fastiso_geo_write_output_
        procedure :: end           => fastiso_geo_end_
    end type fastiso_geo_t

    type, extends(geo_t) :: fastiso_bsl_t
        type(bsl_class) :: bsl                         ! barystatic sea level (global part)
    contains
        procedure :: init          => fastiso_bsl_init_
        procedure :: init_state    => fastiso_bsl_init_state_
        procedure :: update        => fastiso_bsl_update_
        procedure :: reduce        => fastiso_bsl_reduce_
        procedure :: restart_write => fastiso_bsl_restart_write_
        procedure :: restart_read  => fastiso_bsl_restart_read_
        procedure :: write_output  => fastiso_bsl_write_output_
        procedure :: end           => fastiso_bsl_end_
    end type fastiso_bsl_t

contains

    function fastiso_geo_constraints() result(con)
        ! Regional only, any number of domains.
        type(constraints_t) :: con
        con = constraints_t(build_flag="fastisostasy=1")
    end function fastiso_geo_constraints

    function fastiso_bsl_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(build_flag="fastisostasy=1", regional=.false., global=.true., max_instances=1)
    end function fastiso_bsl_constraints

    ! =========================================================================
    ! Regional: bindings
    ! =========================================================================

    subroutine fastiso_geo_init_(self, ctx)
        ! kryos.f90:342-352 (isos_init on grid_isos), bsl method for the local copy.
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        real(sp) :: dy

        call group_read(self, ctx)
        self%prec = sp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = geo_opt_t(A_ocn=.true., dV_bsl=.true., dz_bed_dt=.true., ref=.false.)

        self%dx = real(ctx%grid%G%dx * ctx%grid%cs%xy_conv, sp)
        dy      = real(ctx%grid%G%dy * ctx%grid%cs%xy_conv, sp)
        if (self%dx /= dy) then
            call fesm_error("kryos_geo_fastisostasy::init", "FastIsostasy supports square grid cells only", &
                            "instance "//trim(self%name)//": dx = "//to_str(real(self%dx, dp))//" m, dy = "// &
                            to_str(real(dy, dp))//" m on grid `"//trim(self%grid)//"`")
        end if
        call log_msg(LOG_INFO, "kryos_geo_fastisostasy", "isostasy grid `"//trim(self%grid)//"`, dx = "// &
                     to_str(real(self%dx, dp))//" m", "instance "//trim(self%name))

        call isos_init(self%isos, trim(ctx%par_file), trim(self%group), self%nx, self%ny, self%dx, dy, &
                       cnst=ctx%cnst, domain=trim(ctx%domain), grid_name=trim(ctx%grid_name))
        if (self%isos%par%use_restart) then
            call fesm_error("kryos_geo_fastisostasy::init", "`restart` of group `"//trim(self%group)// &
                            "` must be None", "instance "//trim(self%name)// &
                            ": restarts come from the bundle (`[kryos] restart`)")
        end if

        call bsl_loc_init(self, ctx)
        call reference_narrow(self, ctx)
        self%xc = ctx%grid%G%x
        self%yc = ctx%grid%G%y
    end subroutine fastiso_geo_init_

    subroutine fastiso_geo_init_state_(self, st, pass)
        ! Passes of core-design 3.2; others: nothing.
        class(fastiso_geo_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        integer,              intent(in)    :: pass

        select case (pass)
        case (PASS_INIT)                                ! nothing to put: the state needs the global bsl
        case (PASS_GEOMETRY)                            ! domain_init_isostasy (kryos_startup.f90:236-241)
            call bsl_from_imports(self)
            call isos_init_ref(self%isos, self%z_bed_ref_g, self%H_ice_ref_g)
            call isos_init_state(self%isos, self%z_bed_ref_g, self%H_ice_ref_g, st%t1, self%bsl_loc)
            call exports_state(self, st%t1)
            call exports_volume(self, st%t1)
        case (PASS_RESTART)                             ! the state restart_read rebuilt
            call exports_state(self, st%t1)
        end select
    end subroutine fastiso_geo_init_state_

    subroutine fastiso_geo_update_(self, st)
        ! step_isostasy (kryos_coupling.f90:150-160) on the imports of the step.
        ! `bsl` is exported by geo_scatter, not here.
        class(fastiso_geo_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        real(sp), allocatable :: H_ice(:,:)

        allocate(H_ice(self%nx, self%ny))
        call fields_get(self%imports, self%hi_H_ice, H_ice)
        call bsl_from_imports(self)
        call isos_update(self%isos, H_ice, st%t1, self%bsl_loc)
        call exports_state(self, st%t1, with_bsl=.false.)
        call exports_volume(self, st%t1)
    end subroutine fastiso_geo_update_

    subroutine fastiso_geo_restart_write_(self, dir, t)
        ! kryos_startup.f90:421
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        call isos_restart_write(self%isos, restart_file(self%name, dir), t)
    end subroutine fastiso_geo_restart_write_

    subroutine fastiso_geo_restart_read_(self, dir, t)
        ! domain_restart_read (kryos_startup.f90:500-508): through isos_init_state
        ! with use_restart, not a bare isos_restart_read. The core must have put
        ! the imports H_ice (restored ice) and bsl (restored global part) first.
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t

        character(len=:), allocatable :: file
        real(sp), allocatable :: H_ice(:,:)

        file = restart_file(self%name, dir)
        if (len(file) > len(self%isos%par%restart)) then
            call fesm_error("kryos_geo_fastisostasy::restart_read", "restart file name too long for FastIsostasy", &
                            "instance "//trim(self%name)//": `"//file//"`")
        end if
        call require_published(self, self%hi_H_ice, "H_ice", "restart_read")
        call require_published(self, self%hi_bsl,   "bsl",   "restart_read")

        allocate(H_ice(self%nx, self%ny))
        call fields_get(self%imports, self%hi_H_ice, H_ice)
        call bsl_from_imports(self)
        self%isos%par%use_restart = .true.
        self%isos%par%restart     = file
        ! z_bed is not used on a restart: isos_init_state keeps the restored z_bed.
        call isos_init_state(self%isos, self%z_bed_ref_g, H_ice, t, self%bsl_loc)
    end subroutine fastiso_geo_restart_read_

    subroutine fastiso_geo_write_output_(self, file, t, stream, init)
        ! isos.nc and isos_ts.nc of kryos_output.f90:260-277, 383-399; any other
        ! stream: the default writer.
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init

        select case (trim(stream))
        case ("isos")
            if (init) call dims_create(file, self%xc, self%yc, t)
            call isos_write_2D(self%isos, file, t)
        case ("isos_ts")
            if (init) call time_create(file, t)
            call isos_write_1D(self%isos, file, t)
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine fastiso_geo_write_output_

    subroutine fastiso_geo_end_(self)
        class(fastiso_geo_t), intent(inout) :: self
        call isos_end(self%isos)
    end subroutine fastiso_geo_end_

    ! =========================================================================
    ! Regional: helpers
    ! =========================================================================

    subroutine group_read(self, ctx)
        ! Adapter keys of the instance group; the core has validated them.
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        if (nml_has_param(ctx%par_file, ctx%instance, "isos")) &
            call nml_read(ctx%par_file, ctx%instance, "isos", self%group)
    end subroutine group_read

    subroutine bsl_loc_init(self, ctx)
        ! The local bsl copy: method and ocean area of the global group; its
        ! bsl_now and A_ocean_now come from the imports before every use.
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        call nml_read(ctx%par_file, BSL_GROUP, "method", self%bsl_loc%method)
        if (trim(self%bsl_loc%method) == "mixed") call mixed_refused(self%name)
        call nml_read(ctx%par_file, BSL_GROUP, "A_ocean_pd", self%bsl_loc%A_ocean_pd)
        self%bsl_loc%bsl_now      = 0.0_sp
        self%bsl_loc%A_ocean_now  = self%bsl_loc%A_ocean_pd
        self%bsl_loc%use_restart  = .false.
        self%bsl_loc%time         = ctx%time_rel
    end subroutine bsl_loc_init

    subroutine reference_narrow(self, ctx)
        ! The core's reference on the geo grid, narrowed once (sp).
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        if (.not. (allocated(ctx%ref%z_bed) .and. allocated(ctx%ref%H_ice))) then
            call fesm_error("kryos_geo_fastisostasy::init", "reference topography incomplete", &
                            "instance "//trim(self%name)//": z_bed and H_ice are needed on grid `"// &
                            trim(self%grid)//"`")
        end if
        if (any(shape(ctx%ref%z_bed) /= [self%nx, self%ny]) .or. any(shape(ctx%ref%H_ice) /= [self%nx, self%ny])) then
            call fesm_error("kryos_geo_fastisostasy::init", "reference topography has the wrong shape", &
                            "instance "//trim(self%name)//": expected "//to_str(self%nx)//" x "//to_str(self%ny))
        end if
        self%z_bed_ref_g = real(ctx%ref%z_bed, sp)
        self%H_ice_ref_g = real(ctx%ref%H_ice, sp)
    end subroutine reference_narrow

    subroutine bsl_from_imports(self)
        ! The global part's bsl (and ocean area) as the local bsl_class sees it.
        class(fastiso_geo_t), intent(inout) :: self

        call fields_get(self%imports, self%hi_bsl, self%bsl_loc%bsl_now)
        if (wired(self, self%hi_A_ocn)) then
            call fields_get(self%imports, self%hi_A_ocn, self%bsl_loc%A_ocean_now)
        else if (trim(self%bsl_loc%method) == "fastiso") then
            call fesm_error("kryos_geo_fastisostasy::update", "import A_ocn is not wired", &
                            "instance "//trim(self%name)//": bsl method `fastiso` divides by the ocean area")
        end if
    end subroutine bsl_from_imports

    subroutine require_published(self, h, name, what)
        class(fastiso_geo_t), intent(in) :: self
        integer,              intent(in) :: h
        character(len=*),     intent(in) :: name, what
        if (.not. wired(self, h)) then
            call fesm_error("kryos_geo_fastisostasy::"//what, "import "//name//" has not been published", &
                            "instance "//trim(self%name)//": the core exchanges it (restored ice, global part) "// &
                            "before the restart read")
        end if
    end subroutine require_published

    logical function wired(self, h)
        ! Has the import with handle h been published?
        class(fastiso_geo_t), intent(in) :: self
        integer,              intent(in) :: h
        wired = h > 0
        if (wired) wired = self%imports%f(h)%published
    end function wired

    subroutine mixed_refused(instance)
        character(len=*), intent(in) :: instance
        call fesm_error("kryos_geo_fastisostasy::init", "`[barysealevel] method = mixed` is not supported", &
                        "instance "//trim(instance)//": its accumulator needs per-domain bookkeeping (phase 3)")
    end subroutine mixed_refused

    ! =========================================================================
    ! Regional: exports (adapters-design 3.1): sp values widened, formed in sp
    ! =========================================================================

    subroutine exports_state(self, t1, with_bsl)
        ! INSTANT at t1: kryos_coupling.f90:222-223 and check_isostasy_reference
        ! (kryos_coupling.f90:268-273), sums and differences in sp.
        class(fastiso_geo_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        logical,              intent(in), optional :: with_bsl   ! default .true.

        logical :: put_bsl

        put_bsl = .true.
        if (present(with_bsl)) put_bsl = with_bsl
        associate (out => self%isos%out, now => self%isos%now)
            call fields_put(self%exports, self%h_dz_bed,    out%w + out%we,                t1=t1)
            call fields_put(self%exports, self%h_dz_ss,     out%z_ss - now%bsl,            t1=t1)
            call fields_put(self%exports, self%h_z_bed_ref, out%z_bed - out%w - out%we,    t1=t1)
            call fields_put(self%exports, self%h_dz_bed_dt, out%dwdt,                      t1=t1)
            if (put_bsl) call fields_put(self%exports, self%h_bsl, now%bsl,                t1=t1)
        end associate
    end subroutine exports_state

    subroutine exports_volume(self, t1)
        ! INSTANT: the volume change of this update [m3] (the global part's gather source).
        class(fastiso_geo_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        call fields_put(self%exports, self%h_dV_bsl, self%isos%now%deltaV_bsl, t1=t1)
    end subroutine exports_volume

    ! =========================================================================
    ! Regional: output (kryos_output.f90:196-277, 376-399)
    ! =========================================================================

    subroutine dims_create(file, xc, yc, t)
        ! io_dims_init: xc, yc in km and the time axis.
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: xc(:), yc(:), t
        call nc_create(file)
        call nc_write_dim(file, "xc", x=xc, units="km")
        call nc_write_dim(file, "yc", x=yc, units="km")
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
    end subroutine dims_create

    subroutine time_create(file, t)
        ! isos_write_1D_init
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: t
        call nc_create(file)
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
    end subroutine time_create

    subroutine isos_write_2D(isos, file, t)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: t

        integer :: ncid, n

        call nc_open(file, ncid, writable=.true.)
        n = nc_time_index(file, "time", t, ncid)
        call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)
        call var_2D(file, "z_bed", isos%out%z_bed, n, ncid, "m", "Bedrock elevation")
        call var_2D(file, "z_ss",  isos%out%z_ss,  n, ncid, "m", "Sea-surface height")
        call var_2D(file, "w",     isos%out%w,     n, ncid, "m", "Viscous displacement")
        call var_2D(file, "we",    isos%out%we,    n, ncid, "m", "Elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_2D

    subroutine var_2D(file, vnm, var, n, ncid, units, long_name)
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(sp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                      units=units, long_name=long_name, prec_abs=1e-3_sp)
    end subroutine var_2D

    subroutine isos_write_1D(isos, file, t)
        type(isos_class), intent(in) :: isos
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: t

        integer :: ncid, n, np

        call nc_open(file, ncid, writable=.true.)
        n  = nc_time_index(file, "time", t, ncid)
        np = size(isos%out%z_bed)
        call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)
        call var_ts(file, "bsl",        isos%now%bsl,                      n, ncid, "m", "Barystatic sea level")
        call var_ts(file, "z_bed_mean", sum(isos%out%z_bed)/real(np, sp),  n, ncid, "m", "Mean bedrock elevation")
        call var_ts(file, "z_bed_min",  minval(isos%out%z_bed),            n, ncid, "m", "Min bedrock elevation")
        call var_ts(file, "z_bed_max",  maxval(isos%out%z_bed),            n, ncid, "m", "Max bedrock elevation")
        call var_ts(file, "w_mean",     sum(isos%out%w)/real(np, sp),      n, ncid, "m", "Mean viscous displacement")
        call var_ts(file, "we_mean",    sum(isos%out%we)/real(np, sp),     n, ncid, "m", "Mean elastic displacement")
        call nc_close(ncid)
    end subroutine isos_write_1D

    subroutine var_ts(file, vnm, val, n, ncid, units, long_name)
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(sp),         intent(in) :: val
        integer,          intent(in) :: n, ncid
        call nc_write(file, vnm, val, dim1="time", start=[n], count=[1], ncid=ncid, &
                      units=units, long_name=long_name)
    end subroutine var_ts

    ! =========================================================================
    ! Global part: bindings
    ! =========================================================================

    subroutine fastiso_bsl_init_(self, ctx)
        ! yelmox.f90:92-93 (bsl_init, bsl_update at the start time).
        class(fastiso_bsl_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        character(len=NAME_LEN) :: group
        character(len=512) :: restart

        group = BSL_GROUP
        if (nml_has_param(ctx%par_file, ctx%instance, "bsl")) &
            call nml_read(ctx%par_file, ctx%instance, "bsl", group)
        if (trim(group) /= BSL_GROUP) then
            call fesm_error("kryos_geo_fastisostasy::init", "key `bsl` must name group `"//BSL_GROUP//"`", &
                            "instance "//trim(self%name)//": FastIsostasy reads that group name only (got `"// &
                            trim(group)//"`)")
        end if
        call nml_read(ctx%par_file, BSL_GROUP, "restart", restart)
        select case (trim(restart))
        case ("None", "none", "no")
        case default
            call fesm_error("kryos_geo_fastisostasy::init", "`restart` of group `"//BSL_GROUP//"` must be None", &
                            "instance "//trim(self%name)//": restarts come from the bundle (`[kryos] restart`)")
        end select

        self%prec = sp
        call bsl_init(self%bsl, trim(ctx%par_file), ctx%time_rel)
        if (trim(self%bsl%method) == "mixed") call mixed_refused(self%name)
        call bsl_update(self%bsl, ctx%time_rel)
    end subroutine fastiso_bsl_init_

    subroutine fastiso_bsl_init_state_(self, st, pass)
        ! The state of bsl_init + bsl_update, and after a restart read.
        class(fastiso_bsl_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        integer,              intent(in)    :: pass

        select case (pass)
        case (PASS_INIT, PASS_RESTART)
            call exports_bsl(self, st%t1, with_A=.true.)
        end select
    end subroutine fastiso_bsl_init_state_

    subroutine fastiso_bsl_update_(self, st)
        ! yelmox.f90:180: bsl_update on the forcing time, before any domain.
        class(fastiso_bsl_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        call bsl_update(self%bsl, st%t_eval)
        call exports_bsl(self, st%t1, with_A=.true.)
    end subroutine fastiso_bsl_update_

    subroutine fastiso_bsl_reduce_(self, st)
        ! The gathered volume change of the domains (fastisostasy.f90:838-841):
        ! bsl -= dV / A_ocean in sp, method fastiso only; puts bsl.
        class(fastiso_bsl_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        real(sp) :: dV

        call fields_get(self%imports, self%hi_dV_bsl, dV)
        if (trim(self%bsl%method) == "fastiso") self%bsl%bsl_now = self%bsl%bsl_now - dV / self%bsl%A_ocean_now
        call exports_bsl(self, st%t1, with_A=.false.)
    end subroutine fastiso_bsl_reduce_

    subroutine fastiso_bsl_restart_write_(self, dir, t)
        ! kryos_startup.f90:127
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        call bsl_restart_write(self%bsl, restart_file(self%name, dir), t)
    end subroutine fastiso_bsl_restart_write_

    subroutine fastiso_bsl_restart_read_(self, dir, t)
        ! kryos_startup.f90:54; not updated here, the next step's bsl_update does.
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (unused => t); end associate
        call bsl_restart_read(self%bsl, restart_file(self%name, dir))
    end subroutine fastiso_bsl_restart_read_

    subroutine fastiso_bsl_write_output_(self, file, t, stream, init)
        ! The library's own series (bsl, A_ocean); any other stream: the default writer.
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init

        select case (trim(stream))
        case ("bsl_ts")
            if (init) call bsl_write_init(self%bsl, file, t)
            call bsl_write_step(self%bsl, file, t)
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine fastiso_bsl_write_output_

    subroutine fastiso_bsl_end_(self)
        class(fastiso_bsl_t), intent(inout) :: self
        call bsl_end(self%bsl)
    end subroutine fastiso_bsl_end_

    subroutine exports_bsl(self, t1, with_A)
        ! INSTANT at t1: bsl_now [m] and A_ocean_now [m2], widened.
        class(fastiso_bsl_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1
        logical,              intent(in)    :: with_A

        call fields_put(self%exports, self%h_bsl, self%bsl%bsl_now, t1=t1)
        if (with_A) call fields_put(self%exports, self%h_A_ocn, self%bsl%A_ocean_now, t1=t1)
    end subroutine exports_bsl

    function restart_file(name, dir) result(file)
        character(len=*), intent(in) :: name, dir
        character(len=:), allocatable :: file
        file = trim(dir)//"/"//trim(name)//"_restart.nc"
    end function restart_file

end module kryos_geo_fastisostasy
