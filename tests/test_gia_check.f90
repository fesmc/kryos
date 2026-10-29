module test_gia_toys
    ! Toy geo and ice instances that only exist to be admitted with a backend name and
    ! capabilities: the checks of kryos_gia_check read names, capabilities and options.
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_geo,       only: geo_t
    use kryos_ice,       only: ice_t
    implicit none
    private

    type, extends(geo_t), public :: toy_geo_t
    contains
        procedure :: init => toy_init_geo
        procedure :: update => toy_update_geo
        procedure :: restart_write => toy_rw_geo
        procedure :: restart_read => toy_rr_geo
    end type toy_geo_t

    type, extends(ice_t), public :: toy_ice_t
    contains
        procedure :: init => toy_init_ice
        procedure :: update => toy_update_ice
        procedure :: restart_write => toy_rw_ice
        procedure :: restart_read => toy_rr_ice
    end type toy_ice_t

contains

    subroutine toy_init_geo(self, ctx)
        class(toy_geo_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        associate (u => ctx%time + self%nx); end associate
    end subroutine toy_init_geo

    subroutine toy_update_geo(self, st)
        class(toy_geo_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_update_geo

    subroutine toy_rw_geo(self, dir, t)
        class(toy_geo_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw_geo

    subroutine toy_rr_geo(self, dir, t)
        class(toy_geo_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rr_geo

    subroutine toy_init_ice(self, ctx)
        class(toy_ice_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        associate (u => ctx%time + self%nx); end associate
    end subroutine toy_init_ice

    subroutine toy_update_ice(self, st)
        class(toy_ice_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_update_ice

    subroutine toy_rw_ice(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw_ice

    subroutine toy_rr_ice(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rr_ice

end module test_gia_toys

program test_gia_check
    ! The run-level contract checks of the VILMA2 geo backend (plan 3.2, kryos_gia_check) on
    ! toy instances: a run without vilma2 passes untouched; a consistent vilma2 run passes; each
    ! way to be inconsistent stops with a named message (the program re-runs itself with an
    ! argument per case, as test_wiring). Needs no backend library.

    use precision,       only: dp
    use kryos_config,    only: config_t, config_init
    use kryos_component, only: init_ctx_t, census_t, constraints_t
    use kryos_geo,       only: geo_t, geo_adopt
    use kryos_ice,       only: ice_t, ice_adopt
    use kryos_globe,     only: globe_t
    use kryos_domain,    only: domain_t, C_GEO
    use kryos_gia_check, only: gia_check
    use test_gia_toys,   only: toy_geo_t, toy_ice_t

    implicit none

    character(len=*), parameter :: PAR = "test_gia_check.nml", errfile = "test_gia_check_abort.err"
    integer :: nfail = 0
    type(config_t) :: cfg

    call par_write()
    call config_init(cfg, PAR)

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    call test_passes()
    call expect_abort("mix_global_fastiso", "FastIsostasy together with VILMA2")
    call expect_abort("mix_port_fastiso", "FastIsostasy together with VILMA2")
    call expect_abort("port_no_global", "has no global geo with vilma2")
    call expect_abort("domain_no_geo", "has no geo instance")
    call expect_abort("domain_no_ice", "has no ice instance")
    call expect_abort("no_flotation", "has no flotation binding")
    call expect_abort("no_footprint", "does not export `footprint`")
    call expect_abort("no_h_ice_pd", "does not export `H_ice_pd`")
    call expect_abort("dt_differs", "differs from the global geo")
    call expect_abort("no_spinup", "has no VILMA2 spin-up before its seed")

    call execute_command_line("rm -f "//errfile//" "//PAR)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_gia_check: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_gia_check: OK"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine adopt_geo(slot, name, backend, scope, spinup)
        ! A toy geo admitted under this backend name (the check reads names only); spinup: what the
        ! backend says about its spin-up (default: it has one).
        use kryos_geo, only: geo_slot_t
        type(geo_slot_t), intent(out) :: slot
        character(len=*), intent(in)  :: name, backend, scope
        logical,          intent(in), optional :: spinup
        class(geo_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        ctx%instance = name
        ctx%backend  = backend
        ctx%scope    = scope
        allocate(toy_geo_t :: proto)
        call geo_adopt(slot, proto, constraints_t(regional=.true., global=.true.), ctx, census)
        slot%c%opt%host_spinup = .true.
        if (present(spinup)) slot%c%opt%host_spinup = spinup
    end subroutine adopt_geo

    subroutine make_domain(dom, name, geo_inst, geo_backend, ice, flot, foot, pd)
        ! A domain with a geo of this backend (geo_inst = "": none) and an ice with these capabilities
        ! (ice = .false.: none).
        type(domain_t),   intent(out) :: dom
        character(len=*), intent(in)  :: name, geo_inst, geo_backend
        logical,          intent(in)  :: ice, flot, foot, pd
        class(ice_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        dom%name = name
        dom%inst(C_GEO)    = geo_inst
        dom%backend(C_GEO) = geo_backend
        if (.not. ice) return
        ctx%instance = "ice_"//name
        ctx%backend  = "yelmo"
        allocate(toy_ice_t :: proto)
        call ice_adopt(dom%ice, proto, constraints_t(), ctx, census)
        dom%ice%c%caps%has_flotation = flot
        dom%ice%c%opt%footprint      = foot
        dom%ice%c%opt%H_ice_pd       = pd
    end subroutine make_domain

    subroutine test_passes()
        type(globe_t)  :: globe
        type(domain_t) :: dom(2)

        write(*,"(a)") "runs that pass"
        call gia_check(globe, dom(1:0), cfg)
        call check("no globe, no domain", .true.)

        call adopt_geo(globe%geo, "geo_glob", "fastisostasy", "global")
        call make_domain(dom(1), "NH",  "geo_nh",  "fastisostasy", .true., .false., .false., .false.)
        call make_domain(dom(2), "ANT", "geo_ant", "fastisostasy", .true., .false., .false., .false.)
        call gia_check(globe, dom, cfg)
        call check("fastisostasy everywhere: untouched (the ice needs no vilma2 capabilities)", .true.)

        call adopt_geo(globe%geo, "geo_glob", "vilma2", "global")
        call make_domain(dom(1), "NH",  "geo_nh",  "vilma2", .true., .true., .true., .true.)
        call make_domain(dom(2), "ANT", "geo_ant", "vilma2", .true., .true., .true., .true.)
        call gia_check(globe, dom, cfg)
        call check("vilma2 globe and two vilma2 ports on full ice instances (dt of geo_ant = the global's)", .true.)

        call adopt_geo(globe%geo, "geo_open", "vilma2", "global", spinup=.false.)
        call gia_check(globe, dom, cfg)
        call check("no spin-up with allow_no_spinup = True in the instance group", .true.)
    end subroutine test_passes

    subroutine expect_abort(case_name, text)
        character(len=*), intent(in) :: case_name, text
        character(len=512) :: exe, line
        character(len=:), allocatable :: all_
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > /dev/null 2> "//errfile, exitstat=es)
        all_ = ""
        open(newunit=u, file=errfile, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(u, "(a)", iostat=ios) line
                if (ios /= 0) exit
                all_ = all_//trim(line)//" "
            end do
            close(u)
        end if
        call check(case_name//": aborts (exit 1) with '"//text//"'", es == 1 .and. index(all_, text) > 0)
    end subroutine expect_abort

    subroutine abort_case()
        character(len=64) :: arg
        type(globe_t)  :: globe
        type(domain_t) :: dom(1)

        call get_command_argument(1, arg)
        call adopt_geo(globe%geo, "geo_glob", "vilma2", "global")
        call make_domain(dom(1), "NH", "geo_nh", "vilma2", .true., .true., .true., .true.)
        select case (trim(arg))
        case ("mix_global_fastiso")
            call make_domain(dom(1), "NH", "geo_nh", "fastisostasy", .true., .true., .true., .true.)
        case ("mix_port_fastiso")
            call adopt_geo(globe%geo, "geo_glob", "fastisostasy", "global")
        case ("port_no_global")
            deallocate(globe%geo%c)
        case ("domain_no_geo")
            call make_domain(dom(1), "NH", "", "", .true., .true., .true., .true.)
        case ("domain_no_ice")
            call make_domain(dom(1), "NH", "geo_nh", "vilma2", .false., .false., .false., .false.)
        case ("no_flotation")
            call make_domain(dom(1), "NH", "geo_nh", "vilma2", .true., .false., .true., .true.)
        case ("no_footprint")
            call make_domain(dom(1), "NH", "geo_nh", "vilma2", .true., .true., .false., .true.)
        case ("no_h_ice_pd")
            call make_domain(dom(1), "NH", "geo_nh", "vilma2", .true., .true., .true., .false.)
        case ("no_spinup")
            call adopt_geo(globe%geo, "geo_glob", "vilma2", "global", spinup=.false.)
        case ("dt_differs")
            call make_domain(dom(1), "NH", "geo_bad", "vilma2", .true., .true., .true., .true.)
        end select
        call gia_check(globe, dom, cfg)
    end subroutine abort_case

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&kryos", ' log_level = "warn"', "/", &
                       "&geo_glob", ' backend = "vilma2"', " dt = 10.0", "/", &
                       "&geo_nh",   ' backend = "vilma2"', "/", &
                       "&geo_ant",  ' backend = "vilma2"', " dt = 10.0", "/", &
                       "&geo_bad",  ' backend = "vilma2"', " dt = 5.0", "/", &
                       "&geo_open", ' backend = "vilma2"', " dt = 10.0", " allow_no_spinup = True", "/"
        close(u)
    end subroutine par_write

end program test_gia_check
