program test_climate_paleoemu_stub
    ! With paleoemu=0 the paleoemu climate backend is its stub: the factory admits nothing and stops
    ! with the rebuild hint, for the global part and for the regional port (child runs).
    ! With paleoemu=1 the test is not built (test_skip in config/common.mk).

    use kryos_component, only: init_ctx_t, census_t, constraints_t
    use kryos_climate,   only: climate_slot_t, climate_create
    use kryos_climate_paleoemu, only: paleoemu_global_constraints, paleoemu_port_constraints

    implicit none

    character(len=*), parameter :: errfile = "test_climate_paleoemu_stub_abort.err"
    integer :: nfail = 0
    type(constraints_t) :: ce, cr

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    ce = paleoemu_global_constraints()
    cr = paleoemu_port_constraints()
    call check("constraints unlinked, build flag paleoemu=1", .not. ce%linked .and. .not. cr%linked .and. &
               trim(ce%build_flag) == "paleoemu=1" .and. trim(cr%build_flag) == "paleoemu=1")
    call check("global part: global scope only, one instance; port: regional, no backend restart", &
               ce%global .and. .not. ce%regional .and. ce%max_instances == 1 .and. cr%regional .and. .not. cr%restart)
    call expect_abort("global", "rebuild kryos with `paleoemu=1`")
    call expect_abort("regional", "rebuild kryos with `paleoemu=1`")

    call execute_command_line("rm -f "//errfile)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_climate_paleoemu_stub: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_climate_paleoemu_stub: OK"

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
        type(climate_slot_t) :: g
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        call get_command_argument(1, arg)
        ctx%instance  = "climate_x"
        ctx%component = "climate"
        ctx%backend   = "paleoemu"
        ctx%scope     = trim(arg)
        call climate_create(g, ctx, census)
    end subroutine abort_case

end program test_climate_paleoemu_stub
