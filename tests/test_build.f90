program test_build
    ! Build smoke test (make check): the library links, the run sequence runs
    ! (a clock only, no domain: test_run_clock covers it), and the ice factory
    ! admits the yelmo backend when it is linked (a stub stops in the factory,
    ! so with yelmo=0 that check is skipped). Run from the repo root.

    use precision, only: dp
    use fesm_clock, only: clock_time
    use kryos_run, only: run_t, run_init, run_step, run_done, run_end
    use kryos_component, only: init_ctx_t, census_t, census_count, constraints_t
    use kryos_ice,       only: ice_slot_t, ice_create, ice_present
    use kryos_ice_yelmo, only: yelmo_ice_constraints
    use kryos_geo,       only: geo_slot_t, geo_create, geo_present
    use kryos_geo_vilma2, only: vilma2_earth_constraints

    implicit none

    type(run_t)      :: run
    type(ice_slot_t) :: ice
    type(geo_slot_t) :: geo
    type(init_ctx_t) :: ctx
    type(census_t)   :: census
    type(constraints_t) :: con
    integer :: u, nstep

    call execute_command_line("mkdir -p tmp/test_build")
    open(newunit=u, file="tmp/test_build/par.nml", status="replace", action="write")
    write(u,'(a)') '&kryos'
    write(u,'(a)') '    log_level = "warn"'
    write(u,'(a)') '/'
    write(u,'(a)') '&clock'
    write(u,'(a)') '    time_init = 0.0'
    write(u,'(a)') '    time_end  = 30.0'
    write(u,'(a)') '    dtt       = 10.0'
    write(u,'(a)') '/'
    close(u)

    nstep = 0
    call run_init(run, "tmp/test_build/par.nml")
    do while (.not. run_done(run))
        call run_step(run)
        nstep = nstep + 1
    end do
    call run_end(run)
    if (nstep /= 3 .or. abs(clock_time(run%clk) - 30.0_dp) > 0.0_dp) then
        write(*,'(a)') "test_build: FAILED (run sequence)"
        error stop 1
    end if

    con = yelmo_ice_constraints()
    if (con%linked) then
        ctx%instance = "ice_test"
        ctx%backend  = "yelmo"
        call ice_create(ice, ctx, census)
        if (.not. ice_present(ice) .or. trim(ice%c%backend) /= "yelmo" .or. trim(ice%c%kind) /= "ice" &
            .or. census_count(census, "ice:yelmo:regional") /= 1) then
            write(*,'(a)') "test_build: FAILED (ice factory)"
            error stop 1
        end if
    else
        write(*,'(a)') "test_build: ice factory skipped (yelmo stub)"
    end if

    con = vilma2_earth_constraints()
    if (con%linked) then
        ctx%instance = "geo_test"
        ctx%backend  = "vilma2"
        ctx%scope    = "global"
        census = census_t()
        call geo_create(geo, ctx, census)
        if (.not. geo_present(geo) .or. trim(geo%c%backend) /= "vilma2" .or. trim(geo%c%kind) /= "geo" &
            .or. census_count(census, "geo:vilma2:global") /= 1) then
            write(*,'(a)') "test_build: FAILED (geo factory, vilma2)"
            error stop 1
        end if
    else
        write(*,'(a)') "test_build: geo factory (vilma2) skipped (stub; test_geo_vilma2_stub covers it)"
    end if

    write(*,'(a)') "test_build: OK"

end program test_build
