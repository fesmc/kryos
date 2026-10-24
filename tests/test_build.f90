program test_build
    ! Build smoke test (make check): the library links, the run sequence runs,
    ! and the ice factory admits the yelmo backend when it is linked (a stub
    ! stops in the factory, so with yelmo=0 that check is skipped).

    use precision, only: dp
    use kryos_run, only: run_t, run_init, run_step, run_done, run_end
    use kryos_component, only: init_ctx_t, census_t, census_count, constraints_t
    use kryos_ice,       only: ice_slot_t, ice_create, ice_present
    use kryos_ice_yelmo, only: yelmo_ice_constraints

    implicit none

    type(run_t)      :: run
    type(ice_slot_t) :: ice
    type(init_ctx_t) :: ctx
    type(census_t)   :: census
    type(constraints_t) :: con

    call run_init(run, "none")
    do while (.not. run_done(run))
        call run_step(run)
    end do
    call run_end(run)
    if (run%nstep /= 3 .or. abs(run%time - 3.0_dp) > 0.0_dp) then
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

    write(*,'(a)') "test_build: OK"

end program test_build
