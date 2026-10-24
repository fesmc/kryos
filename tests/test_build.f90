program test_build
    ! Build smoke test (make check): the library links, the run sequence runs,
    ! and the ice factory allocates the yelmo backend (real or stub; the stub
    ! only stops at init, which is not called here).
    ! placeholder: phase 1 fills this in.

    use precision, only: dp
    use kryos_run, only: run_t, run_init, run_step, run_done, run_end
    use kryos_ice, only: ice_t, ice_create

    implicit none

    type(run_t)               :: run
    class(ice_t), allocatable :: ice

    call run_init(run, "none")
    do while (.not. run_done(run))
        call run_step(run)
    end do
    call run_end(run)
    if (run%nstep /= 3 .or. abs(run%time - 3.0_dp) > 0.0_dp) then
        write(*,'(a)') "test_build: FAILED (run sequence)"
        error stop 1
    end if

    call ice_create(ice, "yelmo", "ice_test")
    if (.not. allocated(ice) .or. trim(ice%backend) /= "yelmo" .or. trim(ice%kind) /= "ice") then
        write(*,'(a)') "test_build: FAILED (ice factory)"
        error stop 1
    end if

    write(*,'(a)') "test_build: OK"

end program test_build
