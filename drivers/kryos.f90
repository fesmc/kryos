program kryos
    ! The kryos driver: one program for N domains (design decision 10). The
    ! coupling sequence lives in the library (kryos_run), not here.

    use kryos_run, only: run_t, run_init, run_step, run_done, run_end

    implicit none

    type(run_t)        :: run
    character(len=512) :: par_file

    par_file = "par/kryos.nml"
    if (command_argument_count() >= 1) call get_command_argument(1, par_file)

    call run_init(run, par_file)
    do while (.not. run_done(run))
        call run_step(run)
    end do
    call run_end(run)

end program kryos
