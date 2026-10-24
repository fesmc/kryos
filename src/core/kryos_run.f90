module kryos_run
    ! The run: clock, globe, domains, output schedules, restart bundle, and the
    ! coupling sequence (design decision 10). The driver only calls run_init,
    ! loops run_step until run_done, then run_end; tests and other hosts call
    ! the same sequence.
    !
    ! placeholder: phase 1 fills this in.

    use precision, only: dp

    implicit none

    private
    public :: run_t, run_init, run_step, run_done, run_end

    type :: run_t
        real(dp) :: time     = 0.0_dp       ! [yr]
        real(dp) :: time_end = 0.0_dp       ! [yr]
        real(dp) :: dt       = 1.0_dp       ! [yr]
        integer  :: nstep    = 0
    end type run_t

contains

    subroutine run_init(run, par_file)
        type(run_t),      intent(out) :: run
        character(len=*), intent(in)  :: par_file
        run%time     = 0.0_dp
        run%time_end = 3.0_dp
        run%dt       = 1.0_dp
        run%nstep    = 0
        write(*,'(a)') "kryos: run_init (placeholder) par_file="//trim(par_file)
    end subroutine run_init

    subroutine run_step(run)
        type(run_t), intent(inout) :: run
        run%time  = run%time + run%dt
        run%nstep = run%nstep + 1
        write(*,'(a,f10.2)') "kryos: run_step (placeholder) time=", run%time
    end subroutine run_step

    logical function run_done(run)
        type(run_t), intent(in) :: run
        run_done = run%time >= run%time_end
    end function run_done

    subroutine run_end(run)
        type(run_t), intent(inout) :: run
        write(*,'(a,i0,a)') "kryos: run_end (placeholder) after ", run%nstep, " steps"
    end subroutine run_end

end module kryos_run
