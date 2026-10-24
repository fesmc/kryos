module kryos_forcing_tsgen
    ! forcing backend tsgen (adapters-design section 8, contracts section 7):
    ! the coefficient registry on fesm-utils tsgen; no stub (fesm-utils is
    ! always linked). Phase 1 placeholder: type, bindings, constraints; stop.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_forcing,   only: forcing_t
    use tsgen,           only: tsgen_class

    implicit none

    private
    public :: tsgen_forcing_t, tsgen_forcing_constraints

    type, extends(forcing_t) :: tsgen_forcing_t
        type(tsgen_class), allocatable :: tsg(:)   ! one generator per tsgen coefficient
    contains
        procedure :: init          => tsgen_init_
        procedure :: update        => tsgen_update_
        procedure :: restart_write => tsgen_restart_write_
        procedure :: restart_read  => tsgen_restart_read_
    end type tsgen_forcing_t

contains

    function tsgen_forcing_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(regional=.false., global=.true., max_instances=1)
    end function tsgen_forcing_constraints

    subroutine pending(instance, what)
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_forcing_tsgen::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(instance))
    end subroutine pending

    subroutine tsgen_init_(self, ctx)
        class(tsgen_forcing_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine tsgen_init_

    subroutine tsgen_update_(self, st)
        class(tsgen_forcing_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine tsgen_update_

    subroutine tsgen_restart_write_(self, dir, t)
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine tsgen_restart_write_

    subroutine tsgen_restart_read_(self, dir, t)
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine tsgen_restart_read_

end module kryos_forcing_tsgen
