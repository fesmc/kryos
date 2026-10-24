module kryos_climate_snapshots
    ! climate backend snapshots (adapters-design section 4): the yelmox
    ! snapshots/climate_terms code moves into kryos as private modules of this
    ! backend. Phase 1 placeholder: type, bindings, constraints; bindings stop.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_climate,   only: climate_t

    implicit none

    private
    public :: snapshots_climate_t, snapshots_climate_constraints

    type, extends(climate_t) :: snapshots_climate_t
        character(len=256) :: group = ""        ! the snapshots group (adapter state follows)
    contains
        procedure :: init          => snapshots_init_
        procedure :: update        => snapshots_update_
        procedure :: restart_write => snapshots_restart_write_
        procedure :: restart_read  => snapshots_restart_read_
    end type snapshots_climate_t

contains

    function snapshots_climate_constraints() result(con)
        ! Regional, any grid, several instances.
        type(constraints_t) :: con
        con = constraints_t()
    end function snapshots_climate_constraints

    subroutine pending(instance, what)
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_climate_snapshots::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(instance))
    end subroutine pending

    subroutine snapshots_init_(self, ctx)
        class(snapshots_climate_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine snapshots_init_

    subroutine snapshots_update_(self, st)
        class(snapshots_climate_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine snapshots_update_

    subroutine snapshots_restart_write_(self, dir, t)
        class(snapshots_climate_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine snapshots_restart_write_

    subroutine snapshots_restart_read_(self, dir, t)
        class(snapshots_climate_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine snapshots_restart_read_

end module kryos_climate_snapshots
