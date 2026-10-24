module kryos_shelf_marine_shelf
    ! shelf backend marine_shelf (adapters-design section 6): the yelmox
    ! marine_shelf/pico code moves into kryos with this adapter. Phase 1
    ! placeholder: type, bindings, constraints; the bindings stop.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_shelf,     only: shelf_t

    implicit none

    private
    public :: marshelf_shelf_t, marshelf_shelf_constraints

    type, extends(shelf_t) :: marshelf_shelf_t
        character(len=256) :: group = ""        ! the marine_shelf group (adapter state follows)
    contains
        procedure :: init          => marshelf_init_
        procedure :: update        => marshelf_update_
        procedure :: restart_write => marshelf_restart_write_
        procedure :: restart_read  => marshelf_restart_read_
    end type marshelf_shelf_t

contains

    function marshelf_shelf_constraints() result(con)
        ! Regional, any grid, several instances.
        type(constraints_t) :: con
        con = constraints_t()
    end function marshelf_shelf_constraints

    subroutine pending(instance, what)
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_shelf_marine_shelf::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(instance))
    end subroutine pending

    subroutine marshelf_init_(self, ctx)
        class(marshelf_shelf_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine marshelf_init_

    subroutine marshelf_update_(self, st)
        class(marshelf_shelf_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine marshelf_update_

    subroutine marshelf_restart_write_(self, dir, t)
        class(marshelf_shelf_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine marshelf_restart_write_

    subroutine marshelf_restart_read_(self, dir, t)
        class(marshelf_shelf_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine marshelf_restart_read_

end module kryos_shelf_marine_shelf
