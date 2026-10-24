module kryos_surface_chion
    ! surface backend chion (adapters-design section 5), built with chion=1,
    ! else the stub. Phase 1 placeholder: type, bindings, constraints; the
    ! bindings stop. It holds a chion_class, so it compiles against chion.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_surface,   only: surface_t
    use chion,           only: chion_class

    implicit none

    private
    public :: chion_surface_t, chion_surface_constraints

    type, extends(surface_t) :: chion_surface_t
        type(chion_class) :: chn                ! the snowpack model
    contains
        procedure :: init          => chion_init_
        procedure :: update        => chion_update_
        procedure :: restart_write => chion_restart_write_
        procedure :: restart_read  => chion_restart_read_
    end type chion_surface_t

contains

    function chion_surface_constraints() result(con)
        ! Regional, any grid, several instances.
        type(constraints_t) :: con
        con = constraints_t(build_flag="chion=1")
    end function chion_surface_constraints

    subroutine pending(instance, what)
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_surface_chion::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(instance))
    end subroutine pending

    subroutine chion_init_(self, ctx)
        class(chion_surface_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine chion_init_

    subroutine chion_update_(self, st)
        class(chion_surface_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine chion_update_

    subroutine chion_restart_write_(self, dir, t)
        class(chion_surface_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine chion_restart_write_

    subroutine chion_restart_read_(self, dir, t)
        class(chion_surface_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine chion_restart_read_

end module kryos_surface_chion
