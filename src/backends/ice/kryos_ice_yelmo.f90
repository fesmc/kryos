module kryos_ice_yelmo
    ! ice backend Yelmo (adapters-design section 2), built with yelmo=1;
    ! otherwise kryos_ice_yelmo_stub.f90 provides this module.
    ! Phase 1 placeholder: type, bindings and constraints; every binding stops
    ! until the adapter lands. It holds a yelmo_class, so it compiles against Yelmo.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_ice,       only: ice_t
    use yelmo,           only: yelmo_class

    implicit none

    private
    public :: yelmo_ice_t, yelmo_ice_constraints

    type, extends(ice_t) :: yelmo_ice_t
        type(yelmo_class) :: ylmo
    contains
        procedure :: init          => yelmo_init_
        procedure :: update        => yelmo_update_
        procedure :: restart_write => yelmo_restart_write_
        procedure :: restart_read  => yelmo_restart_read_
    end type yelmo_ice_t

contains

    function yelmo_ice_constraints() result(con)
        ! Any grid, several instances, regional, with a restart file.
        type(constraints_t) :: con
        con = constraints_t(build_flag="yelmo=1")
    end function yelmo_ice_constraints

    subroutine pending(self, what)
        class(yelmo_ice_t), intent(in) :: self
        character(len=*),   intent(in) :: what
        call fesm_error("kryos_ice_yelmo::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(self%name))
    end subroutine pending

    subroutine yelmo_init_(self, ctx)
        class(yelmo_ice_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        call pending(self, "init")              ! placeholder: adapters-design 2.2
    end subroutine yelmo_init_

    subroutine yelmo_update_(self, st)
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        associate (unused => st%t1); end associate
        call pending(self, "update")            ! placeholder: adapters-design 2.4
    end subroutine yelmo_update_

    subroutine yelmo_restart_write_(self, dir, t)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call pending(self, "restart_write")
    end subroutine yelmo_restart_write_

    subroutine yelmo_restart_read_(self, dir, t)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call pending(self, "restart_read")
    end subroutine yelmo_restart_read_

end module kryos_ice_yelmo
