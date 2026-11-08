module kryos_ice_plastic
    ! Stub of the plastic ice backend, compiled when kryos is built with plastic=0: same module and
    ! type as kryos_ice_plastic.f90, no libplastic. Its constraints are unlinked, so the factory stops
    ! with the rebuild hint. (The margin and tau_b modules, src/backends/ice/plastic/, need no
    ! library and are compiled either way.)

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_ice,       only: ice_t

    implicit none

    private
    public :: plastic_ice_t, plastic_ice_constraints

    type, extends(ice_t) :: plastic_ice_t
    contains
        procedure :: init          => plastic_init_
        procedure :: update        => plastic_update_
        procedure :: restart_write => plastic_restart_write_
        procedure :: restart_read  => plastic_restart_read_
    end type plastic_ice_t

contains

    function plastic_ice_constraints() result(con)
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="plastic=1")
    end function plastic_ice_constraints

    subroutine unavailable(self)               ! safety net: backend_admit stops first
        class(plastic_ice_t), intent(in) :: self
        call fesm_error("kryos_ice_plastic", "instance `"//trim(self%name)// &
                        "` needs backend `plastic`: rebuild kryos with `plastic=1`")
    end subroutine unavailable

    subroutine plastic_init_(self, ctx)
        class(plastic_ice_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self)
    end subroutine plastic_init_

    subroutine plastic_update_(self, st)
        class(plastic_ice_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        associate (unused => st%t1); end associate
        call unavailable(self)
    end subroutine plastic_update_

    subroutine plastic_restart_write_(self, dir, t)
        class(plastic_ice_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self)
    end subroutine plastic_restart_write_

    subroutine plastic_restart_read_(self, dir, t)
        class(plastic_ice_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self)
    end subroutine plastic_restart_read_

end module kryos_ice_plastic
