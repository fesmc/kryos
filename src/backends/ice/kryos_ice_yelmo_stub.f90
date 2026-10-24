module kryos_ice_yelmo
    ! Stub of the Yelmo ice backend, compiled when kryos is built with
    ! yelmo=0: same module and type as kryos_ice_yelmo.f90, no libyelmo. Its
    ! constraints are unlinked, so the factory stops with the rebuild hint.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_ice,       only: ice_t

    implicit none

    private
    public :: yelmo_ice_t, yelmo_ice_constraints

    type, extends(ice_t) :: yelmo_ice_t
    contains
        procedure :: init          => yelmo_init_
        procedure :: init_state    => yelmo_init_state_
        procedure :: update        => yelmo_update_
        procedure :: restart_write => yelmo_restart_write_
        procedure :: restart_read  => yelmo_restart_read_
        procedure :: write_output  => yelmo_write_output_
        procedure :: end           => yelmo_end_
    end type yelmo_ice_t

contains

    function yelmo_ice_constraints() result(con)
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="yelmo=1")
    end function yelmo_ice_constraints

    subroutine unavailable(self)               ! safety net: backend_admit stops first
        class(yelmo_ice_t), intent(in) :: self
        call fesm_error("kryos_ice_yelmo", "instance `"//trim(self%name)// &
                        "` needs backend `yelmo`: rebuild kryos with `yelmo=1`")
    end subroutine unavailable

    subroutine yelmo_init_(self, ctx)
        class(yelmo_ice_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self)
    end subroutine yelmo_init_

    subroutine yelmo_init_state_(self, st, pass)
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        integer,            intent(in)    :: pass
        associate (unused => st%t1 + pass); end associate
        call unavailable(self)
    end subroutine yelmo_init_state_

    subroutine yelmo_update_(self, st)
        class(yelmo_ice_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        associate (unused => st%t1); end associate
        call unavailable(self)
    end subroutine yelmo_update_

    subroutine yelmo_restart_write_(self, dir, t)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self)
    end subroutine yelmo_restart_write_

    subroutine yelmo_restart_read_(self, dir, t)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self)
    end subroutine yelmo_restart_read_

    subroutine yelmo_write_output_(self, file, t, stream, init)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init
        associate (unused => len(file) + len(stream) + t); end associate
        associate (unused => init); end associate
        call unavailable(self)
    end subroutine yelmo_write_output_

    subroutine yelmo_end_(self)
        class(yelmo_ice_t), intent(inout) :: self
        call unavailable(self)
    end subroutine yelmo_end_

end module kryos_ice_yelmo
