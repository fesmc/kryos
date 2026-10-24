module kryos_geo_fastisostasy
    ! Stub of the FastIsostasy geo backend (fastisostasy=0): same module and
    ! types as kryos_geo_fastisostasy.f90, no library; the constraints are
    ! unlinked, so the factory stops with the rebuild hint.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_geo,       only: geo_t

    implicit none

    private
    public :: fastiso_geo_t, fastiso_bsl_t, fastiso_geo_constraints, fastiso_bsl_constraints

    type, extends(geo_t) :: fastiso_geo_t
    contains
        procedure :: init          => fastiso_geo_init_
        procedure :: update        => fastiso_geo_update_
        procedure :: restart_write => fastiso_geo_restart_write_
        procedure :: restart_read  => fastiso_geo_restart_read_
    end type fastiso_geo_t

    type, extends(geo_t) :: fastiso_bsl_t
    contains
        procedure :: init          => fastiso_bsl_init_
        procedure :: update        => fastiso_bsl_update_
        procedure :: restart_write => fastiso_bsl_restart_write_
        procedure :: restart_read  => fastiso_bsl_restart_read_
    end type fastiso_bsl_t

contains

    function fastiso_geo_constraints() result(con)
        ! Regional only, any number of domains.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="fastisostasy=1")
    end function fastiso_geo_constraints

    function fastiso_bsl_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="fastisostasy=1", regional=.false., global=.true., max_instances=1)
    end function fastiso_bsl_constraints

    subroutine unavailable(instance, what)        ! safety net: backend_admit stops first
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_geo_fastisostasy::"//what, "instance `"//trim(instance)// &
                        "` needs backend `fastisostasy`: rebuild kryos with `fastisostasy=1`")
    end subroutine unavailable

    subroutine fastiso_geo_init_(self, ctx)
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine fastiso_geo_init_

    subroutine fastiso_geo_update_(self, st)
        class(fastiso_geo_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine fastiso_geo_update_

    subroutine fastiso_geo_restart_write_(self, dir, t)
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self%name, "restart_write")
    end subroutine fastiso_geo_restart_write_

    subroutine fastiso_geo_restart_read_(self, dir, t)
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self%name, "restart_read")
    end subroutine fastiso_geo_restart_read_

    subroutine fastiso_bsl_init_(self, ctx)
        class(fastiso_bsl_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine fastiso_bsl_init_

    subroutine fastiso_bsl_update_(self, st)
        class(fastiso_bsl_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine fastiso_bsl_update_

    subroutine fastiso_bsl_restart_write_(self, dir, t)
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self%name, "restart_write")
    end subroutine fastiso_bsl_restart_write_

    subroutine fastiso_bsl_restart_read_(self, dir, t)
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call unavailable(self%name, "restart_read")
    end subroutine fastiso_bsl_restart_read_

end module kryos_geo_fastisostasy
