module kryos_geo_fastisostasy
    ! geo backend FastIsostasy (adapters-design section 3), built with
    ! fastisostasy=1, else the stub: fastiso_geo_t (regional, isos_class) and
    ! fastiso_bsl_t (global part, bsl_class); the factory picks by scope.
    ! Phase 1 placeholder: types, bindings, constraints; bindings stop.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_geo,       only: geo_t
    use fastisostasy,    only: isos_class
    use barysealevel,    only: bsl_class

    implicit none

    private
    public :: fastiso_geo_t, fastiso_bsl_t, fastiso_geo_constraints, fastiso_bsl_constraints

    type, extends(geo_t) :: fastiso_geo_t
        type(isos_class) :: isos               ! the domain's FastIsostasy state
        type(bsl_class)  :: bsl_loc            ! copy of the global bsl state
    contains
        procedure :: init          => fastiso_geo_init_
        procedure :: update        => fastiso_geo_update_
        procedure :: restart_write => fastiso_geo_restart_write_
        procedure :: restart_read  => fastiso_geo_restart_read_
    end type fastiso_geo_t

    type, extends(geo_t) :: fastiso_bsl_t
        type(bsl_class) :: bsl                 ! barystatic sea level (global part)
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
        con = constraints_t(build_flag="fastisostasy=1")
    end function fastiso_geo_constraints

    function fastiso_bsl_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(build_flag="fastisostasy=1", regional=.false., global=.true., max_instances=1)
    end function fastiso_bsl_constraints

    subroutine pending(instance, what)
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_geo_fastisostasy::"//what, "not implemented (phase 1 adapter pending)", &
                        "instance "//trim(instance))
    end subroutine pending

    subroutine fastiso_geo_init_(self, ctx)
        class(fastiso_geo_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine fastiso_geo_init_

    subroutine fastiso_geo_update_(self, st)
        class(fastiso_geo_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine fastiso_geo_update_

    subroutine fastiso_geo_restart_write_(self, dir, t)
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine fastiso_geo_restart_write_

    subroutine fastiso_geo_restart_read_(self, dir, t)
        class(fastiso_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine fastiso_geo_restart_read_

    subroutine fastiso_bsl_init_(self, ctx)
        class(fastiso_bsl_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call pending(self%name, "init")
    end subroutine fastiso_bsl_init_

    subroutine fastiso_bsl_update_(self, st)
        class(fastiso_bsl_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call pending(self%name, "update")
    end subroutine fastiso_bsl_update_

    subroutine fastiso_bsl_restart_write_(self, dir, t)
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_write")
    end subroutine fastiso_bsl_restart_write_

    subroutine fastiso_bsl_restart_read_(self, dir, t)
        class(fastiso_bsl_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => len(dir) + t); end associate
        call pending(self%name, "restart_read")
    end subroutine fastiso_bsl_restart_read_

end module kryos_geo_fastisostasy
