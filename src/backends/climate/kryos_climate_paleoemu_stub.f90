module kryos_climate_paleoemu
    ! Stub of the paleoemu climate backend (paleoemu=0): same module and types as
    ! kryos_climate_paleoemu.f90, no library; the constraints are unlinked, so the
    ! factory stops with the rebuild hint.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t
    use kryos_climate,   only: climate_t

    implicit none

    private
    public :: paleoemu_climate_t, paleoemu_port_t, paleoemu_global_constraints, paleoemu_port_constraints

    type, extends(climate_t) :: paleoemu_climate_t
    contains
        procedure :: init          => g_init_
        procedure :: init_state    => g_init_state_
        procedure :: update        => g_update_
        procedure :: restart_write => g_restart_write_
        procedure :: restart_read  => g_restart_read_
    end type paleoemu_climate_t

    type, extends(climate_t) :: paleoemu_port_t
    contains
        procedure :: init          => p_init_
        procedure :: init_state    => p_init_state_
        procedure :: update        => p_update_
        procedure :: restart_write => p_restart_write_
        procedure :: restart_read  => p_restart_read_
    end type paleoemu_port_t

contains

    function paleoemu_global_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="paleoemu=1", regional=.false., global=.true., max_instances=1)
    end function paleoemu_global_constraints

    function paleoemu_port_constraints() result(con)
        ! Regional only (the port), any number of domains; no backend restart file.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="paleoemu=1", restart=.false.)
    end function paleoemu_port_constraints

    subroutine unavailable(instance, what)        ! safety net: backend_admit stops first
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_climate_paleoemu::"//what, "instance `"//trim(instance)// &
                        "` needs backend `paleoemu`: rebuild kryos with `paleoemu=1`")
    end subroutine unavailable

    subroutine g_init_(self, ctx)
        class(paleoemu_climate_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine g_init_

    subroutine g_init_state_(self, st, pass)
        class(paleoemu_climate_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        integer, intent(in) :: pass
        associate (unused => st%t1 + pass); end associate
        call unavailable(self%name, "init_state")
    end subroutine g_init_state_

    subroutine g_update_(self, st)
        class(paleoemu_climate_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine g_update_

    subroutine g_restart_write_(self, dir, t)
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_write")
    end subroutine g_restart_write_

    subroutine g_restart_read_(self, dir, t)
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_read")
    end subroutine g_restart_read_

    subroutine p_init_(self, ctx)
        class(paleoemu_port_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine p_init_

    subroutine p_init_state_(self, st, pass)
        class(paleoemu_port_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        integer, intent(in) :: pass
        associate (unused => st%t1 + pass); end associate
        call unavailable(self%name, "init_state")
    end subroutine p_init_state_

    subroutine p_update_(self, st)
        class(paleoemu_port_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine p_update_

    subroutine p_restart_write_(self, dir, t)
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_write")
    end subroutine p_restart_write_

    subroutine p_restart_read_(self, dir, t)
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_read")
    end subroutine p_restart_read_

end module kryos_climate_paleoemu
