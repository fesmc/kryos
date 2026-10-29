module kryos_geo_vilma2
    ! Stub of the VILMA2 geo backend (vilma=0): same module and types as
    ! kryos_geo_vilma2.f90, no library; the constraints are unlinked, so the
    ! factory stops with the rebuild hint.

    use precision,       only: dp
    use coordinates,     only: grid_class
    use fesm_log,        only: fesm_error
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, NAME_LEN
    use kryos_geo,       only: geo_t

    implicit none

    private
    public :: vilma_earth_t, vilma_geo_t, vilma2_earth_constraints, vilma2_geo_constraints

    type, extends(geo_t) :: vilma_earth_t
        type(grid_class)        :: gauss                 ! as the adapter's: the Gauss grid
        character(len=NAME_LEN) :: group  = "vilma"
        character(len=NAME_LEN) :: ice_bg = "globe_ice"
    contains
        procedure :: init          => earth_init_
        procedure :: init_state    => earth_init_state_
        procedure :: update        => earth_update_
        procedure :: restart_write => earth_restart_write_
        procedure :: restart_read  => earth_restart_read_
        procedure :: end           => earth_end_
    end type vilma_earth_t

    type, extends(geo_t) :: vilma_geo_t
    contains
        procedure :: init          => port_init_
        procedure :: init_state    => port_init_state_
        procedure :: update        => port_update_
        procedure :: restart_write => port_restart_write_
        procedure :: restart_read  => port_restart_read_
    end type vilma_geo_t

contains

    function vilma2_earth_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="vilma=1", regional=.false., global=.true., max_instances=1)
    end function vilma2_earth_constraints

    function vilma2_geo_constraints() result(con)
        ! Regional only (the port), any number of domains; no backend restart file.
        type(constraints_t) :: con
        con = constraints_t(linked=.false., build_flag="vilma=1", restart=.false.)
    end function vilma2_geo_constraints

    subroutine unavailable(instance, what)        ! safety net: backend_admit stops first
        character(len=*), intent(in) :: instance, what
        call fesm_error("kryos_geo_vilma2::"//what, "instance `"//trim(instance)// &
                        "` needs backend `vilma2`: rebuild kryos with `vilma=1`")
    end subroutine unavailable

    subroutine earth_init_(self, ctx)
        class(vilma_earth_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine earth_init_

    subroutine earth_init_state_(self, st, pass)
        class(vilma_earth_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        integer, intent(in) :: pass
        associate (unused => st%t1 + pass); end associate
        call unavailable(self%name, "init_state")
    end subroutine earth_init_state_

    subroutine earth_update_(self, st)
        class(vilma_earth_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine earth_update_

    subroutine earth_restart_write_(self, dir, t)
        class(vilma_earth_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_write")
    end subroutine earth_restart_write_

    subroutine earth_restart_read_(self, dir, t)
        class(vilma_earth_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_read")
    end subroutine earth_restart_read_

    subroutine earth_end_(self)
        class(vilma_earth_t), intent(inout) :: self
        call unavailable(self%name, "end")
    end subroutine earth_end_

    subroutine port_init_(self, ctx)
        class(vilma_geo_t), intent(inout) :: self
        type(init_ctx_t), intent(in) :: ctx
        associate (unused => ctx%time); end associate
        call unavailable(self%name, "init")
    end subroutine port_init_

    subroutine port_init_state_(self, st, pass)
        class(vilma_geo_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        integer, intent(in) :: pass
        associate (unused => st%t1 + pass); end associate
        call unavailable(self%name, "init_state")
    end subroutine port_init_state_

    subroutine port_update_(self, st)
        class(vilma_geo_t), intent(inout) :: self
        type(step_time_t), intent(in) :: st
        associate (unused => st%t1); end associate
        call unavailable(self%name, "update")
    end subroutine port_update_

    subroutine port_restart_write_(self, dir, t)
        class(vilma_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_write")
    end subroutine port_restart_write_

    subroutine port_restart_read_(self, dir, t)
        class(vilma_geo_t), intent(inout) :: self
        character(len=*), intent(in) :: dir
        real(dp), intent(in) :: t
        associate (unused => t + len(dir)); end associate
        call unavailable(self%name, "restart_read")
    end subroutine port_restart_read_

end module kryos_geo_vilma2
