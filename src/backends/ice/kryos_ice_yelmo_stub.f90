module kryos_ice_yelmo
    ! Stub for the Yelmo ice backend, compiled when kryos is built with
    ! yelmo=0: same module and type as kryos_ice_yelmo.f90, no libyelmo.
    ! Selecting backend = "yelmo" then stops at init.

    use precision,      only: dp
    use coordinates,    only: grid_class
    use phys_constants, only: phys_const_class
    use kryos_ice,      only: ice_t

    implicit none

    private
    public :: yelmo_ice_t

    type, extends(ice_t) :: yelmo_ice_t
    contains
        procedure :: init          => yelmo_ice_init
        procedure :: update        => yelmo_ice_update
        procedure :: restart_write => yelmo_ice_restart_write
        procedure :: restart_read  => yelmo_ice_restart_read
    end type yelmo_ice_t

contains

    subroutine unavailable(self)
        class(yelmo_ice_t), intent(in) :: self
        write(*,'(a)') "kryos: instance '"//trim(self%name)//"' needs the yelmo backend; "// &
                       "rebuild kryos with yelmo=1"
        error stop 1
    end subroutine unavailable

    subroutine yelmo_ice_init(self, par_file, group, grid, time, cnst)
        class(yelmo_ice_t),     intent(inout) :: self
        character(len=*),       intent(in)    :: par_file
        character(len=*),       intent(in)    :: group
        type(grid_class),       intent(in)    :: grid
        real(dp),               intent(in)    :: time
        type(phys_const_class), intent(in)    :: cnst
        call unavailable(self)
    end subroutine yelmo_ice_init

    subroutine yelmo_ice_update(self, t_target)
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t_target
        call unavailable(self)
    end subroutine yelmo_ice_update

    subroutine yelmo_ice_restart_write(self, dir, time)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: time
        call unavailable(self)
    end subroutine yelmo_ice_restart_write

    subroutine yelmo_ice_restart_read(self, dir, time)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: time
        call unavailable(self)
    end subroutine yelmo_ice_restart_read

end module kryos_ice_yelmo
