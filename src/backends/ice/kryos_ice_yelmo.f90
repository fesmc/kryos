module kryos_ice_yelmo
    ! ice backend: Yelmo. Built with yelmo=1 (config/common.mk); otherwise
    ! kryos_ice_yelmo_stub.f90 provides this module.
    !
    ! placeholder: phase 1 fills this in. It already holds the Yelmo state, so
    ! the build compiles and links against libyelmo.

    use precision,      only: dp
    use coordinates,    only: grid_class
    use phys_constants, only: phys_const_class
    use kryos_ice,      only: ice_t
    use yelmo,          only: yelmo_class

    implicit none

    private
    public :: yelmo_ice_t

    type, extends(ice_t) :: yelmo_ice_t
        type(yelmo_class) :: ylmo
    contains
        procedure :: init          => yelmo_ice_init
        procedure :: update        => yelmo_ice_update
        procedure :: restart_write => yelmo_ice_restart_write
        procedure :: restart_read  => yelmo_ice_restart_read
    end type yelmo_ice_t

contains

    subroutine yelmo_ice_init(self, par_file, group, grid, time, cnst)
        class(yelmo_ice_t),     intent(inout) :: self
        character(len=*),       intent(in)    :: par_file
        character(len=*),       intent(in)    :: group
        type(grid_class),       intent(in)    :: grid
        real(dp),               intent(in)    :: time
        type(phys_const_class), intent(in)    :: cnst
        write(*,'(a)') "kryos_ice_yelmo: init (placeholder) "//trim(self%name)
    end subroutine yelmo_ice_init

    subroutine yelmo_ice_update(self, t_target)
        class(yelmo_ice_t), intent(inout) :: self
        real(dp),           intent(in)    :: t_target
        write(*,'(a,f10.2)') "kryos_ice_yelmo: update (placeholder) t_target=", t_target
    end subroutine yelmo_ice_update

    subroutine yelmo_ice_restart_write(self, dir, time)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: time
        write(*,'(a)') "kryos_ice_yelmo: restart_write (placeholder) "//trim(dir)
    end subroutine yelmo_ice_restart_write

    subroutine yelmo_ice_restart_read(self, dir, time)
        class(yelmo_ice_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: time
        write(*,'(a)') "kryos_ice_yelmo: restart_read (placeholder) "//trim(dir)
    end subroutine yelmo_ice_restart_read

end module kryos_ice_yelmo
