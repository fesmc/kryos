module kryos_ice
    ! The ice component: abstract ice_t with the canonical field handles
    ! (design section 5), and its factory. Backends extend ice_t, one file each
    ! (src/backends/ice/); the factory lives in the submodule kryos_ice_factory
    ! so that it can `use` the backend modules, which `use` this one.
    !
    ! placeholder: phase 1 fills this in.

    use kryos_component, only: component_t

    implicit none

    private
    public :: ice_t, ice_create

    type, abstract, extends(component_t) :: ice_t
        ! Handles of the canonical fields, resolved at init (section 5).
        integer :: h_H_ice = 0, h_z_srf = 0, h_f_ice = 0, h_f_grnd = 0
        integer :: h_z_bed = 0, h_z_sl  = 0, h_smb   = 0, h_T_srf  = 0, h_bmb_shlf = 0
    end type ice_t

    interface
        ! Allocate the dynamic type of `ice` from the backend name (decision 1).
        module subroutine ice_create(ice, backend, name)
            class(ice_t), allocatable, intent(out) :: ice
            character(len=*),          intent(in)  :: backend
            character(len=*),          intent(in)  :: name
        end subroutine ice_create
    end interface

end module kryos_ice
