submodule (kryos_ice) kryos_ice_factory
    ! The ice factory: the only place that names the ice backends.
    ! placeholder: phase 1 fills this in (backend constraints, decision 12).

    use kryos_ice_yelmo, only: yelmo_ice_t

    implicit none

contains

    module subroutine ice_create(ice, backend, name)
        class(ice_t), allocatable, intent(out) :: ice
        character(len=*),          intent(in)  :: backend
        character(len=*),          intent(in)  :: name

        select case (trim(backend))
        case ("yelmo")
            allocate(yelmo_ice_t :: ice)
        case default
            write(*,'(a)') "kryos_ice: unknown ice backend '"//trim(backend)// &
                           "' for instance '"//trim(name)//"'"
            error stop 1
        end select

        ice%name    = name
        ice%kind    = "ice"
        ice%backend = backend
        ice%scope   = "regional"
    end subroutine ice_create

end submodule kryos_ice_factory
