submodule (kryos_ice) kryos_ice_factory
    ! The ice factory: the only place that names the ice backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use fesm_log,        only: fesm_error
    use kryos_ice_yelmo, only: yelmo_ice_t, yelmo_ice_constraints

    implicit none

contains

    module subroutine ice_create(ice, ctx, census)
        type(ice_slot_t), intent(out)   :: ice
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("yelmo")
            allocate(yelmo_ice_t :: ice%c)
            con = yelmo_ice_constraints()
        ! phase 2: case ("data")
        case default
            call fesm_error("kryos_ice::ice_create", "unknown ice backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: yelmo")
        end select
        call component_identity(ice%c, ctx, "ice", con)
        call backend_admit(ice%c, census)
    end subroutine ice_create

end submodule kryos_ice_factory
