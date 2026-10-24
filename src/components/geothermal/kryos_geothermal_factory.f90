submodule (kryos_geothermal) kryos_geothermal_factory
    ! The geothermal factory: the only place that names the geothermal backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    implicit none

contains

    module subroutine geothermal_create(geothermal, ctx, census)
        type(geothermal_slot_t), intent(out)   :: geothermal
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("data")
            allocate(data_geothermal_t :: geothermal%c)
            con = data_geothermal_constraints()
        case default
            call fesm_error("kryos_geothermal::geothermal_create", "unknown geothermal backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: data")
        end select
        call component_identity(geothermal%c, ctx, "geothermal", con)
        call backend_admit(geothermal%c, census)
    end subroutine geothermal_create

end submodule kryos_geothermal_factory
