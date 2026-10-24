submodule (kryos_climate) kryos_climate_factory
    ! The climate factory: the only place that names the climate backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use kryos_climate_snapshots, only: snapshots_climate_t, snapshots_climate_constraints

    implicit none

contains

    module subroutine climate_create(climate, ctx, census)
        type(climate_slot_t), intent(out)   :: climate
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("snapshots")
            allocate(snapshots_climate_t :: climate%c)
            con = snapshots_climate_constraints()
        case default
            call fesm_error("kryos_climate::climate_create", "unknown climate backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: snapshots")
        end select
        call component_identity(climate%c, ctx, "climate", con)
        call backend_admit(climate%c, census)
    end subroutine climate_create

end submodule kryos_climate_factory
