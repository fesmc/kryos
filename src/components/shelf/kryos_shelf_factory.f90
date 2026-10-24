submodule (kryos_shelf) kryos_shelf_factory
    ! The shelf factory: the only place that names the shelf backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use kryos_shelf_marine_shelf, only: marshelf_shelf_t, marshelf_shelf_constraints

    implicit none

contains

    module subroutine shelf_create(shelf, ctx, census)
        type(shelf_slot_t), intent(out)   :: shelf
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("marine_shelf")
            allocate(marshelf_shelf_t :: shelf%c)
            con = marshelf_shelf_constraints()
        case default
            call fesm_error("kryos_shelf::shelf_create", "unknown shelf backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: marine_shelf")
        end select
        call component_identity(shelf%c, ctx, "shelf", con)
        call backend_admit(shelf%c, census)
    end subroutine shelf_create

end submodule kryos_shelf_factory
