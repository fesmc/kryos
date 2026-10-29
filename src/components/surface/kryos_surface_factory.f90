submodule (kryos_surface) kryos_surface_factory
    ! The surface factory: the only place that names the surface backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use kryos_surface_chion, only: chion_surface_t, chion_surface_constraints
    use kryos_surface_smb_simple, only: smb_simple_surface_t, smb_simple_surface_constraints

    implicit none

contains

    module subroutine surface_create(surface, ctx, census)
        type(surface_slot_t), intent(out)   :: surface
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("chion")
            allocate(chion_surface_t :: surface%c)
            con = chion_surface_constraints()
        case ("smb_simple")
            allocate(smb_simple_surface_t :: surface%c)
            con = smb_simple_surface_constraints()
        case default
            call fesm_error("kryos_surface::surface_create", "unknown surface backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: chion, smb_simple")
        end select
        call component_identity(surface%c, ctx, "surface", con)
        call backend_admit(surface%c, census)
    end subroutine surface_create

end submodule kryos_surface_factory
