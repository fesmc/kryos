submodule (kryos_sediments) kryos_sediments_factory
    ! The sediments factory: the only place that names the sediments backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    implicit none

contains

    module subroutine sediments_create(sediments, ctx, census)
        type(sediments_slot_t), intent(out)   :: sediments
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("data")
            allocate(data_sediments_t :: sediments%c)
            con = data_sediments_constraints()
        case default
            call fesm_error("kryos_sediments::sediments_create", "unknown sediments backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: data")
        end select
        call component_identity(sediments%c, ctx, "sediments", con)
        call backend_admit(sediments%c, census)
    end subroutine sediments_create

end submodule kryos_sediments_factory
