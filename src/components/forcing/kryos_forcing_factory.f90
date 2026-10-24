submodule (kryos_forcing) kryos_forcing_factory
    ! The forcing factory: the only place that names the forcing backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use kryos_forcing_tsgen, only: tsgen_forcing_t, tsgen_forcing_constraints

    implicit none

contains

    module subroutine forcing_create(forcing, ctx, census)
        type(forcing_slot_t), intent(out)   :: forcing
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("tsgen")
            allocate(tsgen_forcing_t :: forcing%c)
            con = tsgen_forcing_constraints()
        case default
            call fesm_error("kryos_forcing::forcing_create", "unknown forcing backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: tsgen")
        end select
        call component_identity(forcing%c, ctx, "forcing", con)
        call backend_admit(forcing%c, census)
    end subroutine forcing_create

end submodule kryos_forcing_factory
