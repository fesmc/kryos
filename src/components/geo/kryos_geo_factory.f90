submodule (kryos_geo) kryos_geo_factory
    ! The geo factory: the only place that names the geo backends. Allocates
    ! the dynamic type from the backend name, sets its constraints and admits
    ! it (decision 12), so a stub stops here, before init.

    use kryos_geo_fastisostasy, only: fastiso_geo_t, fastiso_bsl_t, fastiso_geo_constraints, &
                                      fastiso_bsl_constraints
    use kryos_geo_vilma2,       only: vilma_earth_t, vilma_geo_t, vilma2_earth_constraints, &
                                      vilma2_geo_constraints

    implicit none

contains

    module subroutine geo_create(geo, ctx, census)
        type(geo_slot_t), intent(out)   :: geo
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census
        type(constraints_t) :: con

        select case (trim(ctx%backend))
        case ("fastisostasy")
            if (trim(ctx%scope) == "global") then     ! the global part: bsl
                allocate(fastiso_bsl_t :: geo%c)
                con = fastiso_bsl_constraints()
            else
                allocate(fastiso_geo_t :: geo%c)
                con = fastiso_geo_constraints()
            end if
        case ("vilma2")
            if (trim(ctx%scope) == "global") then     ! the Earth (VILMA2 on its Gauss grid)
                allocate(vilma_earth_t :: geo%c)
                con = vilma2_earth_constraints()
            else                                      ! the port of a domain
                allocate(vilma_geo_t :: geo%c)
                con = vilma2_geo_constraints()
            end if
        case default
            call fesm_error("kryos_geo::geo_create", "unknown geo backend '"//trim(ctx%backend)//"'", &
                            "instance "//trim(ctx%instance)//"; known: fastisostasy, vilma2")
        end select
        call component_identity(geo%c, ctx, "geo", con)
        call backend_admit(geo%c, census)
    end subroutine geo_create

end submodule kryos_geo_factory
