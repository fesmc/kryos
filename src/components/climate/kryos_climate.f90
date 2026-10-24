module kryos_climate
    ! The climate component (contracts-design sections 2.3, 3.3, 7): abstract
    ! climate_t with the handles of the canonical vocabulary (atmosphere on
    ! the month axis, ocean on the producer's depth axis, capabilities), its
    ! slot, the wrappers and the factory submodule kryos_climate_factory.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis

    implicit none

    private

    type, abstract, extends(component_t), public :: climate_t
        integer :: hi_z_srf = 0, hi_H_ice = 0, hi_z_bed = 0, hi_f_grnd = 0, hi_z_sl = 0, hi_z_srf_ref = 0
        integer :: h_tas = 0, h_pr = 0, h_z_srf_clim = 0, h_lapse = 0, h_beta_p = 0
        integer :: h_to_ann = 0, h_so_ann = 0, h_dto_ann = 0                   ! has_ocn
        integer :: h_to_shlf = 0, h_so_shlf = 0, h_dto_shlf = 0, h_dso_shlf = 0 ! has_ocn_shelf
        integer :: h_smb = 0, h_T_srf = 0, h_Qd = 0                            ! has_smb, has_Qd
        ! Forcing scalars (section 7), set by the backend in init: the
        ! coefficients its terms name (imports coef_<name>) and the derived
        ! coefficients it produces (exports <name>), with their units.
        character(len=32), allocatable :: coef_names(:), coef_units(:)
        character(len=32), allocatable :: derived_names(:), derived_units(:)
        integer, allocatable :: hi_coef(:), h_derived(:)
    end type climate_t

    type, public :: climate_slot_t                      ! what domain_t (or globe_t) holds
        class(climate_t), allocatable :: c
    end type climate_slot_t

    interface                                       ! factory: submodule kryos_climate_factory
        module subroutine climate_create(climate, ctx, census)
            type(climate_slot_t), intent(out)   :: climate
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine climate_create
    end interface

    public :: climate_create, climate_adopt, climate_present
    public :: climate_init, climate_init_state, climate_update, climate_restart_write, climate_restart_read
    public :: climate_write_output, climate_end

contains

    subroutine climate_adopt(climate, proto, con, ctx, census)
        ! A host-provided extension of climate_t (design decision 1), admitted as a factory would.
        type(climate_slot_t),     intent(out)   :: climate
        class(climate_t),         intent(in)    :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        allocate(climate%c, mold=proto)
        call component_identity(climate%c, ctx, "climate", con)
        call backend_admit(climate%c, census)
    end subroutine climate_adopt

    logical function climate_present(climate)
        type(climate_slot_t), intent(in) :: climate
        climate_present = allocated(climate%c)
    end function climate_present

    subroutine climate_init(climate, ctx)
        type(climate_slot_t), intent(inout) :: climate
        type(init_ctx_t),     intent(in)    :: ctx
        call component_enter(climate%c)
        call climate%c%init(ctx)
        call climate_declare(climate%c, ctx)
        call component_leave()
    end subroutine climate_init

    subroutine climate_init_state(climate, st, pass)
        type(climate_slot_t),  intent(inout) :: climate
        type(step_time_t),     intent(in)    :: st
        integer,               intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(climate%c)
        call climate%c%init_state(s, pass)
        call component_close(climate%c, s)
        call component_leave()
    end subroutine climate_init_state

    subroutine climate_update(climate, st)
        type(climate_slot_t),  intent(inout) :: climate
        type(step_time_t),     intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(climate%c)
        call climate%c%update(st)
        call component_close(climate%c, st)
        call component_leave()
    end subroutine climate_update

    subroutine climate_restart_write(climate, dir, t)
        type(climate_slot_t), intent(inout) :: climate
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. climate%c%con%restart) return             ! backend keeps no restart file
        call component_enter(climate%c)
        call climate%c%restart_write(dir, t)
        call component_leave()
    end subroutine climate_restart_write

    subroutine climate_restart_read(climate, dir, t)
        type(climate_slot_t), intent(inout) :: climate
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. climate%c%con%restart) return
        call component_enter(climate%c)
        call climate%c%restart_read(dir, t)
        call component_leave()
    end subroutine climate_restart_read

    subroutine climate_write_output(climate, file, t, stream, init)
        type(climate_slot_t), intent(inout) :: climate
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init
        call component_enter(climate%c)
        call climate%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine climate_write_output

    subroutine climate_end(climate)
        type(climate_slot_t), intent(inout) :: climate
        call component_enter(climate%c)
        call climate%c%end()
        call component_leave()
    end subroutine climate_end

    subroutine climate_declare(c, ctx)
        ! Contracts-design section 2.3: geometry from the hub, coefficients
        ! from the forcing; atmosphere always, ocean and surface by capability.
        class(climate_t), intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: d2(2), k
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_add(c%imports, "z_srf", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_srf, grid=g, dims=d2)
        call fields_add(c%imports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice, grid=g, dims=d2)
        call fields_add(c%imports, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed, grid=g, dims=d2)
        call fields_add(c%imports, "f_grnd", "1", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_f_grnd, grid=g, dims=d2)
        call fields_add(c%imports, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_sl, grid=g, dims=d2)
        call fields_add(c%imports, "z_srf_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_srf_ref, grid=g, dims=d2)
        if (.not. allocated(c%coef_names)) allocate(c%coef_names(0), c%coef_units(0))
        allocate(c%hi_coef(size(c%coef_names)))
        do k = 1, size(c%coef_names)
            call fields_add(c%imports, "coef_"//trim(c%coef_names(k)), c%coef_units(k), FIELD_INSTANT, &
                            FIELD_INTENSIVE, c%hi_coef(k))
        end do
        call fields_lock(c%imports)

        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_tas, grid=g, &
                        dims=[d2, component_axis(c, ctx, "month")], zdim="month")
        call fields_add(c%exports, "pr", "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_pr, grid=g, &
                        dims=[d2, component_axis(c, ctx, "month")], zdim="month")
        call fields_add(c%exports, "z_srf_clim", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_srf_clim, grid=g, dims=d2)
        call fields_add(c%exports, "lapse", "K m-1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_lapse, grid=g, &
                        dims=[d2, component_axis(c, ctx, "month")], zdim="month")
        call fields_add(c%exports, "beta_p", "K-1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_beta_p, grid=g, dims=d2)
        if (c%caps%has_ocn) then
            call fields_add(c%exports, "to_ann", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_to_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
            call fields_add(c%exports, "so_ann", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%h_so_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
            call fields_add(c%exports, "dto_ann", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dto_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
        end if
        if (c%caps%has_ocn_shelf) then
            call fields_add(c%exports, "to_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_to_shlf, grid=g, dims=d2)
            call fields_add(c%exports, "so_shlf", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%h_so_shlf, grid=g, dims=d2)
            call fields_add(c%exports, "dto_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dto_shlf, grid=g, dims=d2)
            call fields_add(c%exports, "dso_shlf", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dso_shlf, grid=g, dims=d2)
        end if
        if (c%caps%has_smb) then
            call fields_add(c%exports, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_smb, grid=g, dims=d2)
            call fields_add(c%exports, "T_srf", "K", FIELD_MEAN, FIELD_INTENSIVE, c%h_T_srf, grid=g, dims=d2)
        end if
        if (c%caps%has_Qd) then
            call fields_add(c%exports, "Qd", "m3 s-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_Qd, grid=g, dims=d2)
        end if
        if (.not. allocated(c%derived_names)) allocate(c%derived_names(0), c%derived_units(0))
        allocate(c%h_derived(size(c%derived_names)))
        do k = 1, size(c%derived_names)
            call fields_add(c%exports, c%derived_names(k), c%derived_units(k), FIELD_INSTANT, &
                            FIELD_INTENSIVE, c%h_derived(k))
        end do
        call fields_lock(c%exports)
    end subroutine climate_declare

end module kryos_climate
