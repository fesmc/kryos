module kryos_shelf
    ! The shelf component, sub-shelf melt (contracts-design sections 2.5, 3.3,
    ! decision 6): abstract shelf_t with the hub geometry and two alternative
    ! ocean inputs (profiles or shelf-base values), its slot, the wrappers and
    ! the factory submodule kryos_shelf_factory.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, PATH_LEN, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis, axis_find

    implicit none

    private

    type, public :: shelf_opt_t                     ! ocean inputs the backend accepts (set in init)
        logical :: ocn_prof = .false.               ! set A: profiles on the depth axis (has_ocn)
        logical :: ocn_shelf = .false.              ! set B: shelf-base values (has_ocn_shelf)
    end type shelf_opt_t

    type, abstract, extends(component_t), public :: shelf_t
        type(shelf_opt_t) :: opt
        integer :: hi_H_ice = 0, hi_z_bed = 0, hi_f_grnd = 0, hi_z_sl = 0, hi_z_srf = 0
        integer :: hi_to_ann = 0, hi_so_ann = 0, hi_dto_ann = 0                    ! set A
        integer :: hi_to_shlf = 0, hi_so_shlf = 0, hi_dto_shlf = 0, hi_dso_shlf = 0 ! set B
        integer :: h_bmb_shlf = 0, h_T_shlf = 0, h_tf_shlf = 0
    end type shelf_t

    type, public :: shelf_slot_t                      ! what domain_t (or globe_t) holds
        class(shelf_t), allocatable :: c
    end type shelf_slot_t

    interface                                       ! factory: submodule kryos_shelf_factory
        module subroutine shelf_create(shelf, ctx, census)
            type(shelf_slot_t), intent(out)   :: shelf
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine shelf_create
    end interface

    public :: shelf_create, shelf_adopt, shelf_present
    public :: shelf_init, shelf_init_state, shelf_update, shelf_restart_write, shelf_restart_read
    public :: shelf_restart_files
    public :: shelf_write_output, shelf_end

contains

    subroutine shelf_adopt(shelf, proto, con, ctx, census)
        ! A host-provided extension of shelf_t (design decision 1), admitted as a factory would.
        type(shelf_slot_t),   intent(out)   :: shelf
        class(shelf_t), allocatable, intent(in) :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(shelf%c, mold=proto)
        call component_identity(shelf%c, ctx, "shelf", con)
        call backend_admit(shelf%c, census)
    end subroutine shelf_adopt

    logical function shelf_present(shelf)
        type(shelf_slot_t), intent(in) :: shelf
        shelf_present = allocated(shelf%c)
    end function shelf_present

    subroutine shelf_init(shelf, ctx)
        type(shelf_slot_t), intent(inout) :: shelf
        type(init_ctx_t),   intent(in)    :: ctx
        call component_enter(shelf%c)
        call shelf%c%init(ctx)
        call shelf_declare(shelf%c, ctx)
        call component_leave()
    end subroutine shelf_init

    subroutine shelf_init_state(shelf, st, pass)
        type(shelf_slot_t),  intent(inout) :: shelf
        type(step_time_t),   intent(in)    :: st
        integer,             intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(shelf%c)
        call shelf%c%init_state(s, pass)
        call component_close(shelf%c, s)
        call component_leave()
    end subroutine shelf_init_state

    subroutine shelf_update(shelf, st)
        type(shelf_slot_t),  intent(inout) :: shelf
        type(step_time_t),   intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(shelf%c)
        call shelf%c%update(st)
        call component_close(shelf%c, st)
        call component_leave()
    end subroutine shelf_update

    subroutine shelf_restart_write(shelf, dir, t)
        type(shelf_slot_t), intent(inout) :: shelf
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        if (.not. shelf%c%con%restart) return             ! backend keeps no restart file
        call component_enter(shelf%c)
        call shelf%c%restart_write(dir, t)
        call component_leave()
    end subroutine shelf_restart_write

    subroutine shelf_restart_read(shelf, dir, t)
        type(shelf_slot_t), intent(inout) :: shelf
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        if (.not. shelf%c%con%restart) return
        call component_enter(shelf%c)
        call shelf%c%restart_read(dir, t)
        call component_leave()
    end subroutine shelf_restart_read

    function shelf_restart_files(shelf) result(files)
        ! The base names the backend's restart_write creates in the bundle folder.
        type(shelf_slot_t), intent(in) :: shelf
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => shelf%c)                          ! (gfortran: ICE on a call through shelf%c itself)
            files = c%restart_files()
        end associate
    end function shelf_restart_files

    subroutine shelf_write_output(shelf, file, t, stream, init)
        type(shelf_slot_t), intent(inout) :: shelf
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init
        call component_enter(shelf%c)
        call shelf%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine shelf_write_output

    subroutine shelf_end(shelf)
        type(shelf_slot_t), intent(inout) :: shelf
        call component_enter(shelf%c)
        call shelf%c%end()
        call component_leave()
    end subroutine shelf_end

    subroutine shelf_declare(c, ctx)
        ! Contracts-design section 2.5. Both ocean sets are optional imports;
        ! wire_resolve requires exactly one complete and active. Set A needs the
        ! producer's depth axis (in ctx%axes when the climate has_ocn).
        class(shelf_t),   intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: d2(2)
        character(len=:), allocatable :: g
        logical :: has_depth

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        has_depth = .false.
        if (allocated(ctx%axes)) has_depth = axis_find(ctx%axes, "depth") > 0

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_add(c%imports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice, grid=g, dims=d2)
        call fields_add(c%imports, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed, grid=g, dims=d2)
        call fields_add(c%imports, "f_grnd", "1", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_f_grnd, grid=g, dims=d2)
        call fields_add(c%imports, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_sl, grid=g, dims=d2)
        call fields_add(c%imports, "z_srf", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_srf, grid=g, dims=d2)
        if (c%opt%ocn_prof .and. has_depth) then
            call fields_add(c%imports, "to_ann", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_to_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
            call component_optional(c, "to_ann")
            call fields_add(c%imports, "so_ann", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_so_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
            call component_optional(c, "so_ann")
            call fields_add(c%imports, "dto_ann", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dto_ann, grid=g, &
                            dims=[d2, component_axis(c, ctx, "depth")], zdim="depth")
            call component_optional(c, "dto_ann")
        end if
        if (c%opt%ocn_shelf) then
            call fields_add(c%imports, "to_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_to_shlf, grid=g, dims=d2)
            call component_optional(c, "to_shlf")
            call fields_add(c%imports, "so_shlf", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_so_shlf, grid=g, dims=d2)
            call component_optional(c, "so_shlf")
            call fields_add(c%imports, "dto_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dto_shlf, grid=g, dims=d2)
            call component_optional(c, "dto_shlf")
            call fields_add(c%imports, "dso_shlf", "psu", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dso_shlf, grid=g, dims=d2)
            call component_optional(c, "dso_shlf")
        end if
        call fields_lock(c%imports)

        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "bmb_shlf", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_bmb_shlf, grid=g, dims=d2)
        call fields_add(c%exports, "T_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_T_shlf, grid=g, dims=d2)
        call fields_add(c%exports, "tf_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h_tf_shlf, grid=g, dims=d2)
        call fields_lock(c%exports)
    end subroutine shelf_declare

end module kryos_shelf
