module kryos_surface
    ! The surface component (contracts-design sections 2.4, 3.3): abstract
    ! surface_t (SMB and surface temperature from the atmosphere on the hub
    ! geometry), its slot, the wrappers and the factory kryos_surface_factory.
    ! init_state is used by chion (spin-up in PASS_BOUNDARY, a cycle in PASS_FINAL).

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, PATH_LEN, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis

    implicit none

    private

    type, public :: surface_opt_t                   ! the optional vocabulary the backend uses (set in init)
        logical :: melt = .false., runoff = .false., refrz = .false., albedo = .false.   ! exports
        ! Imports. atm: the monthly atmosphere (tas, pr); geom: the land and ice geometry (H_ice, z_bed, z_sl);
        ! both on for a column model (chion), off for a scheme of the elevation alone (smb_simple). tsl_ann
        ! (climate) and H_ice_ref (ice) are the extra imports of smb_simple, required when set.
        logical :: atm = .true., geom = .true., tsl_ann = .false., H_ice_ref = .false.
    end type surface_opt_t

    type, abstract, extends(component_t), public :: surface_t
        type(surface_opt_t) :: opt
        integer :: hi_z_srf = 0, hi_H_ice = 0, hi_z_bed = 0, hi_z_sl = 0, hi_tas = 0, hi_pr = 0
        integer :: hi_tsl_ann = 0, hi_H_ice_ref = 0
        integer :: h_smb = 0, h_T_srf = 0, h_melt = 0, h_runoff = 0, h_refrz = 0, h_albedo = 0
    end type surface_t

    type, public :: surface_slot_t                      ! what domain_t (or globe_t) holds
        class(surface_t), allocatable :: c
    end type surface_slot_t

    interface                                       ! factory: submodule kryos_surface_factory
        module subroutine surface_create(surface, ctx, census)
            type(surface_slot_t), intent(out)   :: surface
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine surface_create
    end interface

    public :: surface_create, surface_adopt, surface_present
    public :: surface_init, surface_init_state, surface_update, surface_restart_write, surface_restart_read
    public :: surface_restart_files
    public :: surface_write_output, surface_end

contains

    subroutine surface_adopt(surface, proto, con, ctx, census)
        ! A host-provided extension of surface_t (design decision 1), admitted as a factory would.
        type(surface_slot_t),     intent(out)   :: surface
        class(surface_t), allocatable, intent(in) :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(surface%c, mold=proto)
        call component_identity(surface%c, ctx, "surface", con)
        call backend_admit(surface%c, census)
    end subroutine surface_adopt

    logical function surface_present(surface)
        type(surface_slot_t), intent(in) :: surface
        surface_present = allocated(surface%c)
    end function surface_present

    subroutine surface_init(surface, ctx)
        type(surface_slot_t), intent(inout) :: surface
        type(init_ctx_t),     intent(in)    :: ctx
        call component_enter(surface%c)
        call surface%c%init(ctx)
        call surface_declare(surface%c, ctx)
        call component_leave()
    end subroutine surface_init

    subroutine surface_init_state(surface, st, pass)
        type(surface_slot_t),  intent(inout) :: surface
        type(step_time_t),     intent(in)    :: st
        integer,               intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(surface%c)
        call surface%c%init_state(s, pass)
        call component_close(surface%c, s)
        call component_leave()
    end subroutine surface_init_state

    subroutine surface_update(surface, st)
        type(surface_slot_t),  intent(inout) :: surface
        type(step_time_t),     intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(surface%c)
        call surface%c%update(st)
        call component_close(surface%c, st)
        call component_leave()
    end subroutine surface_update

    subroutine surface_restart_write(surface, dir, t)
        type(surface_slot_t), intent(inout) :: surface
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. surface%c%con%restart) return             ! backend keeps no restart file
        call component_enter(surface%c)
        call surface%c%restart_write(dir, t)
        call component_leave()
    end subroutine surface_restart_write

    subroutine surface_restart_read(surface, dir, t)
        type(surface_slot_t), intent(inout) :: surface
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. surface%c%con%restart) return
        call component_enter(surface%c)
        call surface%c%restart_read(dir, t)
        call component_leave()
    end subroutine surface_restart_read

    function surface_restart_files(surface) result(files)
        ! The base names the backend's restart_write creates in the bundle folder.
        type(surface_slot_t), intent(in) :: surface
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => surface%c)                          ! (gfortran: ICE on a call through surface%c itself)
            files = c%restart_files()
        end associate
    end function surface_restart_files

    subroutine surface_write_output(surface, file, t, stream, init)
        type(surface_slot_t), intent(inout) :: surface
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init
        call component_enter(surface%c)
        call surface%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine surface_write_output

    subroutine surface_end(surface)
        type(surface_slot_t), intent(inout) :: surface
        call component_enter(surface%c)
        call surface%c%end()
        call component_leave()
    end subroutine surface_end

    subroutine surface_declare(c, ctx)
        ! Contracts-design section 2.4.
        class(surface_t), intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_add(c%imports, "z_srf", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_srf, grid=g, dims=d2)
        if (c%opt%geom) then
            call fields_add(c%imports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice, grid=g, dims=d2)
            call fields_add(c%imports, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed, grid=g, dims=d2)
            call fields_add(c%imports, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_sl, grid=g, dims=d2)
        end if
        if (c%opt%atm) then
            call fields_add(c%imports, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_tas, grid=g, &
                            dims=[d2, component_axis(c, ctx, "month")], zdim="month")
            call fields_add(c%imports, "pr", "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_pr, grid=g, &
                            dims=[d2, component_axis(c, ctx, "month")], zdim="month")
        end if
        if (c%opt%tsl_ann) then
            call fields_add(c%imports, "tsl_ann", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_tsl_ann, grid=g, dims=d2)
        end if
        if (c%opt%H_ice_ref) then
            call fields_add(c%imports, "H_ice_ref", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice_ref, grid=g, dims=d2)
        end if
        call fields_lock(c%imports)

        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_smb, grid=g, dims=d2)
        call fields_add(c%exports, "T_srf", "K", FIELD_MEAN, FIELD_INTENSIVE, c%h_T_srf, grid=g, dims=d2)
        if (c%opt%melt) then
            call fields_add(c%exports, "melt", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_melt, grid=g, dims=d2)
        end if
        if (c%opt%runoff) then
            call fields_add(c%exports, "runoff", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_runoff, grid=g, dims=d2)
        end if
        if (c%opt%refrz) then
            call fields_add(c%exports, "refrz", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h_refrz, grid=g, dims=d2)
        end if
        if (c%opt%albedo) then
            call fields_add(c%exports, "albedo", "1", FIELD_MEAN, FIELD_INTENSIVE, c%h_albedo, grid=g, dims=d2)
        end if
        call fields_lock(c%exports)
    end subroutine surface_declare

end module kryos_surface
