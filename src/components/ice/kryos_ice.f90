module kryos_ice
    ! The ice component (contracts-design sections 2.1, 3.2): abstract ice_t
    ! with the handles of the canonical vocabulary, the slot the domain holds,
    ! and the free-procedure wrappers, the only call sites of the bindings.
    ! Backends extend ice_t (src/backends/ice/); the factory is the submodule
    ! kryos_ice_factory, which names them; a host adds its own with ice_adopt.

    use precision,       only: dp
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional

    implicit none

    private

    type, public :: ice_opt_t                       ! optional fields the backend supports (set in init)
        logical :: T_shlf = .false., tf_shlf = .false., Q_geo = .false., H_sed = .false., Qd = .false.
        logical :: uxy_s = .false., dz_bed_dt_corr = .false., z_bed_ref = .false., H_ice_ref = .false.
    end type ice_opt_t

    type, abstract, extends(component_t), public :: ice_t
        type(ice_opt_t) :: opt
        integer :: hi_z_bed = 0, hi_z_sl = 0, hi_smb = 0, hi_T_srf = 0, hi_bmb_shlf = 0     ! imports
        integer :: hi_T_shlf = 0, hi_tf_shlf = 0, hi_Q_geo = 0, hi_H_sed = 0, hi_Qd = 0
        integer :: h_H_ice = 0, h_z_srf = 0, h_z_bed = 0, h_z_sl = 0, h_f_ice = 0, h_f_grnd = 0   ! exports
        integer :: h_dHidt = 0, h_smb = 0, h_bmb = 0, h_cmb = 0, h_mb_net = 0, h_M_ice = 0
        integer :: h_ux_s = 0, h_uy_s = 0, h_dz_bed_dt_corr = 0, h_z_bed_ref = 0, h_H_ice_ref = 0
    end type ice_t

    type, public :: ice_slot_t                      ! what domain_t holds
        class(ice_t), allocatable :: c
    end type ice_slot_t

    interface                                       ! factory: submodule kryos_ice_factory
        module subroutine ice_create(ice, ctx, census)
            type(ice_slot_t), intent(out)   :: ice
            type(init_ctx_t), intent(in)    :: ctx  ! instance, backend, grid, domain, scope
            type(census_t),   intent(inout) :: census
        end subroutine ice_create
    end interface

    public :: ice_create, ice_adopt, ice_present
    public :: ice_init, ice_init_state, ice_update, ice_restart_write, ice_restart_read
    public :: ice_write_output, ice_end

contains

    subroutine ice_adopt(ice, proto, con, ctx, census)
        ! A host-provided extension of ice_t (design decision 1), admitted as a factory would.
        type(ice_slot_t),    intent(out)   :: ice
        class(ice_t), allocatable, intent(in) :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(ice%c, mold=proto)
        call component_identity(ice%c, ctx, "ice", con)
        call backend_admit(ice%c, census)
    end subroutine ice_adopt

    logical function ice_present(ice)
        type(ice_slot_t), intent(in) :: ice
        ice_present = allocated(ice%c)
    end function ice_present

    subroutine ice_init(ice, ctx)
        type(ice_slot_t), intent(inout) :: ice
        type(init_ctx_t), intent(in)    :: ctx
        call component_enter(ice%c)
        call ice%c%init(ctx)
        call ice_declare(ice%c)
        call component_leave()
    end subroutine ice_init

    subroutine ice_init_state(ice, st, pass)
        type(ice_slot_t),  intent(inout) :: ice
        type(step_time_t), intent(in)    :: st          ! t0 = t1 = time_init (core section 3.2)
        integer,           intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(ice%c)
        call ice%c%init_state(s, pass)
        call component_close(ice%c, s)
        call component_leave()
    end subroutine ice_init_state

    subroutine ice_update(ice, st)
        type(ice_slot_t),  intent(inout) :: ice
        type(step_time_t), intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(ice%c)
        call ice%c%update(st)
        call component_close(ice%c, st)
        call component_leave()
    end subroutine ice_update

    subroutine ice_restart_write(ice, dir, t)
        type(ice_slot_t), intent(inout) :: ice
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        if (.not. ice%c%con%restart) return             ! backend keeps no restart file
        call component_enter(ice%c)
        call ice%c%restart_write(dir, t)
        call component_leave()
    end subroutine ice_restart_write

    subroutine ice_restart_read(ice, dir, t)
        type(ice_slot_t), intent(inout) :: ice
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        if (.not. ice%c%con%restart) return
        call component_enter(ice%c)
        call ice%c%restart_read(dir, t)
        call component_leave()
    end subroutine ice_restart_read

    subroutine ice_write_output(ice, file, t, stream, init)
        type(ice_slot_t), intent(inout) :: ice
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        call component_enter(ice%c)
        call ice%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine ice_write_output

    subroutine ice_end(ice)
        type(ice_slot_t), intent(inout) :: ice
        call component_enter(ice%c)
        call ice%c%end()
        call component_leave()
    end subroutine ice_end

    subroutine ice_declare(c)
        ! The vocabulary of contracts-design section 2.1, once for all backends.
        ! H_ice_ref (optional): the reference load geo is initialised with
        ! (design decision 4, refined in the phase 1 review).
        class(ice_t), intent(inout) :: c
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_add(c%imports, "z_bed",    "m",           FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed,    grid=g, dims=d2)
        call fields_add(c%imports, "z_sl",     "m",           FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_sl,     grid=g, dims=d2)
        call fields_add(c%imports, "smb",      "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_smb,      grid=g, dims=d2)
        call fields_add(c%imports, "T_srf",    "K",           FIELD_INSTANT, FIELD_INTENSIVE, c%hi_T_srf,    grid=g, dims=d2)
        call fields_add(c%imports, "bmb_shlf", "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_bmb_shlf, grid=g, dims=d2)
        if (c%opt%T_shlf) then
            call fields_add(c%imports, "T_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_T_shlf, grid=g, dims=d2)
            call component_optional(c, "T_shlf")
        end if
        if (c%opt%tf_shlf) then
            call fields_add(c%imports, "tf_shlf", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_tf_shlf, grid=g, dims=d2)
            call component_optional(c, "tf_shlf")
        end if
        if (c%opt%Q_geo) then
            call fields_add(c%imports, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_Q_geo, grid=g, dims=d2)
            call component_optional(c, "Q_geo")
        end if
        if (c%opt%H_sed) then
            call fields_add(c%imports, "H_sed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_H_sed, grid=g, dims=d2)
            call component_optional(c, "H_sed")
        end if
        if (c%opt%Qd) then
            call fields_add(c%imports, "Qd", "m3 s-1", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_Qd, grid=g, dims=d2)
            call component_optional(c, "Qd")
        end if
        call fields_lock(c%imports)

        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "H_ice",  "m",           FIELD_INSTANT, FIELD_EXTENSIVE, c%h_H_ice,  grid=g, dims=d2)
        call fields_add(c%exports, "z_srf",  "m",           FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_srf,  grid=g, dims=d2)
        call fields_add(c%exports, "z_bed",  "m",           FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_bed,  grid=g, dims=d2)
        call fields_add(c%exports, "z_sl",   "m",           FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_sl,   grid=g, dims=d2)
        call fields_add(c%exports, "f_ice",  "1",           FIELD_INSTANT, FIELD_INTENSIVE, c%h_f_ice,  grid=g, dims=d2)
        call fields_add(c%exports, "f_grnd", "1",           FIELD_INSTANT, FIELD_INTENSIVE, c%h_f_grnd, grid=g, dims=d2)
        call fields_add(c%exports, "dHidt",  "m yr-1",      FIELD_MEAN,    FIELD_EXTENSIVE, c%h_dHidt,  grid=g, dims=d2)
        call fields_add(c%exports, "smb",    "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, c%h_smb,    grid=g, dims=d2)
        call fields_add(c%exports, "bmb",    "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, c%h_bmb,    grid=g, dims=d2)
        call fields_add(c%exports, "cmb",    "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, c%h_cmb,    grid=g, dims=d2)
        call fields_add(c%exports, "mb_net", "kg m-2 yr-1", FIELD_MEAN,    FIELD_EXTENSIVE, c%h_mb_net, grid=g, dims=d2)
        call fields_add(c%exports, "M_ice",  "Gt",          FIELD_INSTANT, FIELD_EXTENSIVE, c%h_M_ice)      ! 0-D
        if (c%opt%uxy_s) then
            call fields_add(c%exports, "ux_s", "m yr-1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_ux_s, grid=g, dims=d2)
            call fields_add(c%exports, "uy_s", "m yr-1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_uy_s, grid=g, dims=d2)
        end if
        if (c%opt%dz_bed_dt_corr) then
            call fields_add(c%exports, "dz_bed_dt_corr", "m yr-1", FIELD_INSTANT, FIELD_INTENSIVE, &
                            c%h_dz_bed_dt_corr, grid=g, dims=d2)
        end if
        if (c%opt%z_bed_ref) then
            call fields_add(c%exports, "z_bed_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_bed_ref, grid=g, dims=d2)
        end if
        if (c%opt%H_ice_ref) then
            call fields_add(c%exports, "H_ice_ref", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_H_ice_ref, grid=g, dims=d2)
        end if
        call fields_lock(c%exports)
    end subroutine ice_declare

end module kryos_ice
