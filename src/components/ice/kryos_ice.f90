module kryos_ice
    ! The ice component (contracts-design sections 2.1, 3.2): abstract ice_t
    ! with the handles of the canonical vocabulary, the slot the domain holds,
    ! and the free-procedure wrappers, the only call sites of the bindings.
    ! Backends extend ice_t (src/backends/ice/); the factory is the submodule
    ! kryos_ice_factory, which names them; a host adds its own with ice_adopt.

    use precision,       only: dp
    use fesm_log,        only: fesm_error, to_str
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, PATH_LEN, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional

    implicit none

    private

    type, public :: ice_opt_t                       ! optional fields the backend supports (set in init)
        logical :: T_shlf = .false., tf_shlf = .false., Q_geo = .false., H_sed = .false., Qd = .false.
        logical :: uxy_s = .false., dz_bed_dt_corr = .false., z_bed_ref = .false., H_ice_ref = .false.
        logical :: footprint = .false., H_ice_pd = .false.     ! static: where the ice model holds ice; its present day
        logical :: H_grnd = .false., H_grnd_pd = .false.       ! thickness above flotation, and its present day (static)
    end type ice_opt_t

    type, abstract, extends(component_t), public :: ice_t
        type(ice_opt_t) :: opt
        integer :: hi_z_bed = 0, hi_z_sl = 0, hi_smb = 0, hi_T_srf = 0, hi_bmb_shlf = 0     ! imports
        integer :: hi_T_shlf = 0, hi_tf_shlf = 0, hi_Q_geo = 0, hi_H_sed = 0, hi_Qd = 0
        integer :: h_H_ice = 0, h_z_srf = 0, h_z_bed = 0, h_z_sl = 0, h_f_ice = 0, h_f_grnd = 0   ! exports
        integer :: h_dHidt = 0, h_smb = 0, h_bmb = 0, h_cmb = 0, h_mb_net = 0, h_M_ice = 0
        integer :: h_ux_s = 0, h_uy_s = 0, h_dz_bed_dt_corr = 0, h_z_bed_ref = 0, h_H_ice_ref = 0
        integer :: h_footprint = 0, h_H_ice_pd = 0, h_H_grnd = 0, h_H_grnd_pd = 0
    contains
        procedure :: flotation => ice_flotation_default   ! (H, z_bed, z_sl, f_grnd): capability has_flotation
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
    public :: ice_restart_files
    public :: ice_write_output, ice_end, ice_flotation

contains

    subroutine ice_flotation_default(self, H, z_bed, z_sl, f_grnd)
        ! A backend without the binding (caps%has_flotation false): ice_flotation stops first.
        class(ice_t), intent(in)  :: self
        real(dp),     intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),     intent(out) :: f_grnd(:,:)
        call fesm_error("kryos_ice::flotation", "instance `"//trim(self%name)//"` ("//trim(self%backend)// &
                        ") has no flotation binding")
        associate (unused => size(H) + size(z_bed) + size(z_sl)); end associate
        f_grnd = 0.0_dp
    end subroutine ice_flotation_default

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

    function ice_restart_files(ice) result(files)
        ! The base names the backend's restart_write creates in the bundle folder.
        type(ice_slot_t), intent(in) :: ice
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => ice%c)                          ! (gfortran: ICE on a call through ice%c itself)
            files = c%restart_files()
        end associate
    end function ice_restart_files

    subroutine ice_write_output(ice, file, t, stream, init)
        type(ice_slot_t), intent(inout) :: ice
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        call component_enter(ice%c)
        call ice%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine ice_write_output

    subroutine ice_flotation(ice, H, z_bed, z_sl, f_grnd)
        ! The grounded fraction the ice model's own rule gives for ice thickness H
        ! over bed z_bed and sea level z_sl [m, the instance's grid, dp]
        ! (contracts-design 2.1): what the core uses for ice outside the
        ! footprint and in the reference, never the model's routines itself.
        type(ice_slot_t), intent(in)  :: ice
        real(dp),         intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),         intent(out) :: f_grnd(:,:)

        character(len=*), parameter :: proc = "kryos_ice::ice_flotation"

        if (.not. ice%c%caps%has_flotation) then
            call fesm_error(proc, "instance `"//trim(ice%c%name)//"` ("//trim(ice%c%backend)// &
                            ") does not declare has_flotation", "the backend sets caps%has_flotation in init")
        end if
        if (any([size(H,1), size(z_bed,1), size(z_sl,1), size(f_grnd,1)] /= ice%c%nx) .or. &
            any([size(H,2), size(z_bed,2), size(z_sl,2), size(f_grnd,2)] /= ice%c%ny)) then
            call fesm_error(proc, "arrays are not on the instance grid", "instance `"//trim(ice%c%name)//"`, grid "// &
                            to_str(ice%c%nx)//" x "//to_str(ice%c%ny))
        end if
        call component_enter(ice%c)
        call ice%c%flotation(H, z_bed, z_sl, f_grnd)
        call component_leave()
    end subroutine ice_flotation

    subroutine ice_end(ice)
        type(ice_slot_t), intent(inout) :: ice
        call component_enter(ice%c)
        call ice%c%end()
        call component_leave()
    end subroutine ice_end

    subroutine ice_declare(c)
        ! The vocabulary of contracts-design section 2.1, once for all backends.
        ! H_ice_ref (optional): the reference load geo is initialised with
        ! (design decision 4, refined in the phase 1 review). footprint and
        ! H_ice_pd (optional, static, phase 4 and phase 2 WP3): published once,
        ! at PASS_INIT. H_grnd (optional, with the geometry) and H_grnd_pd (static) are
        ! what the thermal-forcing calibration (experiment opt_tf) reads, with H_ice,
        ! dHidt and H_ice_pd.
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
        if (c%opt%footprint) then
            call fields_add(c%exports, "footprint", "1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_footprint, grid=g, dims=d2)
        end if
        if (c%opt%H_ice_pd) then
            call fields_add(c%exports, "H_ice_pd", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_H_ice_pd, grid=g, dims=d2)
        end if
        if (c%opt%H_grnd) then                      ! per step, with the geometry
            call fields_add(c%exports, "H_grnd", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_H_grnd, grid=g, dims=d2)
        end if
        if (c%opt%H_grnd_pd) then                   ! static, with H_ice_pd
            call fields_add(c%exports, "H_grnd_pd", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_H_grnd_pd, grid=g, dims=d2)
        end if
        call fields_lock(c%exports)
    end subroutine ice_declare

end module kryos_ice
