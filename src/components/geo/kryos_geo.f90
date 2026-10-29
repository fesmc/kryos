module kryos_geo
    ! The geo component, solid earth and sea level (contracts-design sections
    ! 2.2, 3.3, 6): one abstract geo_t for both levels, declared by scope
    ! (regional: bed and sea-surface anomalies; global: bsl and the reduction),
    ! its slot, wrappers incl. geo_reduce and geo_scatter, and the factory.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, fields_get, fields_put, &
                               FIELD_INSTANT, FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, PATH_LEN, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_core_filled

    implicit none

    private

    type, public :: geo_opt_t                       ! optional fields the backend supports (set in init)
        logical :: A_ocn = .false., dV_bsl = .false., dz_bed_dt = .false.
        logical :: ref = .false.                    ! regional imports z_bed_ref, H_ice_ref (decision 4)
        logical :: host = .false.                   ! global part resolved on the host's loads (vilma2): rank 2 on its
                                                    ! own grid (grid, nx, ny set in init), no dV_bsl (contracts 2.2)
        logical :: port = .false.                   ! regional port of such a global part: imports the scattered
                                                    ! dz_bed, dz_ss, dphi instead of H_ice
    end type geo_opt_t

    type, abstract, extends(component_t), public :: geo_t
        type(geo_opt_t) :: opt
        integer :: hi_H_ice = 0, hi_bsl = 0, hi_A_ocn = 0, hi_z_bed_ref = 0, hi_H_ice_ref = 0  ! regional imports
        integer :: hi_dV_bsl = 0                                                               ! global import
        integer :: hi_dz_bed = 0, hi_dz_ss = 0, hi_dphi = 0                                    ! port imports
        integer :: hi_H_ice_ext = 0, hi_host_w = 0, hi_host_Ig = 0, hi_host_C = 0, hi_host_F = 0   ! host global imports
        integer :: hi_host_Ig_ref = 0, hi_host_C_ref = 0
        integer :: h_C_ocn = 0, h_bsl_ive = 0, h_mass_resid = 0                                ! host global exports
        integer :: h_dz_bed = 0, h_dz_ss = 0, h_dV_bsl = 0, h_z_bed_ref = 0, h_dz_bed_dt = 0  ! regional exports
        integer :: h_bsl = 0, h_A_ocn = 0                                                      ! both / global
        integer :: h_dphi = 0                                                                   ! regional export
    contains
        procedure :: reduce => geo_reduce_default   ! (st): the global part's solve; default nothing
    end type geo_t

    type, public :: geo_slot_t
        class(geo_t), allocatable :: c
    end type geo_slot_t

    interface
        module subroutine geo_create(geo, ctx, census)
            type(geo_slot_t), intent(out)   :: geo
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine geo_create
    end interface

    public :: geo_create, geo_adopt, geo_present
    public :: geo_init, geo_init_state, geo_update, geo_restart_write, geo_restart_read
    public :: geo_restart_files
    public :: geo_write_output, geo_end, geo_reduce, geo_scatter

contains

    subroutine geo_reduce_default(self, st)
        class(geo_t),      intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (unused => st%t1 + self%nx); end associate
    end subroutine geo_reduce_default

    subroutine geo_adopt(geo, proto, con, ctx, census)
        type(geo_slot_t),    intent(out)   :: geo
        class(geo_t), allocatable, intent(in) :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(geo%c, mold=proto)
        call component_identity(geo%c, ctx, "geo", con)
        call backend_admit(geo%c, census)
    end subroutine geo_adopt

    logical function geo_present(geo)
        type(geo_slot_t), intent(in) :: geo
        geo_present = allocated(geo%c)
    end function geo_present

    subroutine geo_init(geo, ctx)
        type(geo_slot_t), intent(inout) :: geo
        type(init_ctx_t), intent(in)    :: ctx
        call component_enter(geo%c)
        call geo%c%init(ctx)
        call geo_declare(geo%c)
        call component_leave()
    end subroutine geo_init

    subroutine geo_init_state(geo, st, pass)
        type(geo_slot_t),  intent(inout) :: geo
        type(step_time_t), intent(in)    :: st
        integer,           intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(geo%c)
        call geo%c%init_state(s, pass)
        call component_close(geo%c, s)
        call component_leave()
    end subroutine geo_init_state

    subroutine geo_update(geo, st)
        type(geo_slot_t),  intent(inout) :: geo
        type(step_time_t), intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(geo%c)
        call geo%c%update(st)
        call component_close(geo%c, st)
        call component_leave()
    end subroutine geo_update

    subroutine geo_reduce(geo, st)
        ! The global part after the gather (section 6): fastisostasy
        ! bsl -= sum(dV_bsl)/A_ocn; puts bsl.
        type(geo_slot_t),  intent(inout) :: geo
        type(step_time_t), intent(in)    :: st
        if (trim(geo%c%scope) /= "global") then
            call fesm_error("kryos_geo::geo_reduce", "geo_reduce on a regional instance", &
                            "instance "//trim(geo%c%name))
        end if
        call component_enter(geo%c)
        call geo%c%reduce(st)
        call component_close(geo%c, st)
        call component_leave()
    end subroutine geo_reduce

    subroutine geo_scatter(geo, st)
        ! Regional exports bsl := import bsl (the reduced bsl just exchanged) and
        ! dphi := the same value (the uniform sea-surface shift of a barystatic
        ! geo, contracts-design 2.2), with their validity; component layer only
        ! (section 6). A port (opt%port) is scattered its own dphi: not touched.
        type(geo_slot_t),  intent(inout) :: geo
        type(step_time_t), intent(in)    :: st
        real(dp) :: v, t0, t1
        if (trim(geo%c%scope) /= "regional") then
            call fesm_error("kryos_geo::geo_scatter", "geo_scatter on a global instance", &
                            "instance "//trim(geo%c%name))
        end if
        call fields_get(geo%c%imports, geo%c%hi_bsl, v, t0=t0, t1=t1)
        if (t1 > st%t1) then
            call fesm_error("kryos_geo::geo_scatter", "imported bsl is newer than the step", &
                            "instance "//trim(geo%c%name))
        end if
        call fields_put(geo%c%exports, geo%c%h_bsl,  v, t1=t1, t0=t0)
        if (.not. geo%c%opt%port) call fields_put(geo%c%exports, geo%c%h_dphi, v, t1=t1, t0=t0)
    end subroutine geo_scatter

    subroutine geo_restart_write(geo, dir, t)
        type(geo_slot_t), intent(inout) :: geo
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        if (.not. geo%c%con%restart) return
        call component_enter(geo%c)
        call geo%c%restart_write(dir, t)
        call component_leave()
    end subroutine geo_restart_write

    subroutine geo_restart_read(geo, dir, t)
        type(geo_slot_t), intent(inout) :: geo
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        if (.not. geo%c%con%restart) return
        call component_enter(geo%c)
        call geo%c%restart_read(dir, t)
        call component_leave()
    end subroutine geo_restart_read

    function geo_restart_files(geo) result(files)
        ! The base names the backend's restart_write creates in the bundle folder.
        type(geo_slot_t), intent(in) :: geo
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => geo%c)                          ! (gfortran: ICE on a call through geo%c itself)
            files = c%restart_files()
        end associate
    end function geo_restart_files

    subroutine geo_write_output(geo, file, t, stream, init)
        type(geo_slot_t), intent(inout) :: geo
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        call component_enter(geo%c)
        call geo%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine geo_write_output

    subroutine geo_end(geo)
        type(geo_slot_t), intent(inout) :: geo
        call component_enter(geo%c)
        call geo%c%end()
        call component_leave()
    end subroutine geo_end

    subroutine geo_declare(c)
        ! Section 2.2 by scope. Regional z_bed_ref/H_ice_ref imports: the ice
        ! instance's reference, geo's initial state (decision 4 refined). A
        ! global part with opt%host (vilma2) has its own vocabulary on its own
        ! grid; the regional port of it (opt%port) imports the scattered fields
        ! instead of the ice (contracts-design 2.2).
        class(geo_t), intent(inout) :: c
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_init(c%exports, trim(c%name)//":exports")

        if (trim(c%scope) == "global" .and. c%opt%host) then
            call declare_host(c)
        else if (trim(c%scope) == "global") then
            call fields_add(c%imports, "dV_bsl", "m3", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_dV_bsl)
            call fields_add(c%exports, "bsl",    "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_bsl)
            call fields_add(c%exports, "A_ocn",  "m2", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_A_ocn)
        else
            if (c%opt%port) then
                call fields_add(c%imports, "dz_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dz_bed, grid=g, dims=d2)
                call fields_add(c%imports, "dz_ss",  "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dz_ss,  grid=g, dims=d2)
                call fields_add(c%imports, "dphi",   "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_dphi)
            else
                call fields_add(c%imports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice, grid=g, dims=d2)
            end if
            call fields_add(c%imports, "bsl",   "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_bsl)
            if (c%opt%A_ocn) then
                call fields_add(c%imports, "A_ocn", "m2", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_A_ocn)
                call component_optional(c, "A_ocn")
            end if
            if (c%opt%ref) then
                call fields_add(c%imports, "z_bed_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed_ref, &
                                grid=g, dims=d2)
                call fields_add(c%imports, "H_ice_ref", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice_ref, &
                                grid=g, dims=d2)
            end if

            call fields_add(c%exports, "dz_bed",    "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dz_bed,    grid=g, dims=d2)
            call fields_add(c%exports, "dz_ss",     "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dz_ss,     grid=g, dims=d2)
            call fields_add(c%exports, "bsl",       "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_bsl)
            call fields_add(c%exports, "dphi",      "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dphi)
            call fields_add(c%exports, "z_bed_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_z_bed_ref, grid=g, dims=d2)
            if (c%opt%dV_bsl) call fields_add(c%exports, "dV_bsl", "m3", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_dV_bsl)
            if (c%opt%dz_bed_dt) then
                call fields_add(c%exports, "dz_bed_dt", "m yr-1", FIELD_INSTANT, FIELD_INTENSIVE, c%h_dz_bed_dt, &
                                grid=g, dims=d2)
            end if
        end if

        call fields_lock(c%imports)
        call fields_lock(c%exports)
    end subroutine geo_declare

    subroutine declare_host(c)
        ! The global part of a host-resolved Earth (contracts-design 2.2, vilma2):
        ! rank 2 on the instance's own grid, which the backend sets in init (grid,
        ! nx, ny), the weighted fields per unit cell area. The imports with no
        ! component producer (the external ice and the background reference,
        ! globe_ice) are marked as filled by the core.
        class(geo_t), intent(inout) :: c
        integer :: d2(2)
        character(len=:), allocatable :: g

        if (c%nx < 1 .or. c%ny < 1 .or. len_trim(c%grid) == 0) then
            call fesm_error("kryos_geo::geo_init", "a host-resolved global geo has no grid", &
                            "instance "//trim(c%name)//": its init must set grid, nx and ny (opt%host)")
        end if
        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_add(c%imports, "H_ice_ext",   "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice_ext,   grid=g, dims=d2)
        call fields_add(c%imports, "host_w",      "1", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_host_w,      grid=g, dims=d2)
        call fields_add(c%imports, "host_Ig",     "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_host_Ig,     grid=g, dims=d2)
        call fields_add(c%imports, "host_C",      "1", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_host_C,      grid=g, dims=d2)
        call fields_add(c%imports, "host_F",      "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_host_F,      grid=g, dims=d2)
        call fields_add(c%imports, "z_bed_ref",   "m", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_z_bed_ref,   grid=g, dims=d2)
        call fields_add(c%imports, "H_ice_ref",   "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice_ref,   grid=g, dims=d2)
        call fields_add(c%imports, "host_Ig_ref", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_host_Ig_ref, grid=g, dims=d2)
        call fields_add(c%imports, "host_C_ref",  "1", FIELD_INSTANT, FIELD_INTENSIVE, c%hi_host_C_ref,  grid=g, dims=d2)
        call component_core_filled(c, "H_ice_ext")
        call component_core_filled(c, "z_bed_ref")
        call component_core_filled(c, "H_ice_ref")

        call fields_add(c%exports, "dz_bed",     "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_dz_bed, grid=g, dims=d2)
        call fields_add(c%exports, "dz_ss",      "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_dz_ss,  grid=g, dims=d2)
        call fields_add(c%exports, "C_ocn",      "1",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_C_ocn,  grid=g, dims=d2)
        call fields_add(c%exports, "dphi",       "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_dphi)
        call fields_add(c%exports, "bsl",        "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_bsl)
        call fields_add(c%exports, "bsl_ive",    "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_bsl_ive)
        call fields_add(c%exports, "A_ocn",      "m2", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_A_ocn)
        call fields_add(c%exports, "mass_resid", "1",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_mass_resid)
    end subroutine declare_host

end module kryos_geo
