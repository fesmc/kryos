module kryos_geo
    ! The geo component, solid earth and sea level (contracts-design sections
    ! 2.2, 3.3, 6): one abstract geo_t for both levels, declared by scope
    ! (regional: bed and sea-surface anomalies; global: bsl and the reduction),
    ! its slot, wrappers incl. geo_reduce and geo_scatter, and the factory.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, fields_get, fields_put, &
                               FIELD_INSTANT, FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional

    implicit none

    private

    type, public :: geo_opt_t                       ! optional fields the backend supports (set in init)
        logical :: A_ocn = .false., dV_bsl = .false., dz_bed_dt = .false.
        logical :: ref = .false.                    ! regional imports z_bed_ref, H_ice_ref (decision 4)
    end type geo_opt_t

    type, abstract, extends(component_t), public :: geo_t
        type(geo_opt_t) :: opt
        integer :: hi_H_ice = 0, hi_bsl = 0, hi_A_ocn = 0, hi_z_bed_ref = 0, hi_H_ice_ref = 0  ! regional imports
        integer :: hi_dV_bsl = 0                                                               ! global import
        integer :: h_dz_bed = 0, h_dz_ss = 0, h_dV_bsl = 0, h_z_bed_ref = 0, h_dz_bed_dt = 0  ! regional exports
        integer :: h_bsl = 0, h_A_ocn = 0                                                      ! both / global
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
        ! Regional export bsl := import bsl (the reduced bsl just exchanged),
        ! with its validity; component layer only (section 6).
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
        call fields_put(geo%c%exports, geo%c%h_bsl, v, t1=t1, t0=t0)
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
        ! instance's reference, geo's initial state (decision 4 refined).
        class(geo_t), intent(inout) :: c
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_init(c%exports, trim(c%name)//":exports")

        if (trim(c%scope) == "global") then
            call fields_add(c%imports, "dV_bsl", "m3", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_dV_bsl)
            call fields_add(c%exports, "bsl",    "m",  FIELD_INSTANT, FIELD_INTENSIVE, c%h_bsl)
            call fields_add(c%exports, "A_ocn",  "m2", FIELD_INSTANT, FIELD_EXTENSIVE, c%h_A_ocn)
        else
            call fields_add(c%imports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_H_ice, grid=g, dims=d2)
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

end module kryos_geo
