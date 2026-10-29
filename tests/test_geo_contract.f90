module test_geo_toys
    ! Three toy geo instances that only declare their fields (the test puts the
    ! exports): a barystatic regional one, the regional port of a host-resolved
    ! Earth, and the global part of such an Earth, which sets its own grid.
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_geo,       only: geo_t, geo_opt_t
    implicit none
    private

    integer, parameter, public :: GNX = 8, GNY = 4         ! the toy Earth's grid

    type, extends(geo_t), public :: toy_reg_t
    contains
        procedure :: init => reg_init
        procedure :: update => toy_update
        procedure :: restart_write => toy_rw
        procedure :: restart_read => toy_rr
    end type toy_reg_t

    type, extends(toy_reg_t), public :: toy_port_t
    contains
        procedure :: init => port_init
    end type toy_port_t

    type, extends(toy_reg_t), public :: toy_host_t
    contains
        procedure :: init => host_init
    end type toy_host_t

    type, extends(toy_reg_t), public :: toy_nogrid_t
    contains
        procedure :: init => nogrid_init
    end type toy_nogrid_t

contains

    subroutine reg_init(self, ctx)
        class(toy_reg_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = geo_opt_t(A_ocn=.true., dV_bsl=.true., dz_bed_dt=.true., ref=.true.)
    end subroutine reg_init

    subroutine port_init(self, ctx)
        class(toy_port_t), intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = geo_opt_t(ref=.true., port=.true.)
    end subroutine port_init

    subroutine host_init(self, ctx)
        class(toy_host_t), intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        self%grid = "GAUSS-TOY"
        self%nx   = GNX
        self%ny   = GNY
        self%opt  = geo_opt_t(host=.true.)
    end subroutine host_init

    subroutine nogrid_init(self, ctx)
        class(toy_nogrid_t), intent(inout) :: self
        type(init_ctx_t),    intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        self%opt = geo_opt_t(host=.true.)
    end subroutine nogrid_init

    subroutine toy_update(self, st)
        class(toy_reg_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_update

    subroutine toy_rw(self, dir, t)
        class(toy_reg_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw

    subroutine toy_rr(self, dir, t)
        class(toy_reg_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rr

end module test_geo_toys

program test_geo_contract
    ! The geo contract of phase 4 (contracts-design 2.2) on toy instances,
    ! through the real wrappers (geo_adopt, geo_init, geo_scatter):
    !   * the regional vocabulary with dphi, unchanged otherwise (H_ice import,
    !     optional fields by option);
    !   * the regional port of a host-resolved Earth: scattered dz_bed, dz_ss,
    !     dphi, bsl in, no H_ice import, the same exports;
    !   * the global part with opt%host: its rank-2 imports and exports on its
    !     own grid, the 0-D ones, no dV_bsl, the imports the core fills marked;
    !   * geo_scatter: bsl and dphi of a barystatic geo, bsl only for a port;
    !   * a host part without a grid stops (child run, as test_wiring).

    use precision,       only: dp
    use coords,          only: grid_init
    use fesm_fields,     only: field_t, fields_has, fields_put, FIELD_INSTANT
    use kryos_component, only: init_ctx_t, step_time_t, census_t, constraints_t
    use kryos_geo,       only: geo_t, geo_slot_t, geo_adopt, geo_init, geo_scatter
    use test_geo_toys,   only: toy_reg_t, toy_port_t, toy_host_t, toy_nogrid_t, GNX, GNY

    implicit none

    integer, parameter :: NX = 5, NY = 4
    character(len=*), parameter :: errfile = "test_geo_contract_abort.err"
    integer :: nfail = 0

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    call test_regional()
    call test_port()
    call test_host()
    call test_scatter()
    call expect_abort("host_without_grid", "has no grid")

    call execute_command_line("rm -f "//errfile)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_geo_contract: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_geo_contract: OK"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine make_geo(slot, kind, scope)
        ! A toy geo adopted and initialised through the wrappers (regional: on a small cartesian grid).
        type(geo_slot_t), intent(out) :: slot
        character(len=*), intent(in)  :: kind                ! reg | port | host | nogrid
        character(len=*), intent(in)  :: scope
        class(geo_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        ctx%instance  = "geo_toy"
        ctx%component = "geo"
        ctx%backend   = "toy"
        ctx%scope     = scope
        if (trim(scope) == "regional") then
            ctx%domain    = "Toy"
            ctx%grid_name = "TOY"
            call grid_init(ctx%grid, name="TOY", mtype="cartesian", units="km", dx=16.0_dp, nx=NX, dy=16.0_dp, ny=NY, &
                           x0=0.0_dp, y0=0.0_dp)
        end if
        select case (trim(kind))
        case ("reg");    allocate(toy_reg_t :: proto)
        case ("port");   allocate(toy_port_t :: proto)
        case ("host");   allocate(toy_host_t :: proto)
        case ("nogrid"); allocate(toy_nogrid_t :: proto)
        end select
        call geo_adopt(slot, proto, constraints_t(regional=.true., global=.true.), ctx, census)
        call geo_init(slot, ctx)
    end subroutine make_geo

    logical function declared(st, name, rank, units, grid)
        ! `name` is declared in the field state with this rank and units, on `grid` (rank 2).
        use fesm_fields, only: field_state_t
        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name, units
        integer,             intent(in) :: rank
        character(len=*),    intent(in), optional :: grid
        integer :: h

        declared = fields_has(st, name)
        if (.not. declared) return
        h = 0
        do h = 1, st%n
            if (trim(st%f(h)%name) == trim(name)) exit
        end do
        declared = st%f(h)%rank == rank .and. trim(st%f(h)%units) == trim(units) .and. st%f(h)%tsem == FIELD_INSTANT
        if (present(grid) .and. declared) declared = trim(st%f(h)%grid) == trim(grid)
    end function declared

    subroutine test_regional()
        type(geo_slot_t) :: g

        write(*,"(a)") "regional geo (barystatic)"
        call make_geo(g, "reg", "regional")
        associate (i => g%c%imports, e => g%c%exports)
            call check("imports: H_ice (rank 2), bsl, A_ocn, z_bed_ref, H_ice_ref", &
                       declared(i, "H_ice", 2, "m", "TOY") .and. declared(i, "bsl", 0, "m") .and. &
                       declared(i, "A_ocn", 0, "m2") .and. declared(i, "z_bed_ref", 2, "m", "TOY") .and. &
                       declared(i, "H_ice_ref", 2, "m", "TOY"))
            call check("no port imports", .not. fields_has(i, "dz_bed") .and. .not. fields_has(i, "dz_ss") .and. &
                       .not. fields_has(i, "dphi"))
            call check("exports: dz_bed, dz_ss, z_bed_ref (rank 2), bsl, dphi, dV_bsl, dz_bed_dt", &
                       declared(e, "dz_bed", 2, "m", "TOY") .and. declared(e, "dz_ss", 2, "m", "TOY") .and. &
                       declared(e, "z_bed_ref", 2, "m", "TOY") .and. declared(e, "bsl", 0, "m") .and. &
                       declared(e, "dphi", 0, "m") .and. declared(e, "dV_bsl", 0, "m3") .and. &
                       declared(e, "dz_bed_dt", 2, "m yr-1", "TOY"))
            call check("handles set", g%c%h_dphi > 0 .and. g%c%h_bsl > 0 .and. g%c%hi_H_ice > 0)
            call check("nothing of the host vocabulary", .not. fields_has(e, "C_ocn") .and. &
                       .not. fields_has(i, "host_Ig") .and. .not. allocated(g%c%core_imports))
        end associate
    end subroutine test_regional

    subroutine test_port()
        type(geo_slot_t) :: g

        write(*,"(a)") "regional port of a host-resolved Earth"
        call make_geo(g, "port", "regional")
        associate (i => g%c%imports, e => g%c%exports)
            call check("imports: dz_bed, dz_ss (rank 2 on the geo grid), dphi, bsl, z_bed_ref, H_ice_ref", &
                       declared(i, "dz_bed", 2, "m", "TOY") .and. declared(i, "dz_ss", 2, "m", "TOY") .and. &
                       declared(i, "dphi", 0, "m") .and. declared(i, "bsl", 0, "m") .and. &
                       declared(i, "z_bed_ref", 2, "m", "TOY") .and. declared(i, "H_ice_ref", 2, "m", "TOY"))
            call check("no H_ice import (the ice -> geo H_ice row is then inactive)", .not. fields_has(i, "H_ice"))
            call check("exports: dz_bed, dz_ss, z_bed_ref, bsl, dphi; no dV_bsl", &
                       declared(e, "dz_bed", 2, "m", "TOY") .and. declared(e, "dz_ss", 2, "m", "TOY") .and. &
                       declared(e, "z_bed_ref", 2, "m", "TOY") .and. declared(e, "bsl", 0, "m") .and. &
                       declared(e, "dphi", 0, "m") .and. .not. fields_has(e, "dV_bsl"))
            call check("import handles set", g%c%hi_dz_bed > 0 .and. g%c%hi_dz_ss > 0 .and. g%c%hi_dphi > 0)
        end associate
    end subroutine test_port

    subroutine test_host()
        type(geo_slot_t) :: g
        character(len=*), parameter :: G_ = "GAUSS-TOY"

        write(*,"(a)") "global part with opt%host"
        call make_geo(g, "host", "global")
        associate (i => g%c%imports, e => g%c%exports)
            call check("rank-2 imports on the instance's own grid", &
                       declared(i, "H_ice_ext", 2, "m", G_) .and. declared(i, "host_w", 2, "1", G_) .and. &
                       declared(i, "host_Ig", 2, "m", G_) .and. declared(i, "host_C", 2, "1", G_) .and. &
                       declared(i, "host_F", 2, "m", G_) .and. declared(i, "z_bed_ref", 2, "m", G_) .and. &
                       declared(i, "H_ice_ref", 2, "m", G_) .and. declared(i, "host_Ig_ref", 2, "m", G_) .and. &
                       declared(i, "host_C_ref", 2, "1", G_))
            call check("shape of the grid", i%f(g%c%hi_host_Ig)%nx == GNX .and. i%f(g%c%hi_host_Ig)%ny == GNY)
            call check("rank-2 exports: dz_bed, dz_ss, C_ocn", &
                       declared(e, "dz_bed", 2, "m", G_) .and. declared(e, "dz_ss", 2, "m", G_) .and. &
                       declared(e, "C_ocn", 2, "1", G_))
            call check("0-D exports: dphi, bsl, bsl_ive, A_ocn, mass_resid", &
                       declared(e, "dphi", 0, "m") .and. declared(e, "bsl", 0, "m") .and. &
                       declared(e, "bsl_ive", 0, "m") .and. declared(e, "A_ocn", 0, "m2") .and. &
                       declared(e, "mass_resid", 0, "1"))
            call check("no dV_bsl import (the host part is not reduced)", .not. fields_has(i, "dV_bsl"))
            call check("the imports with no component producer are marked", allocated(g%c%core_imports))
            if (allocated(g%c%core_imports)) then
                call check("marked: H_ice_ext, z_bed_ref, H_ice_ref (and only these)", &
                           size(g%c%core_imports) == 3 .and. any(g%c%core_imports == "H_ice_ext") .and. &
                           any(g%c%core_imports == "z_bed_ref") .and. any(g%c%core_imports == "H_ice_ref"))
            end if
        end associate
    end subroutine test_host

    subroutine test_scatter()
        ! geo_scatter: bsl := import bsl and, for a barystatic geo, dphi := the same value.
        type(geo_slot_t) :: g, p

        write(*,"(a)") "geo_scatter"
        call make_geo(g, "reg", "regional")
        call make_geo(p, "port", "regional")
        call fields_put(g%c%imports, g%c%hi_bsl, -1.25_dp/3.0_dp, t1=10.0_dp, t0=0.0_dp)
        call fields_put(p%c%imports, p%c%hi_bsl, -1.25_dp/3.0_dp, t1=10.0_dp, t0=0.0_dp)
        call geo_scatter(g, step_time_t(t0=0.0_dp, t1=10.0_dp, t_eval=10.0_dp, dt=10.0_dp, dtt=10.0_dp))
        call geo_scatter(p, step_time_t(t0=0.0_dp, t1=10.0_dp, t_eval=10.0_dp, dt=10.0_dp, dtt=10.0_dp))
        associate (e => g%c%exports)
            call check("barystatic: bsl and dphi exported, equal, valid at the import's time", &
                       e%f(g%c%h_bsl)%v(1,1,1) == -1.25_dp/3.0_dp .and. e%f(g%c%h_dphi)%v(1,1,1) == -1.25_dp/3.0_dp .and. &
                       e%f(g%c%h_bsl)%t1 == 10.0_dp .and. e%f(g%c%h_dphi)%t1 == 10.0_dp .and. e%f(g%c%h_dphi)%published)
        end associate
        call check("port: bsl exported, dphi left to its own scatter", &
                   p%c%exports%f(p%c%h_bsl)%published .and. .not. p%c%exports%f(p%c%h_dphi)%published)
    end subroutine test_scatter

    subroutine expect_abort(case_name, text)
        character(len=*), intent(in) :: case_name, text
        character(len=512) :: exe, line
        character(len=:), allocatable :: all_
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > /dev/null 2> "//errfile, exitstat=es)
        all_ = ""
        open(newunit=u, file=errfile, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(u, "(a)", iostat=ios) line
                if (ios /= 0) exit
                all_ = all_//trim(line)//" "
            end do
            close(u)
        end if
        call check(case_name//": aborts (exit 1) with '"//text//"'", es == 1 .and. index(all_, text) > 0)
    end subroutine expect_abort

    subroutine abort_case()
        character(len=64) :: arg
        type(geo_slot_t) :: g

        call get_command_argument(1, arg)
        select case (trim(arg))
        case ("host_without_grid")
            call make_geo(g, "nogrid", "global")
        end select
    end subroutine abort_case

end program test_geo_contract
