module test_ice_toys
    ! Two toy ice instances that only declare their fields (the test puts the
    ! exports): a bare one, and one with the phase 4 additions (footprint,
    ! H_ice_pd, a flotation binding).
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_ice,       only: ice_t, ice_opt_t
    implicit none
    private

    real(dp), parameter, public :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp

    type, extends(ice_t), public :: toy_bare_t
    contains
        procedure :: init => bare_init
        procedure :: update => toy_update
        procedure :: restart_write => toy_rw
        procedure :: restart_read => toy_rr
    end type toy_bare_t

    type, extends(toy_bare_t), public :: toy_full_t
    contains
        procedure :: init => full_init
        procedure :: flotation => full_flotation
    end type toy_full_t

contains

    subroutine bare_init(self, ctx)
        class(toy_bare_t), intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
    end subroutine bare_init

    subroutine full_init(self, ctx)
        class(toy_full_t), intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = ice_opt_t(z_bed_ref=.true., footprint=.true., H_ice_pd=.true.)
        self%caps%has_flotation = .true.
    end subroutine full_init

    subroutine full_flotation(self, H, z_bed, z_sl, f_grnd)
        ! Grounded where the ice outweighs the water it would displace.
        class(toy_full_t), intent(in)  :: self
        real(dp),          intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),          intent(out) :: f_grnd(:,:)
        associate (unused => self%nx); end associate
        f_grnd = merge(1.0_dp, 0.0_dp, H*RHO_I/RHO_SW + (z_bed - z_sl) >= 0.0_dp)
    end subroutine full_flotation

    subroutine toy_update(self, st)
        class(toy_bare_t), intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_update

    subroutine toy_rw(self, dir, t)
        class(toy_bare_t), intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw

    subroutine toy_rr(self, dir, t)
        class(toy_bare_t), intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rr

end module test_ice_toys

program test_ice_contract
    ! The ice contract of phase 4 (contracts-design 2.1) on toy ice instances,
    ! through the real wrappers (ice_adopt, ice_init, ice_flotation):
    !   * the optional static exports footprint and H_ice_pd are declared
    !     exactly when the backend sets their option, with their units, shape
    !     and nature;
    !   * the flotation binding: the wrapper passes the arrays to the backend's
    !     rule and back unchanged; a backend without the capability, or arrays
    !     off the instance grid, stop with a named message (the program re-runs
    !     itself with an argument per case, as test_wiring).
    ! (The yelmo adapter's values, and its flotation against Yelmo's own f_grnd,
    ! are checked in test_run_sequence, which links Yelmo.)

    use precision,       only: dp
    use coords,          only: grid_init
    use fesm_fields,     only: fields_has, fields_put, FIELD_INSTANT, FIELD_INTENSIVE
    use kryos_component, only: init_ctx_t, census_t, constraints_t
    use kryos_ice,       only: ice_t, ice_slot_t, ice_adopt, ice_init, ice_flotation
    use test_ice_toys,   only: toy_full_t, toy_bare_t, RHO_I, RHO_SW

    implicit none

    integer, parameter :: NX = 5, NY = 4
    character(len=*), parameter :: errfile = "test_ice_contract_abort.err"
    integer :: nfail = 0

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    call test_declaration()
    call test_flotation()
    call test_flotation_stops()

    call execute_command_line("rm -f "//errfile)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_ice_contract: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_ice_contract: OK"

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

    subroutine make_ice(slot, full)
        ! A toy ice on a small cartesian grid, adopted and initialised through the wrappers.
        type(ice_slot_t), intent(out) :: slot
        logical,          intent(in)  :: full
        class(ice_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        ctx%instance  = "ice_toy"
        ctx%component = "ice"
        ctx%backend   = "toy"
        ctx%domain    = "Toy"
        ctx%grid_name = "TOY"
        call grid_init(ctx%grid, name="TOY", mtype="cartesian", units="km", dx=16.0_dp, nx=NX, dy=16.0_dp, ny=NY, &
                       x0=0.0_dp, y0=0.0_dp)
        if (full) then
            allocate(toy_full_t :: proto)
        else
            allocate(toy_bare_t :: proto)
        end if
        call ice_adopt(slot, proto, constraints_t(), ctx, census)
        call ice_init(slot, ctx)
    end subroutine make_ice

    subroutine test_declaration()
        type(ice_slot_t) :: full, bare
        real(dp) :: a(NX,NY), b(NX,NY)
        integer :: i, j

        write(*,"(a)") "declaration of footprint and H_ice_pd"
        call make_ice(full, .true.)
        call make_ice(bare, .false.)
        associate (c => full%c, e => full%c%exports)
            call check("footprint declared on request", fields_has(e, "footprint") .and. c%h_footprint > 0)
            call check("H_ice_pd declared on request", fields_has(e, "H_ice_pd") .and. c%h_H_ice_pd > 0)
            call check("footprint: 1, instant, intensive, rank 2 on the instance grid", &
                       trim(e%f(c%h_footprint)%units) == "1" .and. e%f(c%h_footprint)%tsem == FIELD_INSTANT .and. &
                       e%f(c%h_footprint)%cons == FIELD_INTENSIVE .and. e%f(c%h_footprint)%rank == 2 .and. &
                       trim(e%f(c%h_footprint)%grid) == "TOY" .and. e%f(c%h_footprint)%nx == NX .and. &
                       e%f(c%h_footprint)%ny == NY)
            call check("H_ice_pd: m, instant, intensive, rank 2 on the instance grid", &
                       trim(e%f(c%h_H_ice_pd)%units) == "m" .and. e%f(c%h_H_ice_pd)%tsem == FIELD_INSTANT .and. &
                       e%f(c%h_H_ice_pd)%cons == FIELD_INTENSIVE .and. e%f(c%h_H_ice_pd)%rank == 2 .and. &
                       trim(e%f(c%h_H_ice_pd)%grid) == "TOY")
            call check("neither is an import", .not. fields_has(c%imports, "footprint") .and. &
                       .not. fields_has(c%imports, "H_ice_pd"))

            ! a static put, once: the value and its time come back
            do j = 1, NY
                do i = 1, NX
                    a(i,j) = merge(1.0_dp, 0.0_dp, i + j > 4)
                    b(i,j) = 100.0_dp*real(i, dp) + real(j, dp)
                end do
            end do
            call fields_put(full%c%exports, full%c%h_footprint, a, t1=0.0_dp)
            call fields_put(full%c%exports, full%c%h_H_ice_pd, b, t1=0.0_dp)
            call check("footprint and H_ice_pd published at t0", e%f(c%h_footprint)%published .and. &
                       e%f(c%h_H_ice_pd)%published .and. e%f(c%h_footprint)%t1 == 0.0_dp .and. &
                       e%f(c%h_H_ice_pd)%t1 == 0.0_dp)
        end associate

        call check("a backend that does not set the options declares neither", &
                   .not. fields_has(bare%c%exports, "footprint") .and. .not. fields_has(bare%c%exports, "H_ice_pd") &
                   .and. bare%c%h_footprint == 0 .and. bare%c%h_H_ice_pd == 0)
        call check("the rest of the vocabulary is unchanged", fields_has(bare%c%exports, "H_ice") .and. &
                   fields_has(bare%c%exports, "f_grnd") .and. .not. fields_has(bare%c%exports, "z_bed_ref"))
    end subroutine test_declaration

    subroutine test_flotation()
        type(ice_slot_t) :: full, bare
        real(dp) :: H(NX,NY), zb(NX,NY), zs(NX,NY), f(NX,NY), H0(NX,NY), zb0(NX,NY), zs0(NX,NY), e(NX,NY)
        integer :: i, j

        write(*,"(a)") "flotation binding"
        call make_ice(full, .true.)
        call make_ice(bare, .false.)
        call check("capability: declared by the backend that implements it only", &
                   full%c%caps%has_flotation .and. .not. bare%c%caps%has_flotation)

        do j = 1, NY
            do i = 1, NX
                H(i,j)  = 400.0_dp*real(i, dp) - 150.0_dp*real(j, dp) + 1.0_dp/3.0_dp     ! some cells ice free
                zb(i,j) = -900.0_dp + 250.0_dp*real(i, dp) - 60.0_dp*real(j, dp)
                zs(i,j) = 12.5_dp - 0.5_dp*real(i + j, dp)
            end do
        end do
        H = max(H, 0.0_dp)
        H0 = H;  zb0 = zb;  zs0 = zs
        e = merge(1.0_dp, 0.0_dp, H*RHO_I/RHO_SW + (zb - zs) >= 0.0_dp)

        call ice_flotation(full, H, zb, zs, f)
        call check("the wrapper returns the backend's rule on the arrays it was given", all(f == e))
        call check("grounded and floating cells both occur (the test is not vacuous)", any(f == 1.0_dp) .and. any(f == 0.0_dp))
        call check("the inputs are unchanged", all(H == H0) .and. all(zb == zb0) .and. all(zs == zs0))

        ! the same rule on a different sea level, per call (nothing is cached)
        zs = zs + 500.0_dp
        call ice_flotation(full, H, zb, zs, f)
        e = merge(1.0_dp, 0.0_dp, H*RHO_I/RHO_SW + (zb - zs) >= 0.0_dp)
        call check("a second call with another sea level", all(f == e))
    end subroutine test_flotation

    subroutine test_flotation_stops()
        write(*,"(a)") "flotation: named stops"
        call expect_abort("no_capability", "does not declare has_flotation")
        call expect_abort("off_grid", "arrays are not on the instance grid")
    end subroutine test_flotation_stops

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
        ! Child: build the situation and call the wrapper.
        character(len=64) :: arg
        type(ice_slot_t) :: full, bare
        real(dp) :: H(NX,NY), f(NX,NY), g(NX-1,NY)

        call get_command_argument(1, arg)
        H = 100.0_dp
        select case (trim(arg))
        case ("no_capability")
            call make_ice(bare, .false.)
            call ice_flotation(bare, H, H, H, f)
        case ("off_grid")
            call make_ice(full, .true.)
            call ice_flotation(full, H(1:NX-1,:), H(1:NX-1,:), H(1:NX-1,:), g)
        end select
    end subroutine abort_case

end program test_ice_contract
