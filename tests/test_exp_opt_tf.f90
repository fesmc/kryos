module test_opt_tf_toys
    ! A toy ice and a toy shelf: only their declared fields and the shelf component's corrections
    ! matter to the experiment opt_tf (it uses no backend type).
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_ice,       only: ice_t, ice_opt_t
    use kryos_shelf,     only: shelf_t
    implicit none
    private

    type, extends(ice_t), public :: toy_ice_t
    contains
        procedure :: init => ice_init_
        procedure :: update => ice_update_
        procedure :: restart_write => ice_rw
        procedure :: restart_read => ice_rw
    end type toy_ice_t

    type, extends(toy_ice_t), public :: toy_ice_nogpd_t     ! no H_grnd_pd export
    contains
        procedure :: init => ice_init_nogpd
    end type toy_ice_nogpd_t

    type, extends(shelf_t), public :: toy_shelf_t
    contains
        procedure :: init => shelf_init_
        procedure :: update => shelf_update_
        procedure :: restart_write => shelf_rw
        procedure :: restart_read => shelf_rw
    end type toy_shelf_t

contains

    subroutine ice_init_(self, ctx)
        class(toy_ice_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = ice_opt_t(H_ice_pd=.true., H_grnd=.true., H_grnd_pd=.true.)
    end subroutine ice_init_
    subroutine ice_init_nogpd(self, ctx)
        class(toy_ice_nogpd_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = ice_opt_t(H_ice_pd=.true., H_grnd=.true., H_grnd_pd=.false.)
    end subroutine ice_init_nogpd
    subroutine ice_update_(self, st)
        class(toy_ice_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine ice_update_
    subroutine ice_rw(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine ice_rw

    subroutine shelf_init_(self, ctx)
        class(toy_shelf_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx
        associate (u => ctx%time + self%nx); end associate
    end subroutine shelf_init_
    subroutine shelf_update_(self, st)
        class(toy_shelf_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine shelf_update_
    subroutine shelf_rw(self, dir, t)
        class(toy_shelf_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine shelf_rw

end module test_opt_tf_toys


program test_exp_opt_tf
    ! The experiment opt_tf (core-design 3.9, phase 2 WP3) on a hand-made domain_t: a toy ice that
    ! exports H_ice, H_grnd, dHidt, H_ice_pd, H_grnd_pd, a toy shelf that holds the corrections
    ! (shelf_t%corr). No run, no Yelmo's state: the experiment reads the contract fields only.
    !
    !   * the hook calls Yelmo's optimize_tf_corr_basin (L21) and optimize_tf_corr (L21-points) on
    !     the shelf component's tf_corr, on its basins and grid spacing, with the ice's exports
    !     narrowed to sp and remapped ice -> shelf grid: same grid (a copy), and ANT-32KM -> ANT-16KM
    !     (bilinear, the shelf being finer); the result equals the same calls made by this test,
    !     bit for bit, and differs from the start state;
    !   * the window of the elapsed time: outside it nothing changes;
    !   * stops (child processes, "abort <case>"): the ice lacks an export, no shelf instance;
    !   * the tf_basins list of L21 (a basin list and the default -1 for all).

    use precision,       only: sp, dp
    use coords,          only: grid_init, grid_cdo_read_desc
    use coupler,         only: coupler_init, remap
    use fesm_clock,      only: clock_t, clock_init, clock_start, cadence_init
    use fesm_fields,     only: fields_put, fields_get
    use fesm_log,        only: log_set_level, LOG_WARN
    use phys_constants,  only: phys_const_class, phys_const_set
    use ice_optimization, only: optimize_tf_corr, optimize_tf_corr_basin
    use kryos_component, only: init_ctx_t, census_t, constraints_t, step_time_t, component_close
    use kryos_config,    only: config_t, config_init
    use kryos_domain,    only: domain_t, C_ICE, C_SHELF
    use kryos_ice,       only: ice_t, ice_slot_t, ice_adopt, ice_init
    use kryos_experiment
    use test_opt_tf_toys

    implicit none

    character(len=*), parameter :: dir = "tmp/test_exp_opt_tf/"
    character(len=*), parameter :: par = dir//"par.nml"
    integer :: nfail = 0
    character(len=512) :: arg
    type(phys_const_class) :: cnst
    type(config_t) :: cfg
    type(clock_t), target :: clk

    call log_set_level(LOG_WARN)
    call phys_const_set(cnst, label="Earth", rho_ice=910.0_dp, rho_sw=1028.0_dp, rho_w=1000.0_dp, g=9.81_dp)
    call clock_init(clk, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")
    call clock_start(clk)

    call get_command_argument(1, arg)
    if (trim(arg) == "abort") then
        call get_command_argument(2, arg)
        call config_init(cfg, par)
        call abort_case(trim(arg))
        stop 0
    end if

    call execute_command_line("mkdir -p "//dir)
    call write_par()
    call config_init(cfg, par)

    call test_case("SamePts",    "TOY",      "TOY",      .true.)
    call test_case("SameBasin",  "TOY",      "TOY",      .false.)
    call test_case("RemapPts",   "ANT-32KM", "ANT-16KM", .true.)
    call test_case("RemapBasin", "ANT-32KM", "ANT-16KM", .false.)
    call test_case("BasinList",  "TOY",      "TOY",      .false., tf_basins=[2, 3])
    call test_stops()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_exp_opt_tf: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_exp_opt_tf: OK"

contains

    subroutine make_domain(d, name, gice, gshelf, with_shelf, no_grnd_pd)
        ! A domain with a toy ice on grid gice (every step), a toy shelf on gshelf holding
        ! corrections, and the coupler of the maps/ folder.
        type(domain_t),   intent(out) :: d
        character(len=*), intent(in)  :: name, gice, gshelf
        logical,          intent(in)  :: with_shelf
        logical,          intent(in), optional :: no_grnd_pd

        class(ice_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census
        type(clock_t)    :: c2
        type(init_ctx_t) :: sctx
        integer :: i, j

        d%name = name
        d%grid(C_ICE)   = gice
        d%grid(C_SHELF) = gshelf
        call coupler_init(d%cpl, "maps")

        ctx%instance = "ice_toy";  ctx%component = "ice";  ctx%backend = "toy";  ctx%domain = name
        ctx%grid_name = gice
        call grid_of(ctx%grid, gice)
        if (present(no_grnd_pd)) then
            if (no_grnd_pd) then
                allocate(toy_ice_nogpd_t :: proto)
            else
                allocate(toy_ice_t :: proto)
            end if
        else
            allocate(toy_ice_t :: proto)
        end if
        call ice_adopt(d%ice, proto, constraints_t(), ctx, census)
        call ice_init(d%ice, ctx)
        call clock_init(c2, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")
        call cadence_init(d%ice%c%cadence, c2, "ice_toy", .true.)

        if (.not. with_shelf) return
        allocate(toy_shelf_t :: d%shelf%c)
        call grid_of(sctx%grid, gshelf)
        associate (corr => d%shelf%c%corr, g => sctx%grid%G)
            corr%nx = g%nx;  corr%ny = g%ny
            corr%dx = real(g%dx * sctx%grid%cs%xy_conv, sp)
            corr%xc = real(g%x, sp);  corr%yc = real(g%y, sp)
            allocate(corr%tf_corr(g%nx, g%ny), corr%tf_corr_basin(g%nx, g%ny), corr%basins(g%nx, g%ny))
            corr%tf_corr = 0.0_sp
            corr%tf_corr_basin = 0.0_sp
            do j = 1, g%ny
                do i = 1, g%nx
                    corr%basins(i,j) = real(mod(i*j/17, 4), sp)     ! ids 0 .. 3
                end do
            end do
        end associate
    end subroutine make_domain

    subroutine grid_of(grid, name)
        use coords, only: grid_class
        type(grid_class), intent(out) :: grid
        character(len=*), intent(in)  :: name
        if (trim(name) == "TOY") then
            call grid_init(grid, name="TOY", mtype="cartesian", units="km", dx=16.0_dp, nx=24, dy=16.0_dp, ny=20, &
                           x0=0.0_dp, y0=0.0_dp)
        else
            call grid_cdo_read_desc(grid, trim(name), "maps")
        end if
    end subroutine grid_of

    subroutine put_ice(d)
        ! Exports of the ice at the end of the previous step: smooth, sp-exact values, an error
        ! against the present day that varies over the grid.
        type(domain_t), intent(inout) :: d
        real(dp), allocatable :: a(:,:)
        type(step_time_t) :: s
        integer :: i, j, nx, ny

        nx = d%ice%c%nx;  ny = d%ice%c%ny
        allocate(a(nx, ny))
        do j = 1, ny
            do i = 1, nx
                a(i,j) = 400.0_dp + 3.0_dp*i + 5.0_dp*j
            end do
        end do
        call fields_put(d%ice%c%exports, d%ice%c%h_H_ice, sp_exact(a), t1=0.0_dp)
        do j = 1, ny
            do i = 1, nx
                a(i,j) = 400.0_dp + 3.0_dp*i + 5.0_dp*j - 40.0_dp + 11.0_dp*mod(i + 2*j, 7)
            end do
        end do
        call fields_put(d%ice%c%exports, d%ice%c%h_H_ice_pd, sp_exact(a), t1=0.0_dp)
        do j = 1, ny
            do i = 1, nx
                a(i,j) = 40.0_dp + 25.0_dp*mod(i + j, 11)
            end do
        end do
        call fields_put(d%ice%c%exports, d%ice%c%h_H_grnd, sp_exact(a), t1=0.0_dp)
        do j = 1, ny
            do i = 1, nx
                a(i,j) = 40.0_dp + 25.0_dp*mod(i + j, 11) + 7.0_dp*mod(j, 3) - 7.0_dp
            end do
        end do
        call fields_put(d%ice%c%exports, d%ice%c%h_H_grnd_pd, sp_exact(a), t1=0.0_dp)
        do j = 1, ny
            do i = 1, nx
                a(i,j) = -0.8_dp + 0.2_dp*mod(i, 5) - 0.1_dp*mod(j, 3)
            end do
        end do
        call fields_put(d%ice%c%exports, d%ice%c%h_dHidt, sp_exact(a), t1=0.0_dp, t0=-10.0_dp)
        a = 0.0_dp                                     ! the other MEAN fields of the ice, for the closing below
        call fields_put(d%ice%c%exports, d%ice%c%h_smb,    a, t1=0.0_dp, t0=-10.0_dp)
        call fields_put(d%ice%c%exports, d%ice%c%h_bmb,    a, t1=0.0_dp, t0=-10.0_dp)
        call fields_put(d%ice%c%exports, d%ice%c%h_cmb,    a, t1=0.0_dp, t0=-10.0_dp)
        call fields_put(d%ice%c%exports, d%ice%c%h_mb_net, a, t1=0.0_dp, t0=-10.0_dp)
        s%t0 = -10.0_dp;  s%t1 = 0.0_dp;  s%dt = 10.0_dp;  s%dtt = 10.0_dp      ! closes the MEAN field, as the wrappers do
        call component_close(d%ice%c, s)
    end subroutine put_ice

    function sp_exact(a) result(b)
        real(dp), intent(in) :: a(:,:)
        real(dp) :: b(size(a,1), size(a,2))
        b = real(real(a, sp), dp)
    end function sp_exact

    subroutine test_case(name, gice, gshelf, points, tf_basins)
        character(len=*), intent(in) :: name, gice, gshelf
        logical,          intent(in) :: points
        integer,          intent(in), optional :: tf_basins(:)      ! L21: the list of the experiment group

        type(domain_t) :: d
        type(experiment_list_t) :: list
        type(hook_time_t) :: st
        real(sp), allocatable :: exp_corr(:,:)
        real(sp), allocatable :: Hm(:,:), Hgm(:,:), dHm(:,:), Hpdm(:,:), Hgpdm(:,:), v(:,:)
        character(len=32) :: mth
        integer :: tfb(100)
        character(len=:), allocatable :: label

        label = name//" ("//gice//" -> "//gshelf//", "//merge("L21-points", "L21       ", points)//")"
        write(*,"(a)") label
        call make_domain(d, name, gice, gshelf, .true.)
        call put_ice(d)
        call experiments_create(list, cfg, d, clk, cnst)
        call check(label//": created, wants STEP_BEGIN only", experiments_count(list) == 1 .and. &
                   list%s(1)%c%wants(HOOK_STEP_BEGIN) .and. count(list%s(1)%c%wants) == 1)

        st%t0 = 90.0_dp;  st%t1 = 100.0_dp;  st%dtt = 10.0_dp;  st%clk => clk

        ! outside the window [10, 50] of the elapsed time: nothing changes
        st%t_elapsed = 5.0_dp
        call experiments_hook(list, d, HOOK_STEP_BEGIN, st)
        st%t_elapsed = 60.0_dp
        call experiments_hook(list, d, HOOK_STEP_BEGIN, st)
        call check(label//": outside the window nothing changes", all(d%shelf%c%corr%tf_corr == 0.0_sp))

        ! inside: the same calls made here
        mth = "con"
        if (gshelf /= gice) mth = "bilin"                 ! the shelf grid is the finer
        v = narrow(d, d%ice%c%h_H_ice);   call remap(d%cpl, v, gice, Hm,   gshelf, method=trim(mth))
        v = narrow(d, d%ice%c%h_H_grnd);  call remap(d%cpl, v, gice, Hgm,  gshelf, method=trim(mth))
        v = narrow(d, d%ice%c%h_dHidt);   call remap(d%cpl, v, gice, dHm,  gshelf, method=trim(mth))
        v = narrow(d, d%ice%c%h_H_ice_pd); call remap(d%cpl, v, gice, Hpdm, gshelf, method=trim(mth))
        allocate(exp_corr(size(d%shelf%c%corr%tf_corr,1), size(d%shelf%c%corr%tf_corr,2)))
        exp_corr = 0.0_sp
        if (points) then
            v = narrow(d, d%ice%c%h_H_grnd_pd); call remap(d%cpl, v, gice, Hgpdm, gshelf, method=trim(mth))
            call optimize_tf_corr(exp_corr, Hm, Hgm, dHm, Hpdm, Hgpdm, 500.0_sp, d%shelf%c%corr%basins, .false., &
                                  100.0_sp, 10.0_sp, -1.0_sp, 1.0_sp, d%shelf%c%corr%dx, sigma=50.0e3_sp, dt=10.0_dp)
        else
            tfb = 0
            tfb(1) = -1
            if (present(tf_basins)) tfb(1:size(tf_basins)) = tf_basins
            call optimize_tf_corr_basin(exp_corr, Hm, Hgm, dHm, Hpdm, d%shelf%c%corr%basins, 500.0_sp, &
                                        100.0_sp, 10.0_sp, -1.0_sp, 1.0_sp, tfb, dt=10.0_dp)
        end if
        st%t_elapsed = 30.0_dp
        call experiments_hook(list, d, HOOK_STEP_BEGIN, st)
        call check(label//": tf_corr = the library call on the narrowed, remapped exports, bit for bit", &
                   all(d%shelf%c%corr%tf_corr == exp_corr))
        call check(label//": tf_corr has moved from 0 ("//trim(itoa(count(exp_corr /= 0.0_sp)))//" cells)", any(exp_corr /= 0.0_sp))
        call check(label//": the corrections by region are untouched", all(d%shelf%c%corr%tf_corr_basin == 0.0_sp))
        call experiments_end(list)
    end subroutine test_case

    function narrow(d, h) result(v)
        type(domain_t), intent(in) :: d
        integer,        intent(in) :: h
        real(sp), allocatable :: v(:,:)
        allocate(v(d%ice%c%nx, d%ice%c%ny))
        call fields_get(d%ice%c%exports, h, v)
    end function narrow

    subroutine test_stops()
        write(*,"(a)") "stops"
        call stops("no_export", "needs the ice export", "the ice does not export H_grnd_pd (L21-points)")
        call stops("no_shelf",  "needs a shelf instance", "no shelf instance")
        call stops("bad_method", "is not recognised", "an unknown method")
        call stops("ok",        "", "a valid case does not stop", ok_expected=.true.)
    end subroutine test_stops

    subroutine abort_case(name)
        character(len=*), intent(in) :: name
        type(domain_t) :: d
        type(experiment_list_t) :: list

        select case (name)
        case ("no_shelf")
            call make_domain(d, "SamePts", "TOY", "TOY", .false.)
            call put_ice(d)
            call experiments_create(list, cfg, d, clk, cnst)
        case ("bad_method")
            call make_domain(d, "BadMethod", "TOY", "TOY", .true.)
            call experiments_create(list, cfg, d, clk, cnst)
        case ("no_export")
            call make_domain(d, "SamePts", "TOY", "TOY", .true., no_grnd_pd=.true.)
            call experiments_create(list, cfg, d, clk, cnst)
        case ("ok")
            call make_domain(d, "SamePts", "TOY", "TOY", .true.)
            call experiments_create(list, cfg, d, clk, cnst)
        end select
    end subroutine abort_case

    subroutine stops(name, text, label, expect_text, ok_expected)
        character(len=*), intent(in)           :: name, text, label
        character(len=*), intent(in), optional :: expect_text
        logical,          intent(in), optional :: ok_expected

        character(len=512) :: exe
        integer :: es, es2
        logical :: want_ok, good

        want_ok = .false.
        if (present(ok_expected)) want_ok = ok_expected
        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" abort "//name//" > "//dir//name//".log 2>&1", exitstat=es)
        if (want_ok) then
            call check(label, es == 0)
            return
        end if
        good = es /= 0
        if (len_trim(text) > 0) then
            call execute_command_line("grep -q '"//text//"' "//dir//name//".log", exitstat=es2)
            good = good .and. es2 == 0
        end if
        if (present(expect_text)) then
            call execute_command_line("grep -q '"//expect_text//"' "//dir//name//".log", exitstat=es2)
            good = good .and. es2 == 0
        end if
        call check(label//" ("//dir//name//".log)", good)
    end subroutine stops

    subroutine write_par()
        integer :: u
        character(len=*), parameter :: win = "    tf_time_init = 10.0", win2 = "    tf_time_end = 50.0"

        open(newunit=u, file=par, status="replace", action="write")
        write(u,'(a)') '&kryos', '    domains = "SamePts"', '/'
        write(u,'(a)') '&SamePts',    '    experiments = "opt_pts"', '/'
        write(u,'(a)') '&SameBasin',  '    experiments = "opt_basin"', '/'
        write(u,'(a)') '&RemapPts',   '    experiments = "opt_pts"', '/'
        write(u,'(a)') '&RemapBasin', '    experiments = "opt_basin"', '/'
        write(u,'(a)') '&BasinList',  '    experiments = "opt_list"', '/'
        write(u,'(a)') '&BadMethod',  '    experiments = "opt_bad"', '/'
        write(u,'(a)') '&opt_pts',   '    kind = "opt_tf"', '    method = "L21-points"', win, win2, '/'
        write(u,'(a)') '&opt_basin', '    kind = "opt_tf"', '    method = "L21"', win, win2, '/'
        write(u,'(a)') '&opt_list',  '    kind = "opt_tf"', '    method = "L21"', '    tf_basins = 2 3', win, win2, '/'
        write(u,'(a)') '&opt_bad',   '    kind = "opt_tf"', '    method = "L22"', '/'
        close(u)
    end subroutine write_par

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=16) :: s
        write(s, '(i0)') i
    end function itoa

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

end program test_exp_opt_tf
