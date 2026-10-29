module test_gi_toy
    ! A toy ice instance that only declares its fields and implements the
    ! flotation binding (grounded where the ice outweighs the water it would
    ! displace, with the densities of the run's constants).
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_ice,       only: ice_t, ice_opt_t
    implicit none
    private

    real(dp), parameter, public :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp

    type, extends(ice_t), public :: toy_ice_t
    contains
        procedure :: init => toy_init
        procedure :: update => toy_update
        procedure :: restart_write => toy_rw
        procedure :: restart_read => toy_rr
        procedure :: flotation => toy_flotation
    end type toy_ice_t

contains

    subroutine toy_init(self, ctx)
        class(toy_ice_t),  intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = ice_opt_t(z_bed_ref=.true., footprint=.true., H_ice_pd=.true.)
        self%caps%has_flotation = .true.
    end subroutine toy_init

    subroutine toy_flotation(self, H, z_bed, z_sl, f_grnd)
        class(toy_ice_t), intent(in)  :: self
        real(dp),         intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),         intent(out) :: f_grnd(:,:)
        associate (unused => self%nx); end associate
        f_grnd = merge(1.0_dp, 0.0_dp, H*RHO_I/RHO_SW + (z_bed - z_sl) >= 0.0_dp)
    end subroutine toy_flotation

    subroutine toy_update(self, st)
        class(toy_ice_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_update

    subroutine toy_rw(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw

    subroutine toy_rr(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rr

end module test_gi_toy

program test_globe_ice
    ! kryos_globe_ice (phase4-plan 4, 5.1-5.3, 13.1) on mock ice instances
    ! (ice_adopt), two domains with disjoint footprints on small polar
    ! stereographic grids, a small Gaussian target grid and a coarse lon-lat
    ! background with a synthetic IceHistory file (years ago, 80 ka to 0):
    !   * init: the [globe_ice] group, the history group's guards (ref_time,
    !     out_of_range), the datasrc bound on the background grid;
    !   * add_domain: the true areas, footprint and coverage of each domain on
    !     the target grid; one ice grid per domain.
    ! Works in tmp/test_globe_ice/ below the repository root (make check runs there).

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use precision,       only: dp
    use coords,          only: grid_class, grid_init, grid_init_gaussian
    use coupler,         only: coupler_class, coupler_init
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr
    use phys_constants,  only: phys_const_class, phys_const_set
    use kryos_config,    only: config_t, config_init
    use kryos_component, only: init_ctx_t, census_t, constraints_t
    use kryos_ice,       only: ice_t, ice_slot_t, ice_adopt, ice_init
    use kryos_globe_ice, only: globe_ice_t, globe_ice_init, globe_ice_add_domain, globe_ice_end
    use test_gi_toy,     only: toy_ice_t, RHO_I, RHO_SW

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_globe_ice/", FLD = "tmp/test_globe_ice/maps"
    character(len=*), parameter :: errfile = "tmp/test_globe_ice/abort.err"
    character(len=*), parameter :: BGN = "TBG-4DEG", TGN = "TGS-L8"
    character(len=*), parameter :: HIST = "tmp/test_globe_ice/{grid_name}_IceHistory-TEST.nc"
    integer, parameter :: NXB = 90, NYB = 45, NT = 5
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    integer :: nfail = 0

    type :: scene_t                                   ! the grids, the coupler, the globe, the mock ice of two domains
        type(grid_class)       :: gbg, gtg, ga, gb
        type(coupler_class)    :: cpl
        type(config_t)         :: cfg
        type(phys_const_class) :: cnst
        type(globe_ice_t)      :: gi
        type(ice_slot_t)       :: ice_a, ice_b
        real(dp), allocatable  :: fp_a(:,:), fp_b(:,:), hpd_a(:,:), hpd_b(:,:), zb_a(:,:), zb_b(:,:)
    end type scene_t

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                 ! the case did not abort
    end if

    call execute_command_line("rm -rf "//DIR//" && mkdir -p "//FLD)
    call write_history()

    call test_coverage()
    call test_stops()

    call execute_command_line("rm -f "//errfile)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_globe_ice: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_globe_ice: OK"

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

    ! =========================================================================
    ! The data
    ! =========================================================================

    real(dp) function gc_deg(lon1, lat1, lon2, lat2)
        ! great-circle distance [deg]
        real(dp), intent(in) :: lon1, lat1, lon2, lat2
        real(dp) :: a, b, dl
        a = lat1*PI/180.0_dp; b = lat2*PI/180.0_dp; dl = (lon1 - lon2)*PI/180.0_dp
        gc_deg = acos(max(-1.0_dp, min(1.0_dp, sin(a)*sin(b) + cos(a)*cos(b)*cos(dl))))*180.0_dp/PI
    end function gc_deg

    real(dp) function cap(lon, lat, lon0, lat0, rad, h0)
        real(dp), intent(in) :: lon, lat, lon0, lat0, rad, h0
        cap = max(0.0_dp, h0*(1.0_dp - gc_deg(lon, lat, lon0, lat0)/rad))
    end function cap

    real(dp) function eur(lon, lat)
        ! the external ice sheet (outside both domains)
        real(dp), intent(in) :: lon, lat
        eur = cap(lon, lat, 20.0_dp, 62.0_dp, 14.0_dp, 1800.0_dp)
    end function eur

    real(dp) function h_now(i, j)
        ! ice on the background grid at 0 ka: a Greenland-like cap, an Antarctic cap, the external sheet
        integer, intent(in) :: i, j
        real(dp) :: lon, lat
        lon = 2.0_dp + 4.0_dp*(i - 1); lat = -88.0_dp + 4.0_dp*(j - 1)
        h_now = cap(lon, lat, 320.0_dp, 72.0_dp, 18.0_dp, 2500.0_dp) + cap(lon, lat, 0.0_dp, -80.0_dp, 25.0_dp, 3000.0_dp) &
                + eur(lon, lat)
    end function h_now

    real(dp) function h_alt(i, j)
        ! the background X of the clip scenario: no external sheet
        integer, intent(in) :: i, j
        h_alt = h_now(i, j) - eur(2.0_dp + 4.0_dp*(i - 1), -88.0_dp + 4.0_dp*(j - 1))
    end function h_alt

    real(dp) function zb_bg(i, j)
        real(dp), intent(in) :: i, j
        zb_bg = -1500.0_dp*sin(0.07_dp*i)*cos(0.11_dp*j) + 200.0_dp
    end function zb_bg

    subroutine write_history()
        ! GLOBAL-like file: ages 80000 ... 0 (years ago), z_bed (constant), H_ice
        ! (more ice in the past: H0 (1 + a/2), a = age/80 ka), H_ice_dec (less
        ! ice in the past: H0 (1 - 0.8 a)), H_alt (constant, without the external sheet).
        character(len=*), parameter :: f = DIR//BGN//"_IceHistory-TEST.nc"
        real(dp) :: ages(NT), zb(NXB,NYB,NT), h(NXB,NYB,NT), hd(NXB,NYB,NT), ha(NXB,NYB,NT), a
        integer :: i, j, k

        ages = [80000.0_dp, 60000.0_dp, 40000.0_dp, 20000.0_dp, 0.0_dp]
        do k = 1, NT
            a = ages(k)/80000.0_dp
            do j = 1, NYB
                do i = 1, NXB
                    zb(i,j,k) = zb_bg(real(i, dp), real(j, dp))
                    h(i,j,k)  = h_now(i, j)*(1.0_dp + 0.5_dp*a)
                    hd(i,j,k) = h_now(i, j)*(1.0_dp - 0.8_dp*a)
                    ha(i,j,k) = h_alt(i, j)
                end do
            end do
        end do
        call nc_create(f)
        call nc_write_dim(f, "lon", x=2.0_dp, dx=4.0_dp, nx=NXB, units="degrees_east")
        call nc_write_dim(f, "lat", x=-88.0_dp, dx=4.0_dp, nx=NYB, units="degrees_north")
        call nc_write_dim(f, "time", x=ages, units="years ago", unlimited=.false.)
        call nc_write(f, "z_bed", zb, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write(f, "H_ice", h, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write(f, "H_ice_dec", hd, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write(f, "H_alt", ha, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write_attr(f, "grid_name", BGN)
    end subroutine write_history

    subroutine write_par(file, scen)
        ! The par file of a scenario: def (X = the history's 0 ka, anomaly), clip (X without the
        ! external sheet, a history with less ice in the past), ref (ext = "ref"), and the faulty
        ! variants (bad_ref_time, bad_hold, tight).
        character(len=*), intent(in) :: file, scen
        character(len=16) :: ext
        character(len=32) :: hvar, xvar
        character(len=16) :: reft, oor, ctol
        integer :: u

        ext = "file"; hvar = "H_ice"; xvar = "H_ice"; reft = "0.0"; oor = "error"; ctol = "1.0e-9"
        select case (trim(scen))
        case ("clip");          hvar = "H_ice_dec"; xvar = "H_alt"
        case ("ref");           ext = "ref"
        case ("bad_ref_time");  reft = "-20000.0"
        case ("bad_hold");      oor = "hold"
        case ("tight");         ctol = "1.0e-30"
        end select

        open(newunit=u, file=file, status="replace", action="write")
        write(u,"(a)") "&kryos", '    map_fldr = "'//FLD//'"', "/"
        write(u,"(a)") "&globe_ice"
        write(u,"(a)") '    grid = "'//BGN//'"', '    z_bed_ref = "gi_X_zb"', '    H_ice_ref = "gi_X_H"'
        write(u,"(a)") '    ext = "'//trim(ext)//'"', '    ext_field = "gi_hist"', "    ext_tref = 0.0", &
                       "    cons_tol = "//trim(ctol), "/"
        write(u,"(a)") "&gi_X_zb", '    filename = "'//HIST//'"', '    name = "z_bed"', '    mode = "transient"', &
                       '    time_units = "years BP"', "/"
        write(u,"(a)") "&gi_X_H", '    filename = "'//HIST//'"', '    name = "'//trim(xvar)//'"', '    mode = "transient"', &
                       '    time_units = "years BP"', "/"
        write(u,"(a)") "&gi_hist", '    filename = "'//HIST//'"', '    name = "'//trim(hvar)//'"', '    mode = "transient"', &
                       '    time_units = "years BP"', '    anomaly_ref = "slice"', "    ref_time = "//trim(reft), &
                       '    out_of_range = "'//trim(oor)//'"', "/"
        close(u)
    end subroutine write_par

    ! =========================================================================
    ! The scene
    ! =========================================================================

    subroutine grids(sc, planet)
        type(scene_t),    intent(inout) :: sc
        character(len=*), intent(in)    :: planet        ! of the two ice grids
        character(len=16) :: tag
        call grid_init(sc%gbg, name=BGN, mtype="latlon", units="degrees", planet="Spherical Earth", &
                       x0=2.0_dp, dx=4.0_dp, nx=NXB, y0=-88.0_dp, dy=4.0_dp, ny=NYB)
        call grid_init_gaussian(sc%gtg, TGN, nlon=32, nlat=16)
        tag = "-SPH"
        if (trim(planet) /= "Spherical Earth") tag = "-"//trim(planet)
        call grid_init(sc%ga, name="TNH-100KM"//trim(tag), mtype="polar_stereographic", units="kilometers", &
                       planet=trim(planet), lon180=.false., x0=-1500.0_dp, dx=100.0_dp, nx=31, y0=0.0_dp, dy=100.0_dp, &
                       ny=31, lambda=-40.0_dp, phi=70.0_dp)
        call grid_init(sc%gb, name="TANT-100KM"//trim(tag), mtype="polar_stereographic", units="kilometers", &
                       planet=trim(planet), lon180=.false., x0=-2500.0_dp, dx=100.0_dp, nx=51, y0=-2500.0_dp, dy=100.0_dp, &
                       ny=51, lambda=0.0_dp, phi=-71.0_dp)
    end subroutine grids

    subroutine statics(g, xc, yc, rad, h0, fp, hpd, zb)
        ! Footprint (a disc of radius rad [km] about (xc, yc)), present-day cap and a bed that
        ! drops off with the distance: grounded inside, floating towards the rim.
        type(grid_class), intent(in) :: g
        real(dp),         intent(in) :: xc, yc, rad, h0
        real(dp), allocatable, intent(out) :: fp(:,:), hpd(:,:), zb(:,:)
        real(dp) :: r
        integer :: i, j
        allocate(fp(g%G%nx, g%G%ny), hpd(g%G%nx, g%G%ny), zb(g%G%nx, g%G%ny))
        do j = 1, g%G%ny
            do i = 1, g%G%nx
                r = sqrt((g%G%x(i) - xc)**2 + (g%G%y(j) - yc)**2)
                fp(i,j)  = merge(1.0_dp, 0.0_dp, r <= rad)
                hpd(i,j) = fp(i,j)*h0*(1.0_dp - r/rad)
                zb(i,j)  = 300.0_dp - 0.5_dp*r
            end do
        end do
    end subroutine statics

    subroutine make_ice(slot, g, name)
        type(ice_slot_t), intent(out) :: slot
        type(grid_class), intent(in)  :: g
        character(len=*), intent(in)  :: name
        class(ice_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        type(census_t)   :: census
        ctx%instance  = "ice_"//name
        ctx%component = "ice"
        ctx%backend   = "toy"
        ctx%domain    = name
        ctx%grid_name = g%name
        ctx%grid      = g
        allocate(toy_ice_t :: proto)
        call ice_adopt(slot, proto, constraints_t(), ctx, census)
        call ice_init(slot, ctx)
    end subroutine make_ice

    subroutine build(sc, scen, planet, order_ab, uniform)
        ! The scene of scenario `scen`: grids, coupler, globe_ice initialised, both domains
        ! registered (A, B or B, A), and the reference not yet assembled. uniform: the footprint of A is
        ! its whole grid, 1000 m of ice over a bed at -500 m (the direct host-term checks).
        type(scene_t),    intent(inout) :: sc
        character(len=*), intent(in)    :: scen
        character(len=*), intent(in), optional :: planet
        logical,          intent(in), optional :: order_ab, uniform
        character(len=*), parameter :: par = DIR//"par.nml"
        character(len=32) :: pl
        logical :: ab, uni

        pl = "Spherical Earth"
        if (present(planet)) pl = planet
        ab = .true.
        if (present(order_ab)) ab = order_ab
        uni = .false.
        if (present(uniform)) uni = uniform

        call grids(sc, trim(pl))
        call statics(sc%ga, 0.0_dp, 1500.0_dp, 900.0_dp, 2000.0_dp, sc%fp_a, sc%hpd_a, sc%zb_a)
        call statics(sc%gb, 0.0_dp, 0.0_dp, 1800.0_dp, 3000.0_dp, sc%fp_b, sc%hpd_b, sc%zb_b)
        if (uni) then
            sc%fp_a = 1.0_dp
            sc%hpd_a = 1000.0_dp
            sc%zb_a = -500.0_dp
        end if
        call phys_const_set(sc%cnst, label="Earth", rho_ice=RHO_I, rho_sw=RHO_SW, rho_w=1000.0_dp, g=9.81_dp, &
                            area_seasurf=3.618e8_dp)
        call write_par(par, scen)
        call config_init(sc%cfg, par)
        call coupler_init(sc%cpl, FLD)
        call globe_ice_init(sc%gi, sc%cfg, sc%cpl, sc%cnst, sc%gtg, grid_bg=sc%gbg, time_ref=1950.0_dp, &
                            time_span=[-60000.0_dp, 0.0_dp])
        call make_ice(sc%ice_a, sc%ga, "A")
        call make_ice(sc%ice_b, sc%gb, "B")
        if (ab) then
            call globe_ice_add_domain(sc%gi, sc%cpl, "A", sc%ga, sc%fp_a, sc%hpd_a, sc%zb_a)
            call globe_ice_add_domain(sc%gi, sc%cpl, "B", sc%gb, sc%fp_b, sc%hpd_b, sc%zb_b)
        else
            call globe_ice_add_domain(sc%gi, sc%cpl, "B", sc%gb, sc%fp_b, sc%hpd_b, sc%zb_b)
            call globe_ice_add_domain(sc%gi, sc%cpl, "A", sc%ga, sc%fp_a, sc%hpd_a, sc%zb_a)
        end if
    end subroutine build

    function areas(g) result(a)
        type(grid_class), intent(in) :: g
        real(dp), allocatable :: a(:,:)
        a = g%area*g%cs%xy_conv**2
    end function areas

    real(dp) function relerr(a, b)
        real(dp), intent(in) :: a, b
        relerr = abs(a - b)/max(abs(b), tiny(1.0_dp))
    end function relerr

    ! =========================================================================
    ! The tests
    ! =========================================================================

    subroutine test_coverage()
        type(scene_t) :: sc
        real(dp), allocatable :: at(:,:)
        real(dp) :: fp_area, fp_cov
        integer :: d

        write(*,"(a)") "footprint and coverage of each domain"
        call build(sc, "def")
        at = areas(sc%gtg)
        do d = 1, size(sc%gi%dom)
            call check("domain "//trim(sc%gi%dom(d)%name)//": footprint fraction <= coverage <= 1", &
                       all(sc%gi%dom(d)%fp_tgt <= sc%gi%dom(d)%f_tgt + 1.0e-12_dp) .and. &
                       all(sc%gi%dom(d)%f_tgt <= 1.0_dp + 1.0e-9_dp) .and. all(sc%gi%dom(d)%fp_bg <= sc%gi%dom(d)%f_bg + 1.0e-12_dp))
            fp_cov = sum(at*sc%gi%dom(d)%fp_tgt)
            if (trim(sc%gi%dom(d)%name) == "A") then
                fp_area = sum(areas(sc%ga)*sc%fp_a)
            else
                fp_area = sum(areas(sc%gb)*sc%fp_b)
            end if
            call check("domain "//trim(sc%gi%dom(d)%name)//": the footprint fraction covers the footprint area (1e-3)", &
                       relerr(fp_cov, fp_area) < 1.0e-3_dp)
        end do
        call check("the domains do not overlap on the target grid", &
                   maxval(sc%gi%dom(1)%f_tgt + sc%gi%dom(2)%f_tgt) <= 1.0_dp + 1.0e-9_dp .and. &
                   any(sc%gi%dom(1)%f_tgt > 0.0_dp .and. sc%gi%dom(1)%f_tgt < 1.0_dp))
        call globe_ice_end(sc%gi)
    end subroutine test_coverage

    ! =========================================================================
    ! Named stops
    ! =========================================================================

    subroutine test_stops()
        write(*,"(a)") "the checks as named stops"
        call expect_abort("bad_ref_time", "has ref_time")
        call expect_abort("bad_hold", "out_of_range")
        call expect_abort("same_grid", "are on one ice grid")
    end subroutine test_stops

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
        ! Child: build the situation and make the call that must stop.
        character(len=64) :: arg
        type(scene_t) :: sc

        call get_command_argument(1, arg)
        call execute_command_line("mkdir -p "//FLD)
        select case (trim(arg))
        case ("bad_ref_time")
            call build(sc, "bad_ref_time")
        case ("bad_hold")
            call build(sc, "bad_hold")
        case ("same_grid")
            call build(sc, "def")
            call globe_ice_add_domain(sc%gi, sc%cpl, "A3", sc%ga, sc%fp_a, sc%hpd_a, sc%zb_a)
        end select
    end subroutine abort_case

end program test_globe_ice
