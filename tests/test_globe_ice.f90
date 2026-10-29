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
    !   * the background datasrc: the 0 ka slice, the history as an anomaly
    !     against its own 0 ka, hold after ext_tref, stop before the first
    !     record, clip at zero thickness, ext = "ref";
    !   * the host terms against a direct computation (uniform fields, partial
    !     cells); footprint and coverage w; mass conservation of the merge
    !     against plain sums;
    !   * the reference: bit identical when assembled twice, the self-test
    !     (a gather at the reference state returns no anomaly), the spliced
    !     reference;
    !   * domain-order invariance (registration and add order swapped);
    !   * the checks of 4.5 as stops (the program re-runs itself with an argument
    !     per case, as test_ice_contract) and the diagnostics.
    ! Works in tmp/test_globe_ice/ below the repository root (make check runs there).

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use precision,       only: dp
    use coords,          only: grid_class, grid_init, grid_init_gaussian
    use coupler,         only: coupler_class, coupler_init
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr
    use phys_constants,  only: phys_const_class, phys_const_set
    use kryos_config,    only: config_t, config_init
    use kryos_component, only: init_ctx_t, census_t, constraints_t
    use kryos_ice,       only: ice_t, ice_slot_t, ice_adopt, ice_init, ice_flotation
    use kryos_globe_ice, only: globe_ice_t, globe_ice_init, globe_ice_add_domain, globe_ice_end, &
                               globe_ice_reference_begin, globe_ice_reference_add, globe_ice_reference_end, &
                               globe_ice_gather_begin, globe_ice_gather_add, globe_ice_gather_end
    use test_gi_toy,     only: toy_ice_t, RHO_I, RHO_SW

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_globe_ice/", FLD = "tmp/test_globe_ice/maps"
    character(len=*), parameter :: errfile = "tmp/test_globe_ice/abort.err"
    character(len=*), parameter :: BGN = "TBG-4DEG", TGN = "TGS-L8"
    character(len=*), parameter :: HIST = "tmp/test_globe_ice/{grid_name}_IceHistory-TEST.nc"
    integer, parameter :: NXB = 90, NYB = 45, NT = 5
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: SP_TOL = 2.0e-6_dp       ! the datasrc reads in single precision (wp)
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

    call test_background()
    call test_host_terms()
    call test_coverage_and_conservation()
    call test_reference()
    call test_order()
    call test_no_domains()
    call test_wgs84()
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

    subroutine reference(sc, order_ab)
        type(scene_t), intent(inout) :: sc
        logical,       intent(in), optional :: order_ab
        logical :: ab
        ab = .true.
        if (present(order_ab)) ab = order_ab
        call globe_ice_reference_begin(sc%gi, sc%cpl)
        if (ab) then
            call globe_ice_reference_add(sc%gi, sc%cpl, "A", sc%ice_a)
            call globe_ice_reference_add(sc%gi, sc%cpl, "B", sc%ice_b)
        else
            call globe_ice_reference_add(sc%gi, sc%cpl, "B", sc%ice_b)
            call globe_ice_reference_add(sc%gi, sc%cpl, "A", sc%ice_a)
        end if
        call globe_ice_reference_end(sc%gi, sc%cpl)
    end subroutine reference

    subroutine gather(sc, t, dHa, order_ab)
        ! A gather at time t with the reference state of both domains (H = H_ice_pd, bed = z_bed_ref,
        ! sea level 0; f_grnd = the ice's own rule inside the footprint, 0.5 outside, where it is
        ! ignored); dHa: ice added to the footprint of A.
        type(scene_t), intent(inout) :: sc
        real(dp),      intent(in)    :: t
        real(dp),      intent(in), optional :: dHa
        logical,       intent(in), optional :: order_ab
        logical :: ab
        ab = .true.
        if (present(order_ab)) ab = order_ab
        call globe_ice_gather_begin(sc%gi, sc%cpl, t)
        if (ab) then
            call add_one(sc, "A", dHa)
            call add_one(sc, "B")
        else
            call add_one(sc, "B")
            call add_one(sc, "A", dHa)
        end if
        call globe_ice_gather_end(sc%gi, sc%cpl)
    end subroutine gather

    subroutine add_one(sc, name, dH)
        type(scene_t), intent(inout) :: sc
        character(len=*), intent(in) :: name
        real(dp), intent(in), optional :: dH
        real(dp), allocatable :: H(:,:), fY(:,:), zs(:,:)
        if (trim(name) == "A") then
            H = sc%hpd_a
            if (present(dH)) H = H + dH*sc%fp_a
            allocate(fY(size(H,1), size(H,2)), zs(size(H,1), size(H,2)))
            zs = 0.0_dp
            call ice_flotation(sc%ice_a, H, sc%zb_a, zs, fY)
            fY = merge(fY, 0.5_dp, sc%fp_a > 0.5_dp)
            call globe_ice_gather_add(sc%gi, sc%cpl, "A", sc%ice_a, H, fY, sc%zb_a, zs, sc%fp_a)
        else
            H = sc%hpd_b
            allocate(fY(size(H,1), size(H,2)), zs(size(H,1), size(H,2)))
            zs = 0.0_dp
            call ice_flotation(sc%ice_b, H, sc%zb_b, zs, fY)
            fY = merge(fY, 0.5_dp, sc%fp_b > 0.5_dp)
            call globe_ice_gather_add(sc%gi, sc%cpl, "B", sc%ice_b, H, fY, sc%zb_b, zs, sc%fp_b)
        end if
    end subroutine add_one

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

    subroutine test_background()
        type(scene_t) :: sc
        real(dp), allocatable :: x(:,:), ab(:,:), at(:,:)
        real(dp) :: tru, ext_t, vbg, vtg
        integer :: i, j

        write(*,"(a)") "background datasrc: the 0 ka slice, the history as an anomaly"
        call build(sc, "def")
        call reference(sc)
        allocate(x(NXB,NYB))
        do j = 1, NYB
            do i = 1, NXB
                x(i,j) = h_now(i, j)
            end do
        end do
        call check("X ice is the 0 ka slice (last record, years ago -> model time 0; read in single precision)", &
                   maxval(abs(sc%gi%H_X - x)) <= SP_TOL*maxval(x))
        do j = 1, NYB
            do i = 1, NXB
                x(i,j) = zb_bg(real(i, dp), real(j, dp))
            end do
        end do
        call check("X bed is the 0 ka slice", maxval(abs(sc%gi%z_bed_X - x)) <= SP_TOL*maxval(abs(x)))

        ! the history at 40 ka and between two records, with X = the history's own 0 ka: H_ext(t) = H_hist(t)
        call gather(sc, -40000.0_dp)
        call check("t = -40 ka: H_ext = the history at 40 ka (X = its 0 ka)", maxdiff_h(sc%gi%H_ext, 1.0_dp + 0.5_dp*0.5_dp) < SP_TOL)
        call gather(sc, -30000.0_dp)
        call check("t = -30 ka (between records, linear): H_ext = H0 (1 + a/2)", &
                   maxdiff_h(sc%gi%H_ext, 1.0_dp + 0.5_dp*30000.0_dp/80000.0_dp) < SP_TOL)
        call gather(sc, 0.0_dp)
        call check("t = 0: the anomaly is zero, H_ext = X exactly", all(sc%gi%H_ext == sc%gi%H_X))
        call gather(sc, 500.0_dp)
        call check("t = +500 yr: held at its present day (H_ext = X)", all(sc%gi%H_ext == sc%gi%H_X))

        ! the target grid carries the same ice (an amount per target area)
        ab = areas(sc%gbg); at = areas(sc%gtg)
        call gather(sc, -40000.0_dp)
        vbg = sum(ab*sc%gi%H_ext); vtg = sum(at*sc%gi%H_ice_ext)
        call check("H_ice_ext conserves the ice of the background grid (plain sums, 1e-12)", relerr(vtg, vbg) < 1.0e-12_dp)
        vbg = sum(ab*sc%gi%H_X); vtg = sum(at*sc%gi%H_ice_ref)
        call check("H_ice_ref (X on the target grid) conserves the ice of X (1e-12)", relerr(vtg, vbg) < 1.0e-12_dp)
        call globe_ice_end(sc%gi)

        ! ext = "ref": the history is not read
        call build(sc, "ref")
        call reference(sc)
        call gather(sc, -40000.0_dp)
        call check("ext = ref: H_ext = X at any time", all(sc%gi%H_ext == sc%gi%H_X) .and. sc%gi%diag%clip_Gt == 0.0_dp)
        call globe_ice_end(sc%gi)

        ! clip: X without the external sheet, a history with less ice in the past: the anomaly is negative
        ! where X has no ice, and the sum is clipped at zero
        call build(sc, "clip")
        call reference(sc)
        call gather(sc, -40000.0_dp)
        ab = areas(sc%gbg)
        ! direct: raw = H_alt + (H_dec(t) - H_dec(0)) = H_alt - 0.8 a H0; the clipped part is max(0, 0.8 a H0 - H_alt)
        tru = 0.0_dp
        do j = 1, NYB
            do i = 1, NXB
                tru = tru + ab(i,j)*max(0.0_dp, 0.8_dp*0.5_dp*h_now(i, j) - h_alt(i, j))
            end do
        end do
        ext_t = tru*RHO_I*1.0e-12_dp
        call check("clip: H_ext >= 0 everywhere, some cells clipped", minval(sc%gi%H_ext) >= 0.0_dp .and. ext_t > 0.0_dp)
        call check("clip: the clipped mass [Gt] of the diagnostics is the direct one (single-precision data)", &
                   relerr(sc%gi%diag%clip_Gt, ext_t) < 10.0_dp*SP_TOL)
        call globe_ice_end(sc%gi)
    end subroutine test_background

    real(dp) function maxdiff_h(H, fac)
        ! max |H - fac H0| / max(H0) on the background grid
        real(dp), intent(in) :: H(:,:), fac
        real(dp) :: m
        integer :: i, j
        maxdiff_h = 0.0_dp
        m = 0.0_dp
        do j = 1, NYB
            do i = 1, NXB
                m = max(m, h_now(i, j))
            end do
        end do
        do j = 1, NYB
            do i = 1, NXB
                maxdiff_h = max(maxdiff_h, abs(H(i,j) - fac*h_now(i, j))/m)
            end do
        end do
    end function maxdiff_h

    subroutine test_host_terms()
        ! A: the footprint is the whole grid, 1000 m of ice over a bed at -500 m (grounded), so the host
        ! terms of a cell are its weight w times a constant: partial cells included. The weights of the
        ! intensive terms are the map's coverage exactly; the amounts (per true area) differ from w times
        ! the constant by the ratio of the true cell area to the map's polygon area (tolerance below, 5e-4: the measured 8e-5).
        type(scene_t) :: sc
        real(dp), allocatable :: H(:,:), fY(:,:), zs(:,:), z(:,:), w(:,:), a(:,:)
        real(dp), parameter :: AMT_TOL = 5.0e-4_dp
        logical :: partial

        write(*,"(a)") "host terms against a direct computation (uniform A, partial cells)"
        call build(sc, "def", uniform=.true.)
        call reference(sc)
        w = sc%gi%host_w
        call mask_a_only(sc, a)
        partial = any(w*a > 0.0_dp .and. w*a < 0.95_dp)
        call check("coverage: partial cells at the rim and cells fully covered", partial .and. maxval(w*a) > 0.999_dp)
        call check("reference: host_Ig_ref = 1000 w on the cells of A (amounts: tolerance 5e-4)", &
                   maxval(abs((sc%gi%host_Ig_ref - 1000.0_dp*w)*a)) <= AMT_TOL*1000.0_dp)
        write(*,"(a,es10.3)") "  note  max |host_Ig_ref - 1000 w| / 1000 on the cells of A: ", &
            maxval(abs((sc%gi%host_Ig_ref - 1000.0_dp*w)*a))/1000.0_dp
        call check("reference: host_C_ref = 0 on the cells of A (grounded)", maxval(abs(sc%gi%host_C_ref*a)) <= 1.0e-12_dp)

        ! the reference state: no anomaly
        call gather(sc, 0.0_dp)
        call check("self-test: a gather at the reference state returns dI_g = 0", &
                   maxval(abs(sc%gi%host_Ig - sc%gi%host_Ig_ref)) <= 1.0e-12_dp*maxval(sc%gi%host_Ig_ref))
        call check("self-test: dC = 0 and F = 0", maxval(abs(sc%gi%host_C - sc%gi%host_C_ref)) <= 1.0e-12_dp .and. &
                   maxval(abs(sc%gi%host_F)) <= 1.0e-9_dp)

        ! a step state: 300 m of ice over the bed at -500 m floats; F = w zeta0 (C - C_ref) = -500 w
        allocate(H(sc%ga%G%nx, sc%ga%G%ny), fY(sc%ga%G%nx, sc%ga%G%ny), zs(sc%ga%G%nx, sc%ga%G%ny), z(sc%ga%G%nx, sc%ga%G%ny))
        H = 300.0_dp; zs = 0.0_dp; z = -500.0_dp
        call ice_flotation(sc%ice_a, H, z, zs, fY)
        call check("(the toy rule: 300 m over -500 m floats)", all(fY == 0.0_dp))
        call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
        call globe_ice_gather_add(sc%gi, sc%cpl, "A", sc%ice_a, H, fY, z, zs, sc%fp_a)
        call add_one(sc, "B")
        call globe_ice_gather_end(sc%gi, sc%cpl)
        call check("floating A: host_Ig = 0 on the cells of A", maxval(abs(sc%gi%host_Ig*a)) <= 1.0e-12_dp)
        call check("floating A: host_C = w (exactly the coverage)", maxval(abs((sc%gi%host_C - w)*a)) <= 1.0e-12_dp)
        call check("floating A: host_F = zeta0 (C - C_ref) w = -500 w (amounts: 5e-4)", &
                   maxval(abs((sc%gi%host_F + 500.0_dp*w)*a)) <= AMT_TOL*500.0_dp)
        call globe_ice_end(sc%gi)
    end subroutine test_host_terms

    subroutine mask_a_only(sc, a)
        ! 1 on the target cells that A covers and B does not.
        type(scene_t), intent(in) :: sc
        real(dp), allocatable, intent(out) :: a(:,:)
        integer :: i
        allocate(a(sc%gi%nxg, sc%gi%nyg))
        a = 0.0_dp
        do i = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(i)%name) == "A") a = merge(1.0_dp, 0.0_dp, sc%gi%dom(i)%f_tgt > 0.0_dp)
        end do
        do i = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(i)%name) == "B") a = a*merge(0.0_dp, 1.0_dp, sc%gi%dom(i)%f_tgt > 0.0_dp)
        end do
    end subroutine mask_a_only

    subroutine test_coverage_and_conservation()
        type(scene_t) :: sc
        real(dp), allocatable :: at(:,:), ext_d(:,:)
        real(dp) :: v_src, v_tgt, fp_area, fp_cov, dva, dvb, ext_dv
        integer :: d, ia, ib

        write(*,"(a)") "footprint and coverage w; conservation of the merge; diagnostics"
        call build(sc, "def")
        call reference(sc)
        at = areas(sc%gtg)
        call check("host_w in [0, 1], some cells covered, some not", minval(sc%gi%host_w) >= 0.0_dp .and. &
                   maxval(sc%gi%host_w) <= 1.0_dp + 1.0e-9_dp .and. any(sc%gi%host_w > 0.0_dp) .and. any(sc%gi%host_w == 0.0_dp))
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
        call check("w is the sum of the domains' coverage (1e-12)", &
                   maxval(abs(sc%gi%host_w - sc%gi%dom(1)%f_tgt - sc%gi%dom(2)%f_tgt)) < 1.0e-12_dp)

        ! mass: the merged grounded ice against the plain sum over the source grids, with more ice in A
        ! (at t = 0 the external ice is X: B has no anomaly)
        call gather(sc, 0.0_dp, dHa=150.0_dp)
        v_tgt = sum(at*sc%gi%host_Ig)
        v_src = 0.0_dp
        do d = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(d)%name) == "A") v_src = v_src + sum(areas(sc%ga)*sc%gi%dom(d)%I_last)
            if (trim(sc%gi%dom(d)%name) == "B") v_src = v_src + sum(areas(sc%gb)*sc%gi%dom(d)%I_last)
        end do
        call check("merge: sum over the target grid = sum over the source grids (plain sums, 1e-12)", relerr(v_tgt, v_src) < 1.0e-12_dp)
        write(*,"(a,2es10.3)") "  note  closure of the merges (grounded ice, external ice): ", sc%gi%diag%cons_I, sc%gi%diag%cons_ext
        call check("merge: the closure checks of the gather are at round-off", sc%gi%diag%cons_I < 1.0e-12_dp .and. &
                   sc%gi%diag%cons_ext < 1.0e-12_dp)
        call check("host fields finite, 0 <= host_C <= host_w, host_Ig >= 0", all(sc%gi%host_C >= -1.0e-9_dp) .and. &
                   all(sc%gi%host_C <= sc%gi%host_w + 1.0e-9_dp) .and. all(sc%gi%host_Ig >= 0.0_dp))

        ! diagnostics of 4.5 against plain sums
        ia = index_of(sc, "A"); ib = index_of(sc, "B")
        dva = 0.0_dp; dvb = 0.0_dp
        do d = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(d)%name) == "A") dva = sum(areas(sc%ga)*(sc%gi%dom(d)%I_last - sc%gi%dom(d)%I_ref))
            if (trim(sc%gi%dom(d)%name) == "B") dvb = sum(areas(sc%gb)*(sc%gi%dom(d)%I_last - sc%gi%dom(d)%I_ref))
        end do
        call check("diag: the grounded-ice change of A against the reference is the plain-sum one (1e-9)", &
                   dva > 0.0_dp .and. relerr(sc%gi%diag%dV_g(ia), dva) < 1.0e-9_dp)
        call check("diag: Gt and ESL of A follow from the volume (rho_ice, rho_sw, sea surface)", &
                   relerr(sc%gi%diag%dM(ia), sc%gi%diag%dV_g(ia)*RHO_I*1.0e-12_dp) < 1.0e-12_dp .and. &
                   relerr(sc%gi%diag%esl(ia), -(RHO_I/RHO_SW)*sc%gi%diag%dV_g(ia)/sc%gi%area_ocn) < 1.0e-12_dp .and. &
                   sc%gi%diag%esl(ia) < 0.0_dp)
        call check("diag: B has no anomaly at t = 0", abs(dvb) <= 1.0e-9_dp*abs(sc%gi%diag%V_g(ib)) .and. &
                   abs(sc%gi%diag%dV_g(ib)) <= 1.0e-9_dp*abs(sc%gi%diag%V_g(ib)))
        call check("diag: the external part is zero at t = 0", sc%gi%diag%ext_dV == 0.0_dp .and. sc%gi%diag%clip_Gt == 0.0_dp)
        call gather(sc, -40000.0_dp, dHa=150.0_dp)
        ext_d = max(0.0_dp, 1.0_dp - sc%gi%host_w)*(sc%gi%H_ice_ext - sc%gi%H_ice_ref)
        ext_dv = sum(at*ext_d)
        call check("diag: the external part is the (1 - w) change of the history's ice (plain sum, 1e-9; more ice: ESL < 0)", &
                   ext_dv > 0.0_dp .and. relerr(sc%gi%diag%ext_dV, ext_dv) < 1.0e-9_dp .and. sc%gi%diag%ext_esl < 0.0_dp)
        call check("diag: the ocean area of the reference is a fraction of the sphere", &
                   sc%gi%diag%a_ocn_ref > 0.2_dp*sum(at) .and. sc%gi%diag%a_ocn_ref < sum(at))
        call globe_ice_end(sc%gi)
    end subroutine test_coverage_and_conservation

    integer function index_of(sc, name)
        ! position of a domain in the diagnostics (merge order)
        type(scene_t), intent(in) :: sc
        character(len=*), intent(in) :: name
        integer :: k
        index_of = 0
        do k = 1, sc%gi%diag%nd
            if (trim(sc%gi%diag%name(k)) == trim(name)) index_of = k
        end do
    end function index_of

    subroutine test_reference()
        type(scene_t) :: sc
        real(dp), allocatable :: w0(:,:), ig0(:,:), c0(:,:), zs0(:,:), hs0(:,:), cov0(:,:), zr0(:,:), hr0(:,:), ab(:,:)
        real(dp), allocatable :: expect(:,:)
        logical, allocatable :: onlya(:,:)
        real(dp) :: lhs, rhs
        integer :: d

        write(*,"(a)") "reference: assembly twice, the spliced reference"
        call build(sc, "def")
        call reference(sc)
        w0 = sc%gi%host_w; ig0 = sc%gi%host_Ig_ref; c0 = sc%gi%host_C_ref
        zs0 = sc%gi%z_bed_spl; hs0 = sc%gi%H_spl; cov0 = sc%gi%cov_bg; zr0 = sc%gi%z_bed_ref; hr0 = sc%gi%H_ice_ref
        call reference(sc)
        call check("the reference assembled twice is bit identical (host terms, X on the target grid)", &
                   all(w0 == sc%gi%host_w) .and. all(ig0 == sc%gi%host_Ig_ref) .and. all(c0 == sc%gi%host_C_ref) .and. &
                   all(zr0 == sc%gi%z_bed_ref) .and. all(hr0 == sc%gi%H_ice_ref))
        call check("... and the spliced reference", all(zs0 == sc%gi%z_bed_spl) .and. all(hs0 == sc%gi%H_spl) .and. &
                   all(cov0 == sc%gi%cov_bg))

        call check("splice: where no domain covers, the reference is X exactly", &
                   all(merge(.true., sc%gi%z_bed_spl == sc%gi%z_bed_X .and. sc%gi%H_spl == sc%gi%H_X, sc%gi%cov_bg > 0.0_dp)) .and. &
                   any(sc%gi%cov_bg == 0.0_dp) .and. any(sc%gi%cov_bg > 0.0_dp))
        ab = areas(sc%gbg)
        lhs = sum(ab*sc%gi%H_spl)
        rhs = sum(ab*(1.0_dp - min(sc%gi%cov_bg, 1.0_dp))*sc%gi%H_X)
        do d = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(d)%name) == "A") rhs = rhs + sum(areas(sc%ga)*sc%gi%dom(d)%H_ref)
            if (trim(sc%gi%dom(d)%name) == "B") rhs = rhs + sum(areas(sc%gb)*sc%gi%dom(d)%H_ref)
        end do
        call check("splice: the ice is the background outside the domains plus the reference ice of the domains (plain sums, 1e-9)", &
                   relerr(lhs, rhs) < 1.0e-9_dp)
        call check("one reference definition: H_ref = H_pd in the footprint, X on the ice grid elsewhere", &
                   all(merge(sc%gi%dom(1)%H_ref == sc%gi%dom(1)%H_pd, .true., sc%gi%dom(1)%footprint > 0.5_dp)) .and. &
                   all(merge(sc%gi%dom(1)%H_ref == sc%gi%dom(1)%ext_ref, .true., sc%gi%dom(1)%footprint < 0.5_dp)))
        call check("splice: coverage in [0, 1]", minval(sc%gi%cov_bg) >= 0.0_dp .and. maxval(sc%gi%cov_bg) <= 1.0_dp + 1.0e-9_dp)
        call globe_ice_end(sc%gi)

        ! a constant reference bed in A: the intensive merge is (1 - cov) z_X + cov c on the cells only A covers
        call build(sc, "def", uniform=.true.)
        call reference(sc)
        allocate(expect(NXB,NYB), onlya(NXB,NYB))
        onlya = sc%gi%cov_bg > 0.0_dp .and. .not. covered_by_b(sc)
        expect = (1.0_dp - min(sc%gi%cov_bg, 1.0_dp))*sc%gi%z_bed_X - 500.0_dp*sc%gi%cov_bg/max(sc%gi%cov_bg, 1.0_dp)
        call check("splice: a domain with a constant bed gives (1 - cov) z_X + cov c where only A covers (1e-9 m)", &
                   count(onlya) > 0 .and. maxval(abs(merge(sc%gi%z_bed_spl - expect, 0.0_dp, onlya))) < 1.0e-9_dp)
        call globe_ice_end(sc%gi)
    end subroutine test_reference

    function covered_by_b(sc) result(m)
        type(scene_t), intent(in) :: sc
        logical, allocatable :: m(:,:)
        integer :: i
        allocate(m(NXB,NYB))
        m = .false.
        do i = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(i)%name) == "B") m = sc%gi%dom(i)%f_bg > 0.0_dp
        end do
    end function covered_by_b

    subroutine test_order()
        type(scene_t) :: s1, s2
        logical :: same

        write(*,"(a)") "domain-order invariance"
        call build(s1, "def", order_ab=.true.)
        call build(s2, "def", order_ab=.false.)
        call reference(s1, order_ab=.true.)
        call reference(s2, order_ab=.false.)
        same = all(s1%gi%host_w == s2%gi%host_w) .and. all(s1%gi%host_Ig_ref == s2%gi%host_Ig_ref) .and. &
               all(s1%gi%host_C_ref == s2%gi%host_C_ref) .and. all(s1%gi%z_bed_spl == s2%gi%z_bed_spl) .and. &
               all(s1%gi%H_spl == s2%gi%H_spl) .and. all(s1%gi%cov_bg == s2%gi%cov_bg)
        call check("registration and add order swapped: the reference is bit identical", same)
        call gather(s1, -40000.0_dp, dHa=120.0_dp, order_ab=.true.)
        call gather(s2, -40000.0_dp, dHa=120.0_dp, order_ab=.false.)
        same = all(s1%gi%host_Ig == s2%gi%host_Ig) .and. all(s1%gi%host_C == s2%gi%host_C) .and. &
               all(s1%gi%host_F == s2%gi%host_F) .and. all(s1%gi%H_ice_ext == s2%gi%H_ice_ext) .and. &
               all(s1%gi%diag%dV_g == s2%gi%diag%dV_g) .and. s1%gi%diag%ext_dV == s2%gi%diag%ext_dV
        call check("... and the gather (host terms, external ice, diagnostics)", same)
        call globe_ice_end(s1%gi)
        call globe_ice_end(s2%gi)
    end subroutine test_order

    subroutine test_no_domains()
        ! No domain: the host terms are zero, the reference is X, the ocean area is the background's.
        type(scene_t) :: sc
        real(dp), allocatable :: ab(:,:), cbg(:,:)
        real(dp) :: a, fg
        integer :: i, j

        write(*,"(a)") "no domain registered"
        call build_bare(sc)
        call globe_ice_reference_begin(sc%gi, sc%cpl)
        call globe_ice_reference_end(sc%gi, sc%cpl)
        call check("host_w, host_Ig_ref and host_C_ref are zero", all(sc%gi%host_w == 0.0_dp) .and. &
                   all(sc%gi%host_Ig_ref == 0.0_dp) .and. all(sc%gi%host_C_ref == 0.0_dp))
        call check("the spliced reference is X", all(sc%gi%z_bed_spl == sc%gi%z_bed_X) .and. all(sc%gi%H_spl == sc%gi%H_X))
        call globe_ice_gather_begin(sc%gi, sc%cpl, -20000.0_dp)
        call globe_ice_gather_end(sc%gi, sc%cpl)
        call check("a gather without domains: host terms zero, the external ice is the whole ice", &
                   all(sc%gi%host_Ig == 0.0_dp) .and. all(sc%gi%host_C == 0.0_dp) .and. maxval(sc%gi%H_ice_ext) > 0.0_dp)
        ab = areas(sc%gbg)
        a = 0.0_dp
        do j = 1, NYB
            do i = 1, NXB
                fg = merge(1.0_dp, 0.0_dp, RHO_I*sc%gi%H_X(i,j)/RHO_SW + sc%gi%z_bed_X(i,j) >= 0.0_dp)
                a = a + ab(i,j)*(1.0_dp - fg)
            end do
        end do
        call check("the ocean area of the reference is the background's (1e-9)", relerr(sc%gi%diag%a_ocn_ref, a) < 1.0e-9_dp)
        if (allocated(cbg)) deallocate(cbg)
        call globe_ice_end(sc%gi)
    end subroutine test_no_domains

    subroutine test_wgs84()
        ! Ice grids on the ellipsoid with true areas (section 4.4) and a spherical target: the merge still
        ! conserves the ice to round-off; w (the map's coverage) and the amounts per area differ at the rim
        ! (the ellipsoid cell is larger than the sphere's near the pole): reported, bounded at 2 %.
        type(scene_t) :: sc
        real(dp), allocatable :: w(:,:), a(:,:), at(:,:)
        real(dp) :: v_src, v_tgt, dmax
        integer :: d

        write(*,"(a)") "ellipsoidal ice grids, spherical target"
        call build(sc, "def", planet="WGS84", uniform=.true.)
        call reference(sc)
        call gather(sc, 0.0_dp, dHa=100.0_dp)
        at = areas(sc%gtg)
        v_tgt = sum(at*sc%gi%host_Ig)
        v_src = 0.0_dp
        do d = 1, size(sc%gi%dom)
            if (trim(sc%gi%dom(d)%name) == "A") v_src = v_src + sum(areas(sc%ga)*sc%gi%dom(d)%I_last)
            if (trim(sc%gi%dom(d)%name) == "B") v_src = v_src + sum(areas(sc%gb)*sc%gi%dom(d)%I_last)
        end do
        call check("WGS84: the merge conserves the grounded ice (plain sums, 1e-12)", relerr(v_tgt, v_src) < 1.0e-12_dp)
        call mask_a_only(sc, a)
        w = sc%gi%host_w
        dmax = maxval(abs((sc%gi%host_Ig_ref - 1000.0_dp*w)*a))/1000.0_dp
        write(*,"(a,es10.3)") "  note  max |host_Ig_ref - 1000 w| / 1000 on the cells of A, ellipsoidal grids: ", dmax
        call check("WGS84: amounts and w differ by less than 2 % of the constant at the rim", dmax < 2.0e-2_dp)
        call globe_ice_end(sc%gi)
    end subroutine test_wgs84

    subroutine build_bare(sc)
        type(scene_t), intent(inout) :: sc
        character(len=*), parameter :: par = DIR//"par.nml"
        call grids(sc, "Spherical Earth")
        call phys_const_set(sc%cnst, label="Earth", rho_ice=RHO_I, rho_sw=RHO_SW, rho_w=1000.0_dp, g=9.81_dp, &
                            area_seasurf=3.618e8_dp)
        call write_par(par, "def")
        call config_init(sc%cfg, par)
        call coupler_init(sc%cpl, FLD)
        call globe_ice_init(sc%gi, sc%cfg, sc%cpl, sc%cnst, sc%gtg, grid_bg=sc%gbg)
    end subroutine build_bare

    ! =========================================================================
    ! Named stops
    ! =========================================================================

    subroutine test_stops()
        write(*,"(a)") "the checks as named stops"
        call expect_abort("before_first", "time outside the axis of the file")
        call expect_abort("bad_ref_time", "has ref_time")
        call expect_abort("bad_hold", "out_of_range")
        call expect_abort("overlap", "domains overlap")
        call expect_abort("fp_changed", "footprint of domain `A` changed")
        call expect_abort("nonfinite", "is not finite")
        call expect_abort("negative_h", "negative ice thickness")
        call expect_abort("not_added", "was not added")
        call expect_abort("unknown", "is not registered")
        call expect_abort("tight", "does not conserve the ice")
        call expect_abort("late_domain", "registered after the reference")
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
        real(dp), allocatable :: H(:,:), fY(:,:), zs(:,:), fp2(:,:)
        type(grid_class) :: g2

        call get_command_argument(1, arg)
        call execute_command_line("mkdir -p "//FLD)
        select case (trim(arg))
        case ("before_first")
            call build(sc, "def")
            call reference(sc)
            call gather(sc, -90000.0_dp)
        case ("bad_ref_time")
            call build(sc, "bad_ref_time")
        case ("bad_hold")
            call build(sc, "bad_hold")
        case ("overlap")
            call build(sc, "def")
            ! a third domain on a grid shifted by 500 km over A
            call grid_init(g2, name="TNH2-100KM-SPH", mtype="polar_stereographic", units="kilometers", planet="Spherical Earth", &
                           lon180=.false., x0=-1000.0_dp, dx=100.0_dp, nx=31, y0=500.0_dp, dy=100.0_dp, ny=31, &
                           lambda=-40.0_dp, phi=70.0_dp)
            allocate(fp2(31,31))
            fp2 = 1.0_dp
            call globe_ice_add_domain(sc%gi, sc%cpl, "A2", g2, fp2, 0.0_dp*fp2, 0.0_dp*fp2 - 500.0_dp)
            call globe_ice_reference_begin(sc%gi, sc%cpl)
            call globe_ice_reference_add(sc%gi, sc%cpl, "A", sc%ice_a)
            call globe_ice_reference_add(sc%gi, sc%cpl, "B", sc%ice_b)
            call globe_ice_reference_add(sc%gi, sc%cpl, "A2", sc%ice_a)     ! (flotation of the grid of A: shapes differ)
            call globe_ice_reference_end(sc%gi, sc%cpl)
        case ("fp_changed")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
            sc%fp_a(16,16) = 1.0_dp - sc%fp_a(16,16)
            call add_one(sc, "A")
        case ("nonfinite")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
            H = sc%hpd_a
            H(3,3) = ieee_value(1.0_dp, ieee_quiet_nan)
            allocate(fY(31,31), zs(31,31))
            fY = 1.0_dp; zs = 0.0_dp
            call globe_ice_gather_add(sc%gi, sc%cpl, "A", sc%ice_a, H, fY, sc%zb_a, zs, sc%fp_a)
        case ("negative_h")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
            H = sc%hpd_a
            H(3,3) = -5.0_dp
            allocate(fY(31,31), zs(31,31))
            fY = 1.0_dp; zs = 0.0_dp
            call globe_ice_gather_add(sc%gi, sc%cpl, "A", sc%ice_a, H, fY, sc%zb_a, zs, sc%fp_a)
        case ("not_added")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
            call add_one(sc, "A")
            call globe_ice_gather_end(sc%gi, sc%cpl)
        case ("unknown")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_gather_begin(sc%gi, sc%cpl, 0.0_dp)
            call globe_ice_gather_add(sc%gi, sc%cpl, "C", sc%ice_a, sc%hpd_a, sc%hpd_a, sc%zb_a, sc%hpd_a, sc%fp_a)
        case ("tight")
            call build(sc, "tight")
            call reference(sc)
            call gather(sc, 0.0_dp, dHa=100.0_dp)
        case ("late_domain")
            call build(sc, "def")
            call reference(sc)
            call globe_ice_add_domain(sc%gi, sc%cpl, "A3", sc%ga, sc%fp_a, sc%hpd_a, sc%zb_a)
        case ("same_grid")
            call build(sc, "def")
            call globe_ice_add_domain(sc%gi, sc%cpl, "A3", sc%ga, sc%fp_a, sc%hpd_a, sc%zb_a)
        end select
    end subroutine abort_case

end program test_globe_ice
