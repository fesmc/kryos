program test_geo_vilma2
    ! The VILMA2 geo adapter against VILMA2 driven directly (phase 4 plan 13.1), with
    ! prescribed host terms and no ice instance: the global part (vilma_earth_t) through
    ! the component wrappers, as the core drives it (reference, PASS_GEOMETRY, spin-up and
    ! seed at PASS_STATE, updates, restart), and the regional port (vilma_geo_t).
    !
    !   * Gauss grid: GAUSS-L<lmax>, its size and name, areas = quadrature areas, the
    !     description file (written, then compared, also its planet radius), the grid's sphere has
    !     VILMA2's r_earth (strict, to round-off);
    !   * densities and r_earth: from the run's constants, not the keys;
    !   * init at the reference (exports zero, reference ocean function), spin-up (a no-op
    !     and pre_spinup_1d), the dt = 0 seed, three updates: every export bitwise equal to
    !     solid_earth_* called with the same arrays (dz_ss as VILMA2's minus dphi);
    !   * restart round trip: write, rebuild at the reference, read, continue: bitwise;
    !   * the port copies its imports exactly, with their validity;
    !   * stops (child runs): a grid VILMA2 does not take as its own (se%remap), a grid file of
    !     another size, a host term that is not admissible, an update before the build.
    ! With vilma=0 the test is not built (test_skip in config/common.mk).

    use precision,       only: sp, dp
    use coordinates,     only: grid_class, grid_init
    use ncio,            only: nc_read
    use constants,       only: pi
    use phys_constants,  only: phys_const_load
    use fesm_fields,     only: fields_get, fields_put, fields_restart_write, fields_restart_read
    use kryos_component, only: init_ctx_t, step_time_t, census_t, census_count, PASS_INIT, PASS_GEOMETRY, &
                               PASS_STATE, PASS_RESTART
    use kryos_geo,       only: geo_slot_t, geo_create, geo_init, geo_init_state, geo_update, &
                               geo_restart_write, geo_restart_read, geo_write_output, geo_end
    use kryos_geo_vilma2, only: vilma_earth_t, vilma2_gauss_dims
    use vilma_params,    only: vilma_par_load
    use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_spinup, &
                               solid_earth_finalize
    use vilma_sht,       only: sht_grid_surface_integral

    implicit none

    integer, parameter :: LMAX = 16, NPHI = 64, NLAT = 34, NXP = 6, NYP = 5
    real(dp), parameter :: DTT = 10.0_dp
    character(len=*), parameter :: DIR = "tmp/test_geo_vilma2/", PAR = DIR//"par.nml", PAR_SP = DIR//"par_sp.nml", &
                                   PAR_BAD = DIR//"par_bad.nml", RST = DIR//"rst", MAPS = DIR//"maps", errfile = DIR//"abort.err"
    integer :: nfail = 0

    type(init_ctx_t) :: ctxg, ctxs
    type(census_t)   :: census, census2, census3         ! one global part per census (max_instances = 1)
    type(geo_slot_t) :: glob, glob2, spin
    type(solid_earth) :: dse, dsp
    real(dp) :: lat(NPHI,NLAT), lon(NPHI,NLAT)
    real(dp) :: z_ref(NPHI,NLAT), H_ref(NPHI,NLAT), w(NPHI,NLAT), Ig_ref(NPHI,NLAT), C_ref(NPHI,NLAT)
    real(dp) :: t
    integer  :: n

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    call execute_command_line("rm -rf "//DIR//" && mkdir -p "//DIR)
    call par_write()
    call contexts()
    call test_dims()

    ! ---- the adapter: init, reference ------------------------------------------------
    call geo_create(glob, ctxg, census)
    call check(census_count(census, "geo:vilma2:global") == 1, "factory: one vilma2 global part admitted")
    call geo_init(glob, ctxg)
    call test_grid(glob)
    call coordinates_of(glob)
    call reference_terms()
    call put_reference(glob)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_INIT)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)

    ! ---- direct: the same calls --------------------------------------------------------
    call direct_load(dse, PAR)
    call solid_earth_init(dse, z_ref, H_ref, grid=gauss_of(glob), w=w, I_g=Ig_ref, C=C_ref)
    call check(.not. dse%remap, "direct: the grid is VILMA2's own (no remap)")
    call compare(glob, dse, "reference")
    call check(maxval(abs(field(glob, "dz_bed"))) == 0.0_dp .and. maxval(abs(field(glob, "dz_ss"))) == 0.0_dp, &
               "reference: dz_bed = dz_ss = 0")
    call check(scalar(glob, "dphi") == 0.0_dp .and. scalar(glob, "bsl") == 0.0_dp, "reference: dphi = bsl = 0")

    ! ---- cold start: no spin-up (a no-op), seed at the t0 load -----------------------------
    call put_load(glob, 0, 0.0_dp)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_STATE)
    call solid_earth_spinup(dse, H_ext(0), w, Ig(0), Cc(0), Ff(0))
    call solid_earth_update(dse, H_ext(0), 0.0_dp, w, Ig(0), Cc(0), Ff(0))
    call compare(glob, dse, "seed")
    call check(abs(dse%time) == 0.0_dp, "no spin-up: VILMA2's clock has not moved")

    ! ---- steps ----------------------------------------------------------------------------
    do n = 1, 3
        t = real(n, dp) * DTT
        call put_load(glob, n, t)
        call geo_update(glob, st(t - DTT, t))
        call solid_earth_update(dse, H_ext(n), DTT, w, Ig(n), Cc(n), Ff(n))
        call compare(glob, dse, "step "//char(48 + n))
    end do
    call check(minval(field(glob, "dz_bed")) < 0.0_dp, "the grown ice load depresses the bed somewhere")
    call check(abs(scalar(glob, "bsl")) > 0.0_dp, "bsl moved")
    call check(abs(dse%time - 3.0_dp*DTT) < 1.0e-9_dp, "VILMA2's clock advanced by the steps")
    call test_series(glob, t)

    ! ---- restart round trip -----------------------------------------------------------------
    call execute_command_line("rm -rf "//RST//" && mkdir -p "//RST)
    call geo_restart_write(glob, RST, t)
    call fields_write(glob, RST)
    call check(file_exists(RST//"/geo_glob_restart.nc"), "restart file <name>_restart.nc in the bundle folder")

    call geo_create(glob2, ctxg, census2)
    call geo_init(glob2, ctxg)
    call put_reference(glob2)
    call geo_init_state(glob2, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)      ! built at the reference, memory zero
    call geo_restart_read(glob2, RST, t)
    call fields_read(glob2, RST)
    call geo_init_state(glob2, st(t, t), PASS_RESTART)
    call compare(glob2, dse, "restored")
    t = t + DTT
    call put_load(glob, 4, t)
    call put_load(glob2, 4, t)
    call geo_update(glob,  st(t - DTT, t))
    call geo_update(glob2, st(t - DTT, t))
    call solid_earth_update(dse, H_ext(4), DTT, w, Ig(4), Cc(4), Ff(4))
    call compare(glob,  dse, "continued")
    call compare(glob2, dse, "restarted")

    ! ---- spin-up (VILMA2's pre_spinup_1d) then the seed, as the direct calls -------------------
    call geo_create(spin, ctxs, census3)
    call geo_init(spin, ctxs)
    call put_reference(spin)
    call geo_init_state(spin, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
    call direct_load(dsp, PAR_SP)
    call solid_earth_init(dsp, z_ref, H_ref, grid=gauss_of(spin), w=w, I_g=Ig_ref, C=C_ref)
    call put_load(spin, 0, 0.0_dp)
    call geo_init_state(spin, st(0.0_dp, 0.0_dp), PASS_STATE)
    call solid_earth_spinup(dsp, H_ext(0), w, Ig(0), Cc(0), Ff(0))
    call solid_earth_update(dsp, H_ext(0), 0.0_dp, w, Ig(0), Cc(0), Ff(0))
    call compare(spin, dsp, "spin-up and seed")
    call check(dsp%time > 0.0_dp, "pre_spinup_1d advanced VILMA2's clock")

    call test_port()
    call test_densities(glob)

    call geo_end(spin); call geo_end(glob2); call geo_end(glob)
    call solid_earth_finalize(dsp); call solid_earth_finalize(dse)

    call expect_abort("wrong_grid", "does not take grid")
    call expect_abort("grid_file_size", "has another size")
    call expect_abort("host_term", "are not admissible")
    call expect_abort("before_build", "has not been built")
    call expect_abort("no_published_import", "has not been published")
    call expect_abort("grid_file_radius", "has another planet radius")
    call expect_abort("group_key", "must name group `vilma`")

    call execute_command_line("rm -rf "//DIR)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_geo_vilma2: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_geo_vilma2: OK"

contains

    ! ------------------------------------------------------------------------------------------
    ! The load: south cap (host region, w = 1), a partial band (w = 0.5), a northern cap of
    ! external ice (the 1 - w part); load n grows the south cap and shrinks the north one.
    ! ------------------------------------------------------------------------------------------

    subroutine coordinates_of(g)
        type(geo_slot_t), intent(inout) :: g
        type(grid_class) :: gg
        integer :: i, j
        gg = gauss_of(g)
        do j = 1, NLAT
            do i = 1, NPHI
                lon(i,j) = gg%G%x(i)
                lat(i,j) = gg%G%y(j)
            end do
        end do
    end subroutine coordinates_of

    subroutine reference_terms()
        w = 0.0_dp; Ig_ref = 0.0_dp; C_ref = 0.0_dp
        z_ref = -3000.0_dp; H_ref = 0.0_dp
        where (lat > 55.0_dp) z_ref = 200.0_dp
        where (lat > 72.0_dp) H_ref = 2000.0_dp
        where (lat < -60.0_dp)                          ! the host region: all land under ice
            w = 1.0_dp
            Ig_ref = 1500.0_dp
        end where
        where (lat >= -60.0_dp .and. lat < -50.0_dp)    ! a partial band: half host, a quarter of that ocean
            w = 0.5_dp
            Ig_ref = 400.0_dp
            C_ref = 0.125_dp
        end where
    end subroutine reference_terms

    function H_ext(k) result(h)
        integer, intent(in) :: k
        real(dp) :: h(NPHI,NLAT)
        h = 0.0_dp
        where (lat > 72.0_dp) h = 2000.0_dp * (1.0_dp - 0.1_dp*real(k, dp))
    end function H_ext

    function Ig(k) result(a)
        integer, intent(in) :: k
        real(dp) :: a(NPHI,NLAT)
        a = Ig_ref * (1.0_dp + 0.2_dp*real(k, dp) + 0.1_dp)
    end function Ig

    function Cc(k) result(a)
        integer, intent(in) :: k
        real(dp) :: a(NPHI,NLAT)
        a = C_ref
        where (w == 0.5_dp) a = 0.5_dp * (0.25_dp + 0.05_dp*real(k, dp))
    end function Cc

    function Ff(k) result(a)
        integer, intent(in) :: k
        real(dp) :: a(NPHI,NLAT)
        a = 0.0_dp
        where (w == 0.5_dp) a = w * (-100.0_dp) * (Cc(k) - C_ref)
    end function Ff

    subroutine put_reference(g)
        type(geo_slot_t), intent(inout) :: g
        call fields_put(g%c%imports, g%c%hi_z_bed_ref,   z_ref,  t1=0.0_dp)
        call fields_put(g%c%imports, g%c%hi_H_ice_ref,   H_ref,  t1=0.0_dp)
        call fields_put(g%c%imports, g%c%hi_host_w,      w,      t1=0.0_dp)
        call fields_put(g%c%imports, g%c%hi_host_Ig_ref, Ig_ref, t1=0.0_dp)
        call fields_put(g%c%imports, g%c%hi_host_C_ref,  C_ref,  t1=0.0_dp)
    end subroutine put_reference

    subroutine put_load(g, k, t1)
        type(geo_slot_t), intent(inout) :: g
        integer,          intent(in)    :: k
        real(dp),         intent(in)    :: t1
        call fields_put(g%c%imports, g%c%hi_H_ice_ext, H_ext(k), t1=t1)
        call fields_put(g%c%imports, g%c%hi_host_w,    w,        t1=t1)
        call fields_put(g%c%imports, g%c%hi_host_Ig,   Ig(k),    t1=t1)
        call fields_put(g%c%imports, g%c%hi_host_C,    Cc(k),    t1=t1)
        call fields_put(g%c%imports, g%c%hi_host_F,    Ff(k),    t1=t1)
    end subroutine put_load

    ! ------------------------------------------------------------------------------------------
    ! Comparison with the direct calls
    ! ------------------------------------------------------------------------------------------

    function field(g, name) result(a)
        type(geo_slot_t), intent(in) :: g
        character(len=*), intent(in) :: name
        real(dp) :: a(NPHI,NLAT)
        integer :: h
        h = handle(g, name)
        call fields_get(g%c%exports, h, a)
    end function field

    function scalar(g, name) result(s)
        type(geo_slot_t), intent(in) :: g
        character(len=*), intent(in) :: name
        real(dp) :: s
        call fields_get(g%c%exports, handle(g, name), s)
    end function scalar

    integer function handle(g, name)
        type(geo_slot_t), intent(in) :: g
        character(len=*), intent(in) :: name
        select case (name)
        case ("dz_bed");     handle = g%c%h_dz_bed
        case ("dz_ss");      handle = g%c%h_dz_ss
        case ("C_ocn");      handle = g%c%h_C_ocn
        case ("dphi");       handle = g%c%h_dphi
        case ("bsl");        handle = g%c%h_bsl
        case ("bsl_ive");    handle = g%c%h_bsl_ive
        case ("A_ocn");      handle = g%c%h_A_ocn
        case ("mass_resid"); handle = g%c%h_mass_resid
        case default;        handle = 0
        end select
    end function handle

    subroutine compare(g, d, what)
        type(geo_slot_t), intent(in) :: g
        type(solid_earth), intent(in) :: d
        character(len=*), intent(in)  :: what

        real(dp) :: a(NPHI,NLAT)
        logical :: ok

        a = field(g, "dz_bed");  ok = all(a == d%gg%dz_bed)
        a = field(g, "dz_ss");   ok = ok .and. all(a == d%gg%dz_ss - d%dphi)
        a = field(g, "C_ocn");   ok = ok .and. all(a == d%gg%C)
        call check(ok, what//": dz_bed, dz_ss (= VILMA2's minus dphi), C_ocn bitwise equal to the direct calls")
        call check(scalar(g, "dphi") == d%dphi .and. scalar(g, "bsl") == d%bsl .and. &
                   scalar(g, "bsl_ive") == d%bsl_ive .and. scalar(g, "mass_resid") == d%worst_mass_resid, &
                   what//": dphi, bsl, bsl_ive, mass_resid bitwise equal")
        call check(scalar(g, "A_ocn") == d%par%r_earth**2 * sht_grid_surface_integral(d%sht, d%gg%C), &
                   what//": A_ocn = r_earth^2 times the integral of C")
    end subroutine compare

    ! ------------------------------------------------------------------------------------------
    ! Pieces
    ! ------------------------------------------------------------------------------------------

    function gauss_of(g) result(gg)
        type(geo_slot_t), intent(in) :: g
        type(grid_class) :: gg
        select type (c => g%c)
        type is (vilma_earth_t)
            gg = c%grid_obj
        end select
    end function gauss_of

    subroutine direct_load(d, par_file)
        type(solid_earth), intent(inout) :: d
        character(len=*),  intent(in)    :: par_file
        call vilma_par_load(d%par, par_file, defaults_file="input/vilma_defaults.nml", group="vilma")
        d%par%rho_ice   = ctxg%cnst%rho_ice
        d%par%rho_water = ctxg%cnst%rho_sw
    end subroutine direct_load

    subroutine test_dims()
        integer :: nlat_, nphi_
        character(len=64) :: name

        call vilma2_gauss_dims(170, 0, 0, nlat_, nphi_, name)
        call check(nlat_ == 342 .and. nphi_ == 680 .and. trim(name) == "GAUSS-L170", "lmax 170: 680 x 342, GAUSS-L170")
        call vilma2_gauss_dims(170, 256, 512, nlat_, nphi_, name)
        call check(nlat_ == 256 .and. nphi_ == 512 .and. trim(name) == "GAUSS-L170-512x256", "own nlat, nphi: size in the name")
    end subroutine test_dims

    subroutine test_grid(g)
        type(geo_slot_t), intent(in) :: g
        type(grid_class) :: gg
        integer :: i

        gg = gauss_of(g)
        call check(trim(g%c%grid) == "GAUSS-L16" .and. g%c%nx == NPHI .and. g%c%ny == NLAT, "grid GAUSS-L16, 64 x 34")
        call check(gg%G%nx == NPHI .and. gg%G%ny == NLAT .and. trim(gg%cs%mtype) == "gaussian" .and. &
                   gg%G%x(1) == 0.0_dp .and. gg%G%y(1) < gg%G%y(NLAT), "Gaussian, longitudes from 0, latitudes ascending")
        select type (c => g%c)
        type is (vilma_earth_t)
            call check(abs(gg%cs%planet%R - c%se%par%r_earth) <= 1.0e-12_dp*c%se%par%r_earth .and. &
                       abs(sum(gg%area) - 4.0_dp*pi*c%se%par%r_earth**2) <= 1.0e-12_dp*4.0_dp*pi*c%se%par%r_earth**2, &
                       "the grid is a sphere of VILMA2's r_earth: the areas sum to 4 pi r_earth^2")
        end select
        call check(abs(sum(gg%area) / (4.0_dp*pi*gg%cs%planet%R**2) - 1.0_dp) < 1.0e-12_dp, &
                   "cell areas are the quadrature areas (sum = 4 pi R^2)")
        i = NLAT / 2
        call check(abs(gg%area(1,i) - 2.0_dp*pi/real(NPHI, dp) * gg%cs%planet%R**2 * &
                       (sin(gg%G%y_bnds(i)*pi/180.0_dp) - sin(gg%G%y_bnds(i-1)*pi/180.0_dp))) < 1.0e-3_dp, &
                   "a cell area is R^2 dlon w_j")
        call check(file_exists(MAPS//"/grid_GAUSS-L16.txt"), "grid description written to the map folder")
    end subroutine test_grid

    subroutine test_densities(g)
        type(geo_slot_t), intent(in) :: g
        select type (c => g%c)
        type is (vilma_earth_t)
            call check(c%se%par%rho_ice == ctxg%cnst%rho_ice .and. c%se%par%rho_water == ctxg%cnst%rho_sw .and. &
                       c%se%par%rho_ice /= 931.0_dp, "densities of VILMA2 are the run's rho_ice, rho_sw")
        end select
    end subroutine test_densities

    subroutine test_port()
        ! The regional port: reference at PASS_GEOMETRY, then it copies what is scattered, with its validity.
        type(geo_slot_t) :: p
        type(init_ctx_t) :: ctx
        type(census_t)   :: cs
        real(dp) :: bed(NXP,NYP), a(NXP,NYP), s, t0, t1
        integer :: i, j

        write(*,"(a)") "regional port"
        ctx%instance  = "geo_nh"
        ctx%component = "geo"
        ctx%backend   = "vilma2"
        ctx%domain    = "Toy"
        ctx%grid_name = "TOY"
        ctx%scope     = "regional"
        call grid_init(ctx%grid, name="TOY", mtype="cartesian", units="km", dx=32.0_dp, nx=NXP, dy=32.0_dp, ny=NYP, &
                       x0=0.0_dp, y0=0.0_dp)
        call geo_create(p, ctx, cs)
        call check(census_count(cs, "geo:vilma2:regional") == 1 .and. .not. p%c%con%restart, &
                   "factory: the regional instance is the port, without a backend restart file")
        call geo_init(p, ctx)
        call check(p%c%opt%port .and. p%c%nx == NXP .and. p%c%ny == NYP .and. p%c%hi_dz_bed > 0 .and. &
                   p%c%hi_z_bed_ref > 0 .and. p%c%hi_H_ice == 0, "port vocabulary on the ice grid, no H_ice import")

        do j = 1, NYP
            do i = 1, NXP
                bed(i,j) = -1000.0_dp + 37.0_dp*real(i, dp) - 11.0_dp*real(j, dp)
            end do
        end do
        call fields_put(p%c%imports, p%c%hi_z_bed_ref, bed, t1=0.0_dp)
        call geo_init_state(p, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
        call fields_get(p%c%exports, p%c%h_dz_bed, a)
        call fields_get(p%c%exports, p%c%h_dphi, s)
        call check(all(a == 0.0_dp) .and. s == 0.0_dp, "PASS_GEOMETRY: zero anomalies, dphi = 0")
        call fields_get(p%c%exports, p%c%h_z_bed_ref, a)
        call check(all(a == bed), "PASS_GEOMETRY: z_bed_ref passed through exactly")

        do j = 1, NYP
            do i = 1, NXP
                a(i,j) = 0.1_dp*real(i, dp)/3.0_dp - 0.07_dp*real(j, dp)
            end do
        end do
        call fields_put(p%c%imports, p%c%hi_dz_bed, a,          t1=10.0_dp, t0=0.0_dp)
        call fields_put(p%c%imports, p%c%hi_dz_ss,  2.0_dp*a,   t1=10.0_dp, t0=0.0_dp)
        call fields_put(p%c%imports, p%c%hi_dphi,   -1.0_dp/7.0_dp, t1=10.0_dp, t0=0.0_dp)
        call fields_put(p%c%imports, p%c%hi_bsl,    -1.0_dp/9.0_dp, t1=10.0_dp, t0=0.0_dp)
        call fields_put(p%c%imports, p%c%hi_z_bed_ref, bed, t1=10.0_dp)
        call geo_update(p, st(0.0_dp, 10.0_dp))
        call check(all(exp2(p, p%c%h_dz_bed) == a) .and. all(exp2(p, p%c%h_dz_ss) == 2.0_dp*a), "update: dz_bed, dz_ss copied exactly")
        call fields_get(p%c%exports, p%c%h_dphi, s, t0=t0, t1=t1)
        call check(s == -1.0_dp/7.0_dp .and. t1 == 10.0_dp .and. t0 == 0.0_dp, "update: dphi copied with the import's validity")
        call fields_get(p%c%exports, p%c%h_bsl, s)
        call check(s == -1.0_dp/9.0_dp .and. all(exp2(p, p%c%h_z_bed_ref) == bed), "update: bsl copied, z_bed_ref still the ice's")
        call geo_end(p)
    end subroutine test_port

    subroutine test_series(g, t1)
        ! The default writer's time series (stream 1D: the 0-D exports) works for the global part
        ! (the 2D stream needs the lon-lat writer of WP5).
        type(geo_slot_t), intent(inout) :: g
        real(dp),         intent(in)    :: t1
        real(sp) :: v(1)

        call geo_write_output(g, DIR//"ts.nc", t1, "1D", .true.)
        call nc_read(DIR//"ts.nc", "bsl", v, start=[1], count=[1])
        call check(abs(real(v(1), dp) - scalar(g, "bsl")) <= 1.0e-3_dp, "stream 1D: bsl written (rounded to 1e-3 m)")
        call nc_read(DIR//"ts.nc", "dphi", v, start=[1], count=[1])
        call check(abs(real(v(1), dp) - scalar(g, "dphi")) <= 1.0e-3_dp, "stream 1D: dphi written")
    end subroutine test_series

    function exp2(g, h) result(a)
        type(geo_slot_t), intent(in) :: g
        integer,          intent(in) :: h
        real(dp) :: a(NXP,NYP)
        call fields_get(g%c%exports, h, a)
    end function exp2

    ! ------------------------------------------------------------------------------------------
    ! Bundle files as kryos_restart writes them
    ! ------------------------------------------------------------------------------------------

    subroutine fields_write(g, dir)
        type(geo_slot_t), intent(in) :: g
        character(len=*), intent(in) :: dir
        call fields_restart_write(g%c%imports, dir//"/"//trim(g%c%name)//"_imports.nc")
        call fields_restart_write(g%c%exports, dir//"/"//trim(g%c%name)//"_exports.nc")
    end subroutine fields_write

    subroutine fields_read(g, dir)
        type(geo_slot_t), intent(inout) :: g
        character(len=*), intent(in)    :: dir
        call fields_restart_read(g%c%imports, dir//"/"//trim(g%c%name)//"_imports.nc")
        call fields_restart_read(g%c%exports, dir//"/"//trim(g%c%name)//"_exports.nc")
    end subroutine fields_read

    function st(t0, t1) result(s)
        real(dp), intent(in) :: t0, t1
        type(step_time_t) :: s
        s = step_time_t(t0=t0, t1=t1, t_eval=t1, dt=DTT, dtt=DTT)
    end function st

    logical function file_exists(f)
        character(len=*), intent(in) :: f
        inquire(file=f, exist=file_exists)
    end function file_exists

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (ok) then
            write(*,"(a)") "  ok    "//what
        else
            write(*,"(a)") "  FAIL  "//what
            nfail = nfail + 1
        end if
    end subroutine check

    ! ------------------------------------------------------------------------------------------
    ! Contexts and the par file
    ! ------------------------------------------------------------------------------------------

    subroutine contexts()
        ctxg%par_file  = PAR
        ctxg%instance  = "geo_glob"
        ctxg%component = "geo"
        ctxg%backend   = "vilma2"
        ctxg%scope     = "global"
        ctxg%map_fldr  = MAPS
        ctxg%time      = 0.0_dp
        ctxg%time_rel  = 0.0_dp
        ctxg%dtt       = DTT
        call phys_const_load(ctxg%cnst, "fesm-utils/par/phys_const_earth.nml", "Earth")
        ctxs = ctxg
        ctxs%instance = "geo_sp"
        ctxs%par_file = PAR_SP
    end subroutine contexts

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&geo_glob", ' backend = "vilma2"', ' vilma = "vilma"', "/", &
                       "&vilma", " lmax = 16", "/"
        close(u)
        open(newunit=u, file=PAR_SP, status="replace", action="write")
        write(u,'(a)') "&geo_sp", ' backend = "vilma2"', "/", &
                       "&vilma", " lmax = 16", " pre_spinup_1d = .true.", "/"
        close(u)
        open(newunit=u, file=PAR_BAD, status="replace", action="write")        ! the stop: a group name VILMA2 cannot read
        write(u,'(a)') "&geo_key", ' backend = "vilma2"', ' vilma = "other"', "/"
        close(u)
    end subroutine par_write

    ! ------------------------------------------------------------------------------------------
    ! Stops, in child runs of this program
    ! ------------------------------------------------------------------------------------------

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
        call check(es /= 0 .and. es /= 77 .and. index(all_, text) > 0, case_name//": stops with '"//text//"'")
    end subroutine expect_abort

    subroutine abort_case()
        character(len=64) :: arg
        type(geo_slot_t) :: g

        call get_command_argument(1, arg)
        call contexts()
        call coordinates_for_child()
        call reference_terms()
        if (trim(arg) == "group_key") then
            ctxg%par_file = PAR_BAD
            ctxg%instance = "geo_key"
        end if
        if (trim(arg) /= "grid_file_size" .and. trim(arg) /= "grid_file_radius") then   ! (those build the instance themselves)
            call geo_create(g, ctxg, census)
            call geo_init(g, ctxg)
        end if
        select case (trim(arg))
        case ("wrong_grid")                             ! the grid is shifted: VILMA2 builds a remap
            select type (c => g%c)
            type is (vilma_earth_t)
                c%grid_obj%G%x = c%grid_obj%G%x + 0.01_dp
            end select
            call put_reference(g)
            call geo_init_state(g, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
        case ("grid_file_size")                         ! a description of another size under the same name
            call small_grid_file(DIR//"maps2", 32, 18)
            ctxg%map_fldr = DIR//"maps2"
            call geo_create(g, ctxg, census)
            call geo_init(g, ctxg)
        case ("grid_file_radius")                       ! the right grid, on the default sphere (6371221 m), not r_earth
            call small_grid_file(DIR//"maps3", NPHI, NLAT)
            ctxg%map_fldr = DIR//"maps3"
            call geo_create(g, ctxg, census)
            call geo_init(g, ctxg)
        case ("host_term")                              ! C > w
            Ig_ref = 0.0_dp
            C_ref = w + 0.5_dp
            call put_reference(g)
            call geo_init_state(g, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
        case ("before_build")
            call put_load(g, 1, 10.0_dp)
            call geo_update(g, st(0.0_dp, 10.0_dp))
        case ("no_published_import")
            call geo_init_state(g, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
        end select
    end subroutine abort_case

    subroutine small_grid_file(fldr, nlon, nlat)
        ! A valid description of a Gaussian grid of this size, on the default sphere (fesm-utils'
        ! "Spherical Earth", not VILMA2's radius), under the name of the adapter's grid.
        use coordinates, only: grid_init_gaussian
        use coords,      only: grid_cdo_write_desc_short
        character(len=*), intent(in) :: fldr
        integer,          intent(in) :: nlon, nlat
        type(grid_class) :: g8
        call execute_command_line("rm -rf "//fldr//" && mkdir -p "//fldr)
        call grid_init_gaussian(g8, "GAUSS-L16", nlon=nlon, nlat=nlat, lon180=.false., x0=0.0_dp)
        call grid_cdo_write_desc_short(g8, fldr)
    end subroutine small_grid_file

    subroutine coordinates_for_child()
        ! The latitudes of the grid the adapter builds, before it exists: from a grid of the same size.
        type(grid_class) :: gg
        integer :: i, j
        call grid_init_gaussian_for_test(gg)
        do j = 1, NLAT
            do i = 1, NPHI
                lon(i,j) = gg%G%x(i)
                lat(i,j) = gg%G%y(j)
            end do
        end do
    end subroutine coordinates_for_child

    subroutine grid_init_gaussian_for_test(gg)
        use coordinates, only: grid_init_gaussian
        type(grid_class), intent(out) :: gg
        call grid_init_gaussian(gg, "TEST", nlon=NPHI, nlat=NLAT, lon180=.false., x0=0.0_dp)
    end subroutine grid_init_gaussian_for_test

end program test_geo_vilma2
