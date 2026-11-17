program plastic_coulomb_calib
    ! The LGM calibration of the Coulomb yield law of the ice backend `plastic` (taub_method = "coulomb"), one grid,
    ! standalone solves without an Earth. For each (H_N, lambda_min) of the lists cf. tests/plastic_coulomb_pd.f90 and
    ! tests/plastic_coulomb_scan.f90:
    !
    !   * present day (pd_region > 0): the margin of the history's 0 ka slice, the topography file's bed, z_sl = 0, the
    !     c_f of pd_region that minimises chi2 = sum ((H - H_obs)/sigma)^2 over its grounded ice cells (grid of n_grid
    !     values of ln c_f in [lcf_lo, lcf_hi], widened at a bound, then golden section), with chi2/n, the model's and
    !     the observed volume above flotation of the region;
    !   * at time_ka: the slice's margin and bed (relative to the sea level of the time, z_sl = 0); pd_region keeps its
    !     present-day c_f, each target_region gets the c_f whose volume above flotation is target_sle [m SLE]
    !     (bisection in ln c_f in [ln cf_min, ln cf_max], all targets in the same solves), the other regions the
    !     present-day c_f (or cf_other without pd_region); then per region the volume above flotation and the
    !     maximum thickness.
    !
    ! The regions are the ids of the region file (region_id); pd_region and target_region name ids of it.
    !
    ! mu = c_f(region) lambda (taub_mu; lambda_bed = "relaxed": on the present-day bed of topo_file relaxed by local
    ! isostasy, fixed, the backend's default), smoothing smooth_km. SLE = V_af rho_i / (rho_sw A_ocean), A_ocean =
    ! 3.618e14 m2. Output <out>: per pair the present-day fit, the c_f, volumes and maximum thickness per region at
    ! time_ka (region, HN, lam_min). Analysis: tools/levante/plastic_coulomb_calib.jl. Needs plastic=1:
    !
    !     make libkryos/bin/plastic_coulomb_calib.x plastic=1
    !     plastic_coulomb_calib.x <calib.nml>          (namelist &calib, below)

    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST, model_time_axis
    use kryos_plastic_taub,   only: taub_t, taub_set_coulomb, taub_set_smooth, taub_mu, taub_set_lambda_fixed
    use kryos_plastic_bed,    only: relaxed_bed
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_init, plastic_solve_coulomb, plastic_coulomb_law_t, &
                           plastic_margin_frac_from_sdist
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

    implicit none

    integer, parameter :: NMAX = 32, MAXREG = 16
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp, A_OC = 3.618e14_dp
    ! &calib
    character(len=512) :: hist = "", topo_file = "", region_file = "", out = ""
    character(len=64)  :: margin_var = "f_ice", margin_position = "sdist", bed_var = "z_bed", H_var = "H_ice"
    character(len=64)  :: mask_var = "mask", lambda_bed = "relaxed"
    real(dp) :: margin_min = 0.5_dp, smooth_km = 0.0_dp, time_ka = 20.0_dp
    real(dp) :: HN(NMAX) = -1.0_dp, lam_min(NMAX) = -1.0_dp, p = 1.0_dp, delta = 0.05_dp, zb0 = -300.0_dp, zb1 = 200.0_dp
    integer  :: pd_region = 0, target_region(MAXREG) = 0
    real(dp) :: target_sle(MAXREG) = 0.0_dp, cf_other = 0.01_dp, cf_min = 1.0e-4_dp, cf_max = 3.0_dp
    real(dp) :: lcf_lo = -5.809142990314028_dp, lcf_hi = -3.506557897319982_dp, tol = 1.0e-3_dp
    real(dp) :: sigma_min = 50.0_dp, sigma_rel = 0.1_dp
    integer  :: n_grid = 13
    namelist /calib/ hist, margin_var, margin_min, margin_position, topo_file, bed_var, H_var, mask_var, region_file, &
                     lambda_bed, smooth_km, time_ka, HN, lam_min, p, delta, zb0, zb1, pd_region, target_region, target_sle, &
                     cf_other, cf_min, cf_max, lcf_lo, lcf_hi, tol, sigma_min, sigma_rel, n_grid, out

    type(margin_t)      :: mg
    type(plastic_par_t) :: par
    type(plastic_t), allocatable :: ps(:)
    real(dp), allocatable :: x(:), y(:), t(:), zt(:,:,:), zpd(:,:), zlgm(:,:), zsl(:,:), Hobs(:,:), sig(:,:), zrel(:,:)
    real(dp), allocatable :: sd0(:,:), sdl(:,:), fr0(:,:), frl(:,:)
    real(dp), allocatable :: lcf_pd(:,:), chi2n(:,:), sle_pd(:,:), lcf(:,:,:), sle(:,:,:), hmax(:,:,:)
    integer,  allocatable :: reg(:,:), mobs(:,:), nwide(:,:), nbis(:,:), rid(:), tpos(:)
    logical,  allocatable :: m0(:,:), ml(:,:), isfit(:,:)
    character(len=256) :: arg, units
    real(dp) :: dx, dy, dA, sle_obs
    integer  :: nx, ny, nh, nl, nt, nth, k, it, q, ih, il, th, u_nml, ios, stat, nfit, nreg

    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=calib, iostat=ios)
    if (ios /= 0) error stop "plastic_coulomb_calib: cannot read &calib"
    close(u_nml)
    nh = count(HN > 0.0_dp);  nl = count(lam_min > 0.0_dp);  nt = count(target_region > 0)
    if (trim(lambda_bed) /= "relaxed" .and. trim(lambda_bed) /= "current") error stop "plastic_coulomb_calib: lambda_bed"

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zpd(nx,ny), zsl(nx,ny), Hobs(nx,ny), mobs(nx,ny), reg(nx,ny), sd0(nx,ny), sdl(nx,ny))
    allocate(m0(nx,ny), ml(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp;  dA = dx*dy
    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, 0.0_dp, m0, sd0)
    call margin_at(mg, -1000.0_dp*time_ka, ml, sdl)
    allocate(t(nc_size(hist, "time")))
    allocate(zt(nx, ny, size(t)))
    call nc_read(hist, "time", t)
    units = "years ago"
    t = model_time_axis(t, units)
    it = minloc(abs(t + 1000.0_dp*time_ka), dim=1)
    if (abs(t(it) + 1000.0_dp*time_ka) > 1.0_dp) error stop "plastic_coulomb_calib: no slice at time_ka"
    call nc_read(hist, bed_var, zt)
    zlgm = zt(:, :, it)
    where (zlgm < -9000.0_dp) zlgm = 0.0_dp
    call nc_read(topo_file, bed_var, zpd)
    call nc_read(topo_file, H_var, Hobs)
    call nc_read(topo_file, mask_var, mobs)
    where (zpd < -9000.0_dp) zpd = 0.0_dp
    where (Hobs < 0.0_dp) Hobs = 0.0_dp
    zsl = 0.0_dp
    zrel = relaxed_bed(zpd, Hobs, RHO_I)
    sig = max(sigma_min, sigma_rel*Hobs)
    call nc_read(region_file, "region", reg)
    nreg = nc_size(region_file, "region")
    allocate(rid(nreg), tpos(nt))
    call nc_read(region_file, "region_id", rid)
    do k = 1, nt
        tpos(k) = findloc(rid, target_region(k), dim=1)
        if (tpos(k) == 0) error stop "plastic_coulomb_calib: a target_region is not a region_id of the region file"
    end do
    if (pd_region > 0 .and. .not. any(rid == pd_region)) error stop "plastic_coulomb_calib: pd_region is not a region_id"
    if (count((m0 .or. ml) .and. reg <= 0) > 0) error stop "plastic_coulomb_calib: ice cells without a region"
    isfit = m0 .and. mobs /= 3 .and. reg == pd_region
    nfit = count(isfit)
    sle_obs = 0.0_dp
    if (pd_region > 0) sle_obs = sle_of(Hobs, zpd, reg == pd_region .and. mobs /= 3)

    par%rho_i = RHO_I;  par%rho_w = RHO_SW;  par%g = GRAV;  par%mode = "surface"
    allocate(fr0(nx, ny), frl(nx, ny))
    if (trim(margin_position) == "sdist") then
        call plastic_margin_frac_from_sdist(sd0, m0, dx, dy, fr0)
        call plastic_margin_frac_from_sdist(sdl, ml, dx, dy, frl)
    else
        error stop "plastic_coulomb_calib: margin_position must be sdist"
    end if
    nth = 1
    !$ nth = omp_get_max_threads()
    allocate(ps(0:nth-1))
    do k = 0, nth - 1
        call plastic_init(ps(k), nx, ny, dx, dy, par, stat=stat)
    end do
    write(*,"(a,f6.1,a,i0,a,i0,a,i0,a,i0,a)") "# plastic_coulomb_calib: "//trim(hist)//", 0 and ", time_ka, " ka: ", &
        count(m0), " and ", count(ml), " ice cells; pd_region ", pd_region, ", ", nt, " targets; lambda on the "// &
        trim(lambda_bed)//" bed"
    if (pd_region > 0) write(*,"(a,i0,a,f8.3,a)") "# observed present day, region ", pd_region, ": ", sle_obs, " m SLE"

    allocate(lcf_pd(nh, nl), chi2n(nh, nl), sle_pd(nh, nl), lcf(nreg, nh, nl), sle(nreg, nh, nl), hmax(nreg, nh, nl))
    allocate(nwide(nh, nl), nbis(nh, nl))
    lcf_pd = log(cf_other);  chi2n = 0.0_dp;  sle_pd = 0.0_dp;  nwide = 0
    !$omp parallel do private(q, ih, il, th) schedule(dynamic)
    do q = 1, nh*nl
        ih = mod(q - 1, nh) + 1
        il = (q - 1)/nh + 1
        th = 0
        !$ th = omp_get_thread_num()
        if (pd_region > 0) call fit_pd(ih, il, th)
        call fit_lgm(ih, il, th)
    end do
    !$omp end parallel do

    write(*,"(a)") "#      H_N  lam_min   c_f(PD)  chi2/n  SLE_PD | per region: c_f, SLE, max H at time_ka"
    do il = 1, nl
        do ih = 1, nh
            write(*,"(f10.0,f9.3,f10.4,f8.3,f8.2,' |',*(f8.4,f7.1,f7.0))") HN(ih), lam_min(il), exp(lcf_pd(ih, il)), &
                chi2n(ih, il), sle_pd(ih, il), (exp(lcf(k, ih, il)), sle(k, ih, il), hmax(k, ih, il), k = 1, nreg)
        end do
    end do

    call nc_create(trim(out))
    call nc_write_dim(trim(out), "region", x=real(rid, dp))
    if (nt > 0) call nc_write_dim(trim(out), "target", x=[(real(k, dp), k = 1, nt)])
    call nc_write_dim(trim(out), "HN", x=HN(1:nh), units="m")
    call nc_write_dim(trim(out), "lam_min", x=lam_min(1:nl), units="1")
    call nc_write_dim(trim(out), "one", x=[1.0_dp])
    call nc_write(trim(out), "lcf_pd", lcf_pd, dim1="HN", dim2="lam_min", long_name="ln c_f of pd_region, PD fit")
    call nc_write(trim(out), "chi2n_pd", chi2n, dim1="HN", dim2="lam_min", long_name="chi2/n of the PD fit")
    call nc_write(trim(out), "sle_pd", sle_pd, dim1="HN", dim2="lam_min", units="m", &
                  long_name="PD volume above flotation of pd_region, model")
    call nc_write(trim(out), "sle_pd_obs", [sle_obs], dim1="one", units="m", long_name="observed, the same")
    call nc_write(trim(out), "lcf", lcf, dim1="region", dim2="HN", dim3="lam_min", long_name="ln c_f at time_ka")
    call nc_write(trim(out), "sle", sle, dim1="region", dim2="HN", dim3="lam_min", units="m", &
                  long_name="volume above flotation at time_ka [m SLE]")
    call nc_write(trim(out), "hmax", hmax, dim1="region", dim2="HN", dim3="lam_min", units="m", &
                  long_name="maximum thickness at time_ka")
    call nc_write(trim(out), "pd_widened", nwide, dim1="HN", dim2="lam_min", long_name="widenings of the PD interval")
    call nc_write(trim(out), "n_bisect", nbis, dim1="HN", dim2="lam_min")
    call nc_write(trim(out), "pd_region", [pd_region], dim1="one")
    call nc_write(trim(out), "region_id", rid, dim1="region", long_name="region id of the region file")
    if (nt > 0) then
        call nc_write(trim(out), "target_region", target_region(1:nt), dim1="target", long_name="region id of the target")
        call nc_write(trim(out), "target_sle", target_sle(1:nt), dim1="target", units="m")
    end if
    call nc_write(trim(out), "time_ka", [time_ka], dim1="one")
    write(*,"(a)") "# wrote "//trim(out)

contains

    real(dp) function sle_of(H, zb, m) result(s)
        ! volume above flotation over the cells m [m SLE]
        real(dp), intent(in) :: H(:,:), zb(:,:)
        logical,  intent(in) :: m(:,:)
        s = sum(max(H - (RHO_SW/RHO_I)*max(-zb, 0.0_dp), 0.0_dp), mask=m .and. H > 0.0_dp)*dA*RHO_I/(RHO_SW*A_OC)
    end function sle_of

    subroutine solve(lc, jh, jl, tt, zb, mask, frac)
        ! the Coulomb solve with ln c_f per region lc on bed zb and mask, solver instance tt
        real(dp), intent(in) :: lc(:), zb(:,:), frac(:,:)
        integer,  intent(in) :: jh, jl, tt
        logical,  intent(in) :: mask(:,:)
        type(taub_t) :: tb
        real(dp), allocatable :: mu(:,:)
        integer :: st
        allocate(mu(nx, ny))
        call taub_set_coulomb(tb, reg, rid, lc, p, HN(jh), delta, lam_min(jl), zb0, zb1, trim(region_file), &
                              "calibration c_f")
        if (trim(lambda_bed) == "relaxed") call taub_set_lambda_fixed(tb, zrel)
        call taub_set_smooth(tb, 1000.0_dp*smooth_km, dx, dy)
        call taub_mu(tb, zb, zsl, mu)
        call plastic_solve_coulomb(ps(tt), zb, zsl, mu, mask, law(jh), warm=.false., stat=st, margin_frac=frac)
    end subroutine solve

    real(dp) function chi2_pd(lc, jh, jl, tt) result(c2)
        real(dp), intent(in) :: lc
        integer,  intent(in) :: jh, jl, tt
        real(dp) :: v(nreg)
        v = lc
        call solve(v, jh, jl, tt, zpd, m0, fr0)
        c2 = sum(((ps(tt)%H - Hobs)/sig)**2, mask=isfit)/nfit
    end function chi2_pd

    subroutine fit_pd(jh, jl, tt)
        ! ln c_f of pd_region by a grid, widened at a bound, then golden section
        integer, intent(in) :: jh, jl, tt
        real(dp), parameter :: GR = 0.6180339887498949_dp
        real(dp) :: g(n_grid), f(n_grid), lo, hi, a, b, c, d, fc, fd, w
        integer  :: i, m, nw
        lo = lcf_lo;  hi = lcf_hi;  nw = 0
        do
            do i = 1, n_grid
                g(i) = lo + (hi - lo)*(i - 1)/real(n_grid - 1, dp)
                f(i) = chi2_pd(g(i), jh, jl, tt)
            end do
            m = minloc(f, dim=1)
            if ((m > 1 .and. m < n_grid) .or. nw >= 6) exit
            w = 0.5_dp*(hi - lo)
            if (m == 1) then
                lo = lo - w
            else
                hi = hi + w
            end if
            nw = nw + 1
        end do
        a = g(max(m - 1, 1));  b = g(min(m + 1, n_grid))
        c = b - GR*(b - a);  d = a + GR*(b - a)
        fc = chi2_pd(c, jh, jl, tt);  fd = chi2_pd(d, jh, jl, tt)
        do while (b - a > tol)
            if (fc < fd) then
                b = d;  d = c;  fd = fc;  c = b - GR*(b - a);  fc = chi2_pd(c, jh, jl, tt)
            else
                a = c;  c = d;  fc = fd;  d = a + GR*(b - a);  fd = chi2_pd(d, jh, jl, tt)
            end if
        end do
        lcf_pd(jh, jl) = merge(c, d, fc < fd)
        if (f(m) < min(fc, fd)) lcf_pd(jh, jl) = g(m)
        chi2n(jh, jl) = chi2_pd(lcf_pd(jh, jl), jh, jl, tt)
        sle_pd(jh, jl) = sle_of(ps(tt)%H, zpd, reg == pd_region .and. mobs /= 3)
        nwide(jh, jl) = nw
    end subroutine fit_pd

    subroutine fit_lgm(jh, jl, tt)
        ! the targets' ln c_f by simultaneous bisection, then the volumes and maximum thickness per region
        integer, intent(in) :: jh, jl, tt
        real(dp) :: lc(nreg), lo(nreg), hi(nreg)
        integer  :: r, j, nb
        lc = lcf_pd(jh, jl)
        lo = log(cf_min);  hi = log(cf_max)
        nb = 0
        if (nt > 0) then
            do while (maxval(hi - lo, mask=[(any(tpos == r), r = 1, nreg)]) > tol)
                do j = 1, nt
                    r = tpos(j)
                    lc(r) = 0.5_dp*(lo(r) + hi(r))
                end do
                call solve(lc, jh, jl, tt, zlgm, ml, frl)
                nb = nb + 1
                do j = 1, nt
                    r = tpos(j)
                    if (sle_of(ps(tt)%H, zlgm, reg == rid(r)) < target_sle(j)) then
                        lo(r) = lc(r)
                    else
                        hi(r) = lc(r)
                    end if
                end do
            end do
            do j = 1, nt
                lc(tpos(j)) = 0.5_dp*(lo(tpos(j)) + hi(tpos(j)))
            end do
        end if
        call solve(lc, jh, jl, tt, zlgm, ml, frl)
        do r = 1, nreg
            lcf(r, jh, jl) = lc(r)
            sle(r, jh, jl) = sle_of(ps(tt)%H, zlgm, reg == rid(r))
            hmax(r, jh, jl) = maxval(ps(tt)%H, mask=reg == rid(r) .and. ml)
            if (.not. any(reg == rid(r) .and. ml)) hmax(r, jh, jl) = 0.0_dp
        end do
        nbis(jh, jl) = nb
    end subroutine fit_lgm


    function law(jh) result(lw)
        ! The library's Coulomb law of H_N number jh (p and delta of the namelist).
        integer, intent(in) :: jh
        type(plastic_coulomb_law_t) :: lw
        lw = plastic_coulomb_law_t(p=p, H_N=HN(jh), delta=delta)
    end function law

end program plastic_coulomb_calib
