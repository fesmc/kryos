program plastic_coulomb_pd
    ! The present-day test of the Coulomb yield law of the ice backend `plastic` (taub_method = "coulomb";
    ! docs/backend-plastic.md): tau_b = c_f lambda(z_bed - z_sl) N_eff, solved with the geometry by
    ! plastic_solve_coulomb. Offline and deterministic, with the backend's modules and the library, set up as
    ! tests/plastic_taub_fit.f90:
    !
    !   * the model is the backend at t = 0: the margin of the history's 0 ka slice (margin_var > margin_min,
    !     margin_position as the backend's), the topography file's bed with z_sl = 0, a cold solve; mu = c_f lambda
    !     with one c_f for every region of the region map (taub_mu, with the smoothing smooth_km), lambda on the
    !     relaxed bed z_bed + (rho_i/rho_m) H_obs (lambda_bed = "relaxed", as the backend's default) or on z_bed
    !     ("current");
    !   * fit cells: the ice cells of the margin where the observation is not floating (mask /= 3), in region
    !     fit_region of the region map (0: all); chi2 = sum over them of ((H - H_obs) / sigma)^2, sigma =
    !     max(sigma_min, sigma_rel H_obs);
    !   * a scan over H_N (HN list) x lambda_min (lam_min list): for each pair the c_f that minimises chi2, by a
    !     grid of n_grid values of ln c_f in [lcf_lo, lcf_hi] and golden section between the neighbours of the
    !     grid's minimum (to tol in ln c_f); a minimum at an end of the interval widens it by half its width, up
    !     to 4 times (reported); OpenMP over the pairs;
    !   * the reference: Gowan's tau_b at 0 ka times the domain factors (taub_domain_factor_file), smoothed with
    !     smooth_km, plastic_solve on the same margin, bed and fit cells.
    !
    ! Output <out>: per pair ln c_f, chi2/n, the RMS of H - H_obs and the volume over the fit cells, the
    ! widenings; the best pair's H, mu and tau_b, the reference's H and tau_b, H_obs, sigma, z_bed, the margin,
    ! fit and region maps; the reference's chi2/n, RMS, volume and the observed volume. Needs plastic=1:
    !
    !     make libkryos/bin/plastic_coulomb_pd.x plastic=1
    !     plastic_coulomb_pd.x <pd.nml>          (namelist &pd, below)

    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST
    use kryos_plastic_taub,   only: taub_t, taub_read_domains, taub_read_domain_factor, taub_at, taub_set_smooth, &
                                    taub_set_coulomb, taub_mu, taub_set_lambda_fixed
    use kryos_plastic_bed,    only: relaxed_bed
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_init, plastic_solve, plastic_solve_coulomb, &
                           plastic_margin_frac_from_sdist
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

    implicit none

    integer, parameter :: NMAX = 32, NREG = 5
    ! &pd
    character(len=512) :: hist = "", topo_file = "", region_file = "", taub_file = "", factor_file = "", out = ""
    character(len=64)  :: margin_var = "f_ice", margin_position = "sdist", mode = "surface"
    character(len=64)  :: bed_var = "z_bed", H_var = "H_ice", mask_var = "mask", lambda_bed = "relaxed"
    real(dp) :: margin_min = 0.5_dp, taub_default = 5000.0_dp, smooth_km = 0.0_dp
    integer  :: fit_region = 0
    real(dp) :: HN(NMAX) = -1.0_dp, lam_min(NMAX) = -1.0_dp
    real(dp) :: p = 1.0_dp, zb0 = -300.0_dp, zb1 = 200.0_dp
    real(dp) :: lcf_lo = -5.809142990314028_dp, lcf_hi = -3.506557897319982_dp      ! ln 0.003, ln 0.03
    integer  :: n_grid = 13
    real(dp) :: tol = 1.0e-3_dp, sigma_min = 50.0_dp, sigma_rel = 0.1_dp
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp    ! the kryos Earth
    namelist /pd/ hist, margin_var, margin_min, margin_position, mode, topo_file, bed_var, H_var, mask_var, &
                  region_file, fit_region, taub_file, factor_file, taub_default, smooth_km, HN, lam_min, p, zb0, zb1, &
                  lambda_bed, &
                  lcf_lo, lcf_hi, n_grid, tol, sigma_min, sigma_rel, out

    type(margin_t)      :: mg
    type(taub_t)        :: tr
    type(plastic_par_t) :: par
    type(plastic_t), allocatable :: ps(:)
    real(dp), allocatable :: x(:), y(:), zb(:,:), zsl(:,:), sd(:,:), frac(:,:), Hobs(:,:), sig(:,:), tau_ref(:,:), zrel(:,:)
    real(dp), allocatable :: Href(:,:), Hbest(:,:), mubest(:,:), taubest(:,:)
    real(dp), allocatable :: lcf(:,:), chi2n(:,:), rms(:,:), vol(:,:)
    integer,  allocatable :: mobs(:,:), reg(:,:), nwide(:,:), nsolve(:,:)
    logical,  allocatable :: mask(:,:), isfit(:,:)
    character(len=256) :: arg
    real(dp) :: dx, dy, c_ref, r_ref, v_ref, v_obs, best
    integer  :: nx, ny, nh, nl, nth, k, ih, il, ib, jb, u_nml, ios, stat, nfit, npair, q, t

    ! --- inputs --------------------------------------------------------------------------------------------------------
    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=pd, iostat=ios)
    if (ios /= 0) error stop "plastic_coulomb_pd: cannot read &pd"
    close(u_nml)
    nh = count(HN > 0.0_dp)
    nl = count(lam_min > 0.0_dp)

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zb(nx,ny), zsl(nx,ny), sd(nx,ny), Hobs(nx,ny), sig(nx,ny), mobs(nx,ny), mask(nx,ny), reg(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp

    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, 0.0_dp, mask, sd)
    call nc_read(topo_file, bed_var, zb)
    call nc_read(topo_file, H_var, Hobs)
    call nc_read(topo_file, mask_var, mobs)
    call nc_read(region_file, "region", reg)
    where (zb < -9000.0_dp) zb = 0.0_dp
    where (Hobs < 0.0_dp) Hobs = 0.0_dp
    zsl = 0.0_dp
    zrel = relaxed_bed(zb, Hobs, RHO_I)
    if (trim(lambda_bed) /= "relaxed" .and. trim(lambda_bed) /= "current") error stop "plastic_coulomb_pd: lambda_bed"
    sig = max(sigma_min, sigma_rel*Hobs)
    isfit = mask .and. mobs /= 3
    if (fit_region > 0) isfit = isfit .and. reg == fit_region
    nfit = count(isfit)
    if (count(mask .and. reg <= 0) > 0) error stop "plastic_coulomb_pd: ice cells without a region"

    par%rho_i = RHO_I;  par%rho_w = RHO_SW;  par%g = GRAV;  par%mode = trim(mode)
    if (trim(margin_position) == "sdist") then
        allocate(frac(nx, ny))
        call plastic_margin_frac_from_sdist(sd, mask, dx, dy, frac)
    else
        par%margin = trim(margin_position)
    end if
    nth = 1
    !$ nth = omp_get_max_threads()
    allocate(ps(0:nth-1))
    do k = 0, nth - 1
        call plastic_init(ps(k), nx, ny, dx, dy, par, stat=stat)
    end do

    write(*,"(a)") "# plastic_coulomb_pd: "//trim(hist)//", "//trim(topo_file)//", "//trim(region_file)
    write(*,"(a,i0,a,i0,a,i0,a,f6.1,a,i0,a)") "# margin cells ", count(mask), ", fit cells ", nfit, " (region ", &
        fit_region, "), smooth_km ", smooth_km, ", ", nth, " threads, lambda on the "//trim(lambda_bed)//" bed"

    ! --- the reference: Gowan's tau_b at 0 ka x the domain factors ------------------------------------------------------
    allocate(tau_ref(nx, ny), Href(nx, ny))
    call taub_read_domains(tr, trim(taub_file), "taub_domain", "domain_id", "taub", "taub_t", taub_default, .true., nx, ny)
    if (len_trim(factor_file) > 0) call taub_read_domain_factor(tr, trim(factor_file))
    call taub_set_smooth(tr, 1000.0_dp*smooth_km, dx, dy)
    call taub_at(tr, 0.0_dp, tau_ref)
    call plastic_solve(ps(0), zb, zsl, tau_ref, mask, warm=.false., stat=stat, margin_frac=frac)
    Href = ps(0)%H
    call stats(Href, c_ref, r_ref, v_ref)
    v_obs = sum(Hobs, mask=isfit)*dx*dy
    write(*,"(a,f10.3,a,f8.1,a,f8.4,a,f8.4,a)") "# reference (Gowan x factors): chi2/n ", c_ref, ", rms ", r_ref, &
        " m, volume ", v_ref*1.0e-15_dp, " vs observed ", v_obs*1.0e-15_dp, " 10^6 km3"

    ! --- the scan --------------------------------------------------------------------------------------------------------
    npair = nh*nl
    allocate(lcf(nh, nl), chi2n(nh, nl), rms(nh, nl), vol(nh, nl), nwide(nh, nl), nsolve(nh, nl))
    !$omp parallel do private(q, ih, il, t) schedule(dynamic)
    do q = 1, npair
        ih = mod(q - 1, nh) + 1
        il = (q - 1)/nh + 1
        t = 0
        !$ t = omp_get_thread_num()
        call fit_pair(ih, il, t)
    end do
    !$omp end parallel do

    write(*,"(a)") "#      H_N  lam_min    c_f   ln c_f    chi2/n    rms [m]  volume/obs  widened  solves"
    do il = 1, nl
        do ih = 1, nh
            write(*,"(f10.0,f9.3,f8.4,f9.4,f10.3,f11.1,f12.4,i9,i8)") HN(ih), lam_min(il), exp(lcf(ih, il)), lcf(ih, il), &
                chi2n(ih, il), rms(ih, il), vol(ih, il)/v_obs, nwide(ih, il), nsolve(ih, il)
        end do
    end do
    best = huge(1.0_dp);  ib = 1;  jb = 1
    do il = 1, nl
        do ih = 1, nh
            if (chi2n(ih, il) < best) then
                best = chi2n(ih, il);  ib = ih;  jb = il
            end if
        end do
    end do
    allocate(Hbest(nx, ny), mubest(nx, ny), taubest(nx, ny))
    call forward(ib, jb, lcf(ib, jb), 0, Hbest, mubest, taubest)
    write(*,"(a,f8.0,a,f6.3,a,f8.5,a,f9.3,a,f8.1,a,f8.4)") "# best: H_N ", HN(ib), " m, lam_min ", lam_min(jb), ", c_f ", &
        exp(lcf(ib, jb)), ", chi2/n ", chi2n(ib, jb), ", rms ", rms(ib, jb), " m, volume/obs ", vol(ib, jb)/v_obs

    ! --- output ---------------------------------------------------------------------------------------------------------
    call nc_create(trim(out))
    call nc_write_dim(trim(out), "xc", x=x, units="km")
    call nc_write_dim(trim(out), "yc", x=y, units="km")
    call nc_write_dim(trim(out), "HN", x=HN(1:nh), units="m")
    call nc_write_dim(trim(out), "lam_min", x=lam_min(1:nl), units="1")
    call nc_write_dim(trim(out), "one", x=[1.0_dp])
    call nc_write(trim(out), "lcf", lcf, dim1="HN", dim2="lam_min", units="1", long_name="ln c_f of the minimum")
    call nc_write(trim(out), "chi2n", chi2n, dim1="HN", dim2="lam_min", units="1", long_name="chi2 / n at the minimum")
    call nc_write(trim(out), "rms", rms, dim1="HN", dim2="lam_min", units="m", long_name="RMS of H - H_obs, fit cells")
    call nc_write(trim(out), "volume", vol, dim1="HN", dim2="lam_min", units="m3", long_name="ice volume, fit cells")
    call nc_write(trim(out), "widened", nwide, dim1="HN", dim2="lam_min", units="1", &
                  long_name="widenings of the ln c_f interval (minimum at a bound)")
    call nc_write(trim(out), "H_best", Hbest, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "mu_best", mubest, dim1="xc", dim2="yc", units="1")
    call nc_write(trim(out), "taub_best", taubest, dim1="xc", dim2="yc", units="Pa")
    call nc_write(trim(out), "H_ref", Href, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "taub_ref", tau_ref, dim1="xc", dim2="yc", units="Pa")
    call nc_write(trim(out), "H_obs", Hobs, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "sigma", sig, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "z_bed", zb, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "margin", merge(1, 0, mask), dim1="xc", dim2="yc")
    call nc_write(trim(out), "fit", merge(1, 0, isfit), dim1="xc", dim2="yc")
    call nc_write(trim(out), "region", reg, dim1="xc", dim2="yc")
    call nc_write(trim(out), "best_HN", [HN(ib)], dim1="one", units="m")
    call nc_write(trim(out), "best_lam_min", [lam_min(jb)], dim1="one")
    call nc_write(trim(out), "best_lcf", [lcf(ib, jb)], dim1="one")
    call nc_write(trim(out), "ref_chi2n", [c_ref], dim1="one")
    call nc_write(trim(out), "ref_rms", [r_ref], dim1="one", units="m")
    call nc_write(trim(out), "ref_volume", [v_ref], dim1="one", units="m3")
    call nc_write(trim(out), "obs_volume", [v_obs], dim1="one", units="m3")
    call nc_write(trim(out), "n_fit", [nfit], dim1="one")
    call nc_write(trim(out), "smooth_km", [smooth_km], dim1="one", units="km")
    write(*,"(a)") "# wrote "//trim(out)

contains

    subroutine forward(jh, jl, lc, tt, Hout, muo, tauo)
        ! the Coulomb solve of pair (jh, jl) with ln c_f = lc, solver instance tt
        integer,  intent(in)  :: jh, jl, tt
        real(dp), intent(in)  :: lc
        real(dp), intent(out) :: Hout(:,:), muo(:,:), tauo(:,:)
        type(taub_t) :: tb
        integer :: st
        call taub_set_coulomb(tb, reg, [(k, k = 1, NREG)], [(lc, k = 1, NREG)], p, HN(jh), lam_min(jl), zb0, zb1, &
                              trim(region_file), "uniform c_f")
        if (trim(lambda_bed) == "relaxed") call taub_set_lambda_fixed(tb, zrel)
        call taub_set_smooth(tb, 1000.0_dp*smooth_km, dx, dy)
        call taub_mu(tb, zb, zsl, muo)
        call plastic_solve_coulomb(ps(tt), zb, zsl, muo, mask, p, HN(jh), warm=.false., stat=st, margin_frac=frac)
        Hout = ps(tt)%H
        tauo = ps(tt)%tau_b
    end subroutine forward

    subroutine fit_pair(jh, jl, tt)
        ! the c_f of the minimum chi2 for pair (jh, jl): grid, widen at a bound, golden section
        integer, intent(in) :: jh, jl, tt
        real(dp), parameter :: GR = 0.6180339887498949_dp
        real(dp), allocatable :: g(:), f(:), H(:,:), mu1(:,:), mu(:,:), ta(:,:)
        real(dp) :: lo, hi, a, b, c, d, fc, fd, w, r, v
        integer  :: i, m, nw, ns
        type(taub_t) :: tb

        allocate(H(nx, ny), mu1(nx, ny), mu(nx, ny), ta(nx, ny), g(n_grid), f(n_grid))
        ! lambda alone (c_f = 1, smoothed); mu = c_f mu1 (the smoothing is linear)
        call taub_set_coulomb(tb, reg, [(i, i = 1, NREG)], [(0.0_dp, i = 1, NREG)], p, HN(jh), lam_min(jl), zb0, zb1, &
                              trim(region_file), "uniform c_f")
        if (trim(lambda_bed) == "relaxed") call taub_set_lambda_fixed(tb, zrel)
        call taub_set_smooth(tb, 1000.0_dp*smooth_km, dx, dy)
        call taub_mu(tb, zb, zsl, mu1)
        lo = lcf_lo;  hi = lcf_hi;  nw = 0;  ns = 0
        do
            do i = 1, n_grid
                g(i) = lo + (hi - lo)*(i - 1)/real(n_grid - 1, dp)
                f(i) = chi2_at(g(i), jh, tt, mu1, mu, ns)
            end do
            m = minloc(f, dim=1)
            if ((m > 1 .and. m < n_grid) .or. nw >= 4) exit
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
        fc = chi2_at(c, jh, tt, mu1, mu, ns);  fd = chi2_at(d, jh, tt, mu1, mu, ns)
        do while (b - a > tol)
            if (fc < fd) then
                b = d;  d = c;  fd = fc
                c = b - GR*(b - a);  fc = chi2_at(c, jh, tt, mu1, mu, ns)
            else
                a = c;  c = d;  fc = fd
                d = a + GR*(b - a);  fd = chi2_at(d, jh, tt, mu1, mu, ns)
            end if
        end do
        lcf(jh, jl) = merge(c, d, fc < fd)
        if (f(m) < min(fc, fd)) lcf(jh, jl) = g(m)
        mu = exp(lcf(jh, jl))*mu1
        call plastic_solve_coulomb(ps(tt), zb, zsl, mu, mask, p, HN(jh), warm=.false., stat=i, margin_frac=frac)
        H = ps(tt)%H
        call stats(H, chi2n(jh, jl), r, v)
        rms(jh, jl) = r;  vol(jh, jl) = v;  nwide(jh, jl) = nw;  nsolve(jh, jl) = ns + 1
    end subroutine fit_pair

    real(dp) function chi2_at(lc, jh, tt, mu1, mu, ns) result(c2)
        ! chi2/n of the Coulomb solve with mu = exp(lc) mu1 (H_N of pair jh, solver instance tt)
        real(dp), intent(in)    :: lc, mu1(:,:)
        integer,  intent(in)    :: jh, tt
        real(dp), intent(inout) :: mu(:,:)
        integer,  intent(inout) :: ns
        real(dp) :: rr, vv
        integer  :: st
        mu = exp(lc)*mu1
        call plastic_solve_coulomb(ps(tt), zb, zsl, mu, mask, p, HN(jh), warm=.false., stat=st, margin_frac=frac)
        call stats(ps(tt)%H, c2, rr, vv)
        ns = ns + 1
    end function chi2_at

    subroutine stats(HH, c2n, rm, vv)
        ! chi2/n, RMS of H - H_obs [m] and the volume [m3] over the fit cells
        real(dp), intent(in)  :: HH(:,:)
        real(dp), intent(out) :: c2n, rm, vv
        c2n = sum(((HH - Hobs)/sig)**2, mask=isfit)/nfit
        rm  = sqrt(sum((HH - Hobs)**2, mask=isfit)/nfit)
        vv  = sum(HH, mask=isfit)*dx*dy
    end subroutine stats

end program plastic_coulomb_pd
