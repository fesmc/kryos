program plastic_taub_fit
    ! Stage 1 of the tau_b calibration of the ice backend `plastic` (docs/dev/plastic-backend.md, "Calibration
    ! v2"): one factor exp(u_d) per Gowan tau_b domain d such that the plastic thickness today matches the
    ! observed present-day thickness. Offline and deterministic, with the backend's modules and the library:
    !
    !   * the model is the backend at t = 0: the margin of the history's 0 ka slice (margin_var > margin_min,
    !     margin_position as the backend's), Gowan's tau_b at 0 ka times exp(u_d) on the cells of domain d
    !     (the backend's taub_domain_factor_file), the topography file's bed with z_sl = 0, a cold solve;
    !   * fit cells: the ice cells of the margin where the observation is not floating (mask /= 3);
    !     chi2(u) = sum over them of ((H(u) - H_obs) / sigma)^2, sigma = max(sigma_min, sigma_rel H_obs);
    !   * free domains: those with fit cells; the others keep u = 0;
    !   * penalty: lambda sum over pairs of adjacent free domains of w (u_d - u_d')^2, w their shared boundary
    !     [km] on the grid (the faces between 4-neighbour cells of the two domains, anywhere in the map);
    !   * minimised by Levenberg-Marquardt on u (Gauss-Newton with a damping mu diag(A)), the Jacobian by
    !     forward differences of step du (one solve per free domain, OpenMP over the domains), |u| <= u_max,
    !     for each lambda of the list in turn, each starting from the solution of the previous one.
    !
    ! Output <out>: per lambda the log factors u (free, lambda), chi2, the roughness u' L u, the objective, the
    ! iterations and H (xc, yc, lambda); H_prior (u = 0), H_obs, sigma, the fit and margin masks, the domain map,
    ! tau_b at 0 ka, the free ids and the adjacency (edge list). Needs plastic=1:
    !
    !     make libkryos/bin/plastic_taub_fit.x plastic=1
    !     plastic_taub_fit.x <fit.nml>          (namelist &fit, below)

    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST
    use kryos_plastic_taub,   only: taub_t, taub_read_domains, taub_at
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_init, plastic_solve, plastic_margin_frac_from_sdist
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

    implicit none

    integer, parameter :: NLAM_MAX = 64
    ! &fit
    character(len=512) :: hist = "", taub_file = "", topo_file = "", out = ""
    character(len=64)  :: margin_var = "f_ice", margin_position = "sdist", mode = "surface"
    character(len=64)  :: bed_var = "z_bed", H_var = "H_ice", mask_var = "mask"
    real(dp) :: margin_min = 0.5_dp, taub_default = 5000.0_dp
    real(dp) :: lambda(NLAM_MAX) = -1.0_dp           ! the list (the first entries >= 0)
    real(dp) :: sigma_min = 50.0_dp, sigma_rel = 0.1_dp
    real(dp) :: du = 0.01_dp, u_max = 2.3_dp, tol = 1.0e-5_dp
    integer  :: n_iter = 40
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp    ! the kryos Earth
    namelist /fit/ hist, margin_var, margin_min, margin_position, taub_file, taub_default, topo_file, bed_var, H_var, &
                   mask_var, mode, lambda, sigma_min, sigma_rel, du, u_max, tol, n_iter, out

    type(margin_t)      :: mg
    type(taub_t)        :: tb
    type(plastic_par_t) :: par
    type(plastic_t), allocatable :: ps(:)
    real(dp), allocatable :: x(:), y(:), zb(:,:), zsl(:,:), sd(:,:), frac(:,:), tau0(:,:), Hobs(:,:), sig(:,:)
    real(dp), allocatable :: Hp(:,:), H(:,:), Hl(:,:,:), Jac(:,:), A(:,:), Lw(:,:), u(:), un(:), grad(:), dlt(:), r(:), rn(:)
    real(dp), allocatable :: ul(:,:), chi2l(:), roughl(:), phil(:), lams(:)
    integer,  allocatable :: mobs(:,:), fit_c(:,:), fidx(:), free_id(:), cell(:,:), itl(:), nclipl(:), ei(:), ej(:)
    real(dp), allocatable :: ew(:)
    logical,  allocatable :: mask(:,:), isfit(:,:)
    character(len=256) :: arg
    real(dp) :: dx, dy, phi, phin, chi2, mu, rough
    integer  :: nx, ny, nl, nf, nfit, nth, i, j, k, l, it, q, p, stat, ne, u_nml, ios
    logical  :: accepted

    ! --- inputs --------------------------------------------------------------------------------------------------------
    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=fit, iostat=ios)
    if (ios /= 0) error stop "plastic_taub_fit: cannot read &fit"
    close(u_nml)
    nl = count(lambda >= 0.0_dp)
    lams = lambda(1:nl)

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zb(nx,ny), zsl(nx,ny), sd(nx,ny), tau0(nx,ny), Hobs(nx,ny), sig(nx,ny), mobs(nx,ny), mask(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp

    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, 0.0_dp, mask, sd)
    call taub_read_domains(tb, trim(taub_file), "taub_domain", "domain_id", "taub", "taub_t", taub_default, .true., nx, ny)
    call taub_at(tb, 0.0_dp, tau0)
    call nc_read(topo_file, bed_var, zb)
    call nc_read(topo_file, H_var, Hobs)
    call nc_read(topo_file, mask_var, mobs)
    where (zb < -9000.0_dp) zb = 0.0_dp
    where (Hobs < 0.0_dp) Hobs = 0.0_dp
    zsl = 0.0_dp
    sig = max(sigma_min, sigma_rel*Hobs)

    ! fit cells and free domains (position in tb%ids -> free index)
    isfit = mask .and. mobs /= 3
    nfit = count(isfit)
    allocate(fidx(size(tb%ids)), cell(2, nfit), fit_c(nx, ny))
    fidx = 0;  nf = 0;  k = 0;  fit_c = 0
    do j = 1, ny
        do i = 1, nx
            if (.not. isfit(i,j)) cycle
            k = k + 1;  cell(:, k) = [i, j];  fit_c(i,j) = 1
            p = tb%idx(i,j)
            if (p > 0) then
                if (fidx(p) == 0) then
                    nf = nf + 1;  fidx(p) = nf
                end if
            end if
        end do
    end do
    allocate(free_id(nf))
    do p = 1, size(tb%ids)
        if (fidx(p) > 0) free_id(fidx(p)) = tb%ids(p)
    end do

    ! adjacency: the faces between 4-neighbour cells of two free domains [km]
    allocate(Lw(nf, nf));  Lw = 0.0_dp
    do j = 1, ny
        do i = 1, nx
            if (i < nx) call face(i, j, i + 1, j, dy*1.0e-3_dp)
            if (j < ny) call face(i, j, i, j + 1, dx*1.0e-3_dp)
        end do
    end do
    ne = 0
    do q = 1, nf
        do p = q + 1, nf
            if (Lw(p, q) > 0.0_dp) ne = ne + 1
        end do
    end do
    allocate(ei(ne), ej(ne), ew(ne))
    ne = 0
    do q = 1, nf
        do p = q + 1, nf
            if (Lw(p, q) > 0.0_dp) then
                ne = ne + 1;  ei(ne) = free_id(q);  ej(ne) = free_id(p);  ew(ne) = Lw(p, q)
            end if
        end do
    end do
    ! Lw becomes the weighted graph Laplacian: u' L u = sum over the pairs of w (u_d - u_d')^2
    do q = 1, nf
        Lw(q, q) = -sum(Lw(:, q))
    end do
    Lw = -Lw

    ! solver instances, one per thread
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

    write(*,"(a)") "# plastic_taub_fit: "//trim(hist)//", "//trim(topo_file)
    write(*,"(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)") "# margin cells ", count(mask), ", fit cells ", nfit, " (", &
        count(mask .and. mobs == 3), " floating left out), free domains ", nf, ", adjacent pairs ", ne, ", ", nth, &
        " threads, ", nl, " lambdas"

    ! --- the fit, lambda by lambda ----------------------------------------------------------------------------------------
    allocate(u(nf), un(nf), grad(nf), dlt(nf), r(nfit), rn(nfit), Jac(nfit, nf), A(nf, nf), Hp(nx, ny), H(nx, ny))
    allocate(Hl(nx, ny, nl), ul(nf, nl), chi2l(nl), roughl(nl), phil(nl), itl(nl), nclipl(nl))
    u = 0.0_dp
    call forward(u, Hp, 0)
    r = resid(Hp)
    write(*,"(a,f14.1,a,f10.3)") "# prior: chi2 ", sum(r*r), ", chi2/n ", sum(r*r)/nfit
    write(*,"(a)") "#      lambda  iter         chi2     chi2/n    roughness    objective  n_clip  u range"
    H = Hp
    do l = 1, nl
        r = resid(H)
        phi = sum(r*r) + lams(l)*dot_product(u, matmul(Lw, u))
        mu = 1.0e-3_dp
        do it = 1, n_iter
            call jacobian(u, H, Jac)
            A = matmul(transpose(Jac), Jac) + lams(l)*Lw
            grad = matmul(transpose(Jac), r) + lams(l)*matmul(Lw, u)
            accepted = .false.
            do while (mu < 1.0e8_dp)
                call damped_solve(A, mu, -grad, dlt)
                un = max(-u_max, min(u_max, u + dlt))
                call forward(un, Hp, 0)
                rn = resid(Hp)
                phin = sum(rn*rn) + lams(l)*dot_product(un, matmul(Lw, un))
                if (phin < phi) then
                    accepted = .true.;  exit
                end if
                mu = mu*4.0_dp
            end do
            if (.not. accepted) exit
            u = un;  H = Hp;  r = rn
            mu = max(mu/3.0_dp, 1.0e-6_dp)
            if ((phi - phin) < tol*phi .or. maxval(abs(dlt)) < 1.0e-4_dp) then
                phi = phin;  exit
            end if
            phi = phin
        end do
        chi2 = sum(r*r);  rough = dot_product(u, matmul(Lw, u))
        Hl(:,:,l) = H;  ul(:, l) = u;  chi2l(l) = chi2;  roughl(l) = rough;  phil(l) = phi;  itl(l) = min(it, n_iter)
        nclipl(l) = count(abs(u) >= u_max)
        write(*,"(es12.4,i6,f13.1,f11.4,es13.4,es13.4,i8,2f8.3)") lams(l), itl(l), chi2, chi2/nfit, rough, phi, nclipl(l), &
                                                                minval(u), maxval(u)
    end do

    ! --- output ---------------------------------------------------------------------------------------------------------
    call forward([(0.0_dp, k = 1, nf)], Hp, 0)
    call nc_create(trim(out))
    call nc_write_dim(trim(out), "xc", x=x, units="km")
    call nc_write_dim(trim(out), "yc", x=y, units="km")
    call nc_write_dim(trim(out), "lambda", x=lams, units="1")
    call nc_write_dim(trim(out), "free", x=[(real(k, dp), k = 1, nf)])
    call nc_write_dim(trim(out), "edge", x=[(real(k, dp), k = 1, max(ne, 1))])
    call nc_write(trim(out), "H", Hl, dim1="xc", dim2="yc", dim3="lambda", units="m")
    call nc_write(trim(out), "H_prior", Hp, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "H_obs", Hobs, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "sigma", sig, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "z_bed", zb, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out), "taub0", tau0, dim1="xc", dim2="yc", units="Pa")
    call nc_write(trim(out), "taub_domain", tb%dom, dim1="xc", dim2="yc")
    call nc_write(trim(out), "margin", merge(1, 0, mask), dim1="xc", dim2="yc")
    call nc_write(trim(out), "fit", fit_c, dim1="xc", dim2="yc")
    call nc_write(trim(out), "mask_obs", mobs, dim1="xc", dim2="yc")
    call nc_write(trim(out), "free_id", free_id, dim1="free")
    call nc_write(trim(out), "u", ul, dim1="free", dim2="lambda", units="1")
    call nc_write(trim(out), "chi2", chi2l, dim1="lambda")
    call nc_write(trim(out), "roughness", roughl, dim1="lambda", units="km")
    call nc_write(trim(out), "objective", phil, dim1="lambda")
    call nc_write(trim(out), "iterations", itl, dim1="lambda")
    call nc_write(trim(out), "n_clip", nclipl, dim1="lambda")
    if (ne > 0) then
        call nc_write(trim(out), "edge_a", ei, dim1="edge")
        call nc_write(trim(out), "edge_b", ej, dim1="edge")
        call nc_write(trim(out), "edge_w", ew, dim1="edge", units="km")
    end if
    write(*,"(a)") "# wrote "//trim(out)

contains

    subroutine face(i1, j1, i2, j2, w)
        ! a face between two cells of different free domains adds w to their pair
        integer,  intent(in) :: i1, j1, i2, j2
        real(dp), intent(in) :: w
        integer :: a, b
        if (tb%idx(i1,j1) <= 0 .or. tb%idx(i2,j2) <= 0) return
        a = fidx(tb%idx(i1,j1));  b = fidx(tb%idx(i2,j2))
        if (a == 0 .or. b == 0 .or. a == b) return
        Lw(a, b) = Lw(a, b) + w;  Lw(b, a) = Lw(b, a) + w
    end subroutine face

    subroutine forward(uu, Hout, t)
        ! the plastic thickness with tau_b(0 ka) times exp(uu) of the free domains, solver instance t
        real(dp), intent(in)  :: uu(:)
        real(dp), intent(out) :: Hout(:,:)
        integer,  intent(in)  :: t
        real(dp), allocatable :: tau(:,:)
        integer :: ii, jj, pp, st
        tau = tau0
        do jj = 1, ny
            do ii = 1, nx
                pp = tb%idx(ii,jj)
                if (pp > 0) then
                    if (fidx(pp) > 0) tau(ii,jj) = tau0(ii,jj)*exp(uu(fidx(pp)))
                end if
            end do
        end do
        call plastic_solve(ps(t), zb, zsl, tau, mask, warm=.false., stat=st, margin_frac=frac)
        Hout = ps(t)%H
    end subroutine forward

    function resid(HH) result(rr)
        ! (H - H_obs) / sigma on the fit cells
        real(dp), intent(in) :: HH(:,:)
        real(dp) :: rr(nfit)
        integer :: kk
        do kk = 1, nfit
            rr(kk) = (HH(cell(1,kk), cell(2,kk)) - Hobs(cell(1,kk), cell(2,kk)))/sig(cell(1,kk), cell(2,kk))
        end do
    end function resid

    subroutine jacobian(uu, H0, JJ)
        ! d resid / d u_d by forward differences, one solve per free domain, in parallel
        real(dp), intent(in)  :: uu(:), H0(:,:)
        real(dp), intent(out) :: JJ(:,:)
        real(dp), allocatable :: up(:), Hq(:,:)
        integer :: d, kk, t
        !$omp parallel do private(d, kk, t, up, Hq) schedule(dynamic)
        do d = 1, nf
            t = 0
            !$ t = omp_get_thread_num()
            up = uu;  up(d) = up(d) + du
            allocate(Hq(nx, ny))
            call forward(up, Hq, t)
            do kk = 1, nfit
                JJ(kk, d) = (Hq(cell(1,kk), cell(2,kk)) - H0(cell(1,kk), cell(2,kk)))/(sig(cell(1,kk), cell(2,kk))*du)
            end do
            deallocate(Hq)
        end do
        !$omp end parallel do
    end subroutine jacobian

    subroutine damped_solve(AA, m, b, xx)
        ! (AA + m diag(AA)) xx = b by Cholesky (AA symmetric positive semi-definite, diag > 0 after damping)
        real(dp), intent(in)  :: AA(:,:), m, b(:)
        real(dp), intent(out) :: xx(:)
        real(dp), allocatable :: Lc(:,:)
        integer :: n, ii, jj
        n = size(b)
        Lc = AA
        do ii = 1, n
            Lc(ii,ii) = AA(ii,ii) + m*max(AA(ii,ii), 1.0e-12_dp)
        end do
        do jj = 1, n
            Lc(jj,jj) = sqrt(max(Lc(jj,jj) - sum(Lc(jj,1:jj-1)**2), 1.0e-300_dp))
            do ii = jj + 1, n
                Lc(ii,jj) = (Lc(ii,jj) - sum(Lc(ii,1:jj-1)*Lc(jj,1:jj-1)))/Lc(jj,jj)
            end do
        end do
        do ii = 1, n
            xx(ii) = (b(ii) - sum(Lc(ii,1:ii-1)*xx(1:ii-1)))/Lc(ii,ii)
        end do
        do ii = n, 1, -1
            xx(ii) = (xx(ii) - sum(Lc(ii+1:n,ii)*xx(ii+1:n)))/Lc(ii,ii)
        end do
    end subroutine damped_solve

end program plastic_taub_fit
