program plastic_closure_fit
    ! Stage 1b of the implied mass balance of the ice backend `plastic` (docs/dev/plastic-backend.md, "Implied mass
    ! balance"): the flow closure (Glen's A and the "yield" sliding law, u_b = u_y (tau/tau_b)^m, u_y = uy (tau_b /
    ! 100 kPa)^b) against the observed present-day surface speed, on the present-day plastic geometry. Offline, with the
    ! backend's modules and the library:
    !
    !   * the geometry is that of Stage 1: the backend at t = 0, the margin of the history's 0 ka slice, Gowan's tau_b at
    !     0 ka times the per-domain factors of the present-day fit (factor_file), the topography file's bed, z_sl = 0;
    !   * for each (A, uy, b) of the lists (m fixed): the surface velocity, the divergence of the flux and the outflow
    !     through the margin faces of the steady face fluxes (plastic_smb_flux_divergence) on that geometry;
    !   * the misfit: over the ice cells with an observed speed > speed_min and a model speed > 0, the RMS of
    !     log10(model / observed) and the median of the ratio; the integrals [Gt yr-1] of the implied SMB (div q, the
    !     present day steady), its positive part, the outflow to land and to the ocean.
    !
    ! Output <out>: a table (text, one row per parameter set, sorted by the RMS) and, for the best set, a NetCDF file
    ! with the speeds, the implied SMB and the outflow. Needs plastic=1:
    !
    !     make libkryos/bin/plastic_closure_fit.x plastic=1
    !     plastic_closure_fit.x <fit.nml>          (namelist &fit, below)
    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST
    use kryos_plastic_taub,   only: taub_t, taub_read_domains, taub_read_domain_factor, taub_at
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_glen_t, plastic_sliding_t, plastic_init, plastic_solve, &
                           plastic_margin_frac_from_sdist, plastic_smb_flux_divergence, PLASTIC_OK

    implicit none

    integer, parameter :: NMAX = 32
    ! &fit
    character(len=512) :: hist = "", taub_file = "", factor_file = "", topo_file = "", obs_file = "", out = ""
    character(len=64)  :: margin_var = "f_ice", mode = "surface", bed_var = "z_bed", obs_var = "uxy_srf", mask_var = "mask"
    real(dp) :: margin_min = 0.5_dp, taub_default = 5000.0_dp, speed_min = 0.1_dp
    real(dp) :: A_list(NMAX) = -1.0_dp, uy_list(NMAX) = -1.0_dp, b_list(NMAX) = 99.0_dp
    real(dp) :: n_glen = 3.0_dp, m = 3.0_dp
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp, KG_GT = 1.0e-12_dp
    namelist /fit/ hist, margin_var, margin_min, taub_file, taub_default, factor_file, topo_file, bed_var, mask_var, &
                   obs_file, obs_var, speed_min, mode, A_list, uy_list, b_list, n_glen, m, out

    type(margin_t)      :: mg
    type(taub_t)        :: tb
    type(plastic_par_t) :: par
    type(plastic_t)     :: ps
    type(plastic_glen_t)    :: rheo
    type(plastic_sliding_t) :: sl
    real(dp), allocatable :: x(:), y(:), zb(:,:), zsl(:,:), sd(:,:), frac(:,:), tau(:,:), uobs(:,:)
    real(dp), allocatable :: divq(:,:), ql(:,:), qm(:,:), ux(:,:), uy(:,:), sp(:,:), best(:,:,:), r(:)
    real(dp), allocatable :: rows(:,:)
    integer,  allocatable :: mobs(:,:), ord(:)
    logical,  allocatable :: mask(:,:), use(:,:)
    character(len=256) :: arg
    real(dp) :: dx, dy, cell, rms, med, net, acc, land, marine
    integer  :: nx, ny, na, nu, nb, ia, iu, ib, k, n, stat, u_nml, ios, u_out, kb

    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=fit, iostat=ios)
    if (ios /= 0) error stop "plastic_closure_fit: cannot read &fit"
    close(u_nml)
    na = count(A_list > 0.0_dp);  nu = count(uy_list > 0.0_dp);  nb = count(b_list < 98.0_dp)

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zb(nx,ny), zsl(nx,ny), sd(nx,ny), tau(nx,ny), uobs(nx,ny), mobs(nx,ny), mask(nx,ny), frac(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp;  cell = dx*dy

    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, 0.0_dp, mask, sd)
    call taub_read_domains(tb, trim(taub_file), "taub_domain", "domain_id", "taub", "taub_t", taub_default, .true., nx, ny)
    if (len_trim(factor_file) > 0) call taub_read_domain_factor(tb, trim(factor_file))
    call taub_at(tb, 0.0_dp, tau)
    call nc_read(topo_file, bed_var, zb)
    call nc_read(topo_file, mask_var, mobs)
    where (zb < -9000.0_dp) zb = 0.0_dp
    call nc_read(obs_file, obs_var, uobs)
    where (.not. (uobs > 0.0_dp .and. uobs < 1.0e5_dp)) uobs = -1.0_dp
    zsl = 0.0_dp

    par%rho_i = RHO_I;  par%rho_w = RHO_SW;  par%g = GRAV;  par%mode = trim(mode)
    call plastic_margin_frac_from_sdist(sd, mask, dx, dy, frac)
    call plastic_init(ps, nx, ny, dx, dy, par, stat=stat)
    call plastic_solve(ps, zb, zsl, tau, mask, warm=.false., stat=stat, margin_frac=frac)
    write(*,"(a,i0,a,f8.3,a)") "# plastic_closure_fit: "//trim(hist)//"; ice cells ", count(mask), ", volume ", &
        RHO_I*sum(ps%H)*cell/RHO_SW/3.618e14_dp, " m SLE"

    allocate(divq(nx,ny), ql(nx,ny), qm(nx,ny), ux(nx,ny), uy(nx,ny), sp(nx,ny), best(nx,ny,5), use(nx,ny))
    allocate(sl%beta_map(nx,ny), rows(9, na*nu*nb))
    sl%law = "weertman";  sl%m = m;  rheo%n = n_glen
    k = 0
    do ia = 1, na
        do iu = 1, nu
            do ib = 1, nb
                rheo%A = A_list(ia)
                where (tau > 0.0_dp)
                    sl%beta_map = tau/(uy_list(iu)*(tau/1.0e5_dp)**b_list(ib))**(1.0_dp/m)
                elsewhere
                    sl%beta_map = huge(1.0_dp)
                end where
                call plastic_smb_flux_divergence(ps, zb, zsl, mask, rheo, sl, divq, margin_frac=frac, stat=stat, &
                                                 q_land=ql, q_marine=qm, us_x=ux, us_y=uy)
                if (stat /= PLASTIC_OK) error stop "plastic_closure_fit: plastic_smb_flux_divergence failed"
                sp = hypot(ux, uy)
                use = mask .and. mobs /= 3 .and. uobs > speed_min .and. sp > 0.0_dp
                n = count(use)
                r = log10(pack(sp, use)/pack(uobs, use))
                rms = sqrt(sum(r**2)/max(n, 1))
                med = 10.0_dp**median(r)
                net = RHO_I*sum(divq, mask = mask)*cell*KG_GT
                acc = RHO_I*sum(max(divq - ql - qm, 0.0_dp), mask = mask)*cell*KG_GT
                land = RHO_I*sum(ql)*cell*KG_GT;  marine = RHO_I*sum(qm)*cell*KG_GT
                k = k + 1
                rows(:, k) = [A_list(ia), uy_list(iu), b_list(ib), rms, med, net, acc, land, marine]
                if (k == 1 .or. rms < minval(rows(4, 1:k-1))) then
                    best(:,:,1) = sp;  best(:,:,2) = divq;  best(:,:,3) = ql;  best(:,:,4) = qm;  best(:,:,5) = uobs
                    kb = k
                end if
            end do
        end do
    end do

    ord = sort_index(rows(4, 1:k))
    open(newunit=u_out, file=trim(out)//".txt", status="replace", action="write")
    write(u_out, "(a)") "# A [Pa-3 yr-1], uy [m/yr at 100 kPa], b, RMS log10(model/obs), median ratio, implied SMB net, "// &
                        "accumulation (div q less the margin outflow, > 0), outflow to land, to the ocean [Gt yr-1]; n = "
    do n = 1, k
        write(u_out, "(es10.3, f9.1, f6.2, f8.3, f8.3, 4f10.1)") rows(:, ord(n))
    end do
    close(u_out)
    write(*,"(a,es10.3,a,f7.1,a,f5.2,a,f6.3,a,f6.2)") "# best: A ", rows(1,kb), ", uy ", rows(2,kb), ", b ", rows(3,kb), &
        ", RMS log10 ", rows(4,kb), ", median ratio ", rows(5,kb)

    call nc_create(trim(out)//".nc")
    call nc_write_dim(trim(out)//".nc", "xc", x=x, units="km")
    call nc_write_dim(trim(out)//".nc", "yc", x=y, units="km")
    call nc_write(trim(out)//".nc", "H", ps%H, dim1="xc", dim2="yc", units="m")
    call nc_write(trim(out)//".nc", "speed", best(:,:,1), dim1="xc", dim2="yc", units="m yr-1")
    call nc_write(trim(out)//".nc", "speed_obs", best(:,:,5), dim1="xc", dim2="yc", units="m yr-1")
    call nc_write(trim(out)//".nc", "div_q", best(:,:,2), dim1="xc", dim2="yc", units="m yr-1")
    call nc_write(trim(out)//".nc", "q_land", best(:,:,3), dim1="xc", dim2="yc", units="m yr-1")
    call nc_write(trim(out)//".nc", "q_marine", best(:,:,4), dim1="xc", dim2="yc", units="m yr-1")
    call nc_write(trim(out)//".nc", "tau_b", tau, dim1="xc", dim2="yc", units="Pa")
    call nc_write(trim(out)//".nc", "mask", merge(1, 0, mask), dim1="xc", dim2="yc")
    call nc_write(trim(out)//".nc", "best", rows(:, kb), dim1="par")

contains

    function median(v) result(med)
        ! the median (the lower one for an even size): quickselect on a copy
        real(dp), intent(in) :: v(:)
        real(dp) :: med
        real(dp), allocatable :: a(:)
        real(dp) :: p, t
        integer  :: lo, hi, i, j, kk
        med = 0.0_dp
        if (size(v) == 0) return
        a = v;  kk = (size(v) + 1)/2;  lo = 1;  hi = size(v)
        do while (lo < hi)
            p = a((lo + hi)/2);  i = lo;  j = hi
            do while (i <= j)
                do while (a(i) < p);  i = i + 1;  end do
                do while (a(j) > p);  j = j - 1;  end do
                if (i <= j) then
                    t = a(i);  a(i) = a(j);  a(j) = t;  i = i + 1;  j = j - 1
                end if
            end do
            if (kk <= j) then
                hi = j
            else if (kk >= i) then
                lo = i
            else
                exit
            end if
        end do
        med = a(kk)
    end function median

    function sort_index(v) result(o)
        ! the indices of v in ascending order (insertion sort: the table of parameter sets only)
        real(dp), intent(in) :: v(:)
        integer, allocatable :: o(:)
        integer :: i, j, t
        o = [(i, i = 1, size(v))]
        do i = 2, size(v)
            t = o(i);  j = i - 1
            do while (j >= 1)
                if (v(o(j)) <= v(t)) exit
                o(j + 1) = o(j);  j = j - 1
            end do
            o(j + 1) = t
        end do
    end function sort_index

end program plastic_closure_fit
