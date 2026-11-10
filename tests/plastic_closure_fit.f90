program plastic_closure_fit
    ! Stage 1b of the implied mass balance of the ice backend `plastic` (docs/dev/plastic-backend.md, "Implied mass
    ! balance"): the flow closure (Glen's A and a sliding law of kryos_plastic_friction, its coefficient per cell from
    ! tau_b) against the observed present-day surface speed, on the present-day plastic geometry. Offline, with the
    ! backend's modules and the library:
    !
    !   * the geometry is that of Stage 1: the backend at t = 0, the margin of the history's 0 ka slice, Gowan's tau_b at
    !     0 ka times the per-domain factors of the present-day fit (factor_file), smoothed as the backend's
    !     taub_smooth_km, the topography file's bed, z_sl = 0;
    !   * law = "yield": each (A, uy, b) of A_list, uy_list, b_list (m fixed); law = "coulomb_n": each (A, Hc, lambda_min,
    !     m) of A_list, Hc_list, lam_list, m_list (u0, umax, p, zb0, zb1 fixed): the surface velocity, the divergence of
    !     the flux and the outflow through the margin faces of the steady face fluxes (plastic_smb_flux_divergence, the
    !     driving stress capped at tau_b as the backend's);
    !   * the misfit: over the ice cells with an observed speed > speed_min and a model speed > 0, the RMS of
    !     log10(model / observed) and the median of the ratio; the integrals [Gt yr-1] of the implied SMB (div q, the
    !     present day steady), its positive part, the outflow to land and to the ocean; the ablation area (the fraction
    !     of the ice cells with div q - q_land < 0).
    !
    ! Output <out>: a table (text, one row per parameter set, sorted by the RMS) and, for the best set, a NetCDF file
    ! with the speeds, the implied SMB and the outflow. Needs plastic=1:
    !
    !     make libkryos/bin/plastic_closure_fit.x plastic=1
    !     plastic_closure_fit.x <fit.nml>          (namelist &fit, below)
    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST
    use kryos_plastic_taub,   only: taub_t, taub_read_domains, taub_read_domain_factor, taub_at, taub_set_smooth
    use kryos_plastic_friction, only: friction_t, friction_law, friction_map
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_glen_t, plastic_sliding_t, plastic_init, plastic_solve, &
                           plastic_margin_frac_from_sdist, plastic_smb_flux_divergence, PLASTIC_OK

    implicit none

    integer, parameter :: NMAX = 32
    ! &fit
    character(len=512) :: hist = "", taub_file = "", factor_file = "", topo_file = "", obs_file = "", out = ""
    character(len=64)  :: margin_var = "f_ice", mode = "surface", bed_var = "z_bed", obs_var = "uxy_srf", mask_var = "mask"
    real(dp) :: margin_min = 0.5_dp, taub_default = 5000.0_dp, speed_min = 0.1_dp
    character(len=16)  :: law = "yield"
    real(dp) :: A_list(NMAX) = -1.0_dp, uy_list(NMAX) = -1.0_dp, b_list(NMAX) = 99.0_dp
    real(dp) :: Hc_list(NMAX) = -1.0_dp, lam_list(NMAX) = -1.0_dp, m_list(NMAX) = -1.0_dp
    real(dp) :: n_glen = 3.0_dp, m = 3.0_dp, u0 = 100.0_dp, umax = 3000.0_dp, p = 0.5_dp, zb0 = -500.0_dp, zb1 = 0.0_dp
    real(dp) :: taub_smooth_km = 50.0_dp
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp, KG_GT = 1.0e-12_dp
    integer,  parameter :: NROW = 11
    namelist /fit/ hist, margin_var, margin_min, taub_file, taub_default, factor_file, taub_smooth_km, topo_file, bed_var, &
                   mask_var, obs_file, obs_var, speed_min, mode, law, A_list, uy_list, b_list, Hc_list, lam_list, m_list, &
                   n_glen, m, u0, umax, p, zb0, zb1, out

    type(margin_t)      :: mg
    type(taub_t)        :: tb
    type(plastic_par_t) :: par
    type(plastic_t)     :: ps
    type(plastic_glen_t)    :: rheo
    type(plastic_sliding_t) :: sl
    type(friction_t)        :: fr
    real(dp), allocatable :: x(:), y(:), zb(:,:), zsl(:,:), sd(:,:), frac(:,:), tau(:,:), uobs(:,:)
    real(dp), allocatable :: divq(:,:), ql(:,:), qm(:,:), ux(:,:), uy(:,:), sp(:,:), best(:,:,:), r(:)
    real(dp), allocatable :: rows(:,:), p1(:), p2(:), p3(:)
    integer,  allocatable :: mobs(:,:), ord(:)
    logical,  allocatable :: mask(:,:), use(:,:)
    character(len=256) :: arg
    real(dp) :: dx, dy, cell, rms, med, net, acc, land, marine, abl
    integer  :: nx, ny, na, ia, i1, i2, i3, k, n, stat, u_nml, ios, u_out, kb

    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=fit, iostat=ios)
    if (ios /= 0) error stop "plastic_closure_fit: cannot read &fit"
    close(u_nml)
    na = count(A_list > 0.0_dp)
    fr%kind = trim(law)
    select case (trim(law))
    case ("yield")
        p1 = pack(uy_list, uy_list > 0.0_dp);  p2 = pack(b_list, b_list < 98.0_dp);  p3 = [m]
    case ("coulomb_n")
        p1 = pack(Hc_list, Hc_list > 0.0_dp);  p2 = pack(lam_list, lam_list > 0.0_dp);  p3 = pack(m_list, m_list > 0.0_dp)
        fr%p = p;  fr%zb0 = zb0;  fr%zb1 = zb1
    case default
        error stop "plastic_closure_fit: law must be yield or coulomb_n"
    end select

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zb(nx,ny), zsl(nx,ny), sd(nx,ny), tau(nx,ny), uobs(nx,ny), mobs(nx,ny), mask(nx,ny), frac(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp;  cell = dx*dy

    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, 0.0_dp, mask, sd)
    call taub_read_domains(tb, trim(taub_file), "taub_domain", "domain_id", "taub", "taub_t", taub_default, .true., nx, ny)
    if (len_trim(factor_file) > 0) call taub_read_domain_factor(tb, trim(factor_file))
    call taub_set_smooth(tb, 1000.0_dp*taub_smooth_km, dx, dy)
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
    allocate(sl%fric_map(nx,ny), rows(NROW, na*size(p1)*size(p2)*size(p3)))
    sl%law = friction_law(fr);  sl%u_0 = u0;  sl%u_max = umax;  rheo%n = n_glen
    k = 0
    do ia = 1, na
        do i1 = 1, size(p1)
            do i2 = 1, size(p2)
                do i3 = 1, size(p3)
                    rheo%A = A_list(ia)
                    if (trim(law) == "yield") then
                        fr%uy = p1(i1);  fr%uy_b = p2(i2)
                    else
                        fr%Hc = p1(i1);  fr%lam_min = p2(i2)
                    end if
                    sl%m = p3(i3)
                    call friction_map(fr, sl%m, tau, ps%H, zb, zsl, mask, RHO_I, RHO_SW, GRAV, sl%fric_map)
                    call plastic_smb_flux_divergence(ps, zb, zsl, mask, rheo, sl, divq, margin_frac=frac, stat=stat, &
                                                     q_land=ql, q_marine=qm, us_x=ux, us_y=uy, tau_cap=tau)
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
                    abl = real(count(mask .and. divq - ql < 0.0_dp), dp)/max(count(mask), 1)
                    k = k + 1
                    rows(:, k) = [A_list(ia), p1(i1), p2(i2), p3(i3), rms, med, net, acc, land, marine, abl]
                    if (k == 1 .or. rms < minval(rows(5, 1:k-1))) then
                        best(:,:,1) = sp;  best(:,:,2) = divq;  best(:,:,3) = ql;  best(:,:,4) = qm;  best(:,:,5) = uobs
                        kb = k
                    end if
                end do
            end do
        end do
    end do

    ord = sort_index(rows(5, 1:k))
    open(newunit=u_out, file=trim(out)//".txt", status="replace", action="write")
    if (trim(law) == "yield") then
        write(u_out, "(a)") "# law yield: A [Pa-3 yr-1], uy [m/yr at 100 kPa], b, m,"
    else
        write(u_out, "(a)") "# law coulomb_n: A [Pa-3 yr-1], Hc [m], lambda_min, m,"
    end if
    write(u_out, "(a)") "#   RMS log10(model/obs), median ratio, implied SMB net, accumulation (div q less the margin outflow, "// &
                        "> 0), outflow to land, to the ocean [Gt yr-1], ablation area fraction"
    do n = 1, k
        write(u_out, "(es10.3, f9.1, f6.2, f6.2, f8.3, f8.3, 4f10.1, f7.3)") rows(:, ord(n))
    end do
    close(u_out)
    write(*,"(a,es10.3,a,f7.1,a,f5.2,a,f5.2,a,f6.3,a,f6.2,a,f6.3)") "# best ("//trim(law)//"): A ", rows(1,kb), &
        ", p1 ", rows(2,kb), ", p2 ", rows(3,kb), ", m ", rows(4,kb), ", RMS log10 ", rows(5,kb), ", median ratio ", &
        rows(6,kb), ", ablation area ", rows(11,kb)

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
    call nc_write_dim(trim(out)//".nc", "par", x=[(real(n, dp), n = 1, NROW)])
    call nc_write(trim(out)//".nc", "best", rows(:, kb), dim1="par", &
                  long_name="A, p1 (yield: uy; coulomb_n: Hc), p2 (b; lambda_min), m, RMS log10, median ratio, net, "// &
                            "accumulation, outflow to land, to the ocean, ablation area fraction")

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
