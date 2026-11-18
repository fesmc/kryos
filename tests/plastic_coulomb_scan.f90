program plastic_coulomb_scan
    ! A parameter scan of the Coulomb yield law of the ice backend `plastic` (taub_method = "coulomb") at one time of a
    ! margin history, without an Earth: the margin of the slice at time_ka, the history's bed of that slice (bed_var;
    ! PaleoMIST: relative to the sea level of the time, so z_sl = 0), mu = c_f lambda with one c_f for every region
    ! of the region map (taub_mu, smoothing smooth_km; lambda_bed = "file": lambda fixed on the bed z_lambda of
    ! lambda_file, the backend's default and the same file ({grid_name} is grid_name); "relaxed": on the present-day
    ! bed of topo_file relaxed by local isostasy, fixed; "current": on the slice's bed), a cold
    ! plastic_solve_coulomb for each (c_f, H_N, lambda_min) of
    ! the lists cf x HN x lam_min (OpenMP over the cases). Output <out>: H (xc, yc, case) with case_cf, case_HN,
    ! case_lam_min (case = cf fastest, then HN, then lam_min), the bed, the
    ! margin and the region map (volumes per ice sheet: ice-earth-recon/tools/plastic_coulomb_lgm.jl). Needs plastic=1:
    !
    !     make libkryos/bin/plastic_coulomb_scan.x plastic=1
    !     plastic_coulomb_scan.x <scan.nml>          (namelist &scan, below)

    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, MARGIN_SDIST, model_time_axis
    use kryos_plastic_taub,   only: taub_t, taub_set_coulomb, taub_set_smooth, taub_mu, taub_set_lambda_fixed, &
                                    taub_read_lambda_bed
    use kryos_plastic_bed,    only: relaxed_bed
    use nml,                  only: nml_replace
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_init, plastic_solve_coulomb, plastic_coulomb_law_t, &
                           plastic_margin_frac_from_sdist
    !$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num

    implicit none

    integer, parameter :: NMAX = 32
    ! &scan
    character(len=512) :: hist = "", region_file = "", out = "", topo_file = ""
    character(len=512) :: lambda_file = "kryos_data/plastic/{grid_name}_lambda_bed.nc"
    character(len=64)  :: margin_var = "f_ice", margin_position = "sdist", bed_var = "z_bed", lambda_bed = "file"
    character(len=64)  :: grid_name = ""
    character(len=64)  :: topo_bed_var = "z_bed", topo_H_var = "H_ice"
    real(dp) :: time_ka = 20.0_dp, margin_min = 0.5_dp, smooth_km = 0.0_dp
    real(dp) :: cf(NMAX) = -1.0_dp, HN(NMAX) = -1.0_dp, lam_min(NMAX) = -1.0_dp
    real(dp) :: p = 1.0_dp, delta = 0.05_dp, zb0 = -300.0_dp, zb1 = 200.0_dp
    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, GRAV = 9.81_dp    ! the kryos Earth
    namelist /scan/ hist, margin_var, margin_min, margin_position, bed_var, time_ka, region_file, smooth_km, cf, HN, &
                    lam_min, p, delta, zb0, zb1, lambda_bed, lambda_file, grid_name, topo_file, topo_bed_var, topo_H_var, out

    type(margin_t)      :: mg
    type(plastic_par_t) :: par
    type(plastic_t), allocatable :: ps(:)
    real(dp), allocatable :: x(:), y(:), t(:), zb(:,:), zsl(:,:), sd(:,:), frac(:,:), zbt(:,:,:), H(:,:,:), zfix(:,:)
    real(dp), allocatable :: zpd(:,:), hpd(:,:)
    integer,  allocatable :: reg(:,:), nsw(:,:,:), rid(:)
    logical,  allocatable :: mask(:,:)
    character(len=256) :: arg, units
    real(dp) :: dx, dy
    integer  :: nx, ny, nc, nh, nl, nth, k, it, q, ic, ih, il, th, u_nml, ios, stat

    call get_command_argument(1, arg)
    open(newunit=u_nml, file=trim(arg), status="old", action="read")
    read(u_nml, nml=scan, iostat=ios)
    if (ios /= 0) error stop "plastic_coulomb_scan: cannot read &scan"
    close(u_nml)
    nc = count(cf > 0.0_dp);  nh = count(HN > 0.0_dp);  nl = count(lam_min > 0.0_dp)

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc")
    allocate(x(nx), y(ny), zb(nx,ny), zsl(nx,ny), sd(nx,ny), mask(nx,ny), reg(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp
    call margin_init(mg, trim(hist), trim(margin_var), margin_min, nx, ny, dx, dy, MARGIN_SDIST, .false.)
    call margin_at(mg, -1000.0_dp*time_ka, mask, sd)
    ! the bed of the slice at time_ka (the history's own time axis, "years ago")
    allocate(t(nc_size(hist, "time")))
    allocate(zbt(nx, ny, size(t)))
    call nc_read(hist, "time", t)
    units = "years ago"
    t = model_time_axis(t, units)
    it = minloc(abs(t + 1000.0_dp*time_ka), dim=1)
    if (abs(t(it) + 1000.0_dp*time_ka) > 1.0_dp) error stop "plastic_coulomb_scan: no slice at time_ka"
    call nc_read(hist, bed_var, zbt)
    zb = zbt(:, :, it)
    where (zb < -9000.0_dp) zb = 0.0_dp
    zsl = 0.0_dp
    call nc_read(region_file, "region", reg)
    allocate(rid(nc_size(region_file, "region")))
    call nc_read(region_file, "region_id", rid)             ! the region ids of the map (one uniform c_f for all)
    select case (trim(lambda_bed))
    case ("file")
        zfix = taub_read_lambda_bed(trim(lambda_path()), nx, ny)
    case ("relaxed")
        allocate(zpd(nx, ny), hpd(nx, ny))
        call nc_read(topo_file, topo_bed_var, zpd)
        call nc_read(topo_file, topo_H_var, hpd)
        where (zpd < -9000.0_dp) zpd = 0.0_dp
        zfix = relaxed_bed(zpd, hpd, RHO_I)
    case ("current")
    case default
        error stop "plastic_coulomb_scan: lambda_bed must be file, relaxed or current"
    end select
    if (count(mask .and. reg <= 0) > 0) error stop "plastic_coulomb_scan: ice cells without a region"

    par%rho_i = RHO_I;  par%rho_w = RHO_SW;  par%g = GRAV;  par%mode = "surface"
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
    write(*,"(a,f6.1,a,i0,a,i0,a,i0,a)") "# plastic_coulomb_scan: "//trim(hist)//" at ", time_ka, " ka: ", count(mask), &
        " ice cells, ", nc*nh*nl, " cases, ", nth, " threads"

    allocate(H(nx, ny, nc*nh*nl), nsw(nc, nh, nl))
    !$omp parallel do private(q, ic, ih, il, th) schedule(dynamic)
    do q = 1, nc*nh*nl
        ic = mod(q - 1, nc) + 1
        ih = mod((q - 1)/nc, nh) + 1
        il = (q - 1)/(nc*nh) + 1
        th = 0
        !$ th = omp_get_thread_num()
        call case_solve(ic, ih, il, th)
    end do
    !$omp end parallel do
    write(*,"(a,i0,a,i0)") "# sweeps per solve: ", minval(nsw), " to ", maxval(nsw)

    call nc_create(trim(out))
    call nc_write_dim(trim(out), "xc", x=x, units="km")
    call nc_write_dim(trim(out), "yc", x=y, units="km")
    call nc_write_dim(trim(out), "case", x=[(real(q, dp), q = 1, nc*nh*nl)])
    call nc_write(trim(out), "case_cf", [(((cf(ic), ic = 1, nc), ih = 1, nh), il = 1, nl)], dim1="case", units="1")
    call nc_write(trim(out), "case_HN", [(((HN(ih), ic = 1, nc), ih = 1, nh), il = 1, nl)], dim1="case", units="m")
    call nc_write(trim(out), "case_lam_min", [(((lam_min(il), ic = 1, nc), ih = 1, nh), il = 1, nl)], dim1="case", &
                  units="1")
    call nc_write(trim(out), "H", H, dim1="xc", dim2="yc", dim3="case", units="m")
    call nc_write(trim(out), "z_bed", zb, dim1="xc", dim2="yc", units="m", &
                  long_name="bed of the slice, relative to the sea level of the time")
    call nc_write(trim(out), "margin", merge(1, 0, mask), dim1="xc", dim2="yc")
    call nc_write(trim(out), "region", reg, dim1="xc", dim2="yc")
    write(*,"(a)") "# wrote "//trim(out)

contains

    subroutine case_solve(jc, jh, jl, tt)
        integer, intent(in) :: jc, jh, jl, tt
        type(taub_t) :: tb
        real(dp), allocatable :: mu(:,:)
        integer :: i, st
        allocate(mu(nx, ny))
        call taub_set_coulomb(tb, reg, rid, [(log(cf(jc)), i = 1, size(rid))], p, HN(jh), delta, lam_min(jl), zb0, &
                              zb1, trim(region_file), "uniform c_f")
        if (allocated(zfix)) call taub_set_lambda_fixed(tb, zfix)
        call taub_set_smooth(tb, 1000.0_dp*smooth_km, dx, dy)
        call taub_mu(tb, zb, zsl, mu)
        call plastic_solve_coulomb(ps(tt), zb, zsl, mu, mask, law(jh), warm=.false., stat=st, margin_frac=frac)
        if (st /= 0) write(*,"(a,3i4,a)") "# case ", jc, jh, jl, ": "//trim(ps(tt)%msg)
        H(:, :, jc + nc*(jh - 1) + nc*nh*(jl - 1)) = ps(tt)%H
        nsw(jc, jh, jl) = ps(tt)%n_sweeps
    end subroutine case_solve


    function lambda_path() result(f)
        ! lambda_file with {grid_name} replaced by grid_name
        character(len=512) :: f
        f = lambda_file
        if (index(f, "{grid_name}") > 0) then
            if (len_trim(grid_name) == 0) error stop "plastic_coulomb_scan: lambda_file has {grid_name}: set grid_name"
            call nml_replace(f, "{grid_name}", trim(grid_name))
        end if
    end function lambda_path

    function law(jh) result(lw)
        ! The library's Coulomb law of H_N number jh (p and delta of the namelist).
        integer, intent(in) :: jh
        type(plastic_coulomb_law_t) :: lw
        lw = plastic_coulomb_law_t(p=p, H_N=HN(jh), delta=delta)
    end function law

end program plastic_coulomb_scan
