program plastic_slices
    ! The offline check of the ice backend `plastic` against PaleoMIST (docs/dev/plastic-backend.md): at each
    ! slice of the ice history, the plastic thickness on PaleoMIST's own bed of that slice (z_bed is relative to
    ! the sea level of the time, so z_sl = 0), with its margin (the slice of <margin_var> > <margin_min>, as the
    ! backend's margin_var and margin_H_min: default f_ice > 0.5) and Gowan's tau_b at that time (the product's
    ! series), cold, against the history's H_ice (whatever the margin variable): rms and bias over the margin
    ! (where the plastic ice is), the volume [m SLE] of both (PaleoMIST's: all its H_ice), and the area of
    ! PaleoMIST's ice above 0, 1 and 10 m. The same modules as the
    ! backend (kryos_plastic_margin, kryos_plastic_taub) and the library; no run. Also the stress residual
    ! (plastic_stress_residual: median and max over the ice). With <out.nc>: H and PaleoMIST's H_ice of every
    ! slice (xc, yc, time [ka]), for maps. Needs plastic=1:
    !
    !     make libkryos/bin/plastic_slices.x plastic=1
    !     plastic_slices.x <history.nc> <taub.nc> <mode> [<margin_var> [<margin_min> [<taub_default> [<out.nc>]]]]

    use precision, only: dp
    use ncio,      only: nc_read, nc_size, nc_read_attr, nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_init, margin_at, model_time_axis, MARGIN_NEAREST
    use kryos_plastic_taub,   only: taub_t, taub_read_domains, taub_at
    use plastic_ice, only: plastic_par_t, plastic_t, plastic_init, plastic_solve, plastic_stress_residual

    implicit none

    real(dp), parameter :: RHO_I = 910.0_dp, RHO_SW = 1028.0_dp, G = 9.81_dp, A_OCN = 3.618e14_dp   ! kryos Earth
    character(len=512) :: hist, tfile, arg, outf
    character(len=64)  :: mvar
    character(len=16)  :: mode
    character(len=256) :: units
    type(margin_t)      :: mg
    type(taub_t)        :: tb
    type(plastic_par_t) :: par
    type(plastic_t)     :: ps
    real(dp), allocatable :: x(:), y(:), t(:), zb(:,:), Hp(:,:), tau(:,:), sd(:,:), zsl(:,:), d(:,:)
    logical,  allocatable :: mask(:,:)
    real(dp) :: m_min, tdef, dx, dy, ca, rms, bias, sle, sle_pm, a0, a1, a10, rmax, rmed
    real(dp), allocatable :: r(:,:), rr(:)
    integer  :: nx, ny, nt, k, stat, n

    call get_command_argument(1, hist)
    call get_command_argument(2, tfile)
    call get_command_argument(3, mode)
    mvar = "f_ice";  m_min = 0.5_dp;  tdef = 5000.0_dp
    if (command_argument_count() >= 4) call get_command_argument(4, mvar)
    if (command_argument_count() >= 5) then
        call get_command_argument(5, arg);  read(arg, *) m_min
    end if
    if (command_argument_count() >= 6) then
        call get_command_argument(6, arg);  read(arg, *) tdef
    end if
    outf = ""
    if (command_argument_count() >= 7) call get_command_argument(7, outf)

    nx = nc_size(hist, "xc");  ny = nc_size(hist, "yc");  nt = nc_size(hist, "time")
    allocate(x(nx), y(ny), t(nt), zb(nx,ny), Hp(nx,ny), tau(nx,ny), sd(nx,ny), zsl(nx,ny), d(nx,ny), mask(nx,ny))
    call nc_read(hist, "xc", x);  call nc_read(hist, "yc", y);  call nc_read(hist, "time", t)
    units = "";  call nc_read_attr(hist, "time", "units", units)
    t = model_time_axis(t, units)
    dx = (x(2) - x(1))*1.0e3_dp;  dy = (y(2) - y(1))*1.0e3_dp;  ca = dx*dy

    call margin_init(mg, trim(hist), trim(mvar), m_min, nx, ny, dx, dy, MARGIN_NEAREST, .false.)
    call taub_read_domains(tb, trim(tfile), "taub_domain", "domain_id", "taub", "taub_t", tdef, .false., nx, ny)
    par%rho_i = RHO_I;  par%rho_w = RHO_SW;  par%g = G;  par%mode = trim(mode)
    call plastic_init(ps, nx, ny, dx, dy, par, stat=stat)
    zsl = 0.0_dp
    allocate(r(nx, ny))
    if (len_trim(outf) > 0) then
        call nc_create(trim(outf))
        call nc_write_dim(trim(outf), "xc", x=x, units="km")
        call nc_write_dim(trim(outf), "yc", x=y, units="km")
        call nc_write_dim(trim(outf), "time", x=-t/1.0e3_dp, units="ka")
    end if

    write(*,"(a,a,a,a,a,a,g0,a,f8.1,a)") "# ", trim(hist), ", mode ", trim(mode), ", margin ", trim(mvar)//" > ", m_min, &
                                        ", taub_default ", tdef, " Pa"
    write(*,"(a)") "#  t [ka]  rms [m] bias [m]  SLE [m] SLE_pm [m]  A0     A1    A10 [1e6 km2]  stat n_outer n_bad"// &
                   "  res_med res_max"
    do k = 1, nt
        call nc_read(hist, "z_bed", zb, start=[1,1,k], count=[nx,ny,1])
        call nc_read(hist, "H_ice", Hp, start=[1,1,k], count=[nx,ny,1])
        where (Hp < 0.0_dp) Hp = 0.0_dp                  ! fill value
        call margin_at(mg, t(k), mask, sd)
        n = count(mask .and. zb < -9000.0_dp)            ! bed missing under the ice
        where (zb < -9000.0_dp) zb = 0.0_dp
        call taub_at(tb, t(k), tau)
        call plastic_solve(ps, zb, zsl, tau, mask, warm=.false., stat=stat)
        d = ps%H - Hp
        rms  = sqrt(sum(d*d, mask=mask)/max(count(mask), 1))
        bias = sum(d, mask=mask)/max(count(mask), 1)
        sle    = sum(ps%H)*ca*RHO_I/RHO_SW/A_OCN
        sle_pm = sum(Hp)*ca*RHO_I/RHO_SW/A_OCN
        a0 = count(Hp > 0.0_dp)*ca*1.0e-12_dp;  a1 = count(Hp > 1.0_dp)*ca*1.0e-12_dp;  a10 = count(Hp > 10.0_dp)*ca*1.0e-12_dp
        call plastic_stress_residual(ps, zb, zsl, tau, mask, rmax, r=r)
        rr = pack(r, mask .and. ps%H > 0.0_dp)
        rmed = 0.0_dp
        if (size(rr) > 0) rmed = sorted_mid(rr)
        write(*,"(f9.1,2f9.1,2f10.2,3f7.3,i9,i8,i6,2f9.3)") -t(k)/1.0e3_dp, rms, bias, sle, sle_pm, a0, a1, a10, stat, &
                                                             ps%n_outer_used, n, rmed, rmax
        if (len_trim(outf) > 0) then
            call nc_write(trim(outf), "H", ps%H, dim1="xc", dim2="yc", dim3="time", start=[1,1,k], count=[nx,ny,1], units="m")
            call nc_write(trim(outf), "H_pm", Hp, dim1="xc", dim2="yc", dim3="time", start=[1,1,k], count=[nx,ny,1], units="m")
            call nc_write(trim(outf), "z_bed", zb, dim1="xc", dim2="yc", dim3="time", start=[1,1,k], count=[nx,ny,1], units="m")
        end if
    end do

contains

    function sorted_mid(v) result(m)
        ! The median: shell sort, then the middle.
        real(dp), intent(in) :: v(:)
        real(dp) :: m
        real(dp), allocatable :: a(:)
        real(dp) :: x
        integer :: n, gap, i, j
        a = v;  n = size(a);  gap = n/2
        do while (gap > 0)
            do i = gap + 1, n
                x = a(i);  j = i
                do while (j > gap)
                    if (a(j - gap) <= x) exit
                    a(j) = a(j - gap);  j = j - gap
                end do
                a(j) = x
            end do
            gap = gap/2
        end do
        m = a((n + 1)/2)
        if (mod(n, 2) == 0) m = 0.5_dp*(a(n/2) + a(n/2 + 1))
    end function sorted_mid

end program plastic_slices
