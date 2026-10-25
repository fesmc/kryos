program test_climate_snapshots
    ! climate = snapshots (make check) on a synthetic 4 x 3 grid: a reference
    ! snapshot (tas, pr, zs, to, so) and two snapshots forming the term `lgm`
    ! (coefficient alpha, tas and pr), the uniform ocean term `ocean`
    ! (derived coefficient dTo = 0.25 x the mean sea-level air anomaly). The
    ! exports are compared with an independent evaluation of the blend, the
    ! lapse-rate inflation, the pr unit conversion and the derived coefficient;
    ! the output file is written and the restore path runs.

    use precision,        only: sp, dp
    use ncio
    use coordinates,      only: grid_class, grid_init
    use regions,          only: regions_init
    use phys_constants,   only: sec_year_tropical
    use fesm_fields,      only: fields_get, fields_put
    use kryos_component,  only: init_ctx_t, step_time_t, census_t, PASS_BOUNDARY, PASS_RESTART
    use kryos_climate,    only: climate_slot_t, climate_create, climate_init, climate_init_state, &
                                climate_update, climate_write_output, climate_end

    implicit none

    integer, parameter :: nx = 4, ny = 3
    character(len=*), parameter :: fdat = "test_snap_data.nc", fpar = "test_snap.nml", fout = "test_snap_out.nc"
    real(dp), parameter :: alpha = 0.5_dp, l1 = 0.0080_dp, l2 = 0.0065_dp, beta = 0.05_dp
    real(dp), parameter :: pi = 3.14159265358979323846_dp
    real(dp), parameter :: dz_in = 50.0_dp

    type(climate_slot_t) :: cl
    type(init_ctx_t)     :: ctx
    type(census_t)       :: census
    type(step_time_t)    :: st
    type(grid_class)     :: grid
    real(sp) :: zs(nx,ny), zin(nx,ny), v0
    real(dp) :: lap(12), tas_e, pr_e, dto_e, scale
    real(dp) :: tas(nx,ny,12), pr(nx,ny,12), dto, to(nx,ny,23), l(12)
    real(sp), allocatable :: tmp3(:,:,:)
    real(dp), allocatable :: d3(:,:,:)
    integer :: i, j, m, nfail

    nfail = 0
    call write_data()
    call write_par()
    call grid_init(grid, name="TEST-SNAP", mtype="cartesian", units="kilometers", &
                   x0=0.0_dp, dx=10.0_dp, nx=nx, y0=0.0_dp, dy=10.0_dp, ny=ny)

    ctx%par_file  = fpar
    ctx%instance  = "clim"
    ctx%component = "climate"
    ctx%backend   = "snapshots"
    ctx%scope     = "regional"
    ctx%domain    = "test"
    ctx%grid_name = "TEST-SNAP"
    ctx%grid      = grid
    ctx%time_rel  = 0.0_dp
    ctx%sec_year  = sec_year_tropical          ! the core passes the run calendar's year
    call regions_init(ctx%reg, nx=nx, ny=ny)

    call climate_create(cl, ctx, census)
    call climate_init(cl, ctx)

    call check("caps: ocean profiles", cl%c%caps%has_ocn)
    call check("coefficients named by the terms", size(cl%c%coef_names) == 2 .and. size(cl%c%hi_coef) == 2)
    call check("derived coefficient dTo produced", size(cl%c%derived_names) == 1 .and. trim(cl%c%derived_names(1)) == "dTo")
    call check("depth axis declared", allocated(cl%c%axes))
    do i = 1, size(cl%c%axes)
        if (trim(cl%c%axes(i)%name) == "depth") call check("depth axis: 23 levels to 3000 m", &
            size(cl%c%axes(i)%v) == 23 .and. cl%c%axes(i)%v(23) == 3000.0_dp)
    end do

    ! the hub geometry at the climate grid (z_srf of the reference + 50 m) and the coefficients
    do j = 1, ny
        do i = 1, nx
            zs(i,j)  = real(500.0 + 10.0*i, sp)
            zin(i,j) = zs(i,j) + real(dz_in, sp)
        end do
    end do
    call fields_put(cl%c%imports, cl%c%hi_z_srf, zin, t1=0.0_dp)
    call put_coefs(0.0_sp, 0.0_sp)

    st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call climate_init_state(cl, st, PASS_BOUNDARY)
    call get_dTo(v0)
    call check("alpha = 0, dTo = 0: derived anomaly zero", abs(v0) < 1.0e-5_sp)

    ! step: alpha = 0.5, dTo as last set (0): the blend and the derived coefficient
    call put_coefs(real(alpha, sp), 0.0_sp)
    st = step_time_t(t0=0.0_dp, t1=10.0_dp, t_eval=10.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call climate_update(cl, st)

    scale = sec_year_tropical / 86400.0_dp
    do m = 1, 12
        l(m) = l1 + (l1 - l2)*cos(2.0_dp*pi*(m*30.0_dp - 30.0_dp)/360.0_dp)
    end do
    lap = l
    do j = 1, ny
        do i = 1, nx
            do m = 1, 12
                call expected(i, m, tas(i,j,m), pr(i,j,m))
            end do
        end do
    end do
    dto = 0.25_dp*sum(tas_mean_anom())/real(nx*ny, dp)

    allocate(d3(nx,ny,12))
    call fields_get(cl%c%exports, cl%c%h_tas, d3)
    call check("tas: blend and lapse inflation", maxval(abs(d3 - tas)) < 2.0e-3_dp)
    call fields_get(cl%c%exports, cl%c%h_pr, d3)
    call check("pr [kg m-2 yr-1] = blend x sec_year/86400", maxval(abs(d3/scale - pr)/pr) < 1.0e-4_dp)
    call fields_get(cl%c%exports, cl%c%h_lapse, d3)
    call check("lapse rate, uniform", maxval(abs(d3(1,1,:) - lap)) < 1.0e-8_dp .and. all(d3(:,:,3) == d3(1,1,3)))
    call fields_get(cl%c%exports, cl%c%h_z_srf_clim, zin)
    call check("z_srf_clim = the import", all(zin == zs + real(dz_in, sp)))
    call get_dTo(v0)
    call check("dTo = 0.25 x mean sea-level air anomaly", abs(real(v0, dp) - dto) < 1.0e-3_dp)

    ! ocean: profiles + the dTo set in this very update (terms_eval(derived) then terms_apply_ocn,
    ! as yelmox_climate_update), dto_ann = now - ref
    allocate(tmp3(nx,ny,23))
    call fields_get(cl%c%exports, cl%c%h_to_ann, tmp3)
    call check("ocean: reference profile + the dTo of this update", abs(tmp3(2,2,1) - (275.5_sp + v0)) < 1.0e-4_sp)
    call fields_get(cl%c%exports, cl%c%h_dto_ann, tmp3)
    call check("dto_ann = dTo of this update", maxval(abs(tmp3 - v0)) < 1.0e-4_sp)

    ! the next update reads dTo as last set (the import = the producer's value); the anomaly does not depend on it
    call put_coefs(real(alpha, sp), v0)
    st = step_time_t(t0=10.0_dp, t1=20.0_dp, t_eval=20.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call climate_update(cl, st)
    call fields_get(cl%c%exports, cl%c%h_dto_ann, tmp3)
    call check("dto_ann = dTo of the import", maxval(abs(tmp3 - v0)) < 1.0e-4_sp)
    call fields_get(cl%c%exports, cl%c%h_to_ann, tmp3)
    call check("ocean: profile + dTo", abs(tmp3(2,2,1) - (275.5_sp + v0)) < 1.0e-4_sp)

    ! restore path and output
    st = step_time_t(t0=20.0_dp, t1=20.0_dp, t_eval=20.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call climate_init_state(cl, st, PASS_RESTART)
    call climate_write_output(cl, fout, 0.0_dp, "2D", .true.)
    call climate_write_output(cl, fout, 10.0_dp, "2D", .false.)
    call check("output file written", nc_exists_var(fout, "tas") .and. nc_exists_var(fout, "term_lgm"))

    call climate_end(cl)
    call execute_command_line("rm -f "//fdat//" "//fpar//" "//fout)

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_climate_snapshots: ", nfail, " FAILED"
        error stop 1
    end if
    write(*,'(a)') "test_climate_snapshots: OK"

contains

    function tsl_of(tas0, zsurf, m) result(t)
        real(dp), intent(in) :: tas0, zsurf
        integer,  intent(in) :: m
        real(dp) :: t
        t = tas0 + zsurf*l(m)
    end function tsl_of

    subroutine expected(i, m, tas_x, pr_x)
        ! The blend of the snapshots (data of write_data), inflated to z_srf + dz_in.
        integer,  intent(in)  :: i, m
        real(dp), intent(out) :: tas_x, pr_x
        real(dp) :: z, tpd, tctl, talt, ppd, pctl, palt, tsl, prc
        z    = 500.0_dp + 10.0_dp*i
        tpd  = 260.0_dp + 0.5_dp*m + 0.1_dp*i
        tctl = 262.0_dp + 0.5_dp*m + 0.1_dp*i
        talt = tctl - 8.0_dp + 0.02_dp*i
        ppd  = 2.0_dp + 0.1_dp*m
        pctl = ppd
        palt = 0.8_dp*pctl
        tsl  = tsl_of(tpd, z, m) + alpha*(tsl_of(talt, z, m) - tsl_of(tctl, z, m))
        tas_x = tsl - (z + dz_in)*l(m)
        prc  = (ppd/exp(beta*(-z*l(m)))) * (alpha*((palt/exp(beta*(-z*l(m))))/(pctl/exp(beta*(-z*l(m)))) - 1.0_dp) + 1.0_dp)
        pr_x = prc*exp(beta*(tas_x - tsl))
    end subroutine expected

    function tas_mean_anom() result(a)
        ! Per cell: the mean over the months of tsl_now - tsl_pd = alpha*(tsl_alt - tsl_ctl),
        ! the air anomaly of the update (the reference's annual mean tsl_ann = ta_ann + l1*zs
        ! equals the mean of its monthly tsl, the lapse cycle having zero mean).
        real(dp) :: a(nx*ny)
        integer :: ii, mm, k
        real(dp) :: z, tctl, talt, s
        k = 0
        do j = 1, ny
            do ii = 1, nx
                k = k + 1
                z = 500.0_dp + 10.0_dp*ii
                s = 0.0_dp
                do mm = 1, 12
                    tctl = 262.0_dp + 0.5_dp*mm + 0.1_dp*ii
                    talt = tctl - 8.0_dp + 0.02_dp*ii
                    s = s + alpha*(tsl_of(talt, z, mm) - tsl_of(tctl, z, mm))
                end do
                a(k) = s/12.0_dp
            end do
        end do
    end function tas_mean_anom

    subroutine put_coefs(a, d)
        real(sp), intent(in) :: a, d
        integer :: k
        do k = 1, size(cl%c%coef_names)
            if (trim(cl%c%coef_names(k)) == "alpha") call fields_put(cl%c%imports, cl%c%hi_coef(k), a, t1=0.0_dp)
            if (trim(cl%c%coef_names(k)) == "dTo")   call fields_put(cl%c%imports, cl%c%hi_coef(k), d, t1=0.0_dp)
        end do
    end subroutine put_coefs

    subroutine get_dTo(x)
        real(sp), intent(out) :: x
        call fields_get(cl%c%exports, cl%c%h_derived(1), x)
    end subroutine get_dTo

    subroutine check(what, ok)
        character(len=*), intent(in) :: what
        logical,          intent(in) :: ok
        write(*,'(a,a)') merge("  ok   ", "  FAIL ", ok), what
        if (.not. ok) nfail = nfail + 1
    end subroutine check

    subroutine write_data()
        real(sp) :: t(nx,ny,12), p(nx,ny,12), z(nx,ny), o(nx,ny,4), s(nx,ny,4)
        integer :: ii, jj, mm
        do jj = 1, ny
            do ii = 1, nx
                z(ii,jj) = real(500.0 + 10.0*ii, sp)
                do mm = 1, 12
                    t(ii,jj,mm) = real(260.0 + 0.5*mm + 0.1*ii, sp)
                    p(ii,jj,mm) = real(2.0 + 0.1*mm, sp)
                end do
                o(ii,jj,:) = [275.5_sp, 274.5_sp, 273.5_sp, 272.5_sp]
                s(ii,jj,:) = 34.5_sp
            end do
        end do
        call nc_create(fdat)
        call nc_write_dim(fdat, "xc",    x=0.0_dp, dx=10.0_dp, nx=nx)
        call nc_write_dim(fdat, "yc",    x=0.0_dp, dx=10.0_dp, nx=ny)
        call nc_write_dim(fdat, "month", x=1.0_dp, dx=1.0_dp, nx=12)
        call nc_write_dim(fdat, "depth", x=[0.0_dp, 100.0_dp, 1000.0_dp, 3000.0_dp])
        call nc_write(fdat, "t_pd",  t,        dim1="xc", dim2="yc", dim3="month")
        call nc_write(fdat, "p_pd",  p,        dim1="xc", dim2="yc", dim3="month")
        call nc_write(fdat, "t_ctl", t + 2.0_sp, dim1="xc", dim2="yc", dim3="month")
        call nc_write(fdat, "p_ctl", p,        dim1="xc", dim2="yc", dim3="month")
        do ii = 1, nx
            t(ii,:,:) = t(ii,:,:) + 2.0_sp - 8.0_sp + real(0.02*ii, sp)
        end do
        call nc_write(fdat, "t_alt", t,        dim1="xc", dim2="yc", dim3="month")
        call nc_write(fdat, "p_alt", 0.8_sp*p, dim1="xc", dim2="yc", dim3="month")
        call nc_write(fdat, "zs",    z,        dim1="xc", dim2="yc")
        call nc_write(fdat, "to",    o,        dim1="xc", dim2="yc", dim3="depth")
        call nc_write(fdat, "so",    s,        dim1="xc", dim2="yc", dim3="depth")
    end subroutine write_data

    subroutine write_par()
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,'(a)') "&clim"
        write(u,'(a)') "    backend   = 'snapshots'"
        write(u,'(a)') "    reference = 'pd'"
        write(u,'(a)') "    terms     = 'lgm' 'ocean'"
        write(u,'(a)') "/"
        write(u,'(a)') "&snapshots"
        write(u,'(a)') "/"
        call snap(u, "pd",  "t_pd",  "p_pd",  .true.)
        call snap(u, "ctl", "t_ctl", "p_ctl", .false.)
        call snap(u, "alt", "t_alt", "p_alt", .false.)
        call var(u, "t_pd");  call var(u, "t_ctl");  call var(u, "t_alt")
        call var(u, "p_pd");  call var(u, "p_ctl");  call var(u, "p_alt")
        write(u,'(a)') "&v_zs"
        write(u,'(a)') "    filename = '"//fdat//"'"
        write(u,'(a)') "    name = 'zs'"
        write(u,'(a)') "    units = 'm' 'm'"
        write(u,'(a)') "    scaling = 1.0 0.0"
        write(u,'(a)') "    time = 0 0 0 0 1"
        write(u,'(a)') "/"
        write(u,'(a)') "&v_to"
        write(u,'(a)') "    filename = '"//fdat//"'"
        write(u,'(a)') "    name = 'to'"
        write(u,'(a)') "    units = 'K' 'K'"
        write(u,'(a)') "    scaling = 1.0 0.0"
        write(u,'(a)') "    time = 0 0 0 0 1"
        write(u,'(a)') "/"
        write(u,'(a)') "&v_so"
        write(u,'(a)') "    filename = '"//fdat//"'"
        write(u,'(a)') "    name = 'so'"
        write(u,'(a)') "    units = 'psu' 'psu'"
        write(u,'(a)') "    scaling = 1.0 0.0"
        write(u,'(a)') "    time = 0 0 0 0 1"
        write(u,'(a)') "/"
        write(u,'(a)') "&term_lgm"
        write(u,'(a)') "    coef   = 'alpha'"
        write(u,'(a)') "    snap0  = 'ctl'"
        write(u,'(a)') "    snap1  = 'alt'"
        write(u,'(a)') "    fields = 'tas' 'pr'"
        write(u,'(a)') "/"
        write(u,'(a)') "&term_ocean"
        write(u,'(a)') "    coef    = 'dTo'"
        write(u,'(a)') "    pattern = 'uniform'"
        write(u,'(a)') "    fields  = 'to'"
        write(u,'(a)') "/"
        write(u,'(a)') "&coef_alpha"
        write(u,'(a)') "    method = 'const'"
        write(u,'(a)') "/"
        write(u,'(a)') "&coef_dTo"
        write(u,'(a)') "    method = 'derived'"
        write(u,'(a)') "    from   = 'air_anomaly'"
        write(u,'(a)') "    factor = 0.25"
        write(u,'(a)') "    region = 'all'"
        write(u,'(a)') "    units  = 'K'"
        write(u,'(a)') "/"
        close(u)
    end subroutine write_par

    subroutine snap(u, name, tg, pg, ref)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, tg, pg
        logical,          intent(in) :: ref
        write(u,'(a)') "&snapshot_"//name
        write(u,'(a)') "    tas = 'v_"//tg//"'"
        write(u,'(a)') "    pr  = 'v_"//pg//"'"
        write(u,'(a)') "    zs  = 'v_zs'"
        if (ref) then
            write(u,'(a)') "    to  = 'v_to'"
            write(u,'(a)') "    so  = 'v_so'"
        end if
        write(u,'(a)') "/"
    end subroutine snap

    subroutine var(u, nm)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: nm
        write(u,'(a)') "&v_"//nm
        write(u,'(a)') "    filename = '"//fdat//"'"
        write(u,'(a)') "    name = '"//nm//"'"
        write(u,'(a)') "    units = 'K' 'K'"
        write(u,'(a)') "    scaling = 1.0 0.0"
        write(u,'(a)') "    time = 0 0 0 0 1"
        write(u,'(a)') "/"
    end subroutine var

end program test_climate_snapshots
