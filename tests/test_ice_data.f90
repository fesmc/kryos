program test_ice_data
    ! The ice backend data (make check): a transient thickness, bed and sea surface of a file on the
    ! instance grid (a toy grid), in "years ago" relative to 1950 with the clock's time_ref 1950:
    !   * the fields at the evaluation time, linear in time between the records;
    !   * z_srf (grounded z_bed + H, floating z_sl + (1 - rho_i/rho_w) H, the bed or the sea elsewhere) and
    !     f_grnd by the flotation rule with the run's densities;
    !   * dHidt and mb_net over the step; smb = cmb = 0; the static z_bed_ref, H_ice_pd, footprint;
    !   * without z_bed and z_sl: the reference bed and z_sl = 0.

    use precision,       only: dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr
    use coords,          only: grid_init
    use phys_constants,  only: phys_const_set
    use fesm_fields,     only: fields_get, fields_handle
    use kryos_component, only: init_ctx_t, step_time_t, census_t, PASS_INIT
    use kryos_ice,       only: ice_slot_t, ice_create, ice_init, ice_init_state, ice_update, ice_end

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_ice_data/", PAR = DIR//"par.nml", F = DIR//"hist.nc"
    integer,  parameter :: nx = 41, ny = 31
    real(dp), parameter :: DXKM = 16.0_dp
    integer :: nfail = 0
    type(init_ctx_t) :: ctx
    type(census_t)   :: census, census2
    type(ice_slot_t) :: a, b
    real(dp), allocatable :: H0(:,:), H1(:,:), zb0(:,:), zb1(:,:), x(:,:), zs(:,:), fg(:,:)
    real(dp) :: rho_i, rho_w, m
    integer  :: i, j

    call execute_command_line("rm -rf "//DIR//" && mkdir -p "//DIR)
    allocate(H0(nx, ny), H1(nx, ny), zb0(nx, ny), zb1(nx, ny))
    do j = 1, ny
        do i = 1, nx
            H0(i, j)  = max(0.0_dp, 2000.0_dp - 120.0_dp*abs(real(i - nx/2, dp)))
            H1(i, j)  = 0.5_dp*H0(i, j)
            zb0(i, j) = -400.0_dp + 10.0_dp*real(j, dp)
            zb1(i, j) = zb0(i, j) + 100.0_dp
        end do
    end do
    call write_file()
    call write_par()

    ctx%par_file  = PAR
    ctx%instance  = "ice_x"
    ctx%component = "ice"
    ctx%backend   = "data"
    ctx%domain    = "Greenland"
    ctx%grid_name = "TEST"
    call grid_init(ctx%grid, name="TEST", mtype="cartesian", units="km", dx=DXKM, nx=nx, dy=DXKM, ny=ny, &
                   x0=0.0_dp, y0=0.0_dp)
    ctx%time_ref  = 1950.0_dp
    rho_i = 917.0_dp;  rho_w = 1028.0_dp
    call phys_const_set(ctx%cnst, label="Earth", rho_ice=rho_i, rho_sw=rho_w, rho_w=1000.0_dp, g=9.81_dp, &
                        area_seasurf=3.618e8_dp)
    allocate(ctx%ref%z_bed(nx, ny), ctx%ref%H_ice(nx, ny))
    ctx%ref%z_bed = 50.0_dp;  ctx%ref%H_ice = 7.0_dp

    ! ---- with z_bed and z_sl from the file -----------------------------------------------------
    call ice_create(a, ctx, census)
    call ice_init(a, ctx)
    call ice_init_state(a, st(-10000.0_dp, -10000.0_dp), PASS_INIT)    ! 10 ka: the first record
    call check(maxval(abs(get(a, "H_ice") - H0)) == 0.0_dp .and. maxval(abs(get(a, "z_bed") - zb0)) == 0.0_dp, &
               "init at 10 ka: H_ice and z_bed of the first record")
    call check(all(get(a, "z_sl") == -50.0_dp), "z_sl = dz_sl of the first record (-50 m)")
    call check(all(get(a, "z_bed_ref") == 50.0_dp) .and. all(get(a, "H_ice_pd") == 7.0_dp) .and. &
               all(get(a, "footprint") == 1.0_dp), "static: z_bed_ref, H_ice_pd (the reference), footprint 1")
    call ice_update(a, st(-10000.0_dp, -7500.0_dp))                    ! halfway to the second record (5 ka)
    x = 0.5_dp*(H0 + H1)
    call check(maxval(abs(get(a, "H_ice") - x)) < 1.0e-9_dp .and. &
               maxval(abs(get(a, "z_bed") - 0.5_dp*(zb0 + zb1))) < 1.0e-9_dp .and. &
               all(abs(get(a, "z_sl") + 25.0_dp) < 1.0e-9_dp), "7.5 ka: linear in time between the records")
    call check(maxval(abs(get(a, "dHidt") - (x - H0)/2500.0_dp)) < 1.0e-12_dp .and. &
               maxval(abs(get(a, "mb_net") - rho_i*(x - H0)/2500.0_dp)) < 1.0e-9_dp .and. &
               all(get(a, "smb") == 0.0_dp) .and. all(get(a, "cmb") == 0.0_dp), "dHidt, mb_net over the step; smb = cmb = 0")
    zs = get(a, "z_srf");  fg = get(a, "f_grnd")
    m = 0.0_dp
    do j = 1, ny
        do i = 1, nx
            associate (H => x(i, j), zb => 0.5_dp*(zb0(i, j) + zb1(i, j)), zl => -25.0_dp)
                if (H > 0.0_dp .and. rho_i*H >= rho_w*(zl - zb)) then
                    m = max(m, abs(zs(i, j) - (zb + H)), abs(fg(i, j) - 1.0_dp))
                else if (H > 0.0_dp) then
                    m = max(m, abs(zs(i, j) - (zl + (1.0_dp - rho_i/rho_w)*H)), abs(fg(i, j)))
                else
                    m = max(m, abs(zs(i, j) - max(zb, zl)))
                end if
            end associate
        end do
    end do
    call check(m < 1.0e-9_dp, "z_srf and f_grnd by the flotation rule with the run's densities")
    call check(any(fg == 0.0_dp .and. x > 0.0_dp) .and. any(fg == 1.0_dp) .and. any(x == 0.0_dp), &
               "the test has grounded, floating and ice-free cells")
    call ice_end(a)

    ! ---- without z_bed and z_sl: the reference bed, sea surface 0 -----------------------------------
    ctx%instance = "ice_y"
    call ice_create(b, ctx, census2)
    call ice_init(b, ctx)
    call ice_init_state(b, st(-5000.0_dp, -5000.0_dp), PASS_INIT)
    call check(maxval(abs(get(b, "H_ice") - H1)) == 0.0_dp .and. all(get(b, "z_bed") == 50.0_dp) .and. &
               all(get(b, "z_sl") == 0.0_dp), "only H_ice bound: the reference bed and z_sl = 0")
    call ice_end(b)

    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_ice_data: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_ice_data: OK"

contains

    subroutine write_file()
        ! Two records, 10 ka and 5 ka ("years ago" relative to 1950), on the toy grid.
        real(dp) :: H(nx, ny, 2), zb(nx, ny, 2), dsl(nx, ny, 2)
        integer  :: i
        H(:, :, 1) = H0;  H(:, :, 2) = H1
        zb(:, :, 1) = zb0;  zb(:, :, 2) = zb1
        dsl(:, :, 1) = -50.0_dp;  dsl(:, :, 2) = 0.0_dp
        call nc_create(F)
        call nc_write_attr(F, "grid_name", "TEST")
        call nc_write_dim(F, "xc", x=[(DXKM*(i - 1), i = 1, nx)], units="km")
        call nc_write_dim(F, "yc", x=[(DXKM*(i - 1), i = 1, ny)], units="km")
        call nc_write_dim(F, "time", x=[10000.0_dp, 5000.0_dp], units="years ago")
        call nc_write(F, "H_ice", H, dim1="xc", dim2="yc", dim3="time", units="m")
        call nc_write(F, "z_bed", zb, dim1="xc", dim2="yc", dim3="time", units="m")
        call nc_write(F, "dz_sl", dsl, dim1="xc", dim2="yc", dim3="time", units="m")
    end subroutine write_file

    subroutine write_par()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&ice_x", ' backend = "data"', ' H_ice = "x_H"', ' z_bed = "x_zb"', ' z_sl = "x_sl"', "/", &
                       "&ice_y", ' backend = "data"', ' H_ice = "x_H"', "/", &
                       "&x_H",  ' filename = "'//F//'"', ' name = "H_ice"', ' mode = "transient"', "/", &
                       "&x_zb", ' filename = "'//F//'"', ' name = "z_bed"', ' mode = "transient"', "/", &
                       "&x_sl", ' filename = "'//F//'"', ' name = "dz_sl"', ' mode = "transient"', "/"
        close(u)
    end subroutine write_par

    function get(s, name) result(v)
        type(ice_slot_t), intent(in) :: s
        character(len=*), intent(in) :: name
        real(dp) :: v(nx, ny)
        call fields_get(s%c%exports, fields_handle(s%c%exports, name), v)
    end function get

    function st(t0, t1) result(s)
        real(dp), intent(in) :: t0, t1
        type(step_time_t) :: s
        s = step_time_t(t0=t0, t1=t1, t_eval=t1, dt=t1 - t0, dtt=t1 - t0)
    end function st

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

end program test_ice_data
