program test_climate_paleoemu
    ! The paleoemu climate adapter against the paleoemu module driven directly, as the core drives the
    ! global part (reference and geometry imports, coef_co2, PASS_INIT, PASS_BOUNDARY, updates, restart):
    !
    !   * the instance grid (GLOBAL-0.5DEG of maps/), opt%host, the coefficient co2;
    !   * the ice geometry is put piecewise constant on the emulator's ice-grid cells, so the adapter's
    !     binning returns the ice-grid fields: every export (tas at z_srf_clim, pr in kg m-2 yr-1, lapse,
    !     beta_p) bitwise equal to paleoemu_update + paleoemu_downscale with the same ice, CO2 and the
    !     orbit of [orbit] (kryos_orbit) at the evaluation time;
    !   * restart round trip: write, a new instance at the reference, read, PASS_RESTART, continue: bitwise;
    !   * the port copies its imports exactly, with their validity.
    ! Needs the emulator file in $PALEOEMU_FILE (PaleoEmulator scripts/06_emulator.jl); skipped without.
    ! With paleoemu=0 the test is not built (test_skip in config/common.mk).

    use precision,       only: dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr
    use phys_constants,  only: phys_const_load
    use fesm_fields,     only: fields_get, fields_put, fields_restart_write, fields_restart_read
    use kryos_component, only: init_ctx_t, step_time_t, census_t, census_count, PASS_INIT, PASS_BOUNDARY, PASS_RESTART
    use kryos_climate,   only: climate_slot_t, climate_create, climate_init, climate_init_state, climate_update, &
                               climate_restart_write, climate_restart_read, climate_end
    use kryos_climate_paleoemu, only: paleoemu_climate_t
    use kryos_orbit,     only: orbit_ctl_class, orbit_pars_class, orbit_ctl_load, orbit_pars_now
    use insolation,      only: insol_class, insol_init
    use paleoemu,        only: paleoemu_class, paleoemu_target_class, paleoemu_init, paleoemu_set_reference, &
                               paleoemu_update, paleoemu_target_init, paleoemu_downscale

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_climate_paleoemu/", PAR = DIR//"par.nml", REF = DIR//"ref.nc", &
                                   RST = DIR//"rst", MAPS = "maps", GRID = "GLOBAL-0.5DEG"
    real(dp), parameter :: DTT = 1000.0_dp, SEC_YEAR = 365.0_dp*86400.0_dp
    integer :: nfail = 0

    character(len=1024) :: efile
    type(init_ctx_t)     :: ctx
    type(census_t)       :: census, census2
    type(climate_slot_t) :: glob, glob2, port
    type(paleoemu_class) :: emu
    type(paleoemu_target_class) :: tgt
    type(orbit_ctl_class) :: octl
    type(insol_class)     :: ois
    integer,  allocatable :: bin(:,:)
    real(dp), allocatable :: Hr(:,:), zbr(:,:), zsr(:,:), lon(:,:), lat(:,:), tas_ref(:,:,:), pr_ref(:,:,:), zs_ref(:,:)
    real(dp) :: t
    integer  :: n, st_len

    call get_environment_variable("PALEOEMU_FILE", efile, length=st_len)
    if (st_len == 0) then
        write(*,"(a)") "test_climate_paleoemu: skipped (set PALEOEMU_FILE to the emulator file)"
        stop
    end if

    call execute_command_line("rm -rf "//DIR//" && mkdir -p "//DIR)
    call par_write()
    ctx%par_file  = PAR
    ctx%instance  = "clim_glob"
    ctx%component = "climate"
    ctx%backend   = "paleoemu"
    ctx%scope     = "global"
    ctx%map_fldr  = MAPS
    ctx%time      = -21000.0_dp
    ctx%time_rel  = -21000.0_dp
    ctx%dtt       = DTT
    ctx%sec_year  = SEC_YEAR
    call phys_const_load(ctx%cnst, "fesm-utils/par/phys_const_earth.nml", "Earth")

    ! ---- direct: the emulator, the reference climate on the grid, the orbit --------------------
    call paleoemu_init(emu, trim(efile))
    call reference_write()                         ! needs the grid: from a first instance below

    ! ---- the adapter --------------------------------------------------------------------------
    call climate_create(glob, ctx, census)
    call check(census_count(census, "climate:paleoemu:global") == 1, "factory: one paleoemu global part admitted")
    call climate_init(glob, ctx)
    call check(trim(glob%c%grid) == GRID .and. glob%c%nx == 720 .and. glob%c%ny == 360 .and. glob%c%opt%host, &
               "global part: grid GLOBAL-0.5DEG (720 x 360), opt%host")
    call check(size(glob%c%coef_names) == 1 .and. trim(glob%c%coef_names(1)) == "co2", "coefficient co2")
    call check(trim(glob%c%orbit) == "orbit", "orbit group [orbit] (component_t%orbit)")
    select type (c => glob%c)
    type is (paleoemu_climate_t)
        bin = c%bin
    end select

    call paleoemu_target_init(tgt, emu, lon, lat, tas_ref, pr_ref, zs_ref)
    call orbit_ctl_load(octl, PAR, "orbit")
    call insol_init(ois, fldr=trim(octl%fldr))

    Hr = emu%ice_H_ref;  zbr = emu%ice_z_bed_ref;  zsr = emu%ice_z_srf_ref
    call put_geometry(glob, Hr, zbr, zsr, ref=.true., t1=-21000.0_dp)
    call climate_init_state(glob, st(-21000.0_dp, -21000.0_dp), PASS_INIT)
    call check(maxval(abs(field3(glob, "tas") - tas_ref)) == 0.0_dp, "PASS_INIT: tas = the reference climate")

    call paleoemu_set_reference(emu, Hr, zbr, zsr)
    t = -21000.0_dp
    call drive(glob, 0, t, PASS_BOUNDARY)
    call compare(glob, "PASS_BOUNDARY")
    do n = 1, 3
        t = t + DTT
        call drive(glob, n, t)
        call compare(glob, "step "//char(48 + n))
    end do
    call check(maxval(abs(field3(glob, "tas") - tas_ref)) > 1.0_dp, "the grown ice changes tas somewhere by more than 1 K")

    ! ---- restart round trip -------------------------------------------------------------------
    call execute_command_line("rm -rf "//RST//" && mkdir -p "//RST)
    call climate_restart_write(glob, RST, t)
    call fields_restart_write(glob%c%imports, RST//"/clim_glob_imports.nc")
    call fields_restart_write(glob%c%exports, RST//"/clim_glob_exports.nc")
    call climate_create(glob2, ctx, census2)
    call climate_init(glob2, ctx)
    call climate_restart_read(glob2, RST, t)
    call fields_restart_read(glob2%c%imports, RST//"/clim_glob_imports.nc")
    call fields_restart_read(glob2%c%exports, RST//"/clim_glob_exports.nc")
    call climate_init_state(glob2, st(t, t), PASS_RESTART)
    t = t + DTT
    call drive(glob, 4, t)
    call drive(glob2, 4, t, direct=.false.)
    call compare(glob,  "continued")
    call check(all(field3(glob2, "tas") == field3(glob, "tas")) .and. all(field3(glob2, "pr") == field3(glob, "pr")), &
               "restarted = continued, bitwise")

    ! ---- the port --------------------------------------------------------------------------------
    call port_test()

    call climate_end(glob)
    call climate_end(glob2)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_climate_paleoemu: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_climate_paleoemu: OK"

contains

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&clim_glob", ' backend = "paleoemu"', ' grid = "'//GRID//'"', ' file = "'//trim(efile)//'"', &
                       ' reference = "'//REF//'"', "/", &
                       "&orbit", ' method = "forcing_time"', "/"
        close(u)
    end subroutine par_write

    subroutine reference_write()
        ! A smooth reference climate on GLOBAL-0.5DEG (read from maps/), with its grid_name.
        type(climate_slot_t) :: tmp
        type(census_t) :: cs
        integer :: m

        ! the grid's axes: a throwaway instance reads the description (its init needs the reference file:
        ! write a placeholder first, then the real one on the same axes)
        call placeholder()
        call climate_create(tmp, ctx, cs)
        call climate_init(tmp, ctx)
        lon = tmp%c%grid_obj%lon
        lat = tmp%c%grid_obj%lat
        call climate_end(tmp)
        allocate(tas_ref(size(lon, 1), size(lon, 2), 12), pr_ref(size(lon, 1), size(lon, 2), 12))
        do m = 1, 12
            tas_ref(:, :, m) = 288.0_dp - 40.0_dp*sin(lat*acos(-1.0_dp)/180.0_dp)**2 + 5.0_dp*cos(2.0_dp*acos(-1.0_dp)*m/12.0_dp)
            pr_ref(:, :, m)  = 2.0_dp + cos(lat*acos(-1.0_dp)/90.0_dp)
        end do
        zs_ref = 100.0_dp*cos(lon*acos(-1.0_dp)/180.0_dp)**2
        call write_ref(lon(:, 1), lat(1, :), tas_ref, pr_ref, zs_ref)
    end subroutine reference_write

    subroutine placeholder()
        real(dp), allocatable :: x(:), y(:), a3(:,:,:), p3(:,:,:), a2(:,:)
        integer :: i
        x = [(-179.75_dp + 0.5_dp*(i - 1), i = 1, 720)]
        y = [(-89.75_dp + 0.5_dp*(i - 1), i = 1, 360)]
        allocate(a3(720, 360, 12), p3(720, 360, 12), a2(720, 360))
        a3 = 273.15_dp;  p3 = 1.0_dp;  a2 = 0.0_dp
        call write_ref(x, y, a3, p3, a2)
    end subroutine placeholder

    subroutine write_ref(x, y, tas, pr, zs)
        real(dp), intent(in) :: x(:), y(:), tas(:,:,:), pr(:,:,:), zs(:,:)
        call nc_create(REF)
        call nc_write_attr(REF, "grid_name", GRID)
        call nc_write_dim(REF, "lon", x=x)
        call nc_write_dim(REF, "lat", x=y)
        call nc_write_dim(REF, "time", x=1.0_dp, dx=1.0_dp, nx=12)
        call nc_write(REF, "tas_ref", tas, dim1="lon", dim2="lat", dim3="time")
        call nc_write(REF, "pr_ref", pr, dim1="lon", dim2="lat", dim3="time")
        call nc_write(REF, "z_srf_ref", zs, dim1="lon", dim2="lat")
    end subroutine write_ref

    subroutine drive(g, n, t, pass, direct)
        ! The load of step n at t into the adapter (and, unless direct = .false., the same update of the
        ! direct emulator): a North American ice that grows with n, CO2 falling with n.
        type(climate_slot_t), intent(inout) :: g
        integer,              intent(in)    :: n
        real(dp),             intent(in)    :: t
        integer,  optional,   intent(in)    :: pass
        logical,  optional,   intent(in)    :: direct

        real(dp) :: H(size(Hr, 1), size(Hr, 2)), zs(size(Hr, 1), size(Hr, 2)), co2
        type(orbit_pars_class) :: orb
        logical :: dir_

        dir_ = .true.
        if (present(direct)) dir_ = direct
        H = Hr
        where (emu%ice_region == 1) H = Hr + 300.0_dp*(n + 1)
        zs = max(zbr + H, zsr)
        co2 = 280.0_dp - 20.0_dp*(4 - n)
        call put_geometry(g, H, zbr, zs, ref=.false., t1=t)
        call fields_put(g%c%imports, g%c%hi_coef(1), co2, t1=t)
        if (present(pass)) then
            call climate_init_state(g, st(t, t), pass)
        else
            call climate_update(g, st(t - DTT, t))
        end if
        if (dir_) then
            orb = orbit_pars_now(octl, ois, t)
            call paleoemu_update(emu, t, co2, H, zbr, zs, orbit=[orb%ecc, orb%obl, orb%per])
            call paleoemu_downscale(tgt, emu)
        end if
    end subroutine drive

    subroutine put_geometry(g, H, zb, zs, ref, t1)
        ! Ice-grid fields as piecewise constant fields on the instance grid (each cell the value of the
        ! ice-grid cell holding its centre; 0 outside the ice grid).
        type(climate_slot_t), intent(inout) :: g
        real(dp),             intent(in)    :: H(:,:), zb(:,:), zs(:,:), t1
        logical,              intent(in)    :: ref
        associate (c => g%c)
            if (ref) then
                call fields_put(c%imports, c%hi_H_ice_ref, spread_(H),  t1=t1)
                call fields_put(c%imports, c%hi_z_bed_ref, spread_(zb), t1=t1)
                call fields_put(c%imports, c%hi_z_srf_ref, spread_(zs), t1=t1)
            else
                call fields_put(c%imports, c%hi_H_ice, spread_(H),  t1=t1)
                call fields_put(c%imports, c%hi_z_bed, spread_(zb), t1=t1)
                call fields_put(c%imports, c%hi_z_srf, spread_(zs), t1=t1)
            end if
        end associate
    end subroutine put_geometry

    function spread_(X) result(Y)
        real(dp), intent(in) :: X(:,:)
        real(dp) :: Y(size(bin, 1), size(bin, 2)), Xf(size(X))
        integer  :: i, j
        Xf = reshape(X, [size(X)])
        Y = 0.0_dp
        do j = 1, size(bin, 2)
            do i = 1, size(bin, 1)
                if (bin(i, j) > 0) Y(i, j) = Xf(bin(i, j))
            end do
        end do
    end function spread_

    subroutine compare(g, what)
        type(climate_slot_t), intent(inout) :: g
        character(len=*),     intent(in)    :: what
        real(dp) :: dpr, dtas
        dpr  = maxval(abs(field3(g, "pr") - tgt%pr*SEC_YEAR/86400.0_dp))
        dtas = maxval(abs(field3(g, "tas") - tgt%tas))
        call check(dtas < 1.0e-10_dp .and. dpr < 1.0e-10_dp*maxval(abs(tgt%pr*SEC_YEAR/86400.0_dp)) .and. &
                   maxval(abs(field2(g, "z_srf_clim") - zs_ref)) == 0.0_dp .and. &
                   all(field3(g, "lapse") == emu%gamma) .and. all(field2(g, "beta_p") == 0.0_dp), &
                   what//": tas, pr, z_srf_clim, lapse, beta_p as paleoemu directly (the binning's round-off: "// &
                   "max |d tas| "//trim(r2s(dtas))//" K, |d pr| "//trim(r2s(dpr))//" kg m-2 yr-1)")
    end subroutine compare

    subroutine port_test()
        ! The port copies the imports (as the rows climate_g -> climate put them) with their validity.
        type(init_ctx_t) :: pc
        type(census_t)   :: cs
        real(dp) :: a3(10, 8, 12), a2(10, 8), b3(10, 8, 12), b2(10, 8), t0, t1
        integer  :: i

        pc = ctx
        pc%instance = "clim_port";  pc%scope = "regional";  pc%grid_name = "TEST-10x8"
        pc%grid%G%nx = 10;  pc%grid%G%ny = 8
        call climate_create(port, pc, cs)
        call climate_init(port, pc)
        call check(port%c%opt%port .and. .not. port%c%con%restart, "port: opt%port, no backend restart file")
        a3 = reshape([(real(i, dp), i = 1, size(a3))], shape(a3));  a2 = a3(:, :, 3)
        call fields_put(port%c%imports, port%c%hi_tas, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(port%c%imports, port%c%hi_pr, 2.0_dp*a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(port%c%imports, port%c%hi_lapse, 3.0_dp*a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(port%c%imports, port%c%hi_z_srf_clim, a2, t1=5.0_dp, t0=4.0_dp)
        call fields_put(port%c%imports, port%c%hi_beta_p, -a2, t1=5.0_dp, t0=4.0_dp)
        call climate_update(port, st(4.0_dp, 5.0_dp))
        call fields_get(port%c%exports, port%c%h_tas, b3, t0=t0, t1=t1)
        call fields_get(port%c%exports, port%c%h_beta_p, b2)
        call check(all(b3 == a3) .and. all(b2 == -a2) .and. t0 == 4.0_dp .and. t1 == 5.0_dp, &
                   "port: tas, beta_p copied with their validity")
        call fields_get(port%c%exports, port%c%h_pr, b3)
        call check(all(b3 == 2.0_dp*a3), "port: pr copied")
    end subroutine port_test

    function field3(g, name) result(a)
        type(climate_slot_t), intent(in) :: g
        character(len=*),     intent(in) :: name
        real(dp), allocatable :: a(:,:,:)
        integer :: h
        h = handle(g, name)
        allocate(a(g%c%nx, g%c%ny, 12))
        call fields_get(g%c%exports, h, a)
    end function field3

    function field2(g, name) result(a)
        type(climate_slot_t), intent(in) :: g
        character(len=*),     intent(in) :: name
        real(dp), allocatable :: a(:,:)
        integer :: h
        h = handle(g, name)
        allocate(a(g%c%nx, g%c%ny))
        call fields_get(g%c%exports, h, a)
    end function field2

    integer function handle(g, name) result(h)
        type(climate_slot_t), intent(in) :: g
        character(len=*),     intent(in) :: name
        do h = 1, g%c%exports%n
            if (trim(g%c%exports%f(h)%name) == trim(name)) return
        end do
        h = 0
    end function handle

    function r2s(x) result(str)
        real(dp), intent(in) :: x
        character(len=16) :: str
        write(str, "(es10.2)") x
    end function r2s

    function st(t0, t1) result(s)
        real(dp), intent(in) :: t0, t1
        type(step_time_t) :: s
        s = step_time_t(t0=t0, t1=t1, t_eval=t1, dt=DTT, dtt=DTT)
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

end program test_climate_paleoemu
