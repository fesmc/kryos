program test_climate_paleoemu
    ! The paleoemu climate adapter against the paleoemu module driven directly, as the core drives the
    ! global part (reference and geometry imports, coef_co2, PASS_INIT, PASS_BOUNDARY, updates, restart):
    !
    !   * the instance grid (GLOBAL-0.5DEG of maps/), opt%host, the coefficient co2;
    !   * the ice geometry is put piecewise constant on the emulator's ice-grid cells, so the adapter's
    !     binning returns the ice-grid fields: every export (tas at z_srf_clim, pr in kg m-2 yr-1, lapse,
    !     beta_p) equal to paleoemu_update + paleoemu_downscale with the same ice, CO2 and the orbital
    !     elements of the imports (the forcing's: insol's tables at the evaluation time), to the round-off
    !     of the binning;
    !   * without the orbit imports, the emulator's own orbit table;
    !   * keys scale_dtas, bias, bias_scale, beta_adj: as paleoemu_downscale with scale, bias_scale, beta and the
    !     bias of paleoemu_bias_init (a constant 2 K, dV_ref 50 m SLE: the grown ice gives a weight > 0);
    !   * the ocean dto_lev: 0 at PASS_INIT; then scale_dto (dto - dto0) + dto_nsd dto_sd filled by ocean_fill
    !     (the ocean cells unchanged) and interpolated per level (paleoemu_target_bilinear), on the axis `level`;
    !   * the anchor dto0 (dto_anchor, default on): the emulator's dto at present-day inputs (0 ka, its own orbit,
    !     CO2 280, the reference ice, the filters in equilibrium) as paleoemu directly on a copy, the run's
    !     emulator untouched; at those inputs dto_lev = 0; recomputed bitwise on restart; without it (dto_anchor
    !     = False) dto0 = 0 and dto_lev the filled emulator dto, not 0 at present day (the emulator's reference
    !     inputs, CESM-3Ma's 0 ka, hold the deglaciation in its filters);
    !   * restart round trip: write, a new instance at the reference, read, PASS_RESTART, continue: bitwise;
    !   * the port copies its imports exactly, with their validity;
    !   * the port with a present-day ocean (to_ref, so_ref: groups read from a file on its grid) declares
    !     dto_lev on the global part's axis `level` (from ctx%axes, as domain_create passes it) and exports
    !     set A on the file's depths: to_ann = to_ref + W dto_lev (W linear in depth, constant beyond the
    !     levels), so_ann = so_ref; dto_lev = 0 gives to_ann = to_ref;
    !   * with dto_basins: the basins extended over the domain by nearest neighbour; in each basin the mean
    !     over its sample (default zone:shelf_break_buffer, offshore of the basins) of the regions on the
    !     port's grid; a basin without a sample keeps its cells' values.
    ! Needs the emulator file in $PALEOEMU_FILE (PaleoEmulator scripts/06_emulator.jl); skipped without.
    ! With paleoemu=0 the test is not built (test_skip in config/common.mk).

    use precision,       only: dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr
    use phys_constants,  only: phys_const_load
    use fesm_fields,     only: fields_get, fields_put, fields_restart_write, fields_restart_read
    use kryos_component, only: init_ctx_t, step_time_t, census_t, census_count, PASS_INIT, PASS_BOUNDARY, PASS_RESTART
    use kryos_climate,   only: climate_slot_t, climate_create, climate_init, climate_init_state, climate_update, &
                               climate_restart_write, climate_restart_read, climate_end
    use kryos_climate_paleoemu, only: paleoemu_climate_t, paleoemu_port_t, ocean_fill
    use kryos_component, only: axis_find, axes_equal
    use insolation,      only: insol_class, insol_init, orbit_pars_class, orbit_pars_at
    use regions,         only: regions_init, regions_add_layer
    use paleoemu,        only: paleoemu_class, paleoemu_target_class, paleoemu_init, paleoemu_set_reference, &
                               paleoemu_update, paleoemu_target_init, paleoemu_downscale, &
                               paleoemu_bias_init, paleoemu_bias_weight, paleoemu_target_bilinear

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_climate_paleoemu/", PAR = DIR//"par.nml", REF = DIR//"ref.nc", &
                                   RST = DIR//"rst", MAPS = "maps", GRID = "GLOBAL-0.5DEG", &
                                   BIAS = DIR//"bias.nc", OCN = DIR//"ocean.nc"
    real(dp), parameter :: ZOCN(5) = [-2.0_dp, -100.0_dp, -500.0_dp, -1000.0_dp, -2000.0_dp]   ! as files give z
    real(dp), parameter :: DTT = 1000.0_dp, SEC_YEAR = 365.0_dp*86400.0_dp
    integer :: nfail = 0

    character(len=1024) :: efile
    type(init_ctx_t)     :: ctx
    type(census_t)       :: census, census2
    type(climate_slot_t) :: glob, glob2, port
    type(paleoemu_class) :: emu
    type(paleoemu_target_class) :: tgt
    type(insol_class)     :: ois
    type(climate_slot_t)  :: glob3, glob4, glob5, glob6
    type(census_t)        :: census3, census4, census5, census6
    real(dp), allocatable :: dto0(:,:,:)
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
    call check(glob%c%opt%orbit .and. glob%c%hi_orbit(1) > 0, "orbit imports declared (optional)")
    call check(glob%c%opt%ocn_lev .and. glob%c%h_dto_lev > 0 .and. axis_find(glob%c%axes, "level") > 0, &
               "dto_lev exported on the axis level")
    if (axis_find(glob%c%axes, "level") > 0) &
        call check(all(glob%c%axes(axis_find(glob%c%axes, "level"))%v == emu%ocean_depth) .and. &
                   trim(glob%c%axes(axis_find(glob%c%axes, "level"))%units) == "m", "  level: the emulator's depths [m]")
    select type (c => glob%c)
    type is (paleoemu_climate_t)
        bin = c%bin
    end select

    call paleoemu_target_init(tgt, emu, lon, lat, tas_ref, pr_ref, zs_ref)
    call insol_init(ois, fldr="input")

    Hr = emu%ice_H_ref;  zbr = emu%ice_z_bed_ref;  zsr = emu%ice_z_srf_ref
    call put_geometry(glob, Hr, zbr, zsr, ref=.true., t1=-21000.0_dp)
    call climate_init_state(glob, st(-21000.0_dp, -21000.0_dp), PASS_INIT)
    call check(maxval(abs(field3(glob, "tas") - tas_ref)) == 0.0_dp, "PASS_INIT: tas = the reference climate")
    call check(all(fieldl(glob) == 0.0_dp), "PASS_INIT: dto_lev = 0")

    call paleoemu_set_reference(emu, Hr, zbr, zsr)
    dto0 = anchor_direct()
    t = -21000.0_dp
    call drive(glob, 0, t, PASS_BOUNDARY)
    select type (c => glob%c)
    type is (paleoemu_climate_t)
        call check(maxval(abs(c%dto0 - dto0)) < 1.0e-10_dp .and. maxval(abs(dto0)) > 0.1_dp, &
                   "anchor: dto0 as paleoemu directly at present-day inputs (max |d| "// &
                   trim(r2s(maxval(abs(c%dto0 - dto0))))//" K); not 0 (max |dto0| "//trim(r2s(maxval(abs(dto0))))//" K)")
    end select
    call compare(glob, "PASS_BOUNDARY")
    call compare_ocean(glob, "PASS_BOUNDARY", 1.0_dp, 0.0_dp)
    do n = 1, 3
        t = t + DTT
        call drive(glob, n, t)
        call compare(glob, "step "//char(48 + n))
    end do
    call compare_ocean(glob, "step 3", 1.0_dp, 0.0_dp)
    call check(maxval(abs(fieldl(glob))) > 0.1_dp, "dto_lev differs from 0 by more than 0.1 K somewhere")
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
    call check(all(field3(glob2, "tas") == field3(glob, "tas")) .and. all(field3(glob2, "pr") == field3(glob, "pr")) &
               .and. all(fieldl(glob2) == fieldl(glob)), "restarted = continued, bitwise (tas, pr, dto_lev)")
    select type (c => glob%c)
    type is (paleoemu_climate_t)
        select type (c2 => glob2%c)
        type is (paleoemu_climate_t)
            call check(all(c2%dto0 == c%dto0), "restart: the anchor dto0 recomputed bitwise")
        end select
    end select

    ! ---- without the orbit: the emulator's own table -----------------------------------------------
    call climate_create(glob3, ctx, census3)
    call climate_init(glob3, ctx)
    call put_geometry(glob3, Hr, zbr, zsr, ref=.true., t1=-21000.0_dp)
    call paleoemu_init(emu, trim(efile))
    call paleoemu_set_reference(emu, Hr, zbr, zsr)
    call drive(glob3, 0, -21000.0_dp, PASS_BOUNDARY, orbit=.false.)
    call compare(glob3, "no orbit imports")
    call climate_end(glob3)

    ! ---- the temperature adjustments: scale_dtas, bias, beta_adj --------------------------------------
    call bias_write()
    ctx%instance = "clim_scl"
    call climate_create(glob4, ctx, census4)
    call climate_init(glob4, ctx)
    call put_geometry(glob4, Hr, zbr, zsr, ref=.true., t1=-21000.0_dp)
    call paleoemu_init(emu, trim(efile))
    call paleoemu_set_reference(emu, Hr, zbr, zsr)
    call paleoemu_bias_init(emu, BIAS)
    call drive(glob4, 2, -21000.0_dp, PASS_BOUNDARY, adjust=.true.)
    call check(paleoemu_bias_weight(emu) > 0.0_dp, "the bias weight of the grown ice is positive")
    call compare(glob4, "scale_dtas 1.5, bias x 0.5, beta_adj 0.05")
    call compare_ocean(glob4, "scale_dto 1.5, dto_nsd 1", 1.5_dp, 1.0_dp)
    call climate_end(glob4)

    ! ---- the anchor at present-day inputs ------------------------------------------------------------
    ctx%instance = "clim_glob"
    call climate_create(glob5, ctx, census5)
    call climate_init(glob5, ctx)
    call drive_pd(glob5)
    call check(maxval(abs(fieldl(glob5))) < 1.0e-12_dp, "anchor: dto_lev = 0 at present-day inputs (max |dto_lev| "// &
               trim(r2s(maxval(abs(fieldl(glob5)))))//" K)")
    call compare_ocean(glob5, "anchor at present day", 1.0_dp, 0.0_dp)
    call climate_end(glob5)
    ctx%instance = "clim_noanc"
    call climate_create(glob6, ctx, census6)
    call climate_init(glob6, ctx)
    call drive_pd(glob6)
    select type (c => glob6%c)
    type is (paleoemu_climate_t)
        call check(all(c%dto0 == 0.0_dp), "dto_anchor = False: dto0 = 0")
    end select
    call compare_ocean(glob6, "no anchor at present day", 1.0_dp, 0.0_dp)
    call check(maxval(abs(fieldl(glob6))) > 0.1_dp, "no anchor: dto_lev at present-day inputs is the emulator's offset "// &
               "(max |dto_lev| "//trim(r2s(maxval(abs(fieldl(glob6)))))//" K)")
    call climate_end(glob6)

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

    function anchor_direct() result(d)
        ! dto of a copy of the direct emulator (its reference just set) at present-day inputs.
        real(dp), allocatable :: d(:,:,:)
        type(paleoemu_class), allocatable :: e0
        allocate(e0)
        e0 = emu
        call paleoemu_update(e0, 0.0_dp, 280.0_dp, Hr, zbr, zsr)
        d = e0%now%dto
    end function anchor_direct

    subroutine drive_pd(g)
        ! Present-day inputs at 0 ka (PASS_BOUNDARY: the filters in equilibrium): the reference ice, CO2 280 ppm,
        ! no orbit imports (the emulator's own table); the direct emulator the same.
        type(climate_slot_t), intent(inout) :: g
        call put_geometry(g, Hr, zbr, zsr, ref=.true., t1=0.0_dp)
        call put_geometry(g, Hr, zbr, zsr, ref=.false., t1=0.0_dp)
        call fields_put(g%c%imports, g%c%hi_coef(1), 280.0_dp, t1=0.0_dp)
        call climate_init_state(g, st(0.0_dp, 0.0_dp), PASS_BOUNDARY)
        call paleoemu_init(emu, trim(efile))
        call paleoemu_set_reference(emu, Hr, zbr, zsr)
        call paleoemu_update(emu, 0.0_dp, 280.0_dp, Hr, zbr, zsr)
    end subroutine drive_pd

    subroutine bias_write()
        ! A constant 2 K bias on the CESM grid with dV_ref 50 m SLE.
        real(dp), allocatable :: b(:,:,:)
        allocate(b(emu%nlon, emu%nlat, 12));  b = 2.0_dp
        call nc_create(BIAS)
        call nc_write_dim(BIAS, "lon", x=emu%lon)
        call nc_write_dim(BIAS, "lat", x=emu%lat)
        call nc_write_dim(BIAS, "month", x=1, dx=1, nx=12)
        call nc_write(BIAS, "dtas_b", b, dim1="lon", dim2="lat", dim3="month")
        call nc_write_attr(BIAS, "dV_ref", 50.0_dp)
    end subroutine bias_write

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&clim_glob", ' backend = "paleoemu"', ' grid = "'//GRID//'"', ' file = "'//trim(efile)//'"', &
                       ' reference = "'//REF//'"', ' beta_adj = 0.0', "/", &
                       "&clim_scl", ' backend = "paleoemu"', ' grid = "'//GRID//'"', ' file = "'//trim(efile)//'"', &
                       ' reference = "'//REF//'"', ' scale_dtas = 1.5', ' bias = "'//BIAS//'"', ' bias_scale = 0.5', &
                       ' beta_adj = 0.05', ' scale_dto = 1.5', ' dto_nsd = 1.0', "/", &
                       "&clim_noanc", ' backend = "paleoemu"', ' grid = "'//GRID//'"', ' file = "'//trim(efile)//'"', &
                       ' reference = "'//REF//'"', ' dto_anchor = False', "/", &
                       "&clim_port_ocn", ' backend = "paleoemu"', ' to_ref = "test_to"', ' so_ref = "test_so"', "/", &
                       "&clim_port_bas", ' backend = "paleoemu"', ' to_ref = "test_to"', ' so_ref = "test_so"', &
                       ' dto_basins = "tb"', "/", &
                       "&test_to", ' filename = "'//OCN//'"', ' name = "to"', ' units = "degC" "K"', &
                       ' scaling = 1.0 273.15', ' time = 0 0 0 0 1', "/", &
                       "&test_so", ' filename = "'//OCN//'"', ' name = "so"', ' units = "PSU" "PSU"', &
                       ' scaling = 1.0 0.0', ' time = 0 0 0 0 1', "/"
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

    subroutine drive(g, n, t, pass, direct, orbit, adjust)
        ! The load of step n at t into the adapter (and, unless direct = .false., the same update of the
        ! direct emulator): a North American ice that grows with n, CO2 falling with n, the orbit of the
        ! tables at t as the forcing would export it (orbit = .false.: none, the emulator's own).
        type(climate_slot_t), intent(inout) :: g
        integer,              intent(in)    :: n
        real(dp),             intent(in)    :: t
        integer,  optional,   intent(in)    :: pass
        logical,  optional,   intent(in)    :: direct, orbit
        logical,  optional,   intent(in)    :: adjust

        real(dp), allocatable :: H(:,:), zs(:,:)
        real(dp) :: co2
        type(orbit_pars_class) :: orb
        logical :: dir_, orb_, adj_

        allocate(H(size(Hr, 1), size(Hr, 2)), zs(size(Hr, 1), size(Hr, 2)))
        dir_ = .true.
        if (present(direct)) dir_ = direct
        orb_ = .true.
        if (present(orbit)) orb_ = orbit
        orb = orbit_pars_at(ois, t)
        H = Hr
        where (emu%ice_region == 1) H = Hr + 300.0_dp*(n + 1)
        zs = max(zbr + H, zsr)
        co2 = 280.0_dp - 20.0_dp*(4 - n)
        call put_geometry(g, H, zbr, zs, ref=.false., t1=t)
        call fields_put(g%c%imports, g%c%hi_coef(1), co2, t1=t)
        if (orb_) then
            call fields_put(g%c%imports, g%c%hi_orbit(1), orb%ecc, t1=t)
            call fields_put(g%c%imports, g%c%hi_orbit(2), orb%obl, t1=t)
            call fields_put(g%c%imports, g%c%hi_orbit(3), orb%per, t1=t)
        end if
        if (present(pass)) then
            call climate_init_state(g, st(t, t), pass)
        else
            call climate_update(g, st(t - DTT, t))
        end if
        if (dir_) then
            if (orb_) then
                call paleoemu_update(emu, t, co2, H, zbr, zs, orbit=[orb%ecc, orb%obl, orb%per])
            else
                call paleoemu_update(emu, t, co2, H, zbr, zs)
            end if
            adj_ = .false.
            if (present(adjust)) adj_ = adjust
            if (adj_) then
                call paleoemu_downscale(tgt, emu, scale=[1.5_dp, 1.0_dp], bias_scale=0.5_dp, beta=0.05_dp)
            else
                call paleoemu_downscale(tgt, emu)
            end if
        end if
    end subroutine drive

    subroutine put_geometry(g, H, zb, zs, ref, t1)
        ! Ice-grid fields as piecewise constant fields on the instance grid (each cell the value of the
        ! ice-grid cell holding its centre; 0 outside the ice grid).
        ! (Local arrays, not a function result as the actual argument: ifx 2023.2 -O2 stops with an internal
        ! compiler error on that.)
        type(climate_slot_t), intent(inout) :: g
        real(dp),             intent(in)    :: H(:,:), zb(:,:), zs(:,:), t1
        logical,              intent(in)    :: ref
        real(dp), allocatable :: aH(:,:), ab(:,:), as(:,:)
        allocate(aH(size(bin, 1), size(bin, 2)), ab(size(bin, 1), size(bin, 2)), as(size(bin, 1), size(bin, 2)))
        call spread_(H, aH);  call spread_(zb, ab);  call spread_(zs, as)
        if (ref) then
            call fields_put(g%c%imports, g%c%hi_H_ice_ref, aH, t1=t1)
            call fields_put(g%c%imports, g%c%hi_z_bed_ref, ab, t1=t1)
            call fields_put(g%c%imports, g%c%hi_z_srf_ref, as, t1=t1)
        else
            call fields_put(g%c%imports, g%c%hi_H_ice, aH, t1=t1)
            call fields_put(g%c%imports, g%c%hi_z_bed, ab, t1=t1)
            call fields_put(g%c%imports, g%c%hi_z_srf, as, t1=t1)
        end if
    end subroutine put_geometry

    subroutine spread_(X, Y)
        real(dp), intent(in)  :: X(:,:)
        real(dp), intent(out) :: Y(:,:)
        real(dp) :: Xf(size(X))
        integer  :: i, j
        Xf = reshape(X, [size(X)])
        Y = 0.0_dp
        do j = 1, size(bin, 2)
            do i = 1, size(bin, 1)
                if (bin(i, j) > 0) Y(i, j) = Xf(bin(i, j))
            end do
        end do
    end subroutine spread_

    subroutine compare(g, what, beta)
        type(climate_slot_t), intent(inout) :: g
        character(len=*),     intent(in)    :: what
        real(dp), optional,   intent(in)    :: beta
        real(dp) :: dpr, dtas, bp
        bp = 0.0_dp
        if (present(beta)) bp = beta
        dpr  = maxval(abs(field3(g, "pr") - tgt%pr*SEC_YEAR/86400.0_dp))
        dtas = maxval(abs(field3(g, "tas") - tgt%tas))
        call check(dtas < 1.0e-10_dp .and. dpr < 1.0e-10_dp*maxval(abs(tgt%pr*SEC_YEAR/86400.0_dp)) .and. &
                   maxval(abs(field2(g, "z_srf_clim") - zs_ref)) == 0.0_dp .and. &
                   all(field3(g, "lapse") == emu%gamma) .and. all(field2(g, "beta_p") == bp), &
                   what//": tas, pr, z_srf_clim, lapse, beta_p as paleoemu directly (the binning's round-off: "// &
                   "max |d tas| "//trim(r2s(dtas))//" K, |d pr| "//trim(r2s(dpr))//" kg m-2 yr-1)")
    end subroutine compare

    subroutine compare_ocean(g, what, s, nsd)
        ! dto_lev against the adapter's own emulator state: s (dto - dto0) + nsd dto_sd, filled, each level
        ! bilinearly onto the grid (bitwise); its dto against the direct emulator's (the binning's round-off).
        type(climate_slot_t), intent(inout) :: g
        character(len=*),     intent(in)    :: what
        real(dp),             intent(in)    :: s, nsd
        real(dp), allocatable :: d(:,:,:), d0(:,:,:), x(:,:,:)
        integer :: k

        select type (c => g%c)
        type is (paleoemu_climate_t)
            call check(maxval(abs(c%emu%now%dto - emu%now%dto)) < 1.0e-10_dp, &
                       what//": dto as paleoemu directly (max |d| "//trim(r2s(maxval(abs(c%emu%now%dto - emu%now%dto))))//" K)")
            d0 = s*(c%emu%now%dto - c%dto0) + nsd*c%emu%now%dto_sd
            d  = d0
            call ocean_fill(c%emu%ocean_mask, d)
            call check(all(d == d0 .or. c%emu%ocean_mask /= 1), what//": ocean_fill leaves the ocean cells unchanged")
            call check(count(c%emu%ocean_mask /= 1) > 0 .and. all(abs(d) < 50.0_dp), &
                       what//": the filled land and below-floor cells are finite and bounded")
            allocate(x(g%c%nx, g%c%ny, c%emu%nlev))
            do k = 1, c%emu%nlev
                x(:, :, k) = paleoemu_target_bilinear(c%tgt, d(:, :, k))
            end do
            call check(all(fieldl(g) == x), what//": dto_lev = the bilinear of the filled field per level, bitwise")
        end select
    end subroutine compare_ocean

    function fieldl(g) result(a)
        ! dto_lev [x, y, level]
        type(climate_slot_t), intent(in) :: g
        real(dp), allocatable :: a(:,:,:)
        allocate(a(g%c%nx, g%c%ny, emu%nlev))
        call fields_get(g%c%exports, handle(g, "dto_lev"), a)
    end function fieldl

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
        call check(.not. port%c%opt%ocn_lev .and. .not. port%c%caps%has_ocn .and. port%c%hi_dto_lev == 0, &
                   "port without to_ref, so_ref: no ocean")
        call port_ocean_test(pc)
    end subroutine port_test

    subroutine port_ocean_test(pc0)
        ! The port with the present-day ocean of a file on its grid (to [degC], so; 5 depths, z negative).
        type(init_ctx_t), intent(in) :: pc0
        type(init_ctx_t) :: pc
        type(census_t)   :: cs
        type(climate_slot_t) :: po
        real(dp) :: to(10, 8, 5), so(10, 8, 5), dl(10, 8, 8), a3(10, 8, 12), a2(10, 8), b(10, 8, 5), w, t0, t1
        real(dp) :: ex(10, 8, 5), lev(8)
        integer  :: i, j, k, l, kg, kp

        do k = 1, 5
            do j = 1, 8
                do i = 1, 10
                    to(i, j, k) = 2.0_dp - 0.001_dp*abs(ZOCN(k)) + 0.1_dp*i - 0.05_dp*j
                    so(i, j, k) = 34.0_dp + 0.0005_dp*abs(ZOCN(k)) + 0.01_dp*i
                end do
            end do
        end do
        call nc_create(OCN)
        call nc_write_dim(OCN, "x", x=1.0_dp, dx=1.0_dp, nx=10)
        call nc_write_dim(OCN, "y", x=1.0_dp, dx=1.0_dp, nx=8)
        call nc_write_dim(OCN, "z", x=ZOCN, units="m")
        call nc_write(OCN, "to", to, dim1="x", dim2="y", dim3="z")
        call nc_write(OCN, "so", so, dim1="x", dim2="y", dim3="z")

        pc = pc0
        pc%instance = "clim_port_ocn"
        pc%axes = glob%c%axes                          ! as domain_create seeds them (axes_g)
        call climate_create(po, pc, cs)
        call climate_init(po, pc)
        kg = axis_find(glob%c%axes, "level");  kp = axis_find(po%c%axes, "level")
        call check(po%c%opt%ocn_lev .and. po%c%caps%has_ocn .and. po%c%hi_dto_lev > 0 .and. kp > 0, &
                   "port with to_ref, so_ref: imports dto_lev, exports set A")
        if (kp == 0) return
        call check(axes_equal(po%c%axes(kp), glob%c%axes(kg)), "  on the global part's axis level (the port sees it)")
        k = axis_find(po%c%axes, "depth")
        call check(k > 0, "  axis depth declared")
        if (k == 0) return
        call check(all(po%c%axes(k)%v == abs(ZOCN)), "  depth: the file's |z| [m]")

        lev = glob%c%axes(kg)%v
        do l = 1, 8
            dl(:, :, l) = 0.5_dp*l + 0.01_dp*spread([(real(i, dp), i = 1, 10)], 2, 8)
        end do
        a3 = 270.0_dp;  a2 = 0.0_dp
        call fields_put(po%c%imports, po%c%hi_tas, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_pr, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_lapse, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_z_srf_clim, a2, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_beta_p, a2, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_dto_lev, dl, t1=5.0_dp, t0=4.0_dp)
        call climate_update(po, st(4.0_dp, 5.0_dp))
        ! the expected change: linear in depth between the levels, constant beyond them
        do k = 1, 5
            if (abs(ZOCN(k)) <= lev(1)) then
                ex(:, :, k) = dl(:, :, 1)
            else if (abs(ZOCN(k)) >= lev(8)) then
                ex(:, :, k) = dl(:, :, 8)
            else
                l = count(lev <= abs(ZOCN(k)))
                w = (abs(ZOCN(k)) - lev(l))/(lev(l + 1) - lev(l))
                ex(:, :, k) = (1.0_dp - w)*dl(:, :, l) + w*dl(:, :, l + 1)
            end if
        end do
        call fields_get(po%c%exports, po%c%h_dto_ann, b, t0=t0, t1=t1)
        call check(maxval(abs(b - ex)) < 1.0e-12_dp .and. t0 == 4.0_dp .and. t1 == 5.0_dp, &
                   "  dto_ann = dto_lev linear in depth (constant beyond the levels), with its validity")
        call fields_get(po%c%exports, po%c%h_to_ann, b)
        call check(maxval(abs(b - (to + 273.15_dp + ex))) < 1.0e-4_dp, "  to_ann = to_ref [K] + dto_ann (to_ref in sp)")
        call fields_get(po%c%exports, po%c%h_so_ann, b)
        call check(maxval(abs(b - so)) < 1.0e-4_dp, "  so_ann = so_ref")
        dl = 0.0_dp
        call fields_put(po%c%imports, po%c%hi_dto_lev, dl, t1=6.0_dp, t0=5.0_dp)
        call climate_update(po, st(5.0_dp, 6.0_dp))
        call fields_get(po%c%exports, po%c%h_to_ann, b)
        select type (c => po%c)
        type is (paleoemu_port_t)
            call check(all(b == c%to0), "  dto_lev = 0: to_ann = to_ref")
        end select
        call port_basins_test(pc)
    end subroutine port_ocean_test

    subroutine port_basins_test(pc0)
        ! Basins tb: 1 in x <= 5, 2 in x > 5, none in rows y = 7, 8 (offshore, extended from row 6); the
        ! shelf-break buffer (zone 1) in cells (2:3, 7:8), offshore of basin 1: basin 1 (with its extension)
        ! takes their mean, basin 2 keeps its values.
        type(init_ctx_t), intent(in) :: pc0
        type(init_ctx_t) :: pc
        type(census_t)   :: cs
        type(climate_slot_t) :: po
        integer  :: tb(10, 8), zone(10, 8), i, l
        real(dp) :: dl(10, 8, 8), a3(10, 8, 12), a2(10, 8), b(10, 8, 5), m(8)

        tb = 1;  tb(6:, :) = 2;  tb(:, 7:8) = 0
        zone = 2;  zone(2:3, 7:8) = 1
        pc = pc0
        pc%instance = "clim_port_bas"
        pc%axes = glob%c%axes
        call regions_init(pc%reg, nx=10, ny=8)
        call regions_add_layer(pc%reg, "tb", tb)
        call regions_add_layer(pc%reg, "zone", zone, codes=[0, 1, 2, 3], &
                               names=[character(len=24) :: "open_ocean", "shelf_break_buffer", "continental_shelf", "land"])
        call climate_create(po, pc, cs)
        call climate_init(po, pc)
        do l = 1, 8
            dl(:, :, l) = 0.5_dp*l + 0.01_dp*spread([(real(i, dp), i = 1, 10)], 2, 8) + 0.003_dp*spread([(real(i, dp), i = 1, 8)], 1, 10)
        end do
        a3 = 270.0_dp;  a2 = 0.0_dp
        call fields_put(po%c%imports, po%c%hi_tas, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_pr, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_lapse, a3, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_z_srf_clim, a2, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_beta_p, a2, t1=5.0_dp, t0=4.0_dp)
        call fields_put(po%c%imports, po%c%hi_dto_lev, dl, t1=5.0_dp, t0=4.0_dp)
        call climate_update(po, st(4.0_dp, 5.0_dp))
        call fields_get(po%c%exports, po%c%h_dto_ann, b)
        do l = 1, 8
            m(l) = sum(dl(2:3, 7:8, l))/4.0_dp
        end do
        select type (c => po%c)
        type is (paleoemu_port_t)
            call check(all(abs(b(3, 5, :) - matmul(c%W, m)) < 1.0e-12_dp) .and. &
                       all(abs(b(1, 1, :) - matmul(c%W, m)) < 1.0e-12_dp) .and. &
                       all(abs(b(2, 8, :) - matmul(c%W, m)) < 1.0e-12_dp), &
                       "basins: basin 1 the mean of its offshore shelf-break cells, in all its (extended) cells")
            call check(all(abs(b(7, 3, :) - matmul(c%W, dl(7, 3, :))) < 1.0e-12_dp) .and. &
                       all(abs(b(8, 8, :) - matmul(c%W, dl(8, 8, :))) < 1.0e-12_dp), &
                       "basins: basin 2 (no sample) keeps its values, its extension too")
        end select
    end subroutine port_basins_test

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
