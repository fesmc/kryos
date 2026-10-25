program test_geo_fastisostasy
    ! The FastIsostasy geo adapter against FastIsostasy driven directly, as yelmox
    ! drives it (kryos_startup.f90, kryos_coupling.f90, yelmox.f90): one regional
    ! instance and the global part (bsl method fastiso) through the component
    ! wrappers. Cold start, three updates with a growing ice load, and a restart
    ! round trip must agree bit for bit with the direct calls. With
    ! fastisostasy=0 the test is not built (test_skip in config/common.mk).

    use precision,       only: sp, dp
    use coordinates,     only: grid_init
    use phys_constants,  only: phys_const_load
    use fesm_fields,     only: fields_get, fields_put
    use kryos_component, only: init_ctx_t, step_time_t, census_t, &
                               PASS_INIT, PASS_GEOMETRY, PASS_RESTART
    use kryos_geo,       only: geo_slot_t, geo_create, geo_init, geo_init_state, geo_update, geo_reduce, &
                               geo_scatter, geo_restart_write, geo_restart_read, geo_write_output, geo_end
    use fastisostasy,    only: isos_class, isos_init, isos_init_ref, isos_init_state, isos_update, &
                               bsl_class, bsl_init, bsl_update

    implicit none

    integer, parameter :: NX = 32, NY = 32
    character(len=*), parameter :: PAR = "test_geo_fastisostasy.nml", RST = "test_geo_fastisostasy_rst"
    real(dp), parameter :: DTT = 10.0_dp

    type(init_ctx_t)    :: ctxr, ctxg
    type(census_t)      :: census, census2
    type(geo_slot_t)    :: reg, glob, reg2, glob2
    type(isos_class)    :: isos
    type(bsl_class)     :: bsl
    real(sp) :: z_bed(NX,NY), H0(NX,NY), H(NX,NY)
    real(dp) :: t
    integer  :: n

    call par_write()
    call contexts()
    call cnst_load()

    ! ---- adapter: cold start --------------------------------------------------
    call geo_create(glob, ctxg, census)
    call geo_create(reg,  ctxr, census)
    call geo_init(glob, ctxg)
    call geo_init(reg,  ctxr)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_INIT)
    call geo_init_state(reg,  st(0.0_dp, 0.0_dp), PASS_INIT)
    call to_regional(glob, reg, 0.0_dp)
    call put_ice_exports(reg, 0.0_dp)
    call geo_init_state(reg, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
    ! isos_init_state moves the shared bsl when the initial ice differs from the
    ! reference (yelmox mutates its shared bsl there): the core reduces once after
    ! every PASS_GEOMETRY.
    call reduce(glob, reg, 0.0_dp)
    call geo_scatter(reg, st(0.0_dp, 0.0_dp))

    ! ---- direct: yelmox order ---------------------------------------------------
    call bsl_init(bsl, PAR, 0.0_dp)
    call bsl_update(bsl, 0.0_dp)
    call isos_init(isos, PAR, "isos", NX, NY, 32.0e3_sp, 32.0e3_sp, cnst=ctxr%cnst, domain="test", grid_name="TEST")
    z_bed = real(ctxr%ref%z_bed, sp)
    H0    = real(ctxr%ref%H_ice, sp)
    call isos_init_ref(isos, z_bed, H0)
    call isos_init_state(isos, z_bed, 0.9_sp*H0, 0.0_dp, bsl)
    call compare_regional(reg, isos, bsl, "cold start")

    do n = 1, 3
        t = real(n, dp) * DTT
        H = H0 * (1.0_sp + 0.5_sp * real(n, sp))
        call bsl_update(bsl, t)
        call bsl_update_global(glob, t)
        call to_regional(glob, reg, t)
        call put_H(reg, H, t - DTT)
        call geo_update(reg, st(t - DTT, t))
        call isos_update(isos, H, t, bsl)
        call reduce(glob, reg, t)
        call geo_scatter(reg, st(t - DTT, t))
        call compare_regional(reg, isos, bsl, "step")
    end do

    ! ---- restart round trip -------------------------------------------------------
    call execute_command_line("rm -rf "//RST//" && mkdir -p "//RST)
    call geo_restart_write(glob, RST, t)
    call geo_restart_write(reg,  RST, t)

    call geo_create(glob2, ctxg, census2)
    call geo_create(reg2,  ctxr, census2)
    call geo_init(glob2, ctxg)
    call geo_init(reg2,  ctxr)
    call geo_restart_read(glob2, RST, t)
    call geo_init_state(glob2, st(t, t), PASS_RESTART)
    call to_regional(glob2, reg2, t)
    call put_H(reg2, H, t)
    call geo_restart_read(reg2, RST, t)
    call geo_init_state(reg2, st(t, t), PASS_RESTART)

    t = t + DTT
    H = H * 1.1_sp
    call bsl_update(bsl, t)
    call isos_update(isos, H, t, bsl)
    call advance(glob, reg, H, t)
    call advance(glob2, reg2, H, t)
    call compare_regional(reg,  isos, bsl, "continued run")
    call compare_regional(reg2, isos, bsl, "restarted run")

    ! backend output streams
    call geo_write_output(reg,  RST//"/geo_t_isos.nc",    t, "isos",    .true.)
    call geo_write_output(reg,  RST//"/geo_t_isos_ts.nc", t, "isos_ts", .true.)
    call geo_write_output(glob, RST//"/bsl_t_bsl_ts.nc",  t, "bsl_ts",  .true.)
    call geo_write_output(reg,  RST//"/geo_t_isos.nc",    t + DTT, "isos", .false.)
    call check(file_exists(RST//"/geo_t_isos.nc") .and. file_exists(RST//"/geo_t_isos_ts.nc") .and. &
               file_exists(RST//"/bsl_t_bsl_ts.nc"), "output streams")

    call geo_end(reg2); call geo_end(glob2); call geo_end(reg); call geo_end(glob)
    call execute_command_line("rm -rf "//RST//" "//PAR)
    write(*,'(a)') "test_geo_fastisostasy: OK"

contains

    logical function file_exists(f)
        character(len=*), intent(in) :: f
        inquire(file=f, exist=file_exists)
    end function file_exists

    function st(t0, t1) result(s)
        real(dp), intent(in) :: t0, t1
        type(step_time_t) :: s
        s = step_time_t(t0=t0, t1=t1, t_eval=t1, dt=DTT, dtt=DTT)
    end function st

    subroutine advance(glob, reg, H, t)
        ! One coupling step of the core: global update, regional update, reduce, scatter.
        type(geo_slot_t), intent(inout) :: glob, reg
        real(sp),         intent(in)    :: H(:,:)
        real(dp),         intent(in)    :: t
        call bsl_update_global(glob, t)
        call to_regional(glob, reg, t)
        call put_H(reg, H, t - DTT)
        call geo_update(reg, st(t - DTT, t))
        call reduce(glob, reg, t)
        call geo_scatter(reg, st(t - DTT, t))
    end subroutine advance

    subroutine bsl_update_global(glob, t)
        type(geo_slot_t), intent(inout) :: glob
        real(dp),         intent(in)    :: t
        call geo_update(glob, st(t - DTT, t))
    end subroutine bsl_update_global

    subroutine to_regional(glob, reg, t)
        ! exchange_globe: bsl and A_ocn into the regional imports.
        type(geo_slot_t), intent(inout) :: glob, reg
        real(dp),         intent(in)    :: t
        real(dp) :: v
        call fields_get(glob%c%exports, glob%c%h_bsl, v)
        call fields_put(reg%c%imports, reg%c%hi_bsl, v, t1=t)
        call fields_get(glob%c%exports, glob%c%h_A_ocn, v)
        call fields_put(reg%c%imports, reg%c%hi_A_ocn, v, t1=t)
    end subroutine to_regional

    subroutine put_ice_exports(reg, t)
        ! The ice exports the core exchanges before PASS_GEOMETRY (Yelmo's reference
        ! and initial geometry; here the test's topography, sp-valued).
        type(geo_slot_t), intent(inout) :: reg
        real(dp),         intent(in)    :: t
        call fields_put(reg%c%imports, reg%c%hi_z_bed_ref, ctxr%ref%z_bed, t1=t)
        call fields_put(reg%c%imports, reg%c%hi_H_ice_ref, ctxr%ref%H_ice, t1=t)
        call fields_put(reg%c%imports, reg%c%hi_H_ice,     real(0.9_sp*real(ctxr%ref%H_ice, sp), dp), t1=t)   ! current /= reference
    end subroutine put_ice_exports

    subroutine put_H(reg, H, t)
        type(geo_slot_t), intent(inout) :: reg
        real(sp),         intent(in)    :: H(:,:)
        real(dp),         intent(in)    :: t
        call fields_put(reg%c%imports, reg%c%hi_H_ice, real(H, dp), t1=t)
    end subroutine put_H

    subroutine reduce(glob, reg, t)
        ! gather_into (one domain) and geo_reduce, then the reduced bsl back.
        type(geo_slot_t), intent(inout) :: glob, reg
        real(dp),         intent(in)    :: t
        real(dp) :: v
        call fields_get(reg%c%exports, reg%c%h_dV_bsl, v)
        call fields_put(glob%c%imports, glob%c%hi_dV_bsl, v, t1=t)
        call geo_reduce(glob, st(t - DTT, t))
        call fields_get(glob%c%exports, glob%c%h_bsl, v)
        call fields_put(reg%c%imports, reg%c%hi_bsl, v, t1=t)
    end subroutine reduce

    subroutine compare_regional(reg, isos, bsl, what)
        type(geo_slot_t), intent(inout) :: reg
        type(isos_class), intent(in)    :: isos
        type(bsl_class),  intent(in)    :: bsl
        character(len=*), intent(in)    :: what

        real(dp) :: a(NX,NY), b(NX,NY), s
        real(sp) :: ref(NX,NY)

        call fields_get(reg%c%exports, reg%c%h_dz_bed, a)
        b = real(isos%out%w + isos%out%we, dp)
        call check(all(a == b), what//": dz_bed")
        call fields_get(reg%c%exports, reg%c%h_dz_ss, a)
        b = real(isos%out%z_ss - isos%now%bsl, dp)
        call check(all(a == b), what//": dz_ss")
        call fields_get(reg%c%exports, reg%c%h_z_bed_ref, a)
        b = real(isos%out%z_bed - isos%out%w - isos%out%we, dp)
        call check(all(a == b), what//": z_bed_ref")
        ref = real(ctxr%ref%z_bed, sp)
        if (what == "cold start") then      ! exact at the start, round-off after
            call check(all(real(a, sp) == ref), what//": z_bed_ref is the core's reference")
        else
            call check(maxval(abs(a - real(ref, dp))) < 1.0e-2_dp, what//": z_bed_ref is the core's reference")
        end if
        call fields_get(reg%c%exports, reg%c%h_dz_bed_dt, a)
        b = real(isos%out%dwdt, dp)
        call check(all(a == b), what//": dz_bed_dt")
        call fields_get(reg%c%exports, reg%c%h_bsl, s)
        call check(s == real(isos%now%bsl, dp), what//": bsl (regional)")
        call check(real(s, sp) == bsl%bsl_now, what//": bsl against the shared bsl")
        call fields_get(reg%c%exports, reg%c%h_dV_bsl, s)
        call check(s == real(isos%now%deltaV_bsl, dp), what//": dV_bsl")
        write(*,'(a,a,es12.4)') "test_geo_fastisostasy: ", what//": bitwise equal, bsl = ", real(bsl%bsl_now, dp)
    end subroutine compare_regional

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            write(*,'(a)') "test_geo_fastisostasy: FAILED ("//what//")"
            error stop 1
        end if
    end subroutine check

    subroutine cnst_load()
        call phys_const_load(ctxr%cnst, "fesm-utils/par/phys_const_earth.nml", "Earth")
        ctxg%cnst = ctxr%cnst
    end subroutine cnst_load

    subroutine contexts()
        integer :: i, j

        ctxr%par_file  = PAR
        ctxr%instance  = "geo_t"
        ctxr%component = "geo"
        ctxr%backend   = "fastisostasy"
        ctxr%domain    = "test"
        ctxr%grid_name = "TEST"
        ctxr%scope     = "regional"
        ctxr%time      = 0.0_dp
        ctxr%time_rel  = 0.0_dp
        ctxr%dtt       = DTT
        call grid_init(ctxr%grid, name="TEST", mtype="cartesian", units="km", dx=32.0_dp, nx=NX, &
                       dy=32.0_dp, ny=NY, x0=-512.0_dp, y0=-512.0_dp)
        allocate(ctxr%ref%z_bed(NX,NY), ctxr%ref%H_ice(NX,NY))
        do j = 1, NY
            do i = 1, NX
                ctxr%ref%z_bed(i,j) = real(real(-300.0_sp + 25.0_sp*real(i - 16, sp), sp), dp)
                ctxr%ref%H_ice(i,j) = 0.0_dp
                if (abs(i - 16) < 6 .and. abs(j - 16) < 6) &
                    ctxr%ref%H_ice(i,j) = real(1500.0_sp - 100.0_sp*real(abs(i - 16) + abs(j - 16), sp), dp)
            end do
        end do

        ctxg%par_file  = PAR
        ctxg%instance  = "bsl_t"
        ctxg%component = "geo"
        ctxg%backend   = "fastisostasy"
        ctxg%domain    = ""
        ctxg%grid_name = ""
        ctxg%scope     = "global"
        ctxg%time      = 0.0_dp
        ctxg%time_rel  = 0.0_dp
        ctxg%dtt       = DTT
    end subroutine contexts

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&geo_t", ' backend = "fastisostasy"', ' isos = "isos"', "/", &
                       "&bsl_t", ' backend = "fastisostasy"', ' bsl = "barysealevel"', "/", &
                       "&isos", &
                       " heterogeneous_ssh = .false.", " interactive_sealevel = .false.", " include_elastic = .true.", &
                       " correct_compression = .false.", " correct_distortion = .false.", " method = 3", &
                       " dt_diagnostics = 10.0", " pad = 0.0", ' mantle = "uniform"', ' lithosphere = "uniform"', &
                       ' layering = "equalized"', ' lumping = "freqdomain"', ' viscosity_scaling_method = "scale"', &
                       " viscosity_scaling = 1.0", " nl = 1", " zl = 88.0", " viscosities = 1e21", &
                       " ref_viscosity = 1e21", " tau = 3000.0", " rheo_smoothing_radius = 40.0", &
                       " litho_thickness_scaling = 1.0", " E = 66.0", " nu = 0.28", " rho_litho = 0.0", &
                       " r_earth = 6.378e3", " m_earth = 5.972e24", ' restart = "None"', ' mask_file = "None"', &
                       ' rheology_file = "None"', ' dt_method = "tsit54"', " dt_init = 1.0", " dt_min = 1e-10", &
                       " dt_max = 10.0", " rtol = 1e-4", " atol = 1e-2", " q = 1.0", "/", &
                       "&barysealevel", ' method = "fastiso"', " bsl_init = 0.0", ' sl_path = "None"', &
                       ' sl_name = "none"', " A_ocean_pd = 3.625e14", ' A_ocean_path = "None"', ' restart = "None"', "/"
        close(u)
    end subroutine par_write

end program test_geo_fastisostasy
