program test_ice_plastic
    ! The ice backend `plastic` through the ice wrappers (needs plastic=1; skipped otherwise), on a
    ! toy grid with a flat bed above sea level (no geo: the reference bed, z_sl = 0) and a margin
    ! history of two disks (radius R1 at t = 0, R2 at t = 1000 yr), tau_b constant, mode "flat":
    !   * PASS_INIT: the thickness at the centre is the plastic profile sqrt(2 tau_b/(rho_i g) d) of
    !     the distance d to the margin (within 10 %: first order), the ice is the disk, grounded, f_ice = 1 on it;
    !     the static exports (footprint = the larger disk, H_ice_pd, z_bed_ref) are published;
    !   * an update to t = 500 yr: the margin halfway (area), dHidt = (H - H0)/dt, mb_net = rho_i
    !     dHidt, smb = bmb = cmb = 0, M_ice = rho_i sum(H) dx dy;
    !   * a restart round trip: written at 500, read into a second instance, both stepped to 1000:
    !     the same thickness bit for bit (the warm start comes back);
    !   * the flotation binding is the backend's rule;
    !   * taub_method = "coulomb" (mode "surface", one region, H_N >> H: the wedge): the centre near c_f d (first
    !     order: within 15 %), tau_b of
    !     the library equal to c_f rho_i g N_eff, the same tau_b after a restart (PASS_RESTART recomputes it);
    !   * taub_lambda_bed = "file" (lambda_min 0.5, z_lambda of {grid_name}_lambda_bed.nc): the backend's lambda and mu
    !     bit for bit those of the calibration programs' path (taub_set_coulomb, taub_read_lambda_bed,
    !     taub_set_lambda_fixed, taub_mu: tests/plastic_coulomb_{calib,scan,pd}.f90) on the same file;
    !   * margin_patches: a patch (two disks around another centre) in the west half (region id 1 of a map):
    !     the ice is the patch's there and the base's elsewhere, the footprint the base's and the patch's in
    !     the west.
    ! Run from the repo root (files in tmp/test_ice_plastic).

    use precision,       only: dp
    use coords,          only: grid_init
    use ncio,            only: nc_create, nc_write_dim, nc_write
    use fesm_fields,     only: fields_get
    use kryos_component, only: init_ctx_t, census_t, step_time_t, PASS_INIT, PASS_GEOMETRY, PASS_SPINUP, PASS_RESTART
    use kryos_ice,       only: ice_slot_t, ice_create, ice_init, ice_init_state, ice_update, ice_restart_write, &
                               ice_restart_read, ice_flotation, ice_end
    use kryos_ice_plastic, only: plastic_ice_t
    use kryos_plastic_taub, only: taub_t, taub_set_coulomb, taub_read_lambda_bed, taub_set_lambda_fixed, taub_mu

    implicit none

    character(len=*), parameter :: dir = "tmp/test_ice_plastic"
    integer,  parameter :: N = 61
    real(dp), parameter :: DXKM = 10.0_dp, R1 = 10.0_dp, R2 = 20.0_dp     ! [km], radii [cells]
    real(dp), parameter :: TAU = 5.0e4_dp, RHO_I = 910.0_dp, RHO_W = 1028.0_dp, G = 9.81_dp, CF = 0.01_dp
    integer :: nfail = 0

    call execute_command_line("mkdir -p "//dir)
    call write_inputs()
    call test_backend()
    call test_coulomb()
    call test_lambda_file()
    call test_patch()

    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_ice_plastic: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_ice_plastic: OK"

contains

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    function radius() result(r)
        real(dp) :: r(N, N)
        integer :: i, j
        do j = 1, N
            do i = 1, N
                r(i,j) = sqrt(real((i - 31)**2 + (j - 31)**2, dp))
            end do
        end do
    end function radius

    subroutine write_inputs()
        ! The margin history (H_ice on xc, yc, time in years) and the par file.
        real(dp) :: H(N, N, 2), x(N), r(N, N)
        integer  :: i, j, u

        r = radius()
        x = [((i - 31)*DXKM, i = 1, N)]
        H(:,:,1) = merge(1000.0_dp, 0.0_dp, r <= R1)
        H(:,:,2) = merge(1000.0_dp, 0.0_dp, r <= R2)
        call nc_create(dir//"/margin.nc")
        call nc_write_dim(dir//"/margin.nc", "xc", x=x, units="km")
        call nc_write_dim(dir//"/margin.nc", "yc", x=x, units="km")
        call nc_write_dim(dir//"/margin.nc", "time", x=[0.0_dp, 1000.0_dp], units="years")
        call nc_write(dir//"/margin.nc", "H_ice", H, dim1="xc", dim2="yc", dim3="time")

        open(newunit=u, file=dir//"/par.nml", status="replace", action="write")
        write(u,'(a)') '&ice_test'
        write(u,'(a)') '    backend      = "plastic"'
        write(u,'(a)') '    mode         = "flat"'
        write(u,'(a)') '    margin_file  = "'//dir//'/margin.nc"'
        write(u,'(a)') '    margin_var   = "H_ice"'
        write(u,'(a)') '    margin_H_min = 1.0'
        write(u,'(a)') '    taub_method  = "const"'
        write(u,'(a)') '    taub_const   = 5.0e4'
        write(u,'(a)') '/'
        close(u)

        ! a margin patch: disks of radius 5 and 15 cells around (21, 31) at 0 and 1000 yr, in the west half
        do j = 1, N
            do i = 1, N
                H(i,j,1) = merge(1000.0_dp, 0.0_dp, (i - 21)**2 + (j - 31)**2 <= 25)
                H(i,j,2) = merge(1000.0_dp, 0.0_dp, (i - 21)**2 + (j - 31)**2 <= 225)
            end do
        end do
        call nc_create(dir//"/patch.nc")
        call nc_write_dim(dir//"/patch.nc", "xc", x=x, units="km")
        call nc_write_dim(dir//"/patch.nc", "yc", x=x, units="km")
        call nc_write_dim(dir//"/patch.nc", "time", x=[0.0_dp, 1000.0_dp], units="years")
        call nc_write(dir//"/patch.nc", "H_ice", H, dim1="xc", dim2="yc", dim3="time")
        call nc_create(dir//"/halves.nc")
        call nc_write_dim(dir//"/halves.nc", "xc", x=x, units="km")
        call nc_write_dim(dir//"/halves.nc", "yc", x=x, units="km")
        call nc_write(dir//"/halves.nc", "half", reshape([((merge(1, 2, i <= 31), i = 1, N), j = 1, N)], [N, N]), &
                      dim1="xc", dim2="yc")
        open(newunit=u, file=dir//"/par_patch.nml", status="replace", action="write")
        write(u,'(a)') '&ice_test'
        write(u,'(a)') '    backend        = "plastic"'
        write(u,'(a)') '    mode           = "flat"'
        write(u,'(a)') '    margin_file    = "'//dir//'/margin.nc"'
        write(u,'(a)') '    margin_var     = "H_ice"'
        write(u,'(a)') '    margin_H_min   = 1.0'
        write(u,'(a)') '    margin_patches = "patch_w"'
        write(u,'(a)') '    taub_method    = "const"'
        write(u,'(a)') '    taub_const     = 5.0e4'
        write(u,'(a)') '/'
        write(u,'(a)') '&patch_w'
        write(u,'(a)') '    file        = "'//dir//'/patch.nc"'
        write(u,'(a)') '    var         = "H_ice"'
        write(u,'(a)') '    H_min       = 1.0'
        write(u,'(a)') '    region_file = "'//dir//'/halves.nc"'
        write(u,'(a)') '    region_var  = "half"'
        write(u,'(a)') '    region_ids  = "1"'
        write(u,'(a)') '/'
        close(u)

        ! the Coulomb law: one region everywhere, c_f = 0.01, the wedge (H_N >> H), no bed factor, no smoothing
        call nc_create(dir//"/regions.nc")
        call nc_write_dim(dir//"/regions.nc", "xc", x=x, units="km")
        call nc_write_dim(dir//"/regions.nc", "yc", x=x, units="km")
        call nc_write(dir//"/regions.nc", "region", reshape([(1, i = 1, N*N)], [N, N]), dim1="xc", dim2="yc")
        call nc_create(dir//"/cf.nc")
        call nc_write_dim(dir//"/cf.nc", "region", x=[1.0_dp])
        call nc_write(dir//"/cf.nc", "region_id", [1], dim1="region")
        call nc_write(dir//"/cf.nc", "log_cf", [log(CF)], dim1="region")
        open(newunit=u, file=dir//"/par_coulomb.nml", status="replace", action="write")
        write(u,'(a)') '&ice_test'
        write(u,'(a)') '    backend          = "plastic"'
        write(u,'(a)') '    mode             = "surface"'
        write(u,'(a)') '    margin_file      = "'//dir//'/margin.nc"'
        write(u,'(a)') '    margin_var       = "H_ice"'
        write(u,'(a)') '    margin_H_min     = 1.0'
        write(u,'(a)') '    taub_method      = "coulomb"'
        write(u,'(a)') '    taub_region_file = "'//dir//'/regions.nc"'
        write(u,'(a)') '    taub_cf_file     = "'//dir//'/cf.nc"'
        write(u,'(a)') '    taub_HN          = 1.0e9'
        write(u,'(a)') '    taub_smooth_km   = 0.0'
        write(u,'(a)') '    taub_lambda_file = "'//dir//'/{grid_name}_lambda_bed.nc"'
        write(u,'(a)') '/'
        close(u)
        ! the bed of lambda (taub_lambda_bed = "file", the default): from -485 m in the west to 535 m in the east
        call nc_create(dir//"/TEST_lambda_bed.nc")
        call nc_write_dim(dir//"/TEST_lambda_bed.nc", "xc", x=x, units="km")
        call nc_write_dim(dir//"/TEST_lambda_bed.nc", "yc", x=x, units="km")
        call nc_write(dir//"/TEST_lambda_bed.nc", "z_lambda", lambda_bed_z(), dim1="xc", dim2="yc", units="m")
        open(newunit=u, file=dir//"/par_coulomb_lam.nml", status="replace", action="write")
        write(u,'(a)') '&ice_test'
        write(u,'(a)') '    backend          = "plastic"'
        write(u,'(a)') '    mode             = "surface"'
        write(u,'(a)') '    margin_file      = "'//dir//'/margin.nc"'
        write(u,'(a)') '    margin_var       = "H_ice"'
        write(u,'(a)') '    margin_H_min     = 1.0'
        write(u,'(a)') '    taub_method      = "coulomb"'
        write(u,'(a)') '    taub_region_file = "'//dir//'/regions.nc"'
        write(u,'(a)') '    taub_cf_file     = "'//dir//'/cf.nc"'
        write(u,'(a)') '    taub_HN          = 1.0e9'
        write(u,'(a)') '    taub_smooth_km   = 0.0'
        write(u,'(a)') '    taub_lam_min     = 0.5'
        write(u,'(a)') '    taub_lambda_bed  = "file"'
        write(u,'(a)') '    taub_lambda_file = "'//dir//'/{grid_name}_lambda_bed.nc"'
        write(u,'(a)') '/'
        close(u)
    end subroutine write_inputs

    subroutine test_coulomb()
        type(ice_slot_t)  :: ice, ice2
        type(step_time_t) :: st
        real(dp) :: H(N, N), r(N, N), Hc, e
        integer  :: i, j

        r = radius()
        write(*,"(a)") "taub_method = coulomb"
        call make_ice(ice, "par_coulomb.nml")
        st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=0)
        call ice_init_state(ice, st, PASS_INIT)
        H = export(ice, ice%c%h_H_ice)
        Hc = CF*(R1 + 0.5_dp)*DXKM*1.0e3_dp
        write(*,"(a,f9.2,a,f9.2,a)") "  centre: H = ", H(31,31), " m, the wedge c_f d ", Hc, " m"
        ! the first-order distance on a 10-cell disk is 12 % short (the Nye check above: -6 % in H ~ sqrt(d))
        call check("centre: the wedge c_f d within 15 %", abs(H(31,31) - Hc) < 0.15_dp*Hc)
        call check("the ice is the disk of t = 0", all((H > 0.0_dp) .eqv. (r <= R1)))
        select type (c => ice%c)
        type is (plastic_ice_t)
            e = 0.0_dp
            do j = 1, N
                do i = 1, N
                    if (c%mask(i,j)) e = max(e, abs(c%tau_b(i,j) - c%mu(i,j)*RHO_I*G*H(i,j)*1.0e9_dp/(H(i,j) + 1.0e9_dp)))
                end do
            end do
            call check("tau_b = c_f rho_i g N_eff (the library's)", e <= 1.0e-9_dp*CF*RHO_I*G*maxval(H) .and. &
                       all(c%tau_b == c%ps%tau_b) .and. all(abs(c%mu - CF) < 1.0e-15_dp))
        end select
        st = step_time_t(t0=0.0_dp, t1=500.0_dp, t_eval=500.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=1)
        call ice_update(ice, st)
        call ice_restart_write(ice, dir, 500.0_dp)
        call make_ice(ice2, "par_coulomb.nml")
        call ice_restart_read(ice2, dir, 500.0_dp)
        st = step_time_t(t0=500.0_dp, t1=500.0_dp, t_eval=500.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=1)
        call ice_init_state(ice2, st, PASS_RESTART)
        select type (c => ice%c)
        type is (plastic_ice_t)
            select type (c2 => ice2%c)
            type is (plastic_ice_t)
                call check("restart: the same tau_b (recomputed on the restored H)", all(c%tau_b == c2%tau_b))
            end select
        end select
        call ice_end(ice)
        call ice_end(ice2)
    end subroutine test_coulomb

    function lambda_bed_z() result(z)
        real(dp) :: z(N, N)
        integer :: i, j
        do j = 1, N
            do i = 1, N
                z(i,j) = -500.0_dp + 15.0_dp*i + 2.0_dp*j
            end do
        end do
    end function lambda_bed_z

    subroutine test_lambda_file()
        ! the backend's lambda and mu against the calibration programs' path on the same file
        type(ice_slot_t)  :: ice
        type(step_time_t) :: st
        type(taub_t) :: tb
        real(dp) :: mu(N, N)
        integer  :: k

        write(*,"(a)") "taub_lambda_bed = file"
        call make_ice(ice, "par_coulomb_lam.nml")
        st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=0)
        call ice_init_state(ice, st, PASS_INIT)
        call taub_set_coulomb(tb, reshape([(1, k = 1, N*N)], [N, N]), [1], [log(CF)], 1.0_dp, 1.0e9_dp, 0.05_dp, &
                              0.5_dp, -300.0_dp, 200.0_dp, "regions", "calibration c_f")
        call taub_set_lambda_fixed(tb, taub_read_lambda_bed(dir//"/TEST_lambda_bed.nc", N, N))
        call taub_mu(tb, 0.0_dp*mu + 100.0_dp, 0.0_dp*mu, mu)
        select type (c => ice%c)
        type is (plastic_ice_t)
            call check("lambda bit for bit the calibration programs' (and spans 0.5 to 1)", &
                       all(c%tb%lam == tb%lam) .and. minval(c%tb%lam) == 0.5_dp .and. maxval(c%tb%lam) == 1.0_dp .and. &
                       any(c%tb%lam > 0.5_dp .and. c%tb%lam < 1.0_dp))
            call check("mu = c_f lambda bit for bit the calibration programs'", all(c%mu == mu))
        end select
        call ice_end(ice)
    end subroutine test_lambda_file

    subroutine test_patch()
        type(ice_slot_t)  :: ice
        type(step_time_t) :: st
        real(dp) :: H(N, N), r(N, N), rp(N, N)
        logical  :: west(N, N)
        integer  :: i, j

        write(*,"(a)") "margin_patches"
        r = radius()
        do j = 1, N
            do i = 1, N
                rp(i,j) = sqrt(real((i - 21)**2 + (j - 31)**2, dp))
                west(i,j) = i <= 31
            end do
        end do
        call make_ice(ice, "par_patch.nml")
        st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=0)
        call ice_init_state(ice, st, PASS_INIT)
        H = export(ice, ice%c%h_H_ice)
        call check("patch: the patch's disk in the west, the base's in the east", &
                   all((H > 0.0_dp) .eqv. merge(rp <= 5.0_dp, r <= R1, west)))
        call check("patch: footprint = the base's, or the patch's in the west", &
                   all((export(ice, ice%c%h_footprint) > 0.5_dp) .eqv. (r <= R2 .or. (west .and. rp <= 15.0_dp))))
        st = step_time_t(t0=0.0_dp, t1=500.0_dp, t_eval=500.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=1)
        call ice_update(ice, st)
        H = export(ice, ice%c%h_H_ice)
        call check("patch: at 500 yr the east half of the base's disk of radius 15 (within one cell)", &
                   all(west .or. (H > 0.0_dp) .or. r > 14.0_dp) .and. all(west .or. .not. (H > 0.0_dp) .or. r < 16.0_dp))
        call check("patch: at 500 yr the west the patch's disk of radius 10 (within one cell)", &
                   all(.not. west .or. (H > 0.0_dp) .or. rp > 9.0_dp) .and. all(.not. west .or. .not. (H > 0.0_dp) .or. rp < 11.0_dp))
        call ice_end(ice)
    end subroutine test_patch

    subroutine make_ice(slot, par)
        type(ice_slot_t), intent(out) :: slot
        character(len=*), intent(in), optional :: par
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        ctx%instance  = "ice_test"
        ctx%component = "ice"
        ctx%backend   = "plastic"
        ctx%domain    = "Test"
        ctx%grid_name = "TEST"
        ctx%par_file  = dir//"/par.nml"
        if (present(par)) ctx%par_file = dir//"/"//par
        call grid_init(ctx%grid, name="TEST", mtype="cartesian", units="km", dx=DXKM, nx=N, dy=DXKM, ny=N, &
                       x0=-30.0_dp*DXKM, y0=-30.0_dp*DXKM)
        allocate(ctx%ref%z_bed(N, N), ctx%ref%H_ice(N, N), ctx%ref%z_srf(N, N), ctx%ref%z_bed_sd(N, N))
        ctx%ref%z_bed = 100.0_dp
        ctx%ref%H_ice = 0.0_dp
        ctx%ref%z_srf = 100.0_dp
        ctx%ref%z_bed_sd = 0.0_dp
        call ice_create(slot, ctx, census)
        call ice_init(slot, ctx)
    end subroutine make_ice

    function export(slot, h) result(a)
        type(ice_slot_t), intent(in) :: slot
        integer,          intent(in) :: h
        real(dp) :: a(N, N)
        call fields_get(slot%c%exports, h, a)
    end function export

    subroutine test_backend()
        type(ice_slot_t)  :: ice, ice2
        type(step_time_t) :: st
        real(dp) :: H0(N, N), H(N, N), r(N, N), a(N, N), Hc, m
        real(dp) :: zb(N, N), zs(N, N), f(N, N)

        r = radius()
        write(*,"(a)") "cold start"
        call make_ice(ice)
        st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=0)
        call ice_init_state(ice, st, PASS_INIT)
        call ice_init_state(ice, st, PASS_GEOMETRY)
        call ice_init_state(ice, st, PASS_SPINUP)
        H0 = export(ice, ice%c%h_H_ice)
        Hc = sqrt(2.0_dp*TAU/(RHO_I*G)*(R1 + 0.5_dp)*DXKM*1.0e3_dp)
        write(*,"(a,f9.2,a,f9.2,a)") "  centre: H = ", H0(31,31), " m, the profile ", Hc, " m"
        ! the first-order scheme on a 10-cell disk is below the profile: -6 % with the margin at the faces
        ! (margin_position "sdist" at a slice's time, or "face"), -4 % at the ice-free cell centres ("centre")
        call check("centre: the plastic profile within 10 %", abs(H0(31,31) - Hc) < 0.10_dp*Hc)
        call check("the ice is the disk of t = 0", all((H0 > 0.0_dp) .eqv. (r <= R1)))
        call check("f_ice and f_grnd are 1 on it", all(pack(export(ice, ice%c%h_f_ice), r <= R1) == 1.0_dp) .and. &
                   all(pack(export(ice, ice%c%h_f_grnd), r <= R1) == 1.0_dp))
        call check("footprint: the disk of t = 1000", all((export(ice, ice%c%h_footprint) > 0.5_dp) .eqv. (r <= R2)))
        call check("z_bed_ref, H_ice_pd: the reference", all(export(ice, ice%c%h_z_bed_ref) == 100.0_dp) .and. &
                   all(export(ice, ice%c%h_H_ice_pd) == 0.0_dp))
        call check("z_bed, z_sl: the reference bed and 0 (no geo)", all(export(ice, ice%c%h_z_bed) == 100.0_dp) .and. &
                   all(export(ice, ice%c%h_z_sl) == 0.0_dp))

        write(*,"(a)") "update to 500 yr"
        st = step_time_t(t0=0.0_dp, t1=500.0_dp, t_eval=500.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=1)
        call ice_update(ice, st)
        H = export(ice, ice%c%h_H_ice)
        call check("the margin halfway: area of radius 15 cells within 5 %", &
                   abs(count(H > 0.0_dp) - acos(-1.0_dp)*15.0_dp**2) < 0.05_dp*acos(-1.0_dp)*15.0_dp**2)
        a = export(ice, ice%c%h_dHidt)
        call check("dHidt = (H - H0)/dt", all(abs(a - (H - H0)/500.0_dp) <= 1.0e-12_dp*max(1.0_dp, abs(a))))
        call check("mb_net = rho_i dHidt", all(abs(export(ice, ice%c%h_mb_net) - RHO_I*a) <= 1.0e-9_dp*max(1.0_dp, RHO_I*abs(a))))
        call check("smb = bmb = cmb = 0", all(export(ice, ice%c%h_smb) == 0.0_dp) .and. &
                   all(export(ice, ice%c%h_bmb) == 0.0_dp) .and. all(export(ice, ice%c%h_cmb) == 0.0_dp))
        call fields_get(ice%c%exports, ice%c%h_M_ice, m)
        call check("M_ice = rho_i sum(H) dx dy", abs(m - RHO_I*sum(H)*(DXKM*1.0e3_dp)**2*1.0e-12_dp) <= 1.0e-12_dp*m)

        write(*,"(a)") "restart round trip"
        call ice_restart_write(ice, dir, 500.0_dp)
        call make_ice(ice2)
        call ice_restart_read(ice2, dir, 500.0_dp)
        st = step_time_t(t0=500.0_dp, t1=500.0_dp, t_eval=500.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=1)
        call ice_init_state(ice2, st, PASS_RESTART)
        select type (c => ice%c)
        type is (plastic_ice_t)
            select type (c2 => ice2%c)
            type is (plastic_ice_t)
                call check("H and w restored bit for bit", all(c%ps%H == c2%ps%H) .and. all(c%ps%w == c2%ps%w))
            end select
        end select
        st = step_time_t(t0=500.0_dp, t1=1000.0_dp, t_eval=1000.0_dp, dt=500.0_dp, dtt=500.0_dp, tick=2)
        call ice_update(ice, st)
        call ice_update(ice2, st)
        call check("stepped on: the same thickness", all(export(ice, ice%c%h_H_ice) == export(ice2, ice2%c%h_H_ice)))
        call check("and the same dHidt", all(export(ice, ice%c%h_dHidt) == export(ice2, ice2%c%h_dHidt)))

        write(*,"(a)") "flotation binding"
        zb = -500.0_dp + 10.0_dp*r
        zs = 0.0_dp
        H  = 400.0_dp + 5.0_dp*r
        call ice_flotation(ice, H, zb, zs, f)
        call check("grounded where rho_i H >= rho_w (z_sl - z_bed)", all(f == merge(1.0_dp, 0.0_dp, RHO_I*H >= RHO_W*(zs - zb))))
        call check("both occur", any(f == 1.0_dp) .and. any(f == 0.0_dp))

        call ice_end(ice)
        call ice_end(ice2)
    end subroutine test_backend

end program test_ice_plastic
