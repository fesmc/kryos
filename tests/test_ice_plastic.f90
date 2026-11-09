program test_ice_plastic
    ! The ice backend `plastic` through the ice wrappers (needs plastic=1; skipped otherwise), on a
    ! toy grid with a flat bed above sea level (no geo: the reference bed, z_sl = 0) and a margin
    ! history of two disks (radius R1 at t = 0, R2 at t = 1000 yr), tau_b constant, mode "flat":
    !   * PASS_INIT: the thickness at the centre is the plastic profile sqrt(2 tau_b/(rho_i g) d) of
    !     the distance d to the margin (within 5 %), the ice is the disk, grounded, f_ice = 1 on it;
    !     the static exports (footprint = the larger disk, H_ice_pd, z_bed_ref) are published;
    !   * an update to t = 500 yr: the margin halfway (area), dHidt = (H - H0)/dt, mb_net = rho_i
    !     dHidt, smb = bmb = cmb = 0, M_ice = rho_i sum(H) dx dy;
    !   * a restart round trip: written at 500, read into a second instance, both stepped to 1000:
    !     the same thickness bit for bit (the warm start comes back);
    !   * the flotation binding is the backend's rule.
    ! Run from the repo root (files in tmp/test_ice_plastic).

    use precision,       only: dp
    use coords,          only: grid_init
    use ncio,            only: nc_create, nc_write_dim, nc_write
    use fesm_fields,     only: fields_get
    use kryos_component, only: init_ctx_t, census_t, step_time_t, PASS_INIT, PASS_GEOMETRY, PASS_SPINUP, PASS_RESTART
    use kryos_ice,       only: ice_slot_t, ice_create, ice_init, ice_init_state, ice_update, ice_restart_write, &
                               ice_restart_read, ice_flotation, ice_end
    use kryos_ice_plastic, only: plastic_ice_t

    implicit none

    character(len=*), parameter :: dir = "tmp/test_ice_plastic"
    integer,  parameter :: N = 61
    real(dp), parameter :: DXKM = 10.0_dp, R1 = 10.0_dp, R2 = 20.0_dp     ! [km], radii [cells]
    real(dp), parameter :: TAU = 5.0e4_dp, RHO_I = 910.0_dp, RHO_W = 1028.0_dp, G = 9.81_dp
    integer :: nfail = 0

    call execute_command_line("mkdir -p "//dir)
    call write_inputs()
    call test_backend()

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
        integer  :: i, u

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
    end subroutine write_inputs

    subroutine make_ice(slot)
        type(ice_slot_t), intent(out) :: slot
        type(init_ctx_t) :: ctx
        type(census_t)   :: census

        ctx%instance  = "ice_test"
        ctx%component = "ice"
        ctx%backend   = "plastic"
        ctx%domain    = "Test"
        ctx%grid_name = "TEST"
        ctx%par_file  = dir//"/par.nml"
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
        call check("centre: the plastic profile within 5 %", abs(H0(31,31) - Hc) < 0.05_dp*Hc)
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
