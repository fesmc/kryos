program test_plastic_coulomb
    ! The library-free parts of the Coulomb yield law of the ice backend `plastic` (kryos_plastic_bed,
    ! kryos_plastic_taub "coulomb", kryos_plastic_friction), so this runs with plastic=0:
    !   * lambda: lam_min at a bed zb0 or more below sea level, 1 at zb1 or above, linear between, relative to the
    !     sea level of the time; 1 everywhere with lam_min = 1;
    !   * the sliding law coulomb_n takes lambda from kryos_plastic_bed bit for bit (against its former inline
    !     expression);
    !   * mu = c_f(region) lambda: c_f by region id (not position), 0 off the regions, a region id without cells
    !     allowed; the smoothing keeps a uniform mu; from the files (region map, c_f table);
    !   * a fixed lambda (taub_lambda_bed = "relaxed"): of the relaxed bed z_bed + (rho_i/3300) H, independent of the
    !     bed and sea level of the solve;
    !   * a fixed lambda from a file (taub_lambda_bed = "file", taub_read_lambda_bed): bit for bit the lambda of the
    !     file's z_lambda, independent of the bed and sea level of the solve; a missing file, a file without z_lambda
    !     or off the grid stop;
    !   * a region of the map without c_f, p outside [0, 1], H_N <= 0, delta outside (0, 1] and taub_at stop (the
    !     program re-runs itself).
    ! Run from the repo root (files in tmp/test_plastic_coulomb).

    use precision, only: dp
    use ncio,      only: nc_create, nc_write_dim, nc_write
    use kryos_plastic_bed,      only: bed_lambda, relaxed_bed
    use kryos_plastic_friction, only: friction_t, friction_map
    use kryos_plastic_taub,     only: taub_t, taub_set_coulomb, taub_read_coulomb, taub_mu, taub_set_smooth, taub_at, &
                                      taub_set_lambda_fixed, taub_read_lambda_bed, &
                                      TAUB_M_COULOMB

    implicit none

    character(len=*), parameter :: dir = "tmp/test_plastic_coulomb"
    character(len=*), parameter :: errfile = dir//"/abort.err"
    integer :: nfail = 0

    call execute_command_line("mkdir -p "//dir)
    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                ! the case did not abort
    end if

    call test_lambda()
    call test_friction()
    call test_mu()
    call test_files()
    call test_lambda_file()
    call expect_abort("noregion", "region 3 of the map has no c_f")
    call expect_abort("badp", "taub_p must be in [0, 1]")
    call expect_abort("badHN", "taub_HN must be positive")
    call expect_abort("baddelta0", "taub_delta must be in (0, 1]")
    call expect_abort("baddelta1", "taub_delta must be in (0, 1]")
    call expect_abort("taub_at", "has no tau_b before the solve")
    call expect_abort("nolamfile", "the lambda bed file does not exist")
    call expect_abort("nolamvar", "the lambda bed file has no `z_lambda`")
    call expect_abort("lamgrid", "`z_lambda` is not on the instance grid")

    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_plastic_coulomb: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_plastic_coulomb: OK"

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

    subroutine test_lambda()
        real(dp) :: zb(5)
        write(*,"(a)") "bed factor lambda"
        zb = [-1000.0_dp, -300.0_dp, -50.0_dp, 200.0_dp, 900.0_dp]
        call check("lam_min at and below zb0, 1 at and above zb1, linear between", &
                   all(abs(bed_lambda(zb, 0.0_dp, 0.1_dp, -300.0_dp, 200.0_dp) - [0.1_dp, 0.1_dp, 0.55_dp, 1.0_dp, 1.0_dp]) &
                       < 1.0e-14_dp))
        call check("relative to the sea level of the time", &
                   all(bed_lambda(zb - 120.0_dp, -120.0_dp, 0.1_dp, -300.0_dp, 200.0_dp) == &
                       bed_lambda(zb, 0.0_dp, 0.1_dp, -300.0_dp, 200.0_dp)))
        call check("lam_min = 1: 1 everywhere", all(bed_lambda(zb, 0.0_dp, 1.0_dp, -300.0_dp, 200.0_dp) == 1.0_dp))
    end subroutine test_lambda

    subroutine test_friction()
        ! coulomb_n against the expression it had inline before kryos_plastic_bed
        integer, parameter :: n = 40
        type(friction_t) :: fr
        real(dp) :: tau(n, n), H(n, n), zb(n, n), zsl(n, n), fric(n, n), Hf(n, n), Nn(n, n), lam(n, n), ref(n, n)
        logical  :: mask(n, n)
        integer  :: i, j
        write(*,"(a)") "coulomb_n with kryos_plastic_bed"
        do j = 1, n
            do i = 1, n
                zb(i, j)   = -900.0_dp + 47.3_dp*i + 3.1_dp*j
                H(i, j)    = 30.0_dp*j + 7.0_dp*i
                tau(i, j)  = 2.0e4_dp + 1.3e3_dp*i
                mask(i, j) = mod(i + j, 7) /= 0
            end do
        end do
        zsl = -37.5_dp
        fr%kind = "coulomb_n";  fr%Hc = 120.0_dp;  fr%lam_min = 0.13_dp;  fr%zb0 = -450.0_dp;  fr%zb1 = 25.0_dp
        fr%p = 0.7_dp
        call friction_map(fr, 3.0_dp, tau, H, zb, zsl, mask, 910.0_dp, 1028.0_dp, 9.81_dp, fric)
        Hf = 1028.0_dp/910.0_dp*max(zsl - zb, 0.0_dp)
        Nn = 0.0_dp
        where (H > Hf) Nn = 910.0_dp*9.81_dp*H*(1.0_dp - Hf/H)**fr%p
        lam = 1.0_dp
        if (fr%lam_min < 1.0_dp) lam = fr%lam_min + (1.0_dp - fr%lam_min)* &
                                       min(max((zb - zsl - fr%zb0)/(fr%zb1 - fr%zb0), 0.0_dp), 1.0_dp)
        where (mask .and. tau > 0.0_dp)
            ref = tau*lam*Nn/(910.0_dp*9.81_dp*fr%Hc)
        elsewhere
            ref = huge(1.0_dp)
        end where
        call check("coulomb_n coefficient bit for bit as before", all(fric == ref))
    end subroutine test_friction

    subroutine test_mu()
        integer, parameter :: n = 30
        type(taub_t) :: tb
        integer  :: reg(n, n)
        real(dp) :: zb(n, n), zsl(n, n), mu(n, n)
        write(*,"(a)") "mu = c_f lambda"
        reg = 0;  reg(1:10, :) = 4;  reg(11:20, :) = 2
        zb = 500.0_dp;  zb(:, 1:5) = -1000.0_dp
        zsl = 0.0_dp
        call taub_set_coulomb(tb, reg, [2, 9, 4], log([0.01_dp, 0.5_dp, 0.02_dp]), 1.0_dp, 1000.0_dp, 0.05_dp, 0.1_dp, &
                              -300.0_dp, 200.0_dp, "map", "table")
        call check("method coulomb", tb%method == TAUB_M_COULOMB)
        call taub_mu(tb, zb, zsl, mu)
        call check("c_f by region id; 0 off the regions", abs(mu(5, 10) - 0.02_dp) < 1.0e-15_dp .and. &
                   abs(mu(15, 10) - 0.01_dp) < 1.0e-15_dp .and. all(mu(21:, :) == 0.0_dp))
        call check("times lambda (lam_min on a deep bed)", abs(mu(5, 2) - 0.002_dp) < 1.0e-15_dp)
        zb(:, 1:5) = 500.0_dp
        reg(11:20, :) = 4
        call taub_set_coulomb(tb, reg, [4], [log(0.02_dp)], 1.0_dp, 1000.0_dp, 0.05_dp, 0.1_dp, -300.0_dp, &
                              200.0_dp, "map", "t")
        call taub_set_smooth(tb, 50.0e3_dp, 10.0e3_dp, 10.0e3_dp)
        call taub_mu(tb, zb, zsl, mu)
        call check("smoothing keeps a uniform mu and 0 off the regions", &
                   all(abs(mu(1:20, :) - 0.02_dp) < 1.0e-15_dp) .and. all(mu(21:, :) == 0.0_dp))
        ! fixed lambda on the relaxed bed: a bed 400 m below sea level under 3300 m of ice relaxes to +510 m
        call check("relaxed bed: z_bed + (rho_i/3300) H", &
                   abs(relaxed_bed(-400.0_dp, 3300.0_dp, 910.0_dp) - 510.0_dp) < 1.0e-12_dp .and. &
                   relaxed_bed(-400.0_dp, -5.0_dp, 910.0_dp) == -400.0_dp)
        call taub_set_coulomb(tb, reg, [4], [log(0.02_dp)], 1.0_dp, 1000.0_dp, 0.05_dp, 0.1_dp, -300.0_dp, &
                              200.0_dp, "map", "t")
        zb = -400.0_dp
        call taub_set_lambda_fixed(tb, relaxed_bed(zb, 0.0_dp*zb + 3300.0_dp, 910.0_dp))
        call taub_mu(tb, -2000.0_dp + 0.0_dp*zb, zsl + 100.0_dp, mu)
        call check("fixed lambda: of the relaxed bed (1 here), not of the solve's deep bed", &
                   all(abs(mu(1:20, :) - 0.02_dp) < 1.0e-15_dp))
    end subroutine test_mu

    subroutine test_files()
        integer, parameter :: nx = 6, ny = 4
        character(len=*), parameter :: fr = dir//"/regions.nc", fc = dir//"/cf.nc"
        type(taub_t) :: tb
        integer  :: reg(nx, ny)
        real(dp) :: zb(nx, ny), zsl(nx, ny), mu(nx, ny)
        write(*,"(a)") "from the files"
        reg = 1;  reg(4:, :) = 2;  reg(1, 1) = 0
        call nc_create(fr)
        call nc_write_dim(fr, "xc", x=[10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 50.0_dp, 60.0_dp], units="km")
        call nc_write_dim(fr, "yc", x=[10.0_dp, 20.0_dp, 30.0_dp, 40.0_dp], units="km")
        call nc_write(fr, "region", reg, dim1="xc", dim2="yc")
        call nc_create(fc)
        call nc_write_dim(fc, "region", x=[1.0_dp, 2.0_dp, 3.0_dp])
        call nc_write(fc, "region_id", [3, 2, 1], dim1="region")
        call nc_write(fc, "log_cf", log([0.03_dp, 0.02_dp, 0.01_dp]), dim1="region")
        call taub_read_coulomb(tb, fr, fc, 6, 4, 0.5_dp, 800.0_dp, 0.1_dp, 1.0_dp, -300.0_dp, 200.0_dp)
        zb = 100.0_dp;  zsl = 0.0_dp
        call taub_mu(tb, zb, zsl, mu)
        call check("region map and c_f table by id", abs(mu(2, 2) - 0.01_dp) < 1.0e-15_dp .and. &
                   abs(mu(5, 3) - 0.02_dp) < 1.0e-15_dp .and. mu(1, 1) == 0.0_dp .and. tb%p == 0.5_dp .and. &
                   tb%H_N == 800.0_dp .and. tb%delta == 0.1_dp)
    end subroutine test_files

    subroutine test_lambda_file()
        ! taub_lambda_bed = "file": lambda of the file's z_lambda, as taub_set_lambda_fixed of the same field
        integer, parameter :: nx = 6, ny = 4
        character(len=*), parameter :: fl = dir//"/lambda_bed.nc"
        type(taub_t) :: tb
        integer  :: reg(nx, ny), i, j
        real(dp) :: z(nx, ny), mu(nx, ny), zb(nx, ny), zsl(nx, ny)
        write(*,"(a)") "lambda from a file"
        do j = 1, ny
            do i = 1, nx
                z(i, j) = -437.3_dp + 113.7_dp*i + 41.9_dp*j
            end do
        end do
        call nc_create(fl)
        call nc_write_dim(fl, "xc", x=[(10.0_dp*i, i = 1, nx)], units="km")
        call nc_write_dim(fl, "yc", x=[(10.0_dp*j, j = 1, ny)], units="km")
        call nc_write(fl, "z_lambda", z, dim1="xc", dim2="yc", units="m")
        reg = 1
        call taub_set_coulomb(tb, reg, [1], [log(0.02_dp)], 1.0_dp, 1000.0_dp, 0.05_dp, 0.3_dp, -300.0_dp, 200.0_dp, &
                              "map", "t")
        call taub_set_lambda_fixed(tb, taub_read_lambda_bed(fl, nx, ny))
        call check("the field of the file", all(taub_read_lambda_bed(fl, nx, ny) == z))
        call check("lambda bit for bit of z_lambda, relative to present sea level", &
                   all(tb%lam == bed_lambda(z, 0.0_dp, 0.3_dp, -300.0_dp, 200.0_dp)) .and. any(tb%lam < 1.0_dp) .and. &
                   any(tb%lam == 1.0_dp) .and. any(tb%lam > 0.3_dp .and. tb%lam < 1.0_dp))
        zb = -2000.0_dp;  zsl = 150.0_dp
        call taub_mu(tb, zb, zsl, mu)
        call check("mu = c_f lambda(z_lambda), not of the solve's bed", all(mu == 0.02_dp*tb%lam))
    end subroutine test_lambda_file

    subroutine expect_abort(case_name, text)
        character(len=*), intent(in) :: case_name, text
        character(len=512) :: exe, line
        character(len=:), allocatable :: all_
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > /dev/null 2> "//errfile, exitstat=es)
        all_ = ""
        open(newunit=u, file=errfile, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(u, "(a)", iostat=ios) line
                if (ios /= 0) exit
                all_ = all_//trim(line)//" "
            end do
            close(u)
        end if
        call check(case_name//": aborts (exit 1) with '"//text//"'", es == 1 .and. index(all_, text) > 0)
    end subroutine expect_abort

    subroutine abort_case()
        ! Child: build the situation that must stop.
        character(len=64) :: arg
        type(taub_t) :: tb
        integer  :: reg(2, 2)
        real(dp) :: tau(2, 2)
        real(dp), allocatable :: z(:,:)

        reg = reshape([1, 3, 0, 1], [2, 2])
        call get_command_argument(1, arg)
        select case (trim(arg))
        case ("noregion")
            call taub_set_coulomb(tb, reg, [1, 2], [0.0_dp, 0.0_dp], 1.0_dp, 1000.0_dp, 0.05_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
        case ("badp")
            call taub_set_coulomb(tb, reg, [1, 3], [0.0_dp, 0.0_dp], 1.5_dp, 1000.0_dp, 0.05_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
        case ("badHN")
            call taub_set_coulomb(tb, reg, [1, 3], [0.0_dp, 0.0_dp], 1.0_dp, 0.0_dp, 0.05_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
        case ("baddelta0")
            call taub_set_coulomb(tb, reg, [1, 3], [0.0_dp, 0.0_dp], 1.0_dp, 1000.0_dp, 0.0_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
        case ("baddelta1")
            call taub_set_coulomb(tb, reg, [1, 3], [0.0_dp, 0.0_dp], 1.0_dp, 1000.0_dp, 1.5_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
        case ("nolamfile")
            z = taub_read_lambda_bed(dir//"/none.nc", 6, 4)
        case ("nolamvar")
            z = taub_read_lambda_bed(dir//"/regions.nc", 6, 4)
        case ("lamgrid")
            z = taub_read_lambda_bed(dir//"/lambda_bed.nc", 5, 4)
        case ("taub_at")
            call taub_set_coulomb(tb, reg, [1, 3], [0.0_dp, 0.0_dp], 1.0_dp, 1000.0_dp, 0.05_dp, 1.0_dp, -300.0_dp, &
                                  200.0_dp, "map", "table")
            call taub_at(tb, 0.0_dp, tau)
        end select
    end subroutine abort_case

end program test_plastic_coulomb
