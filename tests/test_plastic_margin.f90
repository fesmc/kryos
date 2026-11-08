program test_plastic_margin
    ! The inputs of the ice backend `plastic` that need no library (kryos_plastic_margin,
    ! kryos_plastic_taub), so this runs with plastic=0:
    !   * the distance transform is exact: against brute force on an irregular mask, dx /= dy;
    !   * the signed distance of a slice: |sd| >= dx/2, sd > 0 gives the slice back;
    !   * two concentric disks (radius R1 at t1, R2 at t2): the margin of the "sdist" mask at t lies at
    !     R1 + (t - t1)/(t2 - t1) (R2 - R1) to within one cell, the slices come back exactly at their
    !     times, the footprint is the union; "nearest" takes the nearest slice;
    !   * a margin file with "years ago" (and a fill value) in model time, increasing;
    !   * outside the time range: hold, or stop with a named message (the program re-runs itself);
    !   * tau_b: domains with a series (either dimension order in the file), an override table, a
    !     default; a domain without a value stops.
    ! Run from the repo root (files in tmp/test_plastic_margin).

    use precision, only: dp
    use ncio,      only: nc_create, nc_write_dim, nc_write
    use kryos_plastic_margin, only: margin_t, margin_set, margin_init, margin_at, margin_footprint, &
                                    margin_signed_distance, margin_edt2, model_time_axis, MARGIN_SDIST, MARGIN_NEAREST
    use kryos_plastic_taub,   only: taub_t, taub_set_domains, taub_read_domains, taub_at, taub_set_const, &
                                    TAUB_M_DOMAINS

    implicit none

    character(len=*), parameter :: dir = "tmp/test_plastic_margin"
    character(len=*), parameter :: errfile = dir//"/abort.err"
    integer, parameter :: N = 61                     ! disk grid: N x N, centre (31, 31)
    real(dp), parameter :: DX = 1000.0_dp, R1 = 8.0_dp, R2 = 22.0_dp     ! [m], radii [cells]
    integer :: nfail = 0

    call execute_command_line("mkdir -p "//dir)
    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                ! the case did not abort
    end if

    call test_edt()
    call test_sdist_slice()
    call test_disks()
    call test_file()
    call test_out_of_range()
    call test_taub()

    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_plastic_margin: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_plastic_margin: OK"

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

    function irregular(nx, ny) result(m)
        ! A deterministic irregular mask: blobs and isolated cells.
        integer, intent(in) :: nx, ny
        logical :: m(nx, ny)
        integer :: i, j
        do j = 1, ny
            do i = 1, nx
                m(i,j) = mod(7*i + 13*j + i*j, 11) == 0 .or. ((i - 9)**2 + (j - 6)**2 < 10)
            end do
        end do
    end function irregular

    function disk(R) result(m)
        ! Cells whose centre lies within R cells of the centre (31, 31).
        real(dp), intent(in) :: R
        logical :: m(N, N)
        m = radius() <= R
    end function disk

    function radius() result(r)
        ! [cells] distance of each cell centre from the centre (31, 31).
        real(dp) :: r(N, N)
        integer :: i, j
        do j = 1, N
            do i = 1, N
                r(i,j) = sqrt(real((i - 31)**2 + (j - 31)**2, dp))
            end do
        end do
    end function radius

    subroutine test_edt()
        integer, parameter :: nx = 23, ny = 17
        real(dp), parameter :: ddx = 2.0_dp, ddy = 3.0_dp
        logical  :: m(nx, ny)
        real(dp) :: d2(nx, ny), b2(nx, ny)
        integer  :: i, j, k, l

        write(*,"(a)") "distance transform"
        m = irregular(nx, ny)
        d2 = margin_edt2(m, ddx, ddy)
        b2 = huge(1.0_dp)
        do j = 1, ny
            do i = 1, nx
                do l = 1, ny
                    do k = 1, nx
                        if (m(k,l)) b2(i,j) = min(b2(i,j), ((i - k)*ddx)**2 + ((j - l)*ddy)**2)
                    end do
                end do
            end do
        end do
        call check("exact against brute force (dx /= dy)", all(abs(d2 - b2) <= 1.0e-9_dp*max(1.0_dp, b2)))
        call check("zero on the cells themselves", all(d2 == 0.0_dp .eqv. m))
        d2 = margin_edt2(.not. (m .or. .not. m), ddx, ddy)
        call check("no cell: huge everywhere", all(d2 == huge(1.0_dp)))
    end subroutine test_edt

    subroutine test_sdist_slice()
        logical  :: m(23, 17)
        real(dp) :: sd(23, 17)

        write(*,"(a)") "signed distance of a slice"
        m = irregular(23, 17)
        sd = margin_signed_distance(m, 2.0_dp, 3.0_dp)
        call check("sd > 0 is the slice", all((sd > 0.0_dp) .eqv. m))
        call check("|sd| >= min(dx, dy)/2", all(abs(sd) >= 1.0_dp - 1.0e-12_dp))
        call check("an ice cell next to the margin along x has sd = h", any(m .and. abs(sd - 1.0_dp) < 1.0e-12_dp))
        sd = margin_signed_distance(spread(spread(.false., 1, 23), 2, 17), 2.0_dp, 3.0_dp)
        call check("no ice: minus the diagonal (less h) everywhere", &
                   all(abs(sd + sqrt(46.0_dp**2 + 51.0_dp**2) - 1.0_dp) < 1.0e-9_dp))
    end subroutine test_sdist_slice

    subroutine test_disks()
        type(margin_t) :: mg
        logical  :: ice(N, N, 2), mask(N, N)
        real(dp) :: sd(N, N), r(N, N), fp(N, N), Rt, t
        integer  :: k
        logical  :: ok

        write(*,"(a)") "two disks, sdist and nearest"
        r = radius()
        ice(:,:,1) = disk(R2)                        ! given in decreasing time order: sorted by margin_set
        ice(:,:,2) = disk(R1)
        call margin_set(mg, [200.0_dp, 100.0_dp], ice, DX, DX, MARGIN_SDIST, .false., "disks")
        call check("slices sorted in time", all(mg%time == [100.0_dp, 200.0_dp]))

        call margin_at(mg, 100.0_dp, mask, sd)
        call check("t1: slice 1 exactly", all(mask .eqv. disk(R1)))
        call margin_at(mg, 200.0_dp, mask, sd)
        call check("t2: slice 2 exactly", all(mask .eqv. disk(R2)))

        ok = .true.
        do k = 1, 9
            t = 100.0_dp + 10.0_dp*k
            Rt = R1 + (t - 100.0_dp)/100.0_dp*(R2 - R1)
            call margin_at(mg, t, mask, sd)
            ok = ok .and. all(mask .or. r > Rt - 1.0_dp) .and. all(.not. mask .or. r < Rt + 1.0_dp)
        end do
        call check("the margin moves linearly (within one cell of R1 + a (R2 - R1))", ok)
        call margin_at(mg, 150.0_dp, mask, sd)
        call check("t = 150: the area of a disk of radius 15 within 5 %", &
                   abs(real(count(mask), dp) - acos(-1.0_dp)*15.0_dp**2) < 0.05_dp*acos(-1.0_dp)*15.0_dp**2)
        fp = margin_footprint(mg)
        call check("footprint: the union of the slices", all((fp > 0.5_dp) .eqv. disk(R2)))

        call margin_set(mg, [100.0_dp, 200.0_dp], ice(:,:,[2,1]), DX, DX, MARGIN_NEAREST, .false., "disks")
        call margin_at(mg, 149.0_dp, mask, sd)
        call check("nearest: before the midpoint, slice 1", all(mask .eqv. disk(R1)))
        call margin_at(mg, 150.0_dp, mask, sd)
        call check("nearest: a tie keeps the earlier slice", all(mask .eqv. disk(R1)))
        call margin_at(mg, 151.0_dp, mask, sd)
        call check("nearest: after it, slice 2", all(mask .eqv. disk(R2)))
    end subroutine test_disks

    subroutine test_file()
        type(margin_t) :: mg
        character(len=*), parameter :: f = dir//"/margin.nc"
        real(dp) :: H(N, N, 2), sd(N, N), x(N)
        logical  :: mask(N, N)
        integer  :: i

        write(*,"(a)") "margin file, years ago"
        x = [(real(i - 1, dp), i = 1, N)]
        H(:,:,1) = merge(1500.0_dp, 0.0_dp, disk(R2))         ! 20 ka: large
        H(:,:,2) = merge(800.0_dp, -9999.0_dp, disk(R1))      ! 0 ka: small, fill value outside
        call nc_create(f)
        call nc_write_dim(f, "xc", x=x, units="km")
        call nc_write_dim(f, "yc", x=x, units="km")
        call nc_write_dim(f, "time", x=[20.0_dp, 0.0_dp], units="kyr BP")
        call nc_write(f, "H_ice", H, dim1="xc", dim2="yc", dim3="time")
        call margin_init(mg, f, "H_ice", 0.0_dp, N, N, DX, DX, MARGIN_SDIST, .false.)
        call check("model time of 'kyr BP', increasing", all(mg%time == [-20000.0_dp, 0.0_dp]))
        call margin_at(mg, -20000.0_dp, mask, sd)
        call check("the 20 ka slice", all(mask .eqv. disk(R2)))
        call margin_at(mg, 0.0_dp, mask, sd)
        call check("the 0 ka slice (fill value: no ice)", all(mask .eqv. disk(R1)))
        call margin_init(mg, f, "H_ice", 1000.0_dp, N, N, DX, DX, MARGIN_SDIST, .false.)
        call margin_at(mg, 0.0_dp, mask, sd)
        call check("margin_H_min: thinner ice is no ice", .not. any(mask))
        call check("model_time_axis: years ago, ka, years", &
                   all(model_time_axis([5.0_dp], "years ago") == [-5.0_dp]) .and. &
                   all(model_time_axis([5.0_dp], "ka") == [5000.0_dp]) .and. &
                   all(model_time_axis([5.0_dp], "years") == [5.0_dp]))
    end subroutine test_file

    subroutine test_out_of_range()
        type(margin_t) :: mg
        logical  :: ice(N, N, 2), mask(N, N)
        real(dp) :: sd(N, N)

        write(*,"(a)") "outside the time range"
        ice(:,:,1) = disk(R1)
        ice(:,:,2) = disk(R2)
        call margin_set(mg, [0.0_dp, 10.0_dp], ice, DX, DX, MARGIN_SDIST, .true., "disks")
        call margin_at(mg, -5.0_dp, mask, sd)
        call check("hold: the first slice before the range", all(mask .eqv. disk(R1)))
        call margin_at(mg, 50.0_dp, mask, sd)
        call check("hold: the last slice after it", all(mask .eqv. disk(R2)))
        call expect_abort("margin_range", "outside the time range of the margin slices")
        call expect_abort("taub_nodomain", "has no tau_b")
        call expect_abort("taub_badcode", "is not a domain of the map")
    end subroutine test_out_of_range

    subroutine test_taub()
        type(taub_t) :: tb
        character(len=*), parameter :: f1 = dir//"/taub1.nc", f2 = dir//"/taub2.nc"
        integer  :: dom(4, 3)
        real(dp) :: tau(4, 3), ser(2, 3), x(4)
        integer  :: i

        write(*,"(a)") "tau_b"
        dom = reshape([0, 1, 1, 2,  0, 2, 2, 5,  1, 1, 0, 0], [4, 3])
        call taub_set_domains(tb, dom, [1, 2, 5], [1.0e4_dp, 2.0e4_dp, 5.0e4_dp], 9.0e3_dp, "test")
        call taub_at(tb, 0.0_dp, tau)
        call check("static: each cell its domain's value, the default elsewhere", &
                   all(tau == merge(9.0e3_dp, merge(1.0e4_dp, merge(2.0e4_dp, 5.0e4_dp, dom == 2), dom == 1), dom == 0)))
        call taub_set_domains(tb, dom, [1, 2, 5], [1.0e4_dp, 2.0e4_dp, 5.0e4_dp], 9.0e3_dp, "test", &
                              codes=[2], values=[2.5e4_dp], table="&tbl")
        call taub_at(tb, 0.0_dp, tau)
        call check("override of id 2", all(pack(tau, dom == 2) == 2.5e4_dp) .and. all(pack(tau, dom == 1) == 1.0e4_dp))

        ! a series in the file (records at 0, 1000, 2000 years ago: id 1 = 50, 30, 10; id 2 = 60, 40, 20),
        ! (id, time) in Fortran order and transposed; id 5 of the map has no value in the file: the table gives it
        ser = reshape([10.0_dp, 20.0_dp,  30.0_dp, 40.0_dp,  50.0_dp, 60.0_dp], [2, 3])   ! (id, time), oldest first
        x = [(real(i, dp), i = 1, 4)]
        call write_domains(f1, .false., dom, ser, x)
        call write_domains(f2, .true., dom, ser, x)
        do i = 1, 2
            call taub_read_domains(tb, merge(f1, f2, i == 1), "taub_domain", "domain_id", "taub", "taub_t", &
                                   0.0_dp, .false., 4, 3, codes=[5], values=[1.0e5_dp], table="&tbl")
            call check("the series is read (file "//merge("1", "2", i == 1)//")", allocated(tb%val_t) .and. &
                       tb%method == TAUB_M_DOMAINS .and. all(tb%time == [-2000.0_dp, -1000.0_dp, 0.0_dp]))
            call taub_at(tb, -1500.0_dp, tau)
            call check("linear in time between the records", all(pack(tau, dom == 1) == 20.0_dp) .and. &
                       all(pack(tau, dom == 2) == 30.0_dp))
            call taub_at(tb, -500.0_dp, tau)
            call check("the table adds an id of the map, constant in time", all(pack(tau, dom == 5) == 1.0e5_dp) .and. &
                       all(pack(tau, dom == 1) == 40.0_dp) .and. all(pack(tau, dom == 0) == 0.0_dp))
        end do

        call taub_set_const(tb, 5.0e4_dp)
        call taub_at(tb, 123.0_dp, tau)
        call check("const", all(tau == 5.0e4_dp))

    end subroutine test_taub

    subroutine write_domains(f, transposed, dom, ser, x)
        ! A tau_b domain file: map, ids 1 and 2, static values and the series (records at 0, 1000, 2000
        ! years ago = ser(:, 3), ser(:, 2), ser(:, 1)), (id, time) in Fortran order or transposed.
        character(len=*), intent(in) :: f
        logical,          intent(in) :: transposed
        integer,          intent(in) :: dom(:,:)
        real(dp),         intent(in) :: ser(:,:), x(:)
        call nc_create(f)
        call nc_write_dim(f, "xc", x=x(1:size(dom, 1)), units="km")
        call nc_write_dim(f, "yc", x=x(1:size(dom, 2)), units="km")
        call nc_write_dim(f, "domain", x=[1.0_dp, 2.0_dp])
        call nc_write_dim(f, "time", x=[0.0_dp, 1000.0_dp, 2000.0_dp], units="years ago")
        call nc_write(f, "taub_domain", dom, dim1="xc", dim2="yc")
        call nc_write(f, "domain_id", [1, 2], dim1="domain")
        call nc_write(f, "taub", [1.0e4_dp, 2.0e4_dp], dim1="domain")
        if (transposed) then
            call nc_write(f, "taub_t", transpose(ser(:, [3, 2, 1])), dim1="time", dim2="domain")
        else
            call nc_write(f, "taub_t", ser(:, [3, 2, 1]), dim1="domain", dim2="time")
        end if
    end subroutine write_domains

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
        type(margin_t) :: mg
        type(taub_t)   :: tb
        logical  :: ice(N, N, 2), mask(N, N)
        real(dp) :: sd(N, N)
        integer  :: dom(2, 2)

        call get_command_argument(1, arg)
        select case (trim(arg))
        case ("margin_range")
            ice(:,:,1) = disk(R1)
            ice(:,:,2) = disk(R2)
            call margin_set(mg, [0.0_dp, 10.0_dp], ice, DX, DX, MARGIN_SDIST, .false., "disks")
            call margin_at(mg, 10.5_dp, mask, sd)
        case ("taub_nodomain")
            dom = reshape([1, 2, 0, 1], [2, 2])
            call taub_set_domains(tb, dom, [1], [1.0e4_dp], 0.0_dp, "test")
        case ("taub_badcode")
            dom = reshape([1, 2, 0, 1], [2, 2])
            call taub_set_domains(tb, dom, [1, 2], [1.0e4_dp, 2.0e4_dp], 0.0_dp, "test", codes=[3], values=[1.0e4_dp], &
                                  table="&tbl")
        end select
    end subroutine abort_case

end program test_plastic_margin
