program test_remap_sp_dp
    ! C4 / core-design R1, adapters-design P4: the precision rule of the exchange.
    !
    ! For real(sp) fields, widening to dp, remapping with the fesm-utils coupler
    ! in dp and narrowing the result gives bit for bit the values of the
    ! coupler's sp remap (map_field_grid_grid_sp: widen, map in dp, narrow).
    !
    ! Every case compares  real(remap_dp(real(x,dp)), sp)  with  remap_sp(x)
    ! as 32-bit patterns (no tolerance; -0 and NaN patterns count), for
    !   grid pairs  GRL-16KM -> GRL-16KM (identity), ANT-16KM -> ANT-32KM,
    !               ANT-32KM -> ANT-16KM
    !   methods     con, bilin, auto
    !   fields      2D and 3D (12 months); topography and monthly climate from
    !               ice_data where present, plus synthetic fields whose values
    !               are not round numbers
    !   variants    1 plain (default missing value), 2 explicit missing_value
    !               with mask2 (mask2 must agree too), 3 explicit missing_value
    !               with mask_pack
    !
    ! Grids come from maps/grid_<name>.txt, weight caches maps/*.nc are
    ! generated on the fly (gitignored). ice_data is read from $KRYOS_ICE_DATA,
    ! default ./ice_data; a missing file is reported and skipped (the
    ! synthetic fields are always run).

    use, intrinsic :: iso_fortran_env, only: int32, int64
    use precision, only: sp, dp
    use coupler,   only: coupler_class, coupler_init, coupler_method, remap
    use ncio,      only: nc_read

    implicit none

    integer, parameter :: NMON = 12
    integer, parameter :: NVAR = 3
    real(sp), parameter :: mv_sp     = -9999.0_sp     ! the coupler's default missing value
    real(sp), parameter :: mv_custom = -1.0e30_sp     ! explicit missing value (variants 2, 3)
    character(len=*), parameter :: map_fldr = "maps"

    type field_t
        character(len=64) :: name = ""
        logical  :: real_data = .false.
        real(sp), allocatable :: v2(:,:)
        real(sp), allocatable :: v3(:,:,:)
    end type

    type grid_data_t
        character(len=16) :: name = ""
        integer :: nx = 0, ny = 0
        type(field_t), allocatable :: f2(:), f3(:)
    end type

    type stat_t
        integer(int64) :: n = 0, nbad = 0, nround = 0, nmask = 0, nmaskbad = 0, maxulp = 0
        integer :: nfld = 0
    end type

    character(len=16), parameter :: pairs(2,3) = reshape( &
        [character(len=16) :: "GRL-16KM", "GRL-16KM", "ANT-16KM", "ANT-32KM", "ANT-32KM", "ANT-16KM"], [2,3])
    character(len=8), parameter :: methods(3) = [character(len=8) :: "con", "bilin", "auto"]
    character(len=16), parameter :: vname(NVAR) = [character(len=16) :: "plain", "missing+mask2", "missing+pack"]

    type(coupler_class) :: cpl
    type(grid_data_t)   :: gd(3)
    type(stat_t)        :: st(NVAR, 2)          ! variant x (2D, 3D)
    type(stat_t)        :: tot
    character(len=256)  :: root
    character(len=32)   :: resolved
    integer :: ip, im, ig, iv, id, nreal, nfailed_files
    integer(int64) :: t0, t1, rate

    call system_clock(t0, rate)
    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = "ice_data"

    call coupler_init(cpl, map_fldr=map_fldr)

    ! --- input fields per grid --------------------------------------------------
    call build_fields(gd(1), "GRL-16KM", 106, 181, root, &
        "Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc", &
        "Greenland/GRL-16KM/GRL-16KM_MARv3.11-ERA_monmean_1961-1990.nc", ["T_srf", "smb  "])
    call build_fields(gd(2), "ANT-16KM", 381, 381, root, &
        "Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine.nc", &
        "Antarctica/ANT-16KM/ANT-16KM_RACMO23-ERAINT-HYBRID_1981-2010.nc", ["T_srf", "smb  "])
    call build_fields(gd(3), "ANT-32KM", 191, 191, root, &
        "Antarctica/ANT-32KM/ANT-32KM_TOPO-BedMachine.nc", &
        "Antarctica/ANT-32KM/ANT-32KM_RACMO23-ERA-INTERIM_monthly_1981-2010.nc", ["pr ", "smb"])

    nreal = 0
    do ig = 1, 3
        nreal = nreal + count([(gd(ig)%f2(id)%real_data, id = 1, size(gd(ig)%f2))]) &
                      + count([(gd(ig)%f3(id)%real_data, id = 1, size(gd(ig)%f3))])
    end do
    write(*,'(a,i0,a)') "test_remap_sp_dp: ", nreal, " fields read from ice_data (rest synthetic)"

    ! --- cases ------------------------------------------------------------------
    write(*,'(a)') ""
    write(*,'(a)') "  pair                  method(resolved)  variant        dim  nfld  values      mismatch  maxulp  rounded  mask2-mismatch"
    tot = stat_t()
    do ip = 1, size(pairs, 2)
        ig = grid_index(trim(pairs(1,ip)))
        do im = 1, size(methods)
            resolved = coupler_method(cpl, trim(pairs(1,ip)), trim(pairs(2,ip)), method=trim(methods(im)))
            st = stat_t()
            do id = 1, size(gd(ig)%f2)
                call run_field_2d(gd(ig)%f2(id), trim(pairs(1,ip)), trim(pairs(2,ip)), trim(methods(im)), st(:,1))
            end do
            do id = 1, size(gd(ig)%f3)
                call run_field_3d(gd(ig)%f3(id), trim(pairs(1,ip)), trim(pairs(2,ip)), trim(methods(im)), st(:,2))
            end do
            do iv = 1, NVAR
                call report(trim(pairs(1,ip))//" -> "//trim(pairs(2,ip)), &
                            trim(methods(im))//"("//trim(resolved)//")", vname(iv), "2D", st(iv,1))
                call report(trim(pairs(1,ip))//" -> "//trim(pairs(2,ip)), &
                            trim(methods(im))//"("//trim(resolved)//")", vname(iv), "3D", st(iv,2))
                call accumulate(tot, st(iv,1))
                call accumulate(tot, st(iv,2))
            end do
        end do
    end do

    call system_clock(t1)
    write(*,'(a)') ""
    write(*,'(a,i0,a,i0,a,i0,a,f7.1,a)') "test_remap_sp_dp: ", tot%n, " values compared, ", tot%nbad, &
        " not bit-identical, ", tot%nmaskbad, " mask2 mismatches; runtime ", real(t1 - t0)/real(rate), " s"
    if (tot%nbad /= 0 .or. tot%nmaskbad /= 0) then
        write(*,'(a)') "test_remap_sp_dp: FAILED"
        error stop 1
    end if
    write(*,'(a)') "test_remap_sp_dp: OK"

contains

    ! ------------------------------------------------------------------------------
    ! cases
    ! ------------------------------------------------------------------------------

    function grid_index(name) result(ig)
        character(len=*), intent(in) :: name
        integer :: ig
        do ig = 1, size(gd)
            if (trim(gd(ig)%name) == name) return
        end do
        write(*,'(a)') "test_remap_sp_dp: no fields for grid "//name
        error stop 1
    end function grid_index

    subroutine run_field_2d(f, src, dst, method, s)
        type(field_t),    intent(in)    :: f
        character(len=*), intent(in)    :: src, dst, method
        type(stat_t),     intent(inout) :: s(NVAR)

        real(sp), allocatable :: a(:,:)
        real(dp), allocatable :: wd(:,:), d(:,:)
        logical,  allocatable :: ma(:,:), mb(:,:), pack(:,:)
        real(sp), allocatable :: x(:,:)
        integer :: i, j

        ! variant 1: default missing value (-9999)
        call remap(cpl, f%v2, src, a, dst, method=method)
        wd = real(f%v2, dp)
        call remap(cpl, wd, src, d, dst, method=method)
        call compare(a, d, s(1))

        ! variant 2: explicit missing value, mask2
        x = f%v2
        where (x == mv_sp) x = mv_custom
        wd = real(x, dp)
        allocate(ma(size(a,1), size(a,2)), mb(size(a,1), size(a,2)))
        call remap(cpl, x, src, a, dst, method=method, missing_value=mv_custom, mask2=ma)
        call remap(cpl, wd, src, d, dst, method=method, missing_value=real(mv_custom, dp), mask2=mb)
        call compare(a, d, s(2))
        call compare_mask(ma, mb, s(2))

        ! variant 3: explicit missing value, mask_pack on the target grid
        allocate(pack(size(a,1), size(a,2)))
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                pack(i,j) = mod(i + 2*j, 3) /= 0
            end do
        end do
        call remap(cpl, x, src, a, dst, method=method, missing_value=mv_custom, mask_pack=pack)
        call remap(cpl, wd, src, d, dst, method=method, missing_value=real(mv_custom, dp), mask_pack=pack)
        call compare(a, d, s(3))
        s%nfld = s%nfld + 1
    end subroutine run_field_2d

    subroutine run_field_3d(f, src, dst, method, s)
        type(field_t),    intent(in)    :: f
        character(len=*), intent(in)    :: src, dst, method
        type(stat_t),     intent(inout) :: s(NVAR)

        real(sp), allocatable :: a(:,:,:), x(:,:,:)
        real(dp), allocatable :: wd(:,:,:), d(:,:,:)
        logical,  allocatable :: ma(:,:,:), mb(:,:,:), pack(:,:)
        integer :: i, j

        call remap(cpl, f%v3, src, a, dst, method=method)
        wd = real(f%v3, dp)
        call remap(cpl, wd, src, d, dst, method=method)
        call compare(a, d, s(1))

        x = f%v3
        where (x == mv_sp) x = mv_custom
        wd = real(x, dp)
        allocate(ma(size(a,1), size(a,2), size(a,3)), mb(size(a,1), size(a,2), size(a,3)))
        call remap(cpl, x, src, a, dst, method=method, missing_value=mv_custom, mask2=ma)
        call remap(cpl, wd, src, d, dst, method=method, missing_value=real(mv_custom, dp), mask2=mb)
        call compare(a, d, s(2))
        call compare_mask(ma, mb, s(2))

        allocate(pack(size(a,1), size(a,2)))
        do j = 1, size(a,2)
            do i = 1, size(a,1)
                pack(i,j) = mod(i + 2*j, 3) /= 0
            end do
        end do
        call remap(cpl, x, src, a, dst, method=method, missing_value=mv_custom, mask_pack=pack)
        call remap(cpl, wd, src, d, dst, method=method, missing_value=real(mv_custom, dp), mask_pack=pack)
        call compare(a, d, s(3))
        s%nfld = s%nfld + 1
    end subroutine run_field_3d

    ! ------------------------------------------------------------------------------
    ! comparison: 32-bit patterns, ulp distance, how many dp results needed rounding
    ! ------------------------------------------------------------------------------

    subroutine compare(a, d, s)
        ! a: the sp remap; d: the dp remap of the widened input (narrowed here).
        real(sp), intent(in) :: a(..)
        real(dp), intent(in) :: d(..)
        type(stat_t), intent(inout) :: s
        select rank (a)
        rank (2)
            select rank (d)
            rank (2)
                if (any(shape(a) /= shape(d))) call shape_fail()
                call compare_flat(reshape(a, [size(a)]), reshape(d, [size(d)]), s)
            end select
        rank (3)
            select rank (d)
            rank (3)
                if (any(shape(a) /= shape(d))) call shape_fail()
                call compare_flat(reshape(a, [size(a)]), reshape(d, [size(d)]), s)
            end select
        end select
    end subroutine compare

    subroutine shape_fail()
        write(*,'(a)') "test_remap_sp_dp: FAILED (sp and dp results differ in shape)"
        error stop 1
    end subroutine shape_fail

    subroutine compare_flat(a, d, s)
        real(sp), intent(in) :: a(:)
        real(dp), intent(in) :: d(:)
        type(stat_t), intent(inout) :: s
        real(sp) :: b
        integer :: k
        integer(int64) :: u
        do k = 1, size(a)
            b = real(d(k), sp)
            s%n = s%n + 1
            if (real(b, dp) /= d(k)) s%nround = s%nround + 1
            if (transfer(a(k), 0_int32) /= transfer(b, 0_int32)) then
                s%nbad = s%nbad + 1
                u = ulp_dist(a(k), b)
                if (u > s%maxulp) s%maxulp = u
            end if
        end do
    end subroutine compare_flat

    elemental function ulp_dist(a, b) result(u)
        real(sp), intent(in) :: a, b
        integer(int64) :: u
        u = abs(ordered(a) - ordered(b))
    end function ulp_dist

    elemental function ordered(a) result(o)
        ! monotonic integer image of an sp value
        real(sp), intent(in) :: a
        integer(int64) :: o
        o = int(transfer(a, 0_int32), int64)
        if (o < 0) o = -2147483648_int64 - o
    end function ordered

    subroutine compare_mask(ma, mb, s)
        logical, intent(in) :: ma(..), mb(..)
        type(stat_t), intent(inout) :: s
        select rank (ma)
        rank (2)
            select rank (mb)
            rank (2)
                s%nmask    = s%nmask + size(ma)
                s%nmaskbad = s%nmaskbad + count(ma .neqv. mb)
            end select
        rank (3)
            select rank (mb)
            rank (3)
                s%nmask    = s%nmask + size(ma)
                s%nmaskbad = s%nmaskbad + count(ma .neqv. mb)
            end select
        end select
    end subroutine compare_mask

    subroutine accumulate(t, s)
        type(stat_t), intent(inout) :: t
        type(stat_t), intent(in)    :: s
        t%n = t%n + s%n;  t%nbad = t%nbad + s%nbad;  t%nmaskbad = t%nmaskbad + s%nmaskbad
        t%nround = t%nround + s%nround
        t%maxulp = max(t%maxulp, s%maxulp)
    end subroutine accumulate

    subroutine report(pair, method, variant, dim, s)
        character(len=*), intent(in) :: pair, method, variant, dim
        type(stat_t),     intent(in) :: s
        write(*,'(2x,a22,1x,a18,1x,a14,1x,a3,2x,i3,2x,i10,2x,i8,2x,i6,2x,i7,2x,i8)') &
            pair, method, variant, dim, s%nfld, s%n, s%nbad, s%maxulp, s%nround, s%nmaskbad
    end subroutine report

    ! ------------------------------------------------------------------------------
    ! input fields
    ! ------------------------------------------------------------------------------

    subroutine build_fields(g, name, nx, ny, root, topo, clim, cvars)
        type(grid_data_t), intent(out) :: g
        character(len=*),  intent(in)  :: name, root, topo, clim
        integer,           intent(in)  :: nx, ny
        character(len=*),  intent(in)  :: cvars(:)

        character(len=*), parameter :: tvars(3) = [character(len=5) :: "z_bed", "H_ice", "z_srf"]
        type(field_t), allocatable :: l2(:), l3(:)
        real(sp), allocatable :: t2(:,:), t3(:,:,:)
        integer :: k, n2, n3
        logical :: ok

        g%name = name
        g%nx = nx
        g%ny = ny
        allocate(l2(3 + 3), l3(1 + size(cvars)))
        n2 = 0
        n3 = 0

        do k = 1, 3
            call add2(l2, n2, name, synthetic_2d(nx, ny, k), "synthetic_"//trim(synth_name(k)), .false.)
        end do
        do k = 1, size(tvars)
            call read_2d(trim(root)//"/"//topo, trim(tvars(k)), nx, ny, t2, ok)
            if (ok) call add2(l2, n2, name, t2, trim(tvars(k)), .true.)
        end do

        call add3(l3, n3, name, synthetic_3d(nx, ny), "synthetic_monthly", .false.)
        do k = 1, size(cvars)
            call read_3d(trim(root)//"/"//clim, trim(cvars(k)), nx, ny, t3, ok)
            if (ok) call add3(l3, n3, name, t3, trim(cvars(k)), .true.)
        end do

        allocate(g%f2(n2), g%f3(n3))
        do k = 1, n2
            call move_field(l2(k), g%f2(k))
        end do
        do k = 1, n3
            call move_field(l3(k), g%f3(k))
        end do

    end subroutine build_fields

    subroutine add2(l, n, gname, v, fname, is_real)
        type(field_t),    intent(inout) :: l(:)
        integer,          intent(inout) :: n
        character(len=*), intent(in)    :: gname, fname
        real(sp),         intent(in)    :: v(:,:)
        logical,          intent(in)    :: is_real
        n = n + 1
        l(n)%v2 = v
        l(n)%name = trim(gname)//" "//trim(fname)
        l(n)%real_data = is_real
    end subroutine add2

    subroutine add3(l, n, gname, v, fname, is_real)
        type(field_t),    intent(inout) :: l(:)
        integer,          intent(inout) :: n
        character(len=*), intent(in)    :: gname, fname
        real(sp),         intent(in)    :: v(:,:,:)
        logical,          intent(in)    :: is_real
        n = n + 1
        l(n)%v3 = v
        l(n)%name = trim(gname)//" "//trim(fname)
        l(n)%real_data = is_real
    end subroutine add3

    subroutine move_field(a, b)
        type(field_t), intent(inout) :: a, b
        b%name = a%name
        b%real_data = a%real_data
        if (allocated(a%v2)) call move_alloc(a%v2, b%v2)
        if (allocated(a%v3)) call move_alloc(a%v3, b%v3)
    end subroutine move_field

    function synth_name(k) result(s)
        integer, intent(in) :: k
        character(len=8) :: s
        select case (k)
        case (1); s = "smooth"
        case (2); s = "wide"
        case default; s = "noisy"
        end select
    end function synth_name

    function synthetic_2d(nx, ny, kind) result(v)
        ! Values that are not round numbers (full sp mantissas, mixed signs and
        ! magnitudes), with holes set to the default missing value.
        integer, intent(in) :: nx, ny, kind
        real(sp), allocatable :: v(:,:)
        integer :: i, j
        real(dp) :: x, r
        integer(int64) :: h

        allocate(v(nx, ny))
        do j = 1, ny
            do i = 1, nx
                select case (kind)
                case (1)    ! smooth, temperature-like
                    x = 251.3719_dp + 23.7_dp*sin(0.0731_dp*i)*cos(0.0417_dp*j) + 0.0123456_dp*i - 0.00987_dp*j
                case (2)    ! wide dynamic range, both signs
                    x = (-1.0_dp)**mod(i + j, 2) * 10.0_dp**(4.0_dp*sin(0.0571_dp*i + 0.0293_dp*j))/3.0_dp
                case default ! pseudo-random, uncorrelated cell to cell
                    h = iand(int(i, int64)*73856093_int64 + int(j, int64)*19349663_int64 + 12345_int64, 2147483647_int64)
                    h = iand(h*1103515245_int64 + 12345_int64, 2147483647_int64)
                    r = real(h, dp)/2147483647.0_dp
                    x = (r - 0.5_dp)*1234.5678_dp + 0.1_dp*i
                end select
                v(i,j) = real(x, sp)
                if (mod(7*i + 13*j, 23) == 0) v(i,j) = mv_sp
                if (i > nx/2 .and. i <= nx/2 + max(1, nx/8) .and. j > ny/3 .and. j <= ny/3 + max(1, ny/8)) &
                    v(i,j) = mv_sp
            end do
        end do
    end function synthetic_2d

    function synthetic_3d(nx, ny) result(v)
        integer, intent(in) :: nx, ny
        real(sp), allocatable :: v(:,:,:)
        real(sp), allocatable :: b(:,:)
        integer :: m, i, j

        allocate(v(nx, ny, NMON))
        b = synthetic_2d(nx, ny, 1)
        do m = 1, NMON
            do j = 1, ny
                do i = 1, nx
                    if (b(i,j) == mv_sp) then
                        v(i,j,m) = mv_sp
                    else
                        v(i,j,m) = real(real(b(i,j), dp)*(1.0_dp + 0.0713_dp*sin(0.5236_dp*m + 0.01_dp*i)) &
                                        + 0.3141592_dp*m, sp)
                    end if
                end do
            end do
        end do
    end function synthetic_3d

    subroutine read_2d(file, var, nx, ny, v, ok)
        character(len=*), intent(in) :: file, var
        integer, intent(in) :: nx, ny
        real(sp), allocatable, intent(inout) :: v(:,:)
        logical, intent(out) :: ok
        inquire(file=file, exist=ok)
        if (.not. ok) then
            write(*,'(a)') "test_remap_sp_dp: missing "//trim(file)//" (skipped)"
            return
        end if
        if (allocated(v)) deallocate(v)
        allocate(v(nx, ny))
        call nc_read(file, var, v)
        call check_finite_2d(v, file, var)
    end subroutine read_2d

    subroutine read_3d(file, var, nx, ny, v, ok)
        character(len=*), intent(in) :: file, var
        integer, intent(in) :: nx, ny
        real(sp), allocatable, intent(inout) :: v(:,:,:)
        logical, intent(out) :: ok
        inquire(file=file, exist=ok)
        if (.not. ok) then
            write(*,'(a)') "test_remap_sp_dp: missing "//trim(file)//" (skipped)"
            return
        end if
        if (allocated(v)) deallocate(v)
        allocate(v(nx, ny, NMON))
        call nc_read(file, var, v)
        call check_finite_2d(reshape(v, [nx*NMON, ny]), file, var)
    end subroutine read_3d

    subroutine check_finite_2d(v, file, var)
        real(sp), intent(in) :: v(:,:)
        character(len=*), intent(in) :: file, var
        if (any(v /= v) .or. any(abs(v) > huge(1.0_sp))) then
            write(*,'(a)') "test_remap_sp_dp: non-finite values in "//trim(var)//" of "//trim(file)
            error stop 1
        end if
    end subroutine check_finite_2d

end program test_remap_sp_dp
