module kryos_plastic_margin
    ! The margin chronology of the ice backend `plastic` (docs/dev/plastic-backend.md): the slices
    ! of an ice-thickness history on the instance grid, each turned once into the signed distance to
    ! its margin, and the ice mask at a model time. No library: compiled whatever `plastic=`, so the
    ! tests run with the stub.
    !
    ! Signed distance of a slice (margin_signed_distance): positive inside the ice, in m, exact Euclidean
    ! distance between cell centres (margin_edt2, the separable transform of Felzenszwalb and
    ! Huttenlocher 2012), shifted by half a cell so that the zero level lies midway between an ice
    ! cell and its ice-free neighbour: an ice cell is >= +h, an ice-free cell <= -h (h = min(dx,dy)/2),
    ! and `sd > 0` gives the slice back exactly. A slice without ice (or without an ice-free cell)
    ! has the distance of the domain diagonal, so that it interpolates like any other.
    !
    ! Mask at time t (margin_at): "sdist" interpolates the signed distance linearly between the
    ! slices that bracket t and takes sd > 0 (the margin moves smoothly; a cell is ice only where it
    ! is ice in one of the two slices); "nearest" takes the nearest slice (Gowan et al. 2021 offline).
    ! Outside the time range of the slices: stop, or hold the first or last slice.

    use precision, only: dp
    use fesm_log,  only: fesm_error, to_str
    use ncio,      only: nc_read, nc_size, nc_read_attr, nc_dims

    implicit none

    private

    integer, parameter, public :: MARGIN_SDIST = 1, MARGIN_NEAREST = 2

    real(dp), parameter :: BIG = huge(1.0_dp)

    type, public :: margin_t
        integer  :: n = 0, nx = 0, ny = 0
        integer  :: interp = MARGIN_SDIST
        logical  :: hold = .false.                   ! outside the slices' time range: hold the first/last (else stop)
        real(dp) :: dx = 0.0_dp, dy = 0.0_dp         ! [m]
        real(dp), allocatable :: time(:)             ! [yr] model time of each slice, increasing
        real(dp), allocatable :: sd(:,:,:)           ! [m] signed distance of each slice (positive inside the ice)
        character(len=512) :: src = ""               ! where the slices came from (messages)
    end type margin_t

    public :: margin_init, margin_set, margin_at, margin_footprint, margin_interp_code
    public :: margin_signed_distance, margin_edt2, model_time_axis, time_order

contains

    ! =========================================================================
    ! Set-up
    ! =========================================================================

    integer function margin_interp_code(name) result(code)
        ! "sdist" | "nearest" -> MARGIN_SDIST | MARGIN_NEAREST; 0 for anything else.
        character(len=*), intent(in) :: name
        select case (trim(name))
        case ("sdist");   code = MARGIN_SDIST
        case ("nearest"); code = MARGIN_NEAREST
        case default;     code = 0
        end select
    end function margin_interp_code

    subroutine margin_init(mg, file, var, H_min, nx, ny, dx, dy, interp, hold)
        ! The slices of `var` (an ice fraction or thickness on (xc, yc, time), the instance grid) in `file`:
        ! ice where var > H_min (a fill value counts as no ice). The time axis is the variable `time`,
        ! converted to model time from its units (model_time_axis).
        type(margin_t),   intent(out) :: mg
        character(len=*), intent(in)  :: file, var
        real(dp),         intent(in)  :: H_min, dx, dy
        integer,          intent(in)  :: nx, ny, interp
        logical,          intent(in)  :: hold

        character(len=*), parameter :: proc = "kryos_plastic_margin::margin_init"
        character(len=64), allocatable :: names(:)
        integer, allocatable :: dims(:)
        real(dp), allocatable :: t(:), H(:,:)
        logical,  allocatable :: ice(:,:,:)
        character(len=256) :: units
        logical :: exists
        integer :: nt, k

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the margin file does not exist", "margin_file = "//trim(file))
        if (H_min < 0.0_dp) call fesm_error(proc, "margin_H_min must be >= 0", "margin_H_min = "//to_str(H_min))

        call nc_dims(file, var, names, dims)
        if (size(dims) /= 3) then
            call fesm_error(proc, "`"//trim(var)//"` is not a field (x, y, time)", "file "//trim(file)// &
                            ", "//to_str(size(dims))//" dimensions")
        end if
        if (dims(1) /= nx .or. dims(2) /= ny) then
            call fesm_error(proc, "`"//trim(var)//"` is not on the instance grid", "file "//trim(file)//": "// &
                            to_str(dims(1))//" x "//to_str(dims(2))//", the grid "//to_str(nx)//" x "//to_str(ny))
        end if
        nt = nc_size(file, "time")
        if (nt /= dims(3)) then
            call fesm_error(proc, "the third dimension of `"//trim(var)//"` is not the time axis", "file "//trim(file))
        end if

        allocate(t(nt), H(nx, ny), ice(nx, ny, nt))
        call nc_read(file, "time", t)
        units = ""
        call nc_read_attr(file, "time", "units", units)
        t = model_time_axis(t, units)
        do k = 1, nt
            call nc_read(file, var, H, start=[1,1,k], count=[nx,ny,1])
            ice(:,:,k) = H > H_min
        end do

        call margin_set(mg, t, ice, dx, dy, interp, hold, trim(file)//" ("//trim(var)//")")
    end subroutine margin_init

    subroutine margin_set(mg, time, ice, dx, dy, interp, hold, src)
        ! The slices `ice(:,:,k)` at model times `time(k)` [yr] (any order, distinct), on a grid of
        ! spacing dx, dy [m]: their signed distances, in increasing time.
        type(margin_t),   intent(out) :: mg
        real(dp),         intent(in)  :: time(:), dx, dy
        logical,          intent(in)  :: ice(:,:,:)
        integer,          intent(in)  :: interp
        logical,          intent(in)  :: hold
        character(len=*), intent(in)  :: src

        character(len=*), parameter :: proc = "kryos_plastic_margin::margin_set"
        integer, allocatable :: p(:)
        integer :: k

        if (size(time) /= size(ice, 3) .or. size(time) == 0) then
            call fesm_error(proc, "need one time per slice and at least one slice", trim(src))
        end if
        if (interp /= MARGIN_SDIST .and. interp /= MARGIN_NEAREST) then
            call fesm_error(proc, "unknown interpolation code "//to_str(interp), trim(src))
        end if
        if (.not. (dx > 0.0_dp .and. dy > 0.0_dp)) call fesm_error(proc, "the grid spacing must be positive", trim(src))

        p = time_order(time)
        do k = 2, size(p)
            if (.not. time(p(k)) > time(p(k-1))) then
                call fesm_error(proc, "two slices have the same time", trim(src)//": t = "//to_str(time(p(k)))//" yr")
            end if
        end do

        mg%n = size(time);  mg%nx = size(ice, 1);  mg%ny = size(ice, 2)
        mg%dx = dx;  mg%dy = dy;  mg%interp = interp;  mg%hold = hold;  mg%src = src
        mg%time = time(p)
        allocate(mg%sd(mg%nx, mg%ny, mg%n))
        do k = 1, mg%n
            mg%sd(:,:,k) = margin_signed_distance(ice(:,:,p(k)), dx, dy)
        end do
    end subroutine margin_set

    function margin_footprint(mg) result(fp)
        ! 1 where some slice has ice, else 0: every mask margin_at can give lies inside it (sdist: the
        ! interpolated distance is positive only where one of the two slices is).
        type(margin_t), intent(in) :: mg
        real(dp), allocatable :: fp(:,:)
        fp = merge(1.0_dp, 0.0_dp, any(mg%sd > 0.0_dp, dim=3))
    end function margin_footprint

    ! =========================================================================
    ! The mask at a time
    ! =========================================================================

    subroutine margin_at(mg, t, mask, sd)
        ! The ice mask at model time t [yr] and the signed distance [m] it is the positive part of.
        type(margin_t), intent(in)  :: mg
        real(dp),       intent(in)  :: t
        logical,        intent(out) :: mask(:,:)
        real(dp),       intent(out) :: sd(:,:)

        integer  :: k
        real(dp) :: a

        call bracket(mg, t, k, a)
        select case (mg%interp)
        case (MARGIN_SDIST)
            if (a == 0.0_dp) then
                sd = mg%sd(:,:,k)
            else
                sd = (1.0_dp - a)*mg%sd(:,:,k) + a*mg%sd(:,:,k+1)
            end if
        case (MARGIN_NEAREST)
            if (a > 0.5_dp) k = k + 1                   ! (a tie keeps the earlier slice)
            sd = mg%sd(:,:,k)
        end select
        mask = sd > 0.0_dp
    end subroutine margin_at

    subroutine bracket(mg, t, k, a)
        ! t = (1 - a) time(k) + a time(k+1), 0 <= a < 1 (a = 0 and k = n at the last slice). Outside
        ! the range: the end slice (a = 0) when holding, else stop.
        type(margin_t), intent(in)  :: mg
        real(dp),       intent(in)  :: t
        integer,        intent(out) :: k
        real(dp),       intent(out) :: a

        if (t < mg%time(1) .or. t > mg%time(mg%n)) then
            if (.not. mg%hold) then
                call fesm_error("kryos_plastic_margin::margin_at", "t = "//to_str(t)// &
                                " yr is outside the time range of the margin slices", trim(mg%src)//": ["// &
                                to_str(mg%time(1))//", "//to_str(mg%time(mg%n))//"] yr; margin_extrapolate = "// &
                                "'hold' holds the first and last slice")
            end if
            k = merge(1, mg%n, t < mg%time(1))
            a = 0.0_dp
            return
        end if
        if (t == mg%time(mg%n)) then
            k = mg%n
            a = 0.0_dp
            return
        end if
        k = 1
        do while (mg%time(k+1) <= t)
            k = k + 1
        end do
        a = (t - mg%time(k)) / (mg%time(k+1) - mg%time(k))
    end subroutine bracket

    ! =========================================================================
    ! Distances
    ! =========================================================================

    function margin_signed_distance(ice, dx, dy) result(sd)
        ! Signed distance to the margin of `ice` [m]: + inside, - outside, |sd| >= h = min(dx, dy)/2
        ! (see the module header); capped at the domain diagonal.
        logical,  intent(in) :: ice(:,:)
        real(dp), intent(in) :: dx, dy
        real(dp), allocatable :: sd(:,:)

        real(dp), allocatable :: d_out(:,:), d_in(:,:)
        real(dp) :: h, L
        integer  :: nx, ny

        nx = size(ice, 1);  ny = size(ice, 2)
        h = 0.5_dp*min(dx, dy)
        L = sqrt((real(nx, dp)*dx)**2 + (real(ny, dp)*dy)**2)
        d_out = sqrt(min(margin_edt2(.not. ice, dx, dy), L*L))     ! to the nearest ice-free cell
        d_in  = sqrt(min(margin_edt2(ice, dx, dy), L*L))           ! to the nearest ice cell
        allocate(sd(nx, ny))
        where (ice)
            sd = d_out - h
        elsewhere
            sd = -(d_in - h)
        end where
    end function margin_signed_distance

    function margin_edt2(fg, dx, dy) result(d2)
        ! Squared Euclidean distance [m2] from each cell centre to the nearest centre of a cell where
        ! fg holds (0 there), on a grid of spacing dx, dy; huge(1.0_dp) where fg holds nowhere. Exact:
        ! the lower envelope of parabolas along x, then along y (Felzenszwalb and Huttenlocher 2012).
        logical,  intent(in) :: fg(:,:)
        real(dp), intent(in) :: dx, dy
        real(dp), allocatable :: d2(:,:)

        integer :: i, j, nx, ny

        nx = size(fg, 1);  ny = size(fg, 2)
        allocate(d2(nx, ny))
        d2 = merge(0.0_dp, BIG, fg)
        do j = 1, ny
            d2(:,j) = envelope_1d(d2(:,j), dx)
        end do
        do i = 1, nx
            d2(i,:) = envelope_1d(d2(i,:), dy)
        end do
    end function margin_edt2

    function envelope_1d(f, s) result(d)
        ! d(q) = min_p ((q - p) s)**2 + f(p) over the p with a finite f (< BIG); BIG when there is none.
        real(dp), intent(in) :: f(:), s
        real(dp) :: d(size(f))

        integer  :: v(size(f))                      ! the parabolas of the envelope (their cell index)
        real(dp) :: z(size(f)+1)                    ! [m] where parabola k takes over from k-1
        integer  :: n, q, k, m
        real(dp) :: x, xs

        n = size(f)
        m = 0
        do q = 1, n
            if (f(q) >= BIG) cycle
            x = real(q - 1, dp)*s
            if (m == 0) then
                m = 1;  v(1) = q;  z(1) = -BIG;  z(2) = BIG
                cycle
            end if
            do                                      ! drop the parabolas the new one hides (z(1) = -BIG ends it)
                xs = ((f(q) + x*x) - (f(v(m)) + (real(v(m) - 1, dp)*s)**2)) / (2.0_dp*(x - real(v(m) - 1, dp)*s))
                if (xs > z(m)) exit
                m = m - 1
            end do
            m = m + 1;  v(m) = q;  z(m) = xs;  z(m+1) = BIG
        end do

        if (m == 0) then
            d = BIG
            return
        end if
        k = 1
        do q = 1, n
            x = real(q - 1, dp)*s
            do while (z(k+1) < x)
                k = k + 1
            end do
            d(q) = (x - real(v(k) - 1, dp)*s)**2 + f(v(k))
        end do
    end function envelope_1d

    ! =========================================================================
    ! Time axes
    ! =========================================================================

    function model_time_axis(t, units) result(tm)
        ! Model time [yr] of a file's time axis from its units (as recon_time_index): "...ago" or
        ! "...BP" count positive into the past; a leading "k" (ka, kyr, kiloyears) means 1000 years.
        real(dp),         intent(in) :: t(:)
        character(len=*), intent(in) :: units
        real(dp) :: tm(size(t))

        character(len=len(units)) :: u
        real(dp) :: scale

        u = adjustl(units)
        scale = 1.0_dp
        if (u(1:1) == "k") scale = 1.0e3_dp
        if (index(u, "ago") > 0 .or. index(u, "BP") > 0) scale = -scale
        tm = t*scale
    end function model_time_axis

    function time_order(t) result(p)
        ! The permutation that sorts t increasingly (insertion sort: a few dozen records).
        real(dp), intent(in) :: t(:)
        integer :: p(size(t))
        integer :: i, j, k
        p = [(i, i = 1, size(t))]
        do i = 2, size(t)
            k = p(i)
            j = i - 1
            do while (j >= 1)
                if (.not. t(p(j)) > t(k)) exit
                p(j+1) = p(j)
                j = j - 1
            end do
            p(j+1) = k
        end do
    end function time_order

end module kryos_plastic_margin
