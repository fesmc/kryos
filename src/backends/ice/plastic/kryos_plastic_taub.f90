module kryos_plastic_taub
    ! The basal shear stress tau_b [Pa] of the ice backend `plastic` (docs/dev/plastic-backend.md), on
    ! the instance grid at a model time. No library: compiled whatever `plastic=`.
    !
    !   "const"    one value everywhere;
    !   "file"     a static gridded field (x, y);
    !   "domains"  Gowan et al. (2021): a map of domain ids (0 = none) and one value per id, static or a
    !              series per id (linear in time), as the FesmData product {GRID}_TauB-Gowan2021.nc; a par
    !              table (codes, values) overrides the value of some ids, constant in time, or gives ids of
    !              the map the file has no value for; every id of the map must then have one. Cells without
    !              a domain take `default` (0: such a cell in the ice mask stops the solve, kryos_ice_plastic).

    use precision, only: dp
    use fesm_log,  only: fesm_error, to_str
    use ncio,      only: nc_read, nc_read_attr, nc_dims, nc_exists_var
    use kryos_plastic_margin, only: model_time_axis, time_order

    implicit none

    private

    integer, parameter, public :: TAUB_M_CONST = 1, TAUB_M_FILE = 2, TAUB_M_DOMAINS = 3

    type, public :: taub_t
        integer  :: method = 0
        real(dp) :: const = 0.0_dp                   ! [Pa] "const"
        real(dp) :: default = 0.0_dp                 ! [Pa] "domains": cells without a domain
        logical  :: hold = .false.                   ! "domains" series: hold the end values outside their range
        real(dp), allocatable :: field(:,:)          ! [Pa] "file"
        integer,  allocatable :: dom(:,:)            ! "domains": the domain id of each cell (<= 0: none)
        integer,  allocatable :: idx(:,:)            ! its position in ids (0: none)
        integer,  allocatable :: ids(:)              ! the domain ids with a value
        real(dp), allocatable :: val(:)              ! [Pa] their static value
        real(dp), allocatable :: val_t(:,:)          ! [Pa] (id, time): their series; unallocated: static
        real(dp), allocatable :: time(:)             ! [yr] model time of the series, increasing
        character(len=512) :: src = ""               ! messages
    end type taub_t

    public :: taub_method_code, taub_set_const, taub_set_field, taub_read_field
    public :: taub_set_domains, taub_read_domains, taub_at

contains

    integer function taub_method_code(name) result(code)
        ! "const" | "file" | "domains" -> TAUB_*; 0 for anything else.
        character(len=*), intent(in) :: name
        select case (trim(name))
        case ("const");   code = TAUB_M_CONST
        case ("file");    code = TAUB_M_FILE
        case ("domains"); code = TAUB_M_DOMAINS
        case default;     code = 0
        end select
    end function taub_method_code

    subroutine taub_set_const(tb, value)
        type(taub_t), intent(out) :: tb
        real(dp),     intent(in)  :: value
        if (.not. value > 0.0_dp) then
            call fesm_error("kryos_plastic_taub::taub_set_const", "taub_const must be positive", &
                            "taub_const = "//to_str(value)//" Pa")
        end if
        tb%method = TAUB_M_CONST
        tb%const  = value
        tb%src    = "taub_const"
    end subroutine taub_set_const

    subroutine taub_set_field(tb, field, src)
        type(taub_t),     intent(out) :: tb
        real(dp),         intent(in)  :: field(:,:)
        character(len=*), intent(in)  :: src
        tb%method = TAUB_M_FILE
        tb%field  = field
        tb%src    = src
    end subroutine taub_set_field

    subroutine taub_read_field(tb, file, var, nx, ny)
        ! The static field `var` (x, y) [Pa] of `file` on the instance grid.
        type(taub_t),     intent(out) :: tb
        character(len=*), intent(in)  :: file, var
        integer,          intent(in)  :: nx, ny

        real(dp), allocatable :: f(:,:)

        call check_xy(file, var, nx, ny)
        allocate(f(nx, ny))
        call nc_read(file, var, f)
        call taub_set_field(tb, f, trim(file)//" ("//trim(var)//")")
    end subroutine taub_read_field

    subroutine taub_set_domains(tb, dom, ids, val, default, src, val_t, time, hold, codes, values, table)
        ! The domain map `dom` (<= 0: none), the values `val` [Pa] of the ids `ids` and, optionally, their
        ! series val_t(id, time) at model times `time` [yr] (any order, distinct); then the override
        ! table `codes`, `values` [Pa] (named `table` in messages).
        type(taub_t),     intent(out) :: tb
        integer,          intent(in)  :: dom(:,:), ids(:)
        real(dp),         intent(in)  :: val(:), default
        character(len=*), intent(in)  :: src
        real(dp),         intent(in), optional :: val_t(:,:), time(:)
        logical,          intent(in), optional :: hold
        integer,          intent(in), optional :: codes(:)
        real(dp),         intent(in), optional :: values(:)
        character(len=*), intent(in), optional :: table

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_domains"
        integer, allocatable :: p(:)
        integer :: k

        if (size(ids) /= size(val)) call fesm_error(proc, "one value per domain id is needed", trim(src))
        do k = 1, size(ids)
            if (ids(k) <= 0) call fesm_error(proc, "domain ids are positive (0 = no domain)", trim(src)// &
                                             ": id "//to_str(ids(k)))
            if (count(ids == ids(k)) > 1) call fesm_error(proc, "domain id "//to_str(ids(k))//" is listed twice", &
                                                         trim(src))
        end do
        if (default < 0.0_dp) call fesm_error(proc, "taub_default must be >= 0", "taub_default = "//to_str(default))

        tb%method  = TAUB_M_DOMAINS
        tb%dom     = dom
        tb%ids     = ids
        tb%val     = val
        tb%default = default
        tb%src     = src
        if (present(hold)) tb%hold = hold
        if (present(val_t)) then
            if (.not. present(time)) call fesm_error(proc, "a series needs its times", trim(src))
            if (size(val_t, 1) /= size(ids) .or. size(val_t, 2) /= size(time) .or. size(time) == 0) then
                call fesm_error(proc, "the series is not (id, time)", trim(src))
            end if
            p = time_order(time)
            do k = 2, size(p)
                if (.not. time(p(k)) > time(p(k-1))) call fesm_error(proc, "two records of the series have the "// &
                                                                     "same time", trim(src))
            end do
            tb%time  = time(p)
            tb%val_t = val_t(:, p)
        end if
        if (present(codes)) then
            if (.not. (present(values) .and. present(table))) call fesm_error(proc, "an override needs its values", src)
            call override(tb, codes, values, table)
        end if
        call index_cells(tb)
    end subroutine taub_set_domains

    subroutine taub_read_domains(tb, file, dom_var, id_var, val_var, ser_var, default, hold, nx, ny, codes, values, table)
        ! The product {GRID}_TauB-Gowan2021.nc: dom_var (x, y) the domain id of each cell (0 = none),
        ! id_var (n) the ids, val_var (n) their tau_b [Pa]; ser_var (n, time) their series when the
        ! file has it (either dimension order; the time coordinate is the variable named like its
        ! dimension, converted from its units). The override table as taub_set_domains.
        type(taub_t),     intent(out) :: tb
        character(len=*), intent(in)  :: file, dom_var, id_var, val_var, ser_var
        real(dp),         intent(in)  :: default
        logical,          intent(in)  :: hold
        integer,          intent(in)  :: nx, ny
        integer,          intent(in), optional :: codes(:)
        real(dp),         intent(in), optional :: values(:)
        character(len=*), intent(in), optional :: table

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_domains"
        character(len=64), allocatable :: names(:)
        integer,  allocatable :: dims(:), dom(:,:), ids(:)
        real(dp), allocatable :: val(:), ser(:,:), t(:)
        character(len=256) :: units
        character(len=:), allocatable :: src
        integer :: n, it
        logical :: exists

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the tau_b domain file does not exist", "taub_domain_file = "//trim(file))
        src = trim(file)

        call check_xy(file, dom_var, nx, ny)
        allocate(dom(nx, ny))
        call nc_read(file, dom_var, dom)

        call nc_dims(file, id_var, names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`"//trim(id_var)//"` is not a list of ids (1-D)", src)
        n = dims(1)
        call nc_dims(file, val_var, names, dims)
        if (size(dims) /= 1 .or. dims(1) /= n) then
            call fesm_error(proc, "`"//trim(val_var)//"` is not one value per id ("//to_str(n)//", 1-D)", &
                            src//"; taub_value_var names the per-domain values")
        end if
        allocate(ids(n), val(n))
        call nc_read(file, id_var, ids)
        call nc_read(file, val_var, val)

        if (.not. nc_exists_var(file, ser_var)) then
            call taub_set_domains(tb, dom, ids, val, default, src, hold=hold, codes=codes, values=values, table=table)
            return
        end if

        call nc_dims(file, ser_var, names, dims)
        if (size(dims) /= 2) call fesm_error(proc, "`"//trim(ser_var)//"` is not a series (id, time)", src)
        if (dims(1) == n) then
            it = 2
        else if (dims(2) == n) then
            it = 1
        else
            call fesm_error(proc, "`"//trim(ser_var)//"` has no dimension of the "//to_str(n)//" ids", src)
        end if
        allocate(t(dims(it)), ser(dims(1), dims(2)))
        call nc_read(file, trim(names(it)), t)
        units = ""
        call nc_read_attr(file, trim(names(it)), "units", units)
        t = model_time_axis(t, units)
        call nc_read(file, ser_var, ser)
        if (it == 1) ser = transpose(ser)
        call taub_set_domains(tb, dom, ids, val, default, src, val_t=ser, time=t, hold=hold, codes=codes, values=values, &
                              table=table)
    end subroutine taub_read_domains

    subroutine override(tb, codes, values, src)
        ! The par table: the value [Pa] of the domain ids `codes`, constant in time; an id the file has
        ! no value for is added (it must occur in the map).
        type(taub_t),     intent(inout) :: tb
        integer,          intent(in)    :: codes(:)
        real(dp),         intent(in)    :: values(:)
        character(len=*), intent(in)    :: src

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_domains"
        real(dp), allocatable :: w(:,:)
        integer :: k, p

        if (size(codes) /= size(values)) call fesm_error(proc, "one value per code is needed", src)
        do k = 1, size(codes)
            if (.not. values(k) > 0.0_dp) call fesm_error(proc, "tau_b must be positive", src//": code "// &
                                                          to_str(codes(k))//", "//to_str(values(k))//" Pa")
            p = findloc(tb%ids, codes(k), dim=1)
            if (p == 0) then
                if (.not. any(tb%dom == codes(k))) then
                    call fesm_error(proc, "code "//to_str(codes(k))//" is not a domain of the map", src// &
                                    "; map "//trim(tb%src))
                end if
                if (allocated(tb%val_t)) then
                    allocate(w(size(tb%ids) + 1, size(tb%time)))
                    w(1:size(tb%ids), :) = tb%val_t
                    w(size(tb%ids) + 1, :) = values(k)
                    call move_alloc(w, tb%val_t)
                end if
                tb%ids = [tb%ids, codes(k)]
                tb%val = [tb%val, values(k)]
            else
                tb%val(p) = values(k)
                if (allocated(tb%val_t)) tb%val_t(p, :) = values(k)
            end if
        end do
    end subroutine override

    subroutine taub_at(tb, t, tau)
        ! tau_b [Pa] at model time t [yr] on every cell (the solver reads it on the ice mask).
        type(taub_t), intent(in)  :: tb
        real(dp),     intent(in)  :: t
        real(dp),     intent(out) :: tau(:,:)

        real(dp), allocatable :: v(:)
        integer :: i, j

        select case (tb%method)
        case (TAUB_M_CONST)
            tau = tb%const
        case (TAUB_M_FILE)
            tau = tb%field
        case (TAUB_M_DOMAINS)
            v = values_at(tb, t)
            do j = 1, size(tau, 2)
                do i = 1, size(tau, 1)
                    if (tb%idx(i,j) > 0) then
                        tau(i,j) = v(tb%idx(i,j))
                    else
                        tau(i,j) = tb%default
                    end if
                end do
            end do
        case default
            call fesm_error("kryos_plastic_taub::taub_at", "tau_b is not set", trim(tb%src))
        end select
    end subroutine taub_at

    function values_at(tb, t) result(v)
        ! The value of each id at t: static, or the series linear in time.
        type(taub_t), intent(in) :: tb
        real(dp),     intent(in) :: t
        real(dp), allocatable :: v(:)

        integer  :: k, n
        real(dp) :: a

        if (.not. allocated(tb%val_t)) then
            v = tb%val
            return
        end if
        n = size(tb%time)
        if (t < tb%time(1) .or. t > tb%time(n)) then
            if (.not. tb%hold) then
                call fesm_error("kryos_plastic_taub::taub_at", "t = "//to_str(t)//" yr is outside the time range "// &
                                "of the tau_b series", trim(tb%src)//": ["//to_str(tb%time(1))//", "// &
                                to_str(tb%time(n))//"] yr; margin_extrapolate = 'hold' holds the end values")
            end if
            v = tb%val_t(:, merge(1, n, t < tb%time(1)))
            return
        end if
        k = 1
        do while (k < n)
            if (tb%time(k+1) > t) exit
            k = k + 1
        end do
        if (k == n) then
            v = tb%val_t(:, n)
        else
            a = (t - tb%time(k)) / (tb%time(k+1) - tb%time(k))
            v = (1.0_dp - a)*tb%val_t(:, k) + a*tb%val_t(:, k+1)
        end if
    end function values_at

    subroutine index_cells(tb)
        ! idx: the position in ids of each cell's domain; a domain of the map without a value stops.
        type(taub_t), intent(inout) :: tb
        integer :: i, j, p

        if (allocated(tb%idx)) deallocate(tb%idx)
        allocate(tb%idx(size(tb%dom, 1), size(tb%dom, 2)))
        tb%idx = 0
        do j = 1, size(tb%dom, 2)
            do i = 1, size(tb%dom, 1)
                if (tb%dom(i,j) <= 0) cycle
                p = findloc(tb%ids, tb%dom(i,j), dim=1)
                if (p == 0) then
                    call fesm_error("kryos_plastic_taub::taub_set_domains", "domain "//to_str(tb%dom(i,j))// &
                                    " of the map has no tau_b", trim(tb%src)//" (cell "//to_str(i)//", "//to_str(j)// &
                                    "); give it one in the override table (taub_override)")
                end if
                tb%idx(i,j) = p
            end do
        end do
    end subroutine index_cells

    subroutine check_xy(file, var, nx, ny)
        ! `var` of `file` is a field (x, y) on the instance grid.
        character(len=*), intent(in) :: file, var
        integer,          intent(in) :: nx, ny
        character(len=64), allocatable :: names(:)
        integer, allocatable :: dims(:)

        call nc_dims(file, var, names, dims)
        if (size(dims) /= 2) then
            call fesm_error("kryos_plastic_taub::read", "`"//trim(var)//"` is not a static field (x, y)", &
                            "file "//trim(file)//", "//to_str(size(dims))//" dimensions")
        end if
        if (dims(1) /= nx .or. dims(2) /= ny) then
            call fesm_error("kryos_plastic_taub::read", "`"//trim(var)//"` is not on the instance grid", "file "// &
                            trim(file)//": "//to_str(dims(1))//" x "//to_str(dims(2))//", the grid "//to_str(nx)// &
                            " x "//to_str(ny))
        end if
    end subroutine check_xy

end module kryos_plastic_taub
