module kryos_outfile
    ! File helpers shared by the default writer (kryos_component), the backends'
    ! own writers and the output module (core-design section 6): the variable_io
    ! tables `input/kryos-variables-<component>.md`, the creation of a 2D or time
    ! series file, the open-a-record step, and one 2D record written with
    ! bit-rounding precision. Free procedures only; no kryos module is used, so
    ! backends and the core may all use it.
    !
    ! Table cache: a table is loaded the first time a component kind needs it and
    ! kept for the run (variable_io prints one line per row on load).

    use precision,   only: sp, dp
    use fesm_log,    only: fesm_error
    use ncio,        only: nc_create, nc_write_dim, nc_open, nc_write, nc_time_index
    use variable_io, only: var_io_type, load_var_io_table

    implicit none

    private

    character(len=*), parameter, public :: OUTFILE_TABLE_DIR = "input"
    integer, parameter :: PATH_LEN = 512, NAME_LEN = 64

    type :: table_t                         ! one loaded table
        character(len=NAME_LEN) :: kind = ""
        character(len=PATH_LEN) :: file = ""
        type(var_io_type), allocatable :: v(:)
    end type table_t

    type(table_t), allocatable, save :: cache(:)
    character(len=PATH_LEN), save :: table_dir = OUTFILE_TABLE_DIR

    interface outfile_var2d
        ! (file, name, var, n, ncid, units, long_name [, prec_abs, prec_rel]): record n of a
        ! (xc, yc, time) variable; var and the precisions have the same kind.
        module procedure outfile_var2d_sp, outfile_var2d_dp
    end interface outfile_var2d

    public :: outfile_set_table_dir, outfile_table_file, outfile_load, outfile_row
    public :: outfile_path, outfile_dims_2d, outfile_dims_ts, outfile_axis, outfile_open, outfile_var2d

contains

    ! =========================================================================
    ! Variable tables
    ! =========================================================================

    subroutine outfile_set_table_dir(dir)
        ! Folder of the tables (default `input`, relative to the working directory).
        ! Drops the cache: tables are reloaded from the new folder.
        character(len=*), intent(in) :: dir
        table_dir = dir
        if (allocated(cache)) deallocate(cache)
    end subroutine outfile_set_table_dir

    function outfile_table_file(kind) result(file)
        ! `<dir>/kryos-variables-<kind>.md`.
        character(len=*), intent(in) :: kind
        character(len=:), allocatable :: file
        file = trim(table_dir)//"/kryos-variables-"//trim(kind)//".md"
    end function outfile_table_file

    subroutine outfile_load(kind, tab, stat)
        ! The table of component kind `kind` (cached, so one load per run). Stops,
        ! or sets stat = 1, when the file is missing.
        character(len=*),               intent(in)            :: kind
        type(var_io_type), allocatable, intent(out)           :: tab(:)
        integer,                        intent(out), optional :: stat

        integer :: k

        if (present(stat)) stat = 0
        k = find(kind)
        if (k == 0) then
            call load(kind, k, stat)
            if (k == 0) return
        end if
        tab = cache(k)%v
    end subroutine outfile_load

    subroutine outfile_row(kind, name, row, found)
        ! The row `name` of the table of `kind`; found = .false. when it has none
        ! (the caller decides whether that stops). A missing table file stops.
        character(len=*),  intent(in)  :: kind, name
        type(var_io_type), intent(out) :: row
        logical,           intent(out) :: found

        integer :: k, i

        found = .false.
        k = find(kind)
        if (k == 0) call load(kind, k)
        do i = 1, size(cache(k)%v)
            if (trim(cache(k)%v(i)%varname) == trim(name)) then
                row = cache(k)%v(i)
                found = .true.
                return
            end if
        end do
    end subroutine outfile_row

    integer function find(kind) result(k)
        character(len=*), intent(in) :: kind
        integer :: i
        k = 0
        if (.not. allocated(cache)) return
        do i = 1, size(cache)
            if (trim(cache(i)%kind) == trim(kind)) then
                k = i
                return
            end if
        end do
    end function find

    subroutine load(kind, k, stat)
        ! Load the table of `kind` into the cache; k = its index (0 on failure).
        character(len=*), intent(in)            :: kind
        integer,          intent(out)           :: k
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_outfile::load"
        type(table_t) :: t
        type(table_t), allocatable :: tmp(:)
        logical :: exists
        integer :: i

        k = 0
        t%kind = kind
        t%file = outfile_table_file(kind)
        inquire(file=trim(t%file), exist=exists)
        if (.not. exists) then
            call fesm_error(proc, "variable table not found: "//trim(t%file), &
                            "component kind `"//trim(kind)//"` needs input/kryos-variables-"//trim(kind)//".md", stat)
            return
        end if
        call load_var_io_table(t%v, trim(t%file))

        ! Grow by one, element by element (ifx stops on an array-section copy of
        ! a derived type with allocatable components from a zero-size array).
        if (.not. allocated(cache)) allocate(cache(0))
        allocate(tmp(size(cache) + 1))
        do i = 1, size(cache)
            tmp(i) = cache(i)
        end do
        tmp(size(tmp)) = t
        call move_alloc(tmp, cache)
        k = size(cache)
    end subroutine load

    ! =========================================================================
    ! Files
    ! =========================================================================

    function outfile_path(dir, base) result(path)
        ! `base` in `dir`: the folder with one separator, and no prefix at all for
        ! the current directory ("", "." and "./"), so a file keeps the plain name
        ! the backends compare it with.
        character(len=*), intent(in) :: dir, base
        character(len=:), allocatable :: path
        integer :: n

        n = len_trim(dir)
        if (n == 0 .or. trim(dir) == "." .or. trim(dir) == "./") then
            path = trim(base)
        else if (dir(n:n) == "/") then
            path = trim(dir)//trim(base)
        else
            path = trim(dir)//"/"//trim(base)
        end if
    end function outfile_path

    subroutine outfile_dims_2d(file, x, y, t)
        ! A new 2D file: xc, yc [km] and the unlimited time axis [year] at `t`
        ! (yelmox io_dims_init).
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: x(:), y(:), t
        call nc_create(file)
        call nc_write_dim(file, "xc", x=x, units="km")
        call nc_write_dim(file, "yc", x=y, units="km")
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
    end subroutine outfile_dims_2d

    subroutine outfile_dims_ts(file, t)
        ! A new time series file: the unlimited time axis [year] at `t`.
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: t
        call nc_create(file)
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
    end subroutine outfile_dims_ts

    subroutine outfile_axis(file, name, v, units)
        ! A further axis (third axis of rank-3 variables: month, depth).
        character(len=*), intent(in) :: file, name, units
        real(dp),         intent(in) :: v(:)
        call nc_write_dim(file, name, x=v, units=units)
    end subroutine outfile_axis

    subroutine outfile_open(file, t, ncid, n)
        ! Open `file` for writing the record of time `t`: n is its index on the
        ! time axis (a new record unless `t` is already there) and the time value
        ! is written. The caller closes with nc_close(ncid).
        character(len=*), intent(in)  :: file
        real(dp),         intent(in)  :: t
        integer,          intent(out) :: ncid, n
        call nc_open(file, ncid, writable=.true.)
        n = nc_time_index(file, "time", t, ncid)
        call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)
    end subroutine outfile_open

    subroutine outfile_var2d_sp(file, vnm, var, n, ncid, units, long_name, prec_abs, prec_rel)
        ! One 2D record (yelmox io_var2D); precision as ncio bit rounding.
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(sp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        real(sp),         intent(in), optional :: prec_abs, prec_rel
        call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                      units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
    end subroutine outfile_var2d_sp

    subroutine outfile_var2d_dp(file, vnm, var, n, ncid, units, long_name, prec_abs, prec_rel)
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(dp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        real(dp),         intent(in), optional :: prec_abs, prec_rel
        call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                      units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
    end subroutine outfile_var2d_dp

end module kryos_outfile
