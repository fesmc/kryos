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
    !
    ! Tables (`input/kryos-variables-<name>.md`): the table of a component kind
    ! (`<kind>`: the exports of the instance, the default writer) and the table of
    ! a backend's own stream (`<kind>-<backend>`, or `hub`: the variables the backend
    ! writes itself). The `prec` column of a row is the only source of rounding;
    ! a variable the table does not list is written lossless. Every row has a
    ! unit (a row without one stops at the load).

    use precision,   only: sp, dp
    use fesm_log,    only: fesm_error
    use ncio,        only: nc_create, nc_write_dim, nc_open, nc_write, nc_time_index
    use variable_io, only: var_io_type, load_var_io_table
    use fesm_fields, only: units_equivalent

    implicit none

    private

    character(len=*), parameter, public :: OUTFILE_TABLE_DIR = "input"
    integer, parameter :: PATH_LEN = 512, NAME_LEN = 64

    ! The cache is a list of nodes, each allocated once and never copied or
    ! moved. A growing allocatable array of table_t fails under ifx 2023.2: its
    ! intrinsic copy stops with severe (189) "incompatible types", and after
    ! move_alloc of an element's v the deallocation of the module array stops
    ! with severe (153) "not allocated" (both standard code; gfortran is right).
    type :: table_t                         ! one loaded table
        character(len=NAME_LEN) :: kind = ""
        character(len=PATH_LEN) :: file = ""
        type(var_io_type), allocatable :: v(:)
        type(table_t), pointer :: next => null()
    end type table_t

    type(table_t), pointer, save :: cache => null()   ! the most recently loaded table first
    character(len=PATH_LEN), save :: table_dir = OUTFILE_TABLE_DIR

    interface outfile_var2d
        ! (file, name, var, n, ncid, units, long_name [, prec_abs, prec_rel, table]): record n of
        ! a (xc, yc, time) variable; var and the precisions have the same kind. With `table`
        ! (a backend table, see outfile_prec) the precision and the long name are the row's
        ! and prec_abs, prec_rel must be absent; a variable the table does not list is
        ! written lossless with the long_name given.
        module procedure outfile_var2d_sp, outfile_var2d_dp
    end interface outfile_var2d

    interface outfile_var3d
        ! (file, name, var, dim3, n, ncid, units, long_name [, table]): record n of a
        ! (xc, yc, dim3, time) variable, the precision from the backend table `table`
        ! (absent: lossless).
        module procedure outfile_var3d_sp, outfile_var3d_dp
    end interface outfile_var3d

    public :: outfile_set_table_dir, outfile_table_file, outfile_load, outfile_row, outfile_prec
    public :: outfile_path, outfile_dims_2d, outfile_dims_lonlat, outfile_dims_ts, outfile_axis, outfile_open, outfile_var2d
    public :: outfile_var3d

contains

    ! =========================================================================
    ! Variable tables
    ! =========================================================================

    subroutine outfile_set_table_dir(dir)
        ! Folder of the tables (default `input`, relative to the working directory).
        ! Drops the cache: tables are reloaded from the new folder.
        character(len=*), intent(in) :: dir
        type(table_t), pointer :: t
        table_dir = dir
        do while (associated(cache))
            t => cache
            cache => cache%next
            deallocate(t)
        end do
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

        type(table_t), pointer :: t

        if (present(stat)) stat = 0
        t => find(kind)
        if (.not. associated(t)) then
            call load(kind, t, stat)
            if (.not. associated(t)) return
        end if
        tab = t%v
    end subroutine outfile_load

    subroutine outfile_row(kind, name, row, found)
        ! The row `name` of the table of `kind`; found = .false. when it has none
        ! (the caller decides whether that stops). A missing table file stops.
        character(len=*),  intent(in)  :: kind, name
        type(var_io_type), intent(out) :: row
        logical,           intent(out) :: found

        type(table_t), pointer :: t
        integer :: i

        found = .false.
        t => find(kind)
        if (.not. associated(t)) call load(kind, t)
        do i = 1, size(t%v)
            if (trim(t%v(i)%varname) == trim(name)) then
                row = t%v(i)
                found = .true.
                return
            end if
        end do
    end subroutine outfile_row

    subroutine outfile_prec(table, name, units, long_name, prec_abs, prec_rel)
        ! The write precision of variable `name` of a backend's own stream, from the
        ! table `table` (`<kind>-<backend>` or `hub`: input/kryos-variables-<table>.md).
        ! A listed variable: the row's prec_abs and prec_rel (0 = none) and its long
        ! name replace the `long_name` given; the row's unit must be the `units` of the
        ! data, or the run stops. A variable the table does not list: lossless
        ! (prec_abs = prec_rel = 0), `long_name` unchanged. A missing table stops.
        character(len=*),  intent(in)    :: table, name, units
        character(len=*),  intent(inout) :: long_name
        real(dp),          intent(out)   :: prec_abs, prec_rel

        type(var_io_type) :: row
        logical :: found

        prec_abs = 0.0_dp
        prec_rel = 0.0_dp
        call outfile_row(table, name, row, found)
        if (.not. found) return
        if (.not. (trim(row%units) == trim(units) .or. units_equivalent(row%units, units))) then
            call fesm_error("kryos_outfile::outfile_prec", "units `"//trim(row%units)//"` of the table differ from "// &
                            "the data's `"//trim(units)//"`", "table "//outfile_table_file(table)// &
                            ", variable `"//trim(name)//"`")
        end if
        prec_abs = row%prec_abs
        prec_rel = row%prec_rel
        if (len_trim(row%long_name) > 0) long_name = row%long_name
    end subroutine outfile_prec

    function find(kind) result(t)
        ! The cached table of `kind`, null when not loaded.
        character(len=*), intent(in) :: kind
        type(table_t), pointer :: t
        t => cache
        do while (associated(t))
            if (trim(t%kind) == trim(kind)) return
            t => t%next
        end do
    end function find

    subroutine load(kind, t, stat)
        ! Load the table of `kind` into the cache; t = its node (null on failure).
        character(len=*),       intent(in)            :: kind
        type(table_t), pointer, intent(out)           :: t
        integer,                intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_outfile::load"
        character(len=PATH_LEN) :: file
        logical :: exists
        integer :: i

        t => null()
        file = outfile_table_file(kind)
        inquire(file=trim(file), exist=exists)
        if (.not. exists) then
            call fesm_error(proc, "variable table not found: "//trim(file), &
                            "component kind `"//trim(kind)//"` needs input/kryos-variables-"//trim(kind)//".md", stat)
            return
        end if

        allocate(t)
        t%kind = kind
        t%file = file
        call load_var_io_table(t%v, trim(file))       ! stops on a bad `prec` entry
        do i = 1, size(t%v)
            if (len_trim(t%v(i)%units) == 0) then
                call fesm_error(proc, "variable `"//trim(t%v(i)%varname)//"` has no unit in "//trim(file), &
                                "every row needs a unit (`1` for a dimensionless variable)", stat)
                deallocate(t)
                t => null()
                return
            end if
        end do
        t%next => cache
        cache => t
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

    subroutine outfile_dims_lonlat(file, lon, lat, t)
        ! A new 2D file of a global lon-lat grid (a Gaussian grid included): lon, lat [degrees] and the
        ! unlimited time axis [year] at `t`.
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: lon(:), lat(:), t
        call nc_create(file)
        call nc_write_dim(file, "lon", x=lon, units="degrees_east")
        call nc_write_dim(file, "lat", x=lat, units="degrees_north")
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
    end subroutine outfile_dims_lonlat

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

    subroutine outfile_var2d_sp(file, vnm, var, n, ncid, units, long_name, prec_abs, prec_rel, table)
        ! One 2D record (yelmox io_var2D); precision as ncio bit rounding.
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(sp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        real(sp),         intent(in), optional :: prec_abs, prec_rel
        character(len=*), intent(in), optional :: table

        character(len=len(long_name)+256) :: lname
        real(dp) :: pa, pr

        lname = long_name
        if (present(table)) then
            call outfile_prec(table, vnm, units, lname, pa, pr)
            call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                          units=units, long_name=trim(lname), prec_abs=real(pa, sp), prec_rel=real(pr, sp))
        else
            call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                          units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end if
    end subroutine outfile_var2d_sp

    subroutine outfile_var2d_dp(file, vnm, var, n, ncid, units, long_name, prec_abs, prec_rel, table)
        character(len=*), intent(in) :: file, vnm, units, long_name
        real(dp),         intent(in) :: var(:,:)
        integer,          intent(in) :: n, ncid
        real(dp),         intent(in), optional :: prec_abs, prec_rel
        character(len=*), intent(in), optional :: table

        character(len=len(long_name)+256) :: lname
        real(dp) :: pa, pr

        lname = long_name
        if (present(table)) then
            call outfile_prec(table, vnm, units, lname, pa, pr)
            call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                          units=units, long_name=trim(lname), prec_abs=pa, prec_rel=pr)
        else
            call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3="time", &
                          start=[1,1,n], count=[size(var,1), size(var,2), 1], ncid=ncid, &
                          units=units, long_name=long_name, prec_abs=prec_abs, prec_rel=prec_rel)
        end if
    end subroutine outfile_var2d_dp

    subroutine outfile_var3d_sp(file, vnm, var, dim3, n, ncid, units, long_name, table)
        ! One record of a (xc, yc, dim3, time) variable, as outfile_var2d.
        character(len=*), intent(in) :: file, vnm, dim3, units, long_name
        real(sp),         intent(in) :: var(:,:,:)
        integer,          intent(in) :: n, ncid
        character(len=*), intent(in), optional :: table

        character(len=len(long_name)+256) :: lname
        real(dp) :: pa, pr

        lname = long_name
        pa = 0.0_dp
        pr = 0.0_dp
        if (present(table)) call outfile_prec(table, vnm, units, lname, pa, pr)
        call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                      start=[1,1,1,n], count=[size(var,1), size(var,2), size(var,3), 1], ncid=ncid, &
                      units=units, long_name=trim(lname), prec_abs=real(pa, sp), prec_rel=real(pr, sp))
    end subroutine outfile_var3d_sp

    subroutine outfile_var3d_dp(file, vnm, var, dim3, n, ncid, units, long_name, table)
        character(len=*), intent(in) :: file, vnm, dim3, units, long_name
        real(dp),         intent(in) :: var(:,:,:)
        integer,          intent(in) :: n, ncid
        character(len=*), intent(in), optional :: table

        character(len=len(long_name)+256) :: lname
        real(dp) :: pa, pr

        lname = long_name
        pa = 0.0_dp
        pr = 0.0_dp
        if (present(table)) call outfile_prec(table, vnm, units, lname, pa, pr)
        call nc_write(file, vnm, var, dim1="xc", dim2="yc", dim3=dim3, dim4="time", &
                      start=[1,1,1,n], count=[size(var,1), size(var,2), size(var,3), 1], ncid=ncid, &
                      units=units, long_name=trim(lname), prec_abs=pa, prec_rel=pr)
    end subroutine outfile_var3d_dp

end module kryos_outfile
