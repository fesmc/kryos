module test_output_toys
    ! Toy components for test_output: exports declared directly (no backend), a
    ! write_output override that logs every call and passes streams 2D and 1D on
    ! to the default writer, as the backends do.

    use precision,       only: sp, dp
    use fesm_fields,     only: fields_init, fields_add, fields_lock, fields_put, fields_finalize_interval, &
                               FIELD_INSTANT, FIELD_MEAN, FIELD_INTENSIVE, FIELD_EXTENSIVE
    use kryos_component, only: component_t, init_ctx_t, step_time_t, constraints_t, axis_t, &
                               component_identity, component_write_output, month_axis

    implicit none

    private

    integer, parameter, public :: NX = 4, NY = 3, NLOG = 32

    type, extends(component_t), public :: toy_t
        integer :: ncall = 0
        character(len=32)  :: log_stream(NLOG) = ""
        real(dp)           :: log_t(NLOG) = 0.0_dp
        logical            :: log_init(NLOG) = .false.
        character(len=256) :: log_file(NLOG) = ""
        integer :: h(8) = 0
    contains
        procedure :: init          => toy_init
        procedure :: update        => toy_update
        procedure :: restart_write => toy_rw
        procedure :: restart_read  => toy_rr
        procedure :: write_output  => toy_write
    end type toy_t

    public :: toy_make, toy_publish, toy_field, toy_calls, toy_call_init

contains

    subroutine toy_make(c, ctx, kind)
        ! Identity from ctx, then the exports of the kind: ice (2D H_ice, z_bed, MEAN smb; 0-D M_ice),
        ! surface (MEAN T_srf, smb), climate (month tas, 2D z_srf_clim, 0-D dTo), geo (0-D bsl, A_ocn).
        type(toy_t),      intent(out) :: c
        type(init_ctx_t), intent(in)  :: ctx
        character(len=*), intent(in)  :: kind

        character(len=:), allocatable :: g
        integer :: d2(2)

        call component_identity(c, ctx, kind, constraints_t())
        c%nx = NX
        c%ny = NY
        g  = trim(c%grid)
        d2 = [NX, NY]
        call fields_init(c%exports, trim(c%name)//":exports")
        select case (kind)
        case ("ice")
            call fields_add(c%exports, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, c%h(1), grid=g, dims=d2)
            call fields_add(c%exports, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h(2), grid=g, dims=d2)
            call fields_add(c%exports, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h(3), grid=g, dims=d2)
            call fields_add(c%exports, "M_ice", "Gt", FIELD_INSTANT, FIELD_EXTENSIVE, c%h(4))
        case ("surface")
            call fields_add(c%exports, "T_srf", "K", FIELD_MEAN, FIELD_INTENSIVE, c%h(1), grid=g, dims=d2)
            call fields_add(c%exports, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, c%h(2), grid=g, dims=d2)
        case ("climate")
            allocate(c%axes(1))
            c%axes(1) = month_axis()
            call fields_add(c%exports, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h(1), grid=g, &
                            dims=[NX, NY, 12], zdim="month")
            call fields_add(c%exports, "z_srf_clim", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h(2), grid=g, dims=d2)
            call fields_add(c%exports, "dTo", "K", FIELD_INSTANT, FIELD_INTENSIVE, c%h(3))
        case ("geo")
            call fields_add(c%exports, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h(1))
            call fields_add(c%exports, "A_ocn", "m2", FIELD_INSTANT, FIELD_EXTENSIVE, c%h(2))
        end select
        call fields_lock(c%exports)
    end subroutine toy_make

    pure function toy_field(t, i, j, k) result(v)
        ! The value the toys put at time t in cell (i, j) and layer k.
        real(dp), intent(in) :: t
        integer,  intent(in) :: i, j, k
        real(dp) :: v
        v = 1000.1234567_dp + real(i, dp) + 10.0_dp*real(j, dp) + 100.0_dp*real(k, dp) + t
    end function toy_field

    subroutine toy_publish(c, t)
        ! Put every export at t (MEAN ones over [t-10, t]) and close the window.
        type(toy_t), intent(inout) :: c
        real(dp),    intent(in)    :: t

        real(dp) :: a(NX, NY), a3(NX, NY, 12)
        integer :: i, j, k

        do j = 1, NY
            do i = 1, NX
                a(i, j) = toy_field(t, i, j, 0)
                do k = 1, 12
                    a3(i, j, k) = toy_field(t, i, j, k)
                end do
            end do
        end do
        select case (trim(c%kind))
        case ("ice")
            call fields_put(c%exports, c%h(1), a, t1=t)
            call fields_put(c%exports, c%h(2), a - 2000.0_dp, t1=t)
            call fields_put(c%exports, c%h(3), 3.0_dp*a, t1=t, t0=t - 10.0_dp)
            call fields_put(c%exports, c%h(4), 2.0_dp*t + 5.0_dp, t1=t)
            call fields_finalize_interval(c%exports, t)
        case ("surface")
            call fields_put(c%exports, c%h(1), a, t1=t, t0=t - 10.0_dp)
            call fields_put(c%exports, c%h(2), 3.0_dp*a, t1=t, t0=t - 10.0_dp)
            call fields_finalize_interval(c%exports, t)
        case ("climate")
            call fields_put(c%exports, c%h(1), a3, t1=t)
            call fields_put(c%exports, c%h(2), a, t1=t)
            call fields_put(c%exports, c%h(3), toy_field(t, 1, 1, 0), t1=t)
        case ("geo")
            call fields_put(c%exports, c%h(1), 0.5_dp*t, t1=t)
            call fields_put(c%exports, c%h(2), 3.6e14_dp, t1=t)
        end select
    end subroutine toy_publish

    integer function toy_calls(c, stream) result(n)
        ! Calls of write_output with `stream`.
        type(toy_t),      intent(in) :: c
        character(len=*), intent(in) :: stream
        n = count(c%log_stream(1:c%ncall) == stream)
    end function toy_calls

    function toy_call_init(c, stream) result(init)
        ! The init flags of the calls with `stream`, in order.
        type(toy_t),      intent(in) :: c
        character(len=*), intent(in) :: stream
        logical, allocatable :: init(:)
        init = pack(c%log_init(1:c%ncall), c%log_stream(1:c%ncall) == stream)
    end function toy_call_init

    subroutine toy_write(self, file, t, stream, init)
        class(toy_t),     intent(inout) :: self
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init

        self%ncall = self%ncall + 1
        self%log_stream(self%ncall) = stream
        self%log_t(self%ncall)      = t
        self%log_init(self%ncall)   = init
        self%log_file(self%ncall)   = file
        if (trim(stream) == "2D" .or. trim(stream) == "1D") call component_write_output(self, file, t, stream, init)
    end subroutine toy_write

    subroutine toy_init(self, ctx)
        class(toy_t),     intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        associate (unused => ctx%time + self%nx); end associate
    end subroutine toy_init

    subroutine toy_update(self, st)
        class(toy_t),      intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (unused => st%t1 + self%nx); end associate
    end subroutine toy_update

    subroutine toy_rw(self, dir, t)
        class(toy_t),     intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (unused => t + self%nx + len(dir)); end associate
    end subroutine toy_rw

    subroutine toy_rr(self, dir, t)
        class(toy_t),     intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (unused => t + self%nx + len(dir)); end associate
    end subroutine toy_rr

end module test_output_toys


program test_output
    ! kryos_output and the default writer (core C8), on toy components and a toy
    ! hub: registration from the out_<stream> keys of the instance groups, one
    ! schedule per distinct group, dispatch over a run, files with the right
    ! names, variables, units, time records and values, out_vars filtering,
    ! a backend stream routed to the backend's write_output, init versus later
    ! records, a writer producer (hub) and the timing table, the registration
    ! errors, [kryos] write_lossless, and the variable tables against the
    ! vocabulary of contracts-design section 2. Output precision (phase 3, WP3):
    ! every row of every table round-trips within its tolerance (instance tables
    ! through the variable_io writer, backend tables through the backend helpers),
    ! the moved hard-coded precisions give the values the literals gave, a variable
    ! a backend table does not list is lossless, write_lossless gives exact
    ! values, time series are lossless, the restart state of a rounded run is a
    ! lossless run's, and a bad prec entry, a row without a unit and a unit that
    ! differs from the data's stop (child processes of this executable).
    !
    ! Files go to tmp/test_output/ (run from the repo root, as make check does).
    ! Abort paths take `stat` (fesm_error returns instead of stopping).

    use precision,        only: sp, dp
    use coordinates,      only: grid_init
    use ncio,             only: nc_read, nc_size, nc_exists_var, nc_read_attr, nc_create, nc_write_dim, nc_write, &
                                nc_close, nc_set_lossless
    use variable_io,      only: var_io_type, var_io_write
    use fesm_fields,      only: fields_restart_write
    use timer,            only: timer_class, timer_init, timer_step
    use fesm_clock,       only: clock_t, clock_init, clock_start, clock_advance, clock_finished, clock_time
    use kryos_component,  only: init_ctx_t
    use kryos_config,     only: config_t, config_init
    use kryos_outfile,    only: outfile_load, outfile_set_table_dir, outfile_path, outfile_row, outfile_open, &
                                outfile_var2d, outfile_var3d, outfile_dims_2d, outfile_axis
    use kryos_output
    use test_output_toys

    implicit none

    character(len=*), parameter :: dir = "tmp/test_output/", rundir = "tmp/test_output/run/"
    character(len=*), parameter :: par = dir//"par.nml"
    integer :: nfail = 0
    character(len=64) :: arg

    type :: hub_ctx_t                           ! context of the toy hub writer
        integer :: ncall = 0
        real(dp) :: t(16) = 0.0_dp
        logical :: init(16) = .false.
        character(len=256) :: file = ""
    end type hub_ctx_t

    call get_command_argument(1, arg)
    if (len_trim(arg) > 0) then                 ! child process: one abort case
        call abort_case(trim(arg))
        stop 0
    end if

    call execute_command_line("rm -rf "//dir//" && mkdir -p "//dir)
    call write_par(par, rundir, lossless=.false.)

    call test_tables()
    call test_run()
    call test_errors()
    call test_lossless()
    call test_precision()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_output: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_output: OK"

contains

    ! -------------------------------------------------------------------------
    ! helpers
    ! -------------------------------------------------------------------------

    subroutine check(ok, label)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: label
        if (ok) then
            write(*,'(a)') "  ok    "//label
        else
            write(*,'(a)') "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    logical function exists(file)
        character(len=*), intent(in) :: file
        inquire(file=file, exist=exists)
    end function exists

    subroutine write_par(file, outdir, lossless)
        character(len=*), intent(in) :: file, outdir
        logical,          intent(in) :: lossless
        integer :: u

        open(newunit=u, file=file, status="replace", action="write")
        write(u,'(a)') "&kryos"
        write(u,'(a)') "    out_dir = '"//outdir//"'"
        write(u,'(a)') "    write_lossless = "//merge(".true. ", ".false.", lossless)
        write(u,'(a)') "/"
        write(u,'(a)') "&out_a"                    ! 0 20 40 60
        write(u,'(a)') "    method = 'const'"
        write(u,'(a)') "    dt = 20.0"
        write(u,'(a)') "/"
        write(u,'(a)') "&out_b"                    ! 10 30 50 and the end, 60
        write(u,'(a)') "    method = 'times'"
        write(u,'(a)') "    times = '10 30 50'"
        write(u,'(a)') "/"
        write(u,'(a)') "&ice_toy"                  ! template [instance] + [ice_yelmo]: backend streams
        write(u,'(a)') "    backend = 'yelmo'"
        write(u,'(a)') "    out_2D = 'out_a'"
        write(u,'(a)') "    out_1D = 'out_b'"
        write(u,'(a)') "    out_yelmo = 'out_a'"
        write(u,'(a)') "    out_yelmo_ts = 'out_b'"
        write(u,'(a)') "    out_vars = 'H_ice' 'smb' 'M_ice'"
        write(u,'(a)') "/"
        write(u,'(a)') "&sf_toy"
        write(u,'(a)') "    backend = 'chion'"
        write(u,'(a)') "    out_2D = 'out_b'"
        write(u,'(a)') "/"
        write(u,'(a)') "&clim_toy"
        write(u,'(a)') "    backend = 'snapshots'"
        write(u,'(a)') "    out_2D = 'out_a'"
        write(u,'(a)') "    out_1D = 'out_a'"
        write(u,'(a)') "/"
        write(u,'(a)') "&geo_glob"
        write(u,'(a)') "    backend = 'fastisostasy'"
        write(u,'(a)') "    out_1D = 'out_a'"
        write(u,'(a)') "/"
        write(u,'(a)') "&Dom"
        write(u,'(a)') "    grid_hub = 'TOY'"
        write(u,'(a)') "    out_2D = 'out_b'"
        write(u,'(a)') "/"
        write(u,'(a)') "&ice_bad"                  ! out_vars names a non-export
        write(u,'(a)') "    backend = 'yelmo'"
        write(u,'(a)') "    out_2D = 'out_a'"
        write(u,'(a)') "    out_vars = 'H_ice' 'bogus'"
        write(u,'(a)') "/"
        write(u,'(a)') "&ice_nosched"              ! schedule group missing
        write(u,'(a)') "    backend = 'yelmo'"
        write(u,'(a)') "    out_2D = 'out_nowhere'"
        write(u,'(a)') "/"
        close(u)
    end subroutine write_par

    subroutine ctx_of(ctx, name, regional)
        ! An init context for instance `name`: backend from the par file is not needed
        ! by component_identity; the grid has 4 x 3 cells at 50 km.
        type(init_ctx_t), intent(out) :: ctx
        character(len=*), intent(in)  :: name
        logical,          intent(in)  :: regional

        ctx%instance = name
        ctx%backend  = "x"
        ctx%domain   = "Dom"
        if (regional) then
            ctx%scope = "regional"
            ctx%grid_name = "TOY"
            call grid_init(ctx%grid, name="TOY", mtype="cartesian", units="kilometers", &
                           x0=25.0_dp, dx=50.0_dp, nx=NX, y0=25.0_dp, dy=50.0_dp, ny=NY)
        else
            ctx%scope = "global"
        end if
    end subroutine ctx_of

    subroutine make(c, name, kind, backend, regional)
        type(toy_t),      intent(out) :: c
        character(len=*), intent(in)  :: name, kind, backend
        logical,          intent(in)  :: regional
        type(init_ctx_t) :: ctx
        call ctx_of(ctx, name, regional)
        ctx%backend = backend
        call toy_make(c, ctx, kind)
    end subroutine make

    function times_of(file) result(t)
        character(len=*), intent(in) :: file
        real(dp), allocatable :: t(:)
        allocate(t(nc_size(file, "time")))
        call nc_read(file, "time", t)
    end function times_of

    function slab(file, name, k, kz) result(a)
        ! Record k of a 2D variable (layer kz of a rank-3 one: kz > 0).
        character(len=*), intent(in) :: file, name
        integer,          intent(in) :: k, kz
        real(sp) :: a(NX, NY)
        real(sp) :: b(NX, NY, 1)
        if (kz > 0) then
            call nc_read(file, name, b, start=[1, 1, kz, k], count=[NX, NY, 1, 1])
        else
            call nc_read(file, name, b, start=[1, 1, k], count=[NX, NY, 1])
        end if
        a = b(:, :, 1)
    end function slab

    function attr(file, name, att) result(s)
        character(len=*), intent(in) :: file, name, att
        character(len=:), allocatable :: s
        character(len=256) :: buf
        buf = ""
        call nc_read_attr(file, name, att, buf)
        s = trim(buf)
    end function attr

    function line_count(file) result(n)
        character(len=*), intent(in) :: file
        integer :: n, u, io
        character(len=1024) :: line
        n = 0
        open(newunit=u, file=file, status="old", action="read")
        do
            read(u, '(a)', iostat=io) line
            if (io /= 0) exit
            n = n + 1
        end do
        close(u)
    end function line_count

    logical function same_times(t, expect)
        real(dp), intent(in) :: t(:), expect(:)
        same_times = size(t) == size(expect)
        if (same_times) same_times = all(abs(t - expect) < 1.0e-9_dp)
    end function same_times

    ! -------------------------------------------------------------------------
    ! variable tables against contracts-design section 2
    ! -------------------------------------------------------------------------

    subroutine test_tables()
        integer :: stat
        type(var_io_type), allocatable :: tab(:)

        call expect("ice", [character(len=16) :: "H_ice", "z_srf", "z_bed", "z_sl", "f_ice", "f_grnd", "dHidt", &
                            "smb", "bmb", "cmb", "mb_net", "M_ice", "ux_s", "uy_s", "dz_bed_dt_corr", &
                            "z_bed_ref", "H_ice_ref", "footprint", "H_ice_pd", "H_grnd", "H_grnd_pd"], &
                    [character(len=16) :: "m", "m", "m", "m", "1", "1", "m yr-1", "kg m-2 yr-1", "kg m-2 yr-1", &
                            "kg m-2 yr-1", "kg m-2 yr-1", "Gt", "m yr-1", "m yr-1", "m yr-1", "m", "m", "1", "m", "m", "m"], &
                    [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2])
        call expect("geo", [character(len=16) :: "dz_bed", "dz_ss", "z_bed_ref", "dz_bed_dt", "bsl", "A_ocn", "dV_bsl", "dphi", &
                            "C_ocn", "bsl_ive", "mass_resid"], &
                    [character(len=16) :: "m", "m", "m", "m yr-1", "m", "m2", "m3", "m", "1", "m", "1"], &
                    [2, 2, 2, 2, 0, 0, 0, 0, 2, 0, 0])
        call expect("climate", [character(len=16) :: "tas", "pr", "z_srf_clim", "lapse", "beta_p", "to_ann", "so_ann", &
                            "dto_ann", "to_shlf", "so_shlf", "dto_shlf", "dso_shlf", "smb", "T_srf", "Qd", "tsl_ann"], &
                    [character(len=16) :: "K", "kg m-2 yr-1", "m", "K m-1", "K-1", "K", "psu", "K", "K", "psu", "K", &
                            "psu", "kg m-2 yr-1", "K", "m3 s-1", "K"], &
                    [3, 3, 2, 3, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2])
        call expect("surface", [character(len=16) :: "smb", "T_srf", "melt", "runoff", "refrz", "albedo"], &
                    [character(len=16) :: "kg m-2 yr-1", "K", "kg m-2 yr-1", "kg m-2 yr-1", "kg m-2 yr-1", "1"], &
                    [2, 2, 2, 2, 2, 2])
        call expect("shelf", [character(len=16) :: "bmb_shlf", "T_shlf", "tf_shlf"], &
                    [character(len=16) :: "kg m-2 yr-1", "K", "K"], [2, 2, 2])
        call expect("geothermal", [character(len=16) :: "Q_geo"], [character(len=16) :: "W m-2"], [2])
        call expect("sediments", [character(len=16) :: "H_sed"], [character(len=16) :: "m"], [2])
        call expect("forcing", [character(len=16) ::], [character(len=16) ::], [integer ::])

        call outfile_load("no_such_kind", tab, stat)
        call check(stat == 1, "tables: a missing table is an error")
    end subroutine test_tables

    subroutine expect(kind, names, units, ranks)
        ! The table of `kind` holds exactly the vocabulary: names, units, rank
        ! (dimensions), a long name and a precision.
        character(len=*),  intent(in) :: kind
        character(len=*),  intent(in) :: names(:), units(:)
        integer,           intent(in) :: ranks(:)

        type(var_io_type), allocatable :: tab(:)
        integer :: i, k, stat, nd
        logical :: ok_u, ok_r, ok_p, ok_l

        call outfile_load(kind, tab, stat)
        call check(stat == 0, "tables: "//kind//" loads")
        if (stat /= 0) return
        call check(size(tab) == size(names), "tables: "//kind//" has "//itoa(size(names))//" rows")
        ok_u = .true.; ok_r = .true.; ok_p = .true.; ok_l = .true.
        do i = 1, size(names)
            k = 0
            do nd = 1, size(tab)
                if (trim(tab(nd)%varname) == trim(names(i))) k = nd
            end do
            if (k == 0) then
                call check(.false., "tables: "//kind//" lists "//trim(names(i)))
                cycle
            end if
            ok_u = ok_u .and. trim(tab(k)%units) == trim(units(i))
            nd = tab(k)%ndims
            if (len_trim(tab(k)%dims(1)) == 0) nd = 0
            ok_r = ok_r .and. nd == ranks(i)
            if (nd >= 2) ok_r = ok_r .and. trim(tab(k)%dims(1)) == "xc" .and. trim(tab(k)%dims(2)) == "yc"
            if (nd >= 2) then
                ok_p = ok_p .and. (tab(k)%prec_abs > 0.0d0 .or. tab(k)%prec_rel > 0.0d0)
            else                                    ! 0-D: time series are lossless
                ok_p = ok_p .and. tab(k)%prec_abs == 0.0d0 .and. tab(k)%prec_rel == 0.0d0
            end if
            ok_l = ok_l .and. len_trim(tab(k)%long_name) > 0
        end do
        call check(ok_u, "tables: "//kind//" units")
        call check(ok_r, "tables: "//kind//" dimensions")
        call check(ok_p, "tables: "//kind//" every 2D, 3D row has a precision, every 0-D row is lossless")
        call check(ok_l, "tables: "//kind//" every row has a long_name")
    end subroutine expect

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: b
        write(b, '(i0)') i
        s = trim(b)
    end function itoa

    ! -------------------------------------------------------------------------
    ! a run: registration, dispatch, files
    ! -------------------------------------------------------------------------

    subroutine hub_writer(ctx, file, t, stream, init)
        class(*),         intent(inout) :: ctx
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        associate (unused => stream); end associate
        select type (ctx)
        type is (hub_ctx_t)
            ctx%ncall = ctx%ncall + 1
            ctx%t(ctx%ncall) = t
            ctx%init(ctx%ncall) = init
            ctx%file = file
        end select
    end subroutine hub_writer

    subroutine test_run()
        type(config_t) :: cfg
        type(clock_t)  :: clk
        type(output_t) :: out
        type(toy_t), target :: ice, sf, clim, geo
        type(hub_ctx_t), target :: hub
        type(timer_class), target :: tmr
        character(len=:), allocatable :: f
        character(len=512) :: files(32)
        real(dp), allocatable :: t(:)
        real(sp) :: a(NX, NY), b(NX, NY)
        real(dp) :: tt, tb(4) = [10.0_dp, 30.0_dp, 50.0_dp, 60.0_dp]
        integer :: stat, i, j, k, nstep
        logical :: ok, differs

        call config_init(cfg, par, stat=stat)
        call check(stat == 0, "run: config_init")
        call clock_init(clk, 0.0_dp, 60.0_dp, 10.0_dp, "365_day", "tropical")

        call make(ice,  "ice_toy",  "ice",     "yelmo",        .true.)
        call make(sf,   "sf_toy",   "surface", "chion",        .true.)
        call make(clim, "clim_toy", "climate", "snapshots",    .true.)
        call make(geo,  "geo_glob", "geo",     "fastisostasy", .false.)
        call check(allocated(ice%xkm) .and. size(ice%xkm) == NX .and. abs(ice%xkm(2) - 75.0_dp) < 1e-9_dp, &
                   "identity: grid axes [km] from the init context")
        call check(.not. allocated(geo%xkm), "identity: a global instance has no grid axes")

        call timer_init(tmr, [character(len=16) :: "output", "step"])
        call output_init(out, cfg)
        call output_add_component(out, cfg, clk, ice)
        call output_add_component(out, cfg, clk, sf)
        call output_add_component(out, cfg, clk, clim)
        call output_add_component(out, cfg, clk, geo)
        call output_add_hub(out, cfg, clk, "Dom", hub_writer, hub)
        call output_add_timing(out, cfg, clk, tmr)
        call output_describe(out)

        ! registration
        call check(size(out%sched) == 2, "registration: one schedule per distinct group (out_a, out_b)")
        call check(size(out%str) == 10, "registration: 10 streams")
        call check(size(out%prod) == 6, "registration: 6 producers (4 components, hub, timing)")
        call check(allocated(ice%out_vars) .and. size(ice%out_vars) == 3, "registration: out_vars of ice_toy")
        call check(.not. allocated(sf%out_vars) .or. size(sf%out_vars) == 0, "registration: out_vars unset = all")
        files(1:size(output_streams(out, "ice_toy"))) = output_streams(out, "ice_toy")
        call check(size(output_streams(out, "ice_toy")) == 4 .and. &
                   trim(files(1)) == rundir//"ice_toy.nc" .and. trim(files(2)) == rundir//"ice_toy_ts.nc" .and. &
                   trim(files(3)) == rundir//"ice_toy_yelmo.nc" .and. trim(files(4)) == rundir//"ice_toy_yelmo_ts.nc", &
                   "registration: file names <instance>[_ts|_<stream>].nc in out_dir")
        call check(size(output_streams(out, "sf_toy")) == 1 .and. size(output_streams(out, "geo_glob")) == 1, &
                   "registration: streams set to none are not registered")
        call check(trim(out%str(size(out%str))%stream) == "timing" .and. &
                   trim(out%str(size(out%str))%group) == "out_b", "registration: timing follows the first 1D group")

        ! the run
        call clock_start(clk)
        nstep = 0
        call publish_all(clock_time(clk), ice, sf, clim, geo)
        call output_write(out, clk)
        do while (.not. clock_finished(clk))
            call clock_advance(clk)
            call publish_all(clock_time(clk), ice, sf, clim, geo)
            call timer_step(tmr, 1)
            call output_write(out, clk)
            nstep = nstep + 1
        end do
        call check(nstep == 6, "run: six steps")

        ! files
        call check(exists(rundir//"ice_toy.nc") .and. exists(rundir//"ice_toy_ts.nc"), "files: ice_toy 2D and 1D")
        call check(exists(rundir//"sf_toy.nc") .and. .not. exists(rundir//"sf_toy_ts.nc"), "files: sf_toy 2D only")
        call check(exists(rundir//"clim_toy.nc") .and. exists(rundir//"clim_toy_ts.nc"), "files: clim_toy 2D and 1D")
        call check(exists(rundir//"geo_glob_ts.nc") .and. .not. exists(rundir//"geo_glob.nc"), "files: geo_glob 1D only")
        call check(.not. exists(rundir//"ice_toy_yelmo.nc") .and. .not. exists(rundir//"ice_toy_yelmo_ts.nc"), &
                   "files: backend streams are the backend's (toy writes none)")

        ! time records: out_a = 0 20 40 60; out_b = 10 30 50 60 (no time_init)
        t = times_of(rundir//"ice_toy.nc")
        call check(same_times(t, [0.0_dp, 20.0_dp, 40.0_dp, 60.0_dp]), "ice_toy.nc: records at out_a times")
        t = times_of(rundir//"ice_toy_ts.nc")
        call check(same_times(t, [10.0_dp, 30.0_dp, 50.0_dp, 60.0_dp]), "ice_toy_ts.nc: records at out_b times")
        t = times_of(rundir//"sf_toy.nc")
        call check(same_times(t, [10.0_dp, 30.0_dp, 50.0_dp, 60.0_dp]), "sf_toy.nc: records at out_b times")
        t = times_of(rundir//"geo_glob_ts.nc")
        call check(same_times(t, [0.0_dp, 20.0_dp, 40.0_dp, 60.0_dp]), "geo_glob_ts.nc: records at out_a times")

        ! ice_toy.nc: out_vars = H_ice smb (2D) / M_ice (1D)
        f = rundir//"ice_toy.nc"
        call check(nc_exists_var(f, "H_ice") .and. nc_exists_var(f, "smb"), "ice_toy.nc: H_ice, smb")
        call check(.not. nc_exists_var(f, "z_bed") .and. .not. nc_exists_var(f, "M_ice"), &
                   "ice_toy.nc: out_vars leaves out z_bed and the 0-D M_ice")
        call check(attr(f, "H_ice", "units") == "m" .and. attr(f, "smb", "units") == "kg m-2 yr-1", &
                   "ice_toy.nc: units attributes")
        call check(attr(f, "H_ice", "long_name") == "Ice thickness", "ice_toy.nc: long_name from the table")
        call check(nc_exists_var(f, "xc") .and. nc_exists_var(f, "yc") .and. attr(f, "xc", "units") == "km", &
                   "ice_toy.nc: xc, yc [km]")
        ok = .true.
        differs = .false.
        do k = 1, 4
            tt = 20.0_dp*real(k - 1, dp)
            a = slab(f, "H_ice", k, 0)
            b = slab(f, "smb", k, 0)
            do j = 1, NY
                do i = 1, NX
                    ok = ok .and. abs(real(a(i, j), dp) - toy_field(tt, i, j, 0)) < 2.0e-3_dp
                    ok = ok .and. abs(real(b(i, j), dp) - 3.0_dp*toy_field(tt, i, j, 0)) < 0.2_dp
                    differs = differs .or. a(i, j) /= real(toy_field(tt, i, j, 0), sp)
                end do
            end do
        end do
        call check(ok, "ice_toy.nc: values within the table precision at every record")
        call check(differs, "ice_toy.nc: values are rounded to the table precision")

        f = rundir//"ice_toy_ts.nc"
        call check(nc_exists_var(f, "M_ice") .and. .not. nc_exists_var(f, "H_ice") .and. .not. nc_exists_var(f, "xc"), &
                   "ice_toy_ts.nc: the 0-D export only, no grid axes")
        call check(attr(f, "M_ice", "units") == "Gt", "ice_toy_ts.nc: units")
        ok = .true.
        call nc_read(f, "M_ice", a(1:4, 1))
        do k = 1, 4
            tt = tb(k)
            ok = ok .and. abs(real(a(k, 1), dp) - (2.0_dp*tt + 5.0_dp)) < 1.0e-4_dp*(2.0_dp*tt + 5.0_dp)
        end do
        call check(ok, "ice_toy_ts.nc: M_ice values")

        ! surface: no out_vars = every export, MEAN exports published at the window end
        f = rundir//"sf_toy.nc"
        call check(nc_exists_var(f, "T_srf") .and. nc_exists_var(f, "smb"), "sf_toy.nc: all exports")
        a = slab(f, "T_srf", 2, 0)
        call check(abs(real(a(2, 3), dp) - toy_field(30.0_dp, 2, 3, 0)) < 2.0e-3_dp, "sf_toy.nc: MEAN export value at its record")

        ! climate: rank 3 with the month axis, and a 0-D export without a table row
        f = rundir//"clim_toy.nc"
        call check(nc_exists_var(f, "tas") .and. nc_exists_var(f, "z_srf_clim") .and. nc_exists_var(f, "month"), &
                   "clim_toy.nc: tas, z_srf_clim, month axis")
        call check(nc_size(f, "month") == 12 .and. nc_size(f, "time") == 4, "clim_toy.nc: 12 months, 4 records")
        a = slab(f, "tas", 3, 7)
        call check(abs(real(a(3, 2), dp) - toy_field(40.0_dp, 3, 2, 7)) < 2.0e-3_dp, "clim_toy.nc: tas(month 7) at record 3")
        f = rundir//"clim_toy_ts.nc"
        call check(nc_exists_var(f, "dTo") .and. attr(f, "dTo", "units") == "K", "clim_toy_ts.nc: derived 0-D export, units")
        call nc_read(f, "dTo", a(1:4, 1))
        call check(all(a(1:4, 1) == real([(toy_field(20.0_dp*real(k - 1, dp), 1, 1, 0), k = 1, 4)], sp)), &
                   "clim_toy_ts.nc: a 0-D export without a row is written lossless (as float)")

        ! global geo
        f = rundir//"geo_glob_ts.nc"
        call check(nc_exists_var(f, "bsl") .and. nc_exists_var(f, "A_ocn"), "geo_glob_ts.nc: bsl, A_ocn")
        call check(attr(f, "A_ocn", "units") == "m2" .and. attr(f, "bsl", "long_name") == "Barystatic sea level", &
                   "geo_glob_ts.nc: units, long_name")
        call nc_read(f, "A_ocn", a(1:4, 1))
        call check(all(a(1:4, 1) == real(3.6e14_dp, sp)), "geo_glob_ts.nc: a 0-D export is lossless in a rounded run (time series)")

        ! backend streams routed to the backend's write_output, init once
        call check(toy_calls(ice, "yelmo") == 4 .and. toy_calls(ice, "yelmo_ts") == 4, &
                   "backend streams: write_output called at each due time")
        call check(all(toy_call_init(ice, "yelmo") .eqv. [.true., .false., .false., .false.]), &
                   "backend streams: init only on the first call (yelmo)")
        call check(all(toy_call_init(ice, "yelmo_ts") .eqv. [.true., .false., .false., .false.]), &
                   "backend streams: init only on the first call (yelmo_ts, first write at t = 10)")
        call check(toy_calls(ice, "2D") == 4 .and. toy_calls(ice, "1D") == 4, "default streams also pass through write_output")
        call check(trim(ice%log_file(1)) == rundir//"ice_toy.nc" .or. trim(ice%log_file(1)) == rundir//"ice_toy_yelmo.nc", &
                   "backend streams: the file handed over is in out_dir")
        call check(abs(ice%log_t(1)) < 1e-12_dp .and. abs(maxval(ice%log_t(1:ice%ncall)) - 60.0_dp) < 1e-9_dp, &
                   "backend streams: times of the calls")

        ! hub writer
        call check(hub%ncall == 4 .and. same_times(hub%t(1:4), [10.0_dp, 30.0_dp, 50.0_dp, 60.0_dp]), &
                   "hub: writer called at the out_b times")
        call check(all(hub%init(1:4) .eqv. [.true., .false., .false., .false.]), "hub: init on the first call")
        call check(trim(hub%file) == rundir//"Dom_hub.nc", "hub: file <domain>_hub.nc")

        ! timing table: header + one line per due time
        call check(exists(rundir//"timing.txt") .and. line_count(rundir//"timing.txt") == 5, &
                   "timing: timing.txt with a header and 4 lines")

        ! stream counters
        ok = .true.
        do k = 1, size(out%str)
            ok = ok .and. out%str(k)%created .and. out%str(k)%nwrite == 4
        end do
        call check(ok, "dispatch: every stream written 4 times")
    end subroutine test_run

    subroutine publish_all(tnow, ice, sf, clim, geo)
        real(dp),    intent(in)    :: tnow
        type(toy_t), intent(inout) :: ice, sf, clim, geo
        call toy_publish(ice, tnow)
        call toy_publish(sf, tnow)
        call toy_publish(clim, tnow)
        call toy_publish(geo, tnow)
    end subroutine publish_all

    ! -------------------------------------------------------------------------
    ! registration errors
    ! -------------------------------------------------------------------------

    subroutine test_errors()
        type(config_t) :: cfg
        type(clock_t)  :: clk
        type(output_t) :: out
        type(toy_t), target :: ice, bad, nosched
        integer :: stat

        call config_init(cfg, par, stat=stat)
        call make(ice, "ice_toy", "ice", "yelmo", .true.)
        call make(bad, "ice_bad", "ice", "yelmo", .true.)
        call make(nosched, "ice_nosched", "ice", "yelmo", .true.)

        call clock_init(clk, 0.0_dp, 60.0_dp, 10.0_dp, "365_day", "tropical")
        call output_init(out, cfg)
        call output_add_component(out, cfg, clk, bad, stat)
        call check(stat == 1, "errors: out_vars naming a non-export")

        call clock_init(clk, 0.0_dp, 60.0_dp, 10.0_dp, "365_day", "tropical")
        call output_init(out, cfg)
        call output_add_component(out, cfg, clk, nosched, stat)
        call check(stat == 1, "errors: a schedule group missing from the par file")

        call clock_init(clk, 0.0_dp, 60.0_dp, 10.0_dp, "365_day", "tropical")
        call output_init(out, cfg)
        call output_add_component(out, cfg, clk, ice, stat)
        call check(stat == 0, "errors: first registration of ice_toy")
        call output_add_component(out, cfg, clk, ice, stat)
        call check(stat == 1, "errors: the same instance twice (one file, two streams)")
    end subroutine test_errors

    ! -------------------------------------------------------------------------
    ! [kryos] write_lossless
    ! -------------------------------------------------------------------------

    subroutine test_lossless()
        character(len=*), parameter :: par2 = dir//"par_lossless.nml", rd = "tmp/test_output/lossless/"
        type(config_t) :: cfg
        type(clock_t)  :: clk
        type(output_t) :: out
        type(toy_t), target :: ice
        real(sp) :: a(NX, NY)
        integer :: stat, i, j

        call write_par(par2, rd, lossless=.true.)
        call config_init(cfg, par2, stat=stat)
        call check(cfg%write_lossless, "lossless: [kryos] write_lossless read")
        call clock_init(clk, 0.0_dp, 60.0_dp, 10.0_dp, "365_day", "tropical")
        call make(ice, "ice_toy", "ice", "yelmo", .true.)
        call output_init(out, cfg)
        call output_add_component(out, cfg, clk, ice)
        call clock_start(clk)
        call toy_publish(ice, 0.0_dp)
        call output_write(out, clk)
        a = slab(rd//"ice_toy.nc", "H_ice", 1, 0)
        call check(all([((a(i, j) == real(toy_field(0.0_dp, i, j, 0), sp), i = 1, NX), j = 1, NY)]), &
                   "lossless: values equal the float of the export, table precision ignored")
    end subroutine test_lossless

    ! -------------------------------------------------------------------------
    ! output precision (phase 3, WP3)
    ! -------------------------------------------------------------------------

    subroutine test_precision()
        character(len=*), parameter :: pdir = dir//"prec/"
        character(len=16), parameter :: KINDS(8) = [character(len=16) :: "ice", "geo", "climate", "surface", "shelf", &
                                                    "geothermal", "sediments", "forcing"]
        character(len=24), parameter :: BACKENDS(5) = [character(len=24) :: "geo-fastisostasy", "climate-snapshots", &
                                                       "surface-chion", "shelf-marine_shelf", "hub"]
        integer :: i

        call execute_command_line("rm -rf "//pdir//" && mkdir -p "//pdir)
        call outfile_set_table_dir("input")
        call nc_set_lossless(.false.)

        do i = 1, size(KINDS)
            call table_roundtrip(pdir, trim(KINDS(i)), .false.)
        end do
        do i = 1, size(BACKENDS)
            call table_roundtrip(pdir, trim(BACKENDS(i)), .true.)
        end do
        call moved_literals(pdir)
        call unlisted_lossless(pdir)
        call lossless_switch(pdir)
        call restart_state(pdir)
        call test_stops(pdir)
    end subroutine test_precision

    function sample(nx, ny, nz, salt) result(a)
        ! Deterministic values of both signs over magnitudes 1e-2 .. 1e4 (and the odd zero).
        integer, intent(in) :: nx, ny, nz, salt
        real(sp) :: a(nx, ny, nz)
        integer :: i, j, k
        real(dp) :: u, mag
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    u = sin(12.9898_dp*real(i + salt, dp) + 78.233_dp*real(j, dp) + 37.719_dp*real(k, dp))*43758.5453_dp
                    u = u - floor(u)
                    mag = 10.0_dp**(real(mod(i + 3*j + k, 7), dp) - 2.0_dp)
                    a(i, j, k) = real((2.0_dp*u - 1.0_dp)*mag, sp)
                end do
            end do
        end do
        a(1, 1, 1) = 0.0_sp
    end function sample

    subroutine make_file(file, dim3, nz)
        ! A 2D file (xc, yc, time) with a third axis `dim3` of nz levels when given.
        character(len=*), intent(in) :: file, dim3
        integer,          intent(in) :: nz
        integer :: i
        call outfile_dims_2d(file, 50.0_dp*real([(i, i = 1, 20)], dp), 50.0_dp*real([(i, i = 1, 15)], dp), 0.0_dp)
        if (len_trim(dim3) > 0) call outfile_axis(file, trim(dim3), real([(i, i = 1, nz)], dp), "1")
    end subroutine make_file

    logical function within(dat, got, abs_p, rel_p)
        ! |got - dat| <= abs_p, or <= rel_p*|dat| (both set: the larger of the two bounds, as ncio).
        real(sp), intent(in) :: dat(:), got(:)
        real(dp), intent(in) :: abs_p, rel_p
        real(dp) :: err, tol
        integer :: i
        within = .true.
        do i = 1, size(dat)
            err = abs(real(got(i), dp) - real(dat(i), dp))
            tol = max(abs_p, rel_p*abs(real(dat(i), dp)))
            if (abs_p <= 0.0_dp .and. rel_p <= 0.0_dp) tol = 0.0_dp
            if (err > tol*1.0001_dp) within = .false.
        end do
    end function within

    subroutine table_roundtrip(pdir, table, backend)
        ! Every row of `table` (the writer's way: variable_io for an instance
        ! table, the backend helpers for a backend table) round-trips within the
        ! row's tolerance, exact when it has none and under nc_set_lossless; a
        ! rank-0 row (a time series) carries no precision.
        character(len=*), intent(in) :: pdir, table
        logical,          intent(in) :: backend

        type(var_io_type), allocatable :: tab(:)
        type(var_io_type) :: row
        character(len=:), allocatable :: file
        character(len=32) :: dim3
        real(sp), allocatable :: dat(:,:,:), got(:,:,:)
        integer :: i, ncid, n, rank, nz, nbad, nlossy
        logical :: ok_r, ok_l

        call outfile_load(table, tab)
        call check(size(tab) > 0 .or. trim(table) == "forcing", "precision: table "//table//" has rows")
        nbad = 0; nlossy = 0
        do i = 1, size(tab)
            row = tab(i)
            rank = row%ndims
            if (len_trim(row%dims(1)) == 0) rank = 0
            if (rank == 0) then                 ! time series stay lossless
                if (row%prec_abs > 0.0_dp .or. row%prec_rel > 0.0_dp) then
                    nbad = nbad + 1
                    write(*,'(a)') "  FAIL  precision: "//table//" 0-D row `"//trim(row%varname)//"` has a prec (time series are lossless)"
                end if
                cycle
            end if
            nz = 1
            dim3 = ""
            if (rank == 3) then
                nz = 3
                dim3 = row%dims(3)
            end if
            dat = sample(20, 15, nz, i)
            file = pdir//"rt_"//table//".nc"
            call make_file(file, trim(dim3), nz)
            call outfile_open(file, 0.0_dp, ncid, n)
            call write_one(file, row, dat, n, ncid, backend, table, dim3)
            call nc_close(ncid)
            allocate(got(20, 15, nz))
            call read_var(file, trim(row%varname), got, rank)
            ok_r = within(reshape(dat, [size(dat)]), reshape(got, [size(got)]), row%prec_abs, row%prec_rel)
            if (.not. ok_r) then
                nbad = nbad + 1
                write(*,'(a)') "  FAIL  precision: "//table//" `"//trim(row%varname)//"` outside its tolerance"
            end if
            if (row%prec_abs > 0.0_dp .or. row%prec_rel > 0.0_dp) then
                if (any(got /= dat)) nlossy = nlossy + 1
            else
                if (any(got /= dat)) then
                    nbad = nbad + 1
                    write(*,'(a)') "  FAIL  precision: "//table//" `"//trim(row%varname)//"` has no prec but is not exact"
                end if
            end if
            deallocate(got)

            ! the same under nc_set_lossless: exact
            call nc_set_lossless(.true.)
            call make_file(file, trim(dim3), nz)
            call outfile_open(file, 0.0_dp, ncid, n)
            call write_one(file, row, dat, n, ncid, backend, table, dim3)
            call nc_close(ncid)
            allocate(got(20, 15, nz))
            call read_var(file, trim(row%varname), got, rank)
            ok_l = all(got == dat)
            call nc_set_lossless(.false.)
            deallocate(got)
            if (.not. ok_l) then
                nbad = nbad + 1
                write(*,'(a)') "  FAIL  precision: "//table//" `"//trim(row%varname)//"` is not exact under write_lossless"
            end if
        end do
        call check(nbad == 0, "precision: "//table//": every row within its tolerance, exact without prec and under write_lossless, 0-D rows lossless")
        call check(nlossy > 0 .or. all(tab%prec_abs == 0.0_dp .and. tab%prec_rel == 0.0_dp), &
                   "precision: "//table//": the rounding rows do round the sample")
    end subroutine table_roundtrip

    subroutine write_one(file, row, dat, n, ncid, backend, table, dim3)
        character(len=*),  intent(in) :: file, table, dim3
        type(var_io_type), intent(in) :: row
        real(sp),          intent(in) :: dat(:,:,:)
        integer,           intent(in) :: n, ncid
        logical,           intent(in) :: backend

        if (backend) then
            if (len_trim(dim3) == 0) then
                call outfile_var2d(file, trim(row%varname), dat(:,:,1), n, ncid, trim(row%units), "long", table=table)
            else
                call outfile_var3d(file, trim(row%varname), dat, trim(dim3), n, ncid, trim(row%units), "long", table=table)
            end if
        else
            if (len_trim(dim3) == 0) then
                call var_io_write(file, row, dat(:,:,1), n, ncid=ncid)
            else
                call var_io_write(file, row, dat, n, ncid=ncid)
            end if
        end if
    end subroutine write_one

    subroutine read_var(file, name, got, rank)
        character(len=*), intent(in)  :: file, name
        real(sp),         intent(out) :: got(:,:,:)
        integer,          intent(in)  :: rank
        real(sp), allocatable :: b2(:,:,:), b3(:,:,:,:)
        if (rank == 3) then
            allocate(b3(size(got,1), size(got,2), size(got,3), 1))
            call nc_read(file, name, b3, start=[1, 1, 1, 1], count=[size(got,1), size(got,2), size(got,3), 1])
            got = b3(:,:,:,1)
        else
            allocate(b2(size(got,1), size(got,2), 1))
            call nc_read(file, name, b2, start=[1, 1, 1], count=[size(got,1), size(got,2), 1])
            got = b2
        end if
    end subroutine read_var

    subroutine moved_literals(pdir)
        ! The precisions that were hard-coded at the call sites (the literals below, as
        ! they were written there) and the table rows give the same file values.
        character(len=*), intent(in) :: pdir

        type :: lit_t
            character(len=24) :: table
            character(len=24) :: name
            real(sp) :: pa, pr
        end type lit_t
        type(lit_t), parameter :: LIT(28) = [ &
            lit_t("geo-fastisostasy", "z_bed", 1e-3_sp, 0.0_sp), lit_t("geo-fastisostasy", "z_ss", 1e-3_sp, 0.0_sp), &
            lit_t("geo-fastisostasy", "w", 1e-3_sp, 0.0_sp), lit_t("geo-fastisostasy", "we", 1e-3_sp, 0.0_sp), &
            lit_t("climate-snapshots", "t2m_ann", 1e-3_sp, 0.0_sp), lit_t("climate-snapshots", "t2m_sum", 1e-3_sp, 0.0_sp), &
            lit_t("climate-snapshots", "pr_ann", 0.0_sp, 1e-4_sp), lit_t("climate-snapshots", "to_ann", 1e-3_sp, 0.0_sp), &
            lit_t("climate-snapshots", "so_ann", 0.0_sp, 1e-4_sp), lit_t("climate-snapshots", "tas", 1e-3_sp, 0.0_sp), &
            lit_t("climate-snapshots", "pr", 0.0_sp, 1e-4_sp), lit_t("climate-snapshots", "tsl_ann", 1e-3_sp, 0.0_sp), &
            lit_t("climate-snapshots", "z_srf", 1e-3_sp, 0.0_sp), lit_t("climate-snapshots", "bias_dtas_ann_now", 1e-3_sp, 0.0_sp), &
            lit_t("surface-chion", "smb", 0.1_sp, 0.0_sp), lit_t("surface-chion", "tsrf", 1e-3_sp, 0.0_sp), &
            lit_t("surface-chion", "melt", 0.1_sp, 0.0_sp), lit_t("surface-chion", "runoff", 0.1_sp, 0.0_sp), &
            lit_t("surface-chion", "refrz", 0.1_sp, 0.0_sp), lit_t("surface-chion", "pr", 0.0_sp, 1e-4_sp), &
            lit_t("surface-chion", "sf", 0.0_sp, 1e-4_sp), lit_t("surface-chion", "t2m", 1e-3_sp, 0.0_sp), &
            lit_t("surface-chion", "S", 0.0_sp, 1e-4_sp), lit_t("surface-chion", "PDDs", 0.0_sp, 1e-4_sp), &
            lit_t("surface-chion", "alb_s", 1e-4_sp, 0.0_sp), lit_t("surface-chion", "H_snow", 0.1_sp, 0.0_sp), &
            lit_t("shelf-marine_shelf", "bmb_shlf", 1e-4_sp, 0.0_sp), lit_t("hub", "f_ice", 1e-4_sp, 0.0_sp)]

        type(var_io_type) :: row
        character(len=:), allocatable :: f1, f2
        real(sp), allocatable :: dat(:,:,:), g1(:,:,:), g2(:,:,:)
        integer :: i, ncid, n, nz, rank, nbad
        logical :: found
        character(len=32) :: dim3

        nbad = 0
        do i = 1, size(LIT)
            call outfile_row(trim(LIT(i)%table), trim(LIT(i)%name), row, found)
            if (.not. found) then
                nbad = nbad + 1
                write(*,'(a)') "  FAIL  moved literals: "//trim(LIT(i)%table)//" lists no `"//trim(LIT(i)%name)//"`"
                cycle
            end if
            rank = row%ndims
            nz = 1; dim3 = ""
            if (rank == 3) then
                nz = 3; dim3 = row%dims(3)
            end if
            dat = sample(20, 15, nz, 100 + i)
            f1 = pdir//"lit_table.nc";  f2 = pdir//"lit_literal.nc"
            ! the table
            call make_file(f1, trim(dim3), nz)
            call outfile_open(f1, 0.0_dp, ncid, n)
            call write_one(f1, row, dat, n, ncid, .true., trim(LIT(i)%table), dim3)
            call nc_close(ncid)
            ! the literal, as the call sites wrote it
            call make_file(f2, trim(dim3), nz)
            call outfile_open(f2, 0.0_dp, ncid, n)
            if (rank == 3) then
                if (LIT(i)%pa > 0.0_sp) then
                    call nc_write(f2, trim(row%varname), dat, dim1="xc", dim2="yc", dim3=trim(dim3), dim4="time", &
                                  start=[1,1,1,n], count=[20,15,nz,1], ncid=ncid, units=trim(row%units), &
                                  long_name="long", prec_abs=LIT(i)%pa)
                else
                    call nc_write(f2, trim(row%varname), dat, dim1="xc", dim2="yc", dim3=trim(dim3), dim4="time", &
                                  start=[1,1,1,n], count=[20,15,nz,1], ncid=ncid, units=trim(row%units), &
                                  long_name="long", prec_rel=LIT(i)%pr)
                end if
            else
                if (LIT(i)%pa > 0.0_sp) then
                    call outfile_var2d(f2, trim(row%varname), dat(:,:,1), n, ncid, trim(row%units), "long", prec_abs=LIT(i)%pa)
                else
                    call outfile_var2d(f2, trim(row%varname), dat(:,:,1), n, ncid, trim(row%units), "long", prec_rel=LIT(i)%pr)
                end if
            end if
            call nc_close(ncid)
            allocate(g1(20, 15, nz), g2(20, 15, nz))
            call read_var(f1, trim(row%varname), g1, rank)
            call read_var(f2, trim(row%varname), g2, rank)
            if (any(g1 /= g2)) then
                nbad = nbad + 1
                write(*,'(a)') "  FAIL  moved literals: "//trim(LIT(i)%table)//" `"//trim(LIT(i)%name)//"` differs from its literal"
            end if
            deallocate(g1, g2)
        end do
        call check(nbad == 0, "precision: the 28 moved literals and their table rows give identical rounded values")
    end subroutine moved_literals

    subroutine unlisted_lossless(pdir)
        ! A backend variable the table does not list is written lossless, with the long_name given.
        character(len=*), intent(in) :: pdir
        character(len=:), allocatable :: f
        real(sp), allocatable :: dat(:,:,:), got(:,:,:)
        integer :: ncid, n

        f = pdir//"unlisted.nc"
        dat = sample(20, 15, 1, 7)
        call make_file(f, "", 0)
        call outfile_open(f, 0.0_dp, ncid, n)
        call outfile_var2d(f, "not_in_the_table", dat(:,:,1), n, ncid, "m", "An unlisted variable", table="geo-fastisostasy")
        call outfile_var2d(f, "w", dat(:,:,1), n, ncid, "m", "Not the table's long name", table="geo-fastisostasy")
        call nc_close(ncid)
        allocate(got(20, 15, 1))
        call read_var(f, "not_in_the_table", got, 2)
        call check(all(got == dat), "precision: a backend variable the table does not list is lossless")
        call check(attr(f, "not_in_the_table", "long_name") == "An unlisted variable", &
                   "precision: an unlisted variable keeps the long_name of the writer")
        call check(attr(f, "w", "long_name") == "Viscous displacement", "precision: a listed variable takes the table's long_name")
        call read_var(f, "w", got, 2)
        call check(any(got /= dat) .and. within(reshape(dat, [size(dat)]), reshape(got, [size(got)]), 1.0d-3, 0.0d0), &
                   "precision: a listed variable is rounded to its row")
    end subroutine unlisted_lossless

    subroutine lossless_switch(pdir)
        ! write_lossless (nc_set_lossless): exact values from a rounding row of the default
        ! writer's table and of a backend table; a rounded 0-D time series row is exact too.
        character(len=*), intent(in) :: pdir
        character(len=:), allocatable :: f
        real(sp), allocatable :: dat(:,:,:), got(:,:,:)
        type(var_io_type) :: row
        logical :: found
        integer :: ncid, n

        dat = sample(20, 15, 1, 11)
        call nc_set_lossless(.true.)
        f = pdir//"lossless_default.nc"
        call outfile_row("ice", "H_ice", row, found)
        call make_file(f, "", 0)
        call outfile_open(f, 0.0_dp, ncid, n)
        call var_io_write(f, row, dat(:,:,1), n, ncid=ncid)
        call outfile_var2d(f, "z_bed", dat(:,:,1), n, ncid, "m", "x", table="geo-fastisostasy")
        call nc_close(ncid)
        allocate(got(20, 15, 1))
        call read_var(f, "H_ice", got, 2)
        call check(found .and. all(got == dat), "precision: write_lossless: a rounding row of an instance table is exact")
        call read_var(f, "z_bed", got, 2)
        call check(all(got == dat), "precision: write_lossless: a rounding row of a backend table is exact")
        call nc_set_lossless(.false.)
    end subroutine lossless_switch

    subroutine restart_state(pdir)
        ! The restart state of a rounded run is a lossless run's: the exports' state file
        ! (fesm_fields restart, double precision) is the same file whatever the switch.
        character(len=*), intent(in) :: pdir
        type(toy_t), target :: ice_r, ice_l
        character(len=*), parameter :: nm = "ice_toy"
        integer :: es

        call make(ice_r, nm, "ice", "yelmo", .true.)
        call make(ice_l, nm, "ice", "yelmo", .true.)
        call toy_publish(ice_r, 30.0_dp)
        call toy_publish(ice_l, 30.0_dp)
        call nc_set_lossless(.false.)
        call fields_restart_write(ice_r%exports, pdir//"state_rounded.nc")
        call nc_set_lossless(.true.)
        call fields_restart_write(ice_l%exports, pdir//"state_lossless.nc")
        call nc_set_lossless(.false.)
        call execute_command_line("cmp -s "//pdir//"state_rounded.nc "//pdir//"state_lossless.nc", exitstat=es)
        call check(es == 0, "precision: the restart state of a rounded run equals that of a lossless run")
    end subroutine restart_state

    ! -------------------------------------------------------------------------
    ! stops: a bad prec entry, a row without a unit, a unit that differs from the data's
    ! -------------------------------------------------------------------------

    subroutine test_stops(pdir)
        character(len=*), intent(in) :: pdir
        character(len=*), parameter :: hdr = "| id | variable | dimensions | units | long_name | prec |"
        integer :: u

        call execute_command_line("mkdir -p "//pdir//"tabs")
        open(newunit=u, file=pdir//"tabs/kryos-variables-badprec.md", status="replace", action="write")
        write(u,'(a)') "# bad prec"
        write(u,'(a)') hdr
        write(u,'(a)') "|----|----------|------------|-------|-----------|------|"
        write(u,'(a)') "|  1 | a        | xc, yc     | m     | A         | 1e-3x |"
        close(u)
        open(newunit=u, file=pdir//"tabs/kryos-variables-nounit.md", status="replace", action="write")
        write(u,'(a)') "# no unit"
        write(u,'(a)') hdr
        write(u,'(a)') "|----|----------|------------|-------|-----------|------|"
        write(u,'(a)') "|  1 | a        | xc, yc     | m     | A         | 1e-3 |"
        write(u,'(a)') "|  2 | b        | xc, yc     |       | B         | 1e-3 |"
        close(u)
        open(newunit=u, file=pdir//"tabs/kryos-variables-unitclash.md", status="replace", action="write")
        write(u,'(a)') "# unit clash"
        write(u,'(a)') hdr
        write(u,'(a)') "|----|----------|------------|-------|-----------|------|"
        write(u,'(a)') "|  1 | a        | xc, yc     | m     | A         | 1e-3 |"
        close(u)

        call expect_stop("bad_prec", "bad prec entry", .true.)
        call expect_stop("no_unit", "has no unit", .false.)
        call expect_stop("unit_clash", "differ from the data", .false.)
        call expect_stop("no_table", "variable table not found", .false.)
    end subroutine test_stops

    subroutine expect_stop(case_name, text, plain_stop)
        ! Run this program as a child on one case: it must print `text` on stderr and not
        ! reach its last line. plain_stop: the stop is variable_io's (a bare `stop`, exit
        ! status 0); otherwise kryos' fesm_error (exit status 1).
        character(len=*), intent(in) :: case_name, text
        logical,          intent(in) :: plain_stop

        character(len=512) :: exe, line
        character(len=:), allocatable :: err, out
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > "//dir//"child.out 2> "//dir//"child.err", exitstat=es)
        err = ""; out = ""
        open(newunit=u, file=dir//"child.err", status="old", action="read", iostat=ios)
        do while (ios == 0)
            read(u, "(a)", iostat=ios) line
            if (ios == 0) err = err//trim(line)//" "
        end do
        close(u)
        open(newunit=u, file=dir//"child.out", status="old", action="read", iostat=ios)
        do while (ios == 0)
            read(u, "(a)", iostat=ios) line
            if (ios == 0) out = out//trim(line)//" "
        end do
        close(u)
        call check(index(err, text) > 0 .and. index(out, "SURVIVED") == 0 .and. (plain_stop .or. es == 1), &
                   "stops: "//case_name//" ('"//text//"')")
    end subroutine expect_stop

    subroutine abort_case(name)
        ! Child: the situation of one case, without `stat`.
        character(len=*), intent(in) :: name
        character(len=*), parameter :: pdir = dir//"prec/"
        type(var_io_type), allocatable :: tab(:)
        real(sp) :: a(NX, NY)
        integer :: ncid, n

        call outfile_set_table_dir(pdir//"tabs")
        select case (name)
        case ("bad_prec")
            call outfile_load("badprec", tab)
        case ("no_unit")
            call outfile_load("nounit", tab)
        case ("no_table")
            call outfile_load("no_such_backend_table", tab)
        case ("unit_clash")                     ! the data are in K, the row says m
            a = 1.0_sp
            call execute_command_line("mkdir -p "//pdir)
            call outfile_dims_2d(pdir//"clash.nc", 50.0_dp*real([(n, n = 1, NX)], dp), 50.0_dp*real([(n, n = 1, NY)], dp), 0.0_dp)
            call outfile_open(pdir//"clash.nc", 0.0_dp, ncid, n)
            call outfile_var2d(pdir//"clash.nc", "a", a, n, ncid, "K", "A", table="unitclash")
            call nc_close(ncid)
        end select
        write(*,'(a)') "SURVIVED"
    end subroutine abort_case

end program test_output
