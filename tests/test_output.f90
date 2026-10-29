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
    ! vocabulary of contracts-design section 2.
    !
    ! Files go to tmp/test_output/ (run from the repo root, as make check does).
    ! Abort paths take `stat` (fesm_error returns instead of stopping).

    use precision,        only: sp, dp
    use coordinates,      only: grid_init
    use ncio,             only: nc_read, nc_size, nc_exists_var, nc_read_attr
    use variable_io,      only: var_io_type
    use timer,            only: timer_class, timer_init, timer_step
    use fesm_clock,       only: clock_t, clock_init, clock_start, clock_advance, clock_finished, clock_time
    use kryos_component,  only: init_ctx_t
    use kryos_config,     only: config_t, config_init
    use kryos_outfile,    only: outfile_load, outfile_set_table_dir, outfile_path
    use kryos_output
    use test_output_toys

    implicit none

    character(len=*), parameter :: dir = "tmp/test_output/", rundir = "tmp/test_output/run/"
    character(len=*), parameter :: par = dir//"par.nml"
    integer :: nfail = 0

    type :: hub_ctx_t                           ! context of the toy hub writer
        integer :: ncall = 0
        real(dp) :: t(16) = 0.0_dp
        logical :: init(16) = .false.
        character(len=256) :: file = ""
    end type hub_ctx_t

    call execute_command_line("rm -rf "//dir//" && mkdir -p "//dir)
    call write_par(par, rundir, lossless=.false.)

    call test_tables()
    call test_run()
    call test_errors()
    call test_lossless()

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
        call expect("geo", [character(len=16) :: "dz_bed", "dz_ss", "z_bed_ref", "dz_bed_dt", "bsl", "A_ocn", "dV_bsl", "dphi"], &
                    [character(len=16) :: "m", "m", "m", "m yr-1", "m", "m2", "m3", "m"], [2, 2, 2, 2, 0, 0, 0, 0])
        call expect("climate", [character(len=16) :: "tas", "pr", "z_srf_clim", "lapse", "beta_p", "to_ann", "so_ann", &
                            "dto_ann", "to_shlf", "so_shlf", "dto_shlf", "dso_shlf", "smb", "T_srf", "Qd"], &
                    [character(len=16) :: "K", "kg m-2 yr-1", "m", "K m-1", "K-1", "K", "psu", "K", "K", "psu", "K", &
                            "psu", "kg m-2 yr-1", "K", "m3 s-1"], &
                    [3, 3, 2, 3, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2])
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
            ok_p = ok_p .and. (tab(k)%prec_abs > 0.0d0 .or. tab(k)%prec_rel > 0.0d0)
            ok_l = ok_l .and. len_trim(tab(k)%long_name) > 0
        end do
        call check(ok_u, "tables: "//kind//" units")
        call check(ok_r, "tables: "//kind//" dimensions")
        call check(ok_p, "tables: "//kind//" every row has a precision")
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

end program test_output
