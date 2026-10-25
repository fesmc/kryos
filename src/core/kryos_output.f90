module kryos_output
    ! Output streams and their dispatch (core-design section 6). A stream is
    ! (producer, file, schedule): the producers are component instances (any
    ! class(component_t), regional or global), writer procedures (the domain's
    ! hub, the run's timing table) and, in C9, whatever else registers one.
    !
    ! Registration (before clock_start: schedules register with the clock):
    !
    !     call output_init(out, cfg)                         ! out_dir, lossless switch
    !     call output_add_component(out, cfg, clk, comp)     ! after <comp>_init: exports declared
    !     call output_add_hub(out, cfg, clk, domain, writer, ctx)
    !     call output_add_timing(out, cfg, clk, tmr)         ! after the components
    !
    ! A component registers one stream per `out_<stream>` key of its group that
    ! names a schedule group (not "none"): the keys of the [instance] template
    ! (out_2D, out_1D; out_vars lists the exports of those two) and of the
    ! backend's template [<kind>_<backend>] (yelmo: out_yelmo, out_yelmo_sm,
    ! out_yelmo_ts). The file is <out_dir>/<instance>.nc (2D),
    ! <instance>_ts.nc (1D) or <instance>_<stream>.nc. Every distinct schedule
    ! group of the run is one schedule_t, initialised once.
    !
    ! Dispatch, at the end of run_init and of run_step:
    !
    !     call output_write(out, clk)
    !
    ! calls schedule_due once per distinct schedule (it moves a cursor), then
    ! write_output(file, t, stream, init) of every stream whose schedule is due,
    ! in registration order; init is .true. on a stream's first write of the run,
    ! which is when its file is created (each run or chunk writes anew).
    !
    ! Plug-in point for the producers that are not components (hub, globe): a
    ! writer procedure of the interface output_writer_i and a context object it
    ! receives back, registered with output_add_writer (output_add_hub is the
    ! domain's form). The context is held by pointer: keep it at a fixed address
    ! (a target that is not reallocated) from registration to the last write.
    ! Components are held by pointer too (a slot's `c`, which needs the target
    ! attribute on the slot).

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_clock,     only: clock_t, schedule_t, schedule_init, schedule_due, clock_time
    use ncio,            only: nc_set_lossless
    use nml,             only: nml_group_keys, nml_key_len
    use timer,           only: timer_class, timer_write_totals
    use fesm_fields,     only: fields_has, FIELD_NAME_LEN
    use kryos_component, only: component_t, component_output, NAME_LEN, PATH_LEN
    use kryos_config,    only: config_t, config_read, config_validate
    use kryos_outfile,   only: outfile_path

    implicit none

    private

    ! The streams of a component the [instance] template declares. Their key
    ! prefix is also how backend streams are found in the backend's template.
    character(len=*), parameter :: KEY_PREFIX = "out_", KEY_VARS = "out_vars"
    integer, parameter :: MAX_VARS = 64

    abstract interface
        subroutine output_writer_i(ctx, file, t, stream, init)
            ! A writer producer: write the record at time t of `stream` to `file`;
            ! init = .true. on the first write of the run (create the file).
            import :: dp
            class(*),         intent(inout) :: ctx
            character(len=*), intent(in)    :: file, stream
            real(dp),         intent(in)    :: t
            logical,          intent(in)    :: init
        end subroutine output_writer_i
    end interface

    type, public :: producer_t
        character(len=NAME_LEN) :: name = ""               ! instance, domain or "run"
        class(component_t), pointer :: comp => null()      ! a component producer, or
        procedure(output_writer_i), pointer, nopass :: write => null()   ! a writer producer
        class(*), pointer :: ctx => null()                 ! passed back to the writer
    end type producer_t

    type, public :: stream_t
        character(len=NAME_LEN) :: producer = "", stream = "", group = ""
        character(len=PATH_LEN) :: file = ""               ! as handed to write_output
        integer :: ip = 0, is = 0                          ! producer and schedule (output_t)
        logical :: created = .false.                       ! first write done
        integer :: nwrite = 0
    end type stream_t

    type, public :: output_t
        character(len=PATH_LEN) :: out_dir = "./"
        type(producer_t), allocatable :: prod(:)
        type(stream_t),   allocatable :: str(:)
        type(schedule_t), allocatable :: sched(:)          ! one per distinct schedule group
        character(len=NAME_LEN), allocatable :: group(:)   ! its name, parallel to sched
    end type output_t

    public :: output_writer_i
    public :: output_init, output_add_component, output_add_writer, output_add_hub, output_add_timing
    public :: output_write, output_describe, output_streams

contains

    ! =========================================================================
    ! Registration
    ! =========================================================================

    subroutine output_init(out, cfg)
        ! An empty output: the folder [kryos] out_dir (created) and the lossless
        ! switch of ncio ([kryos] write_lossless, global for every later write).
        type(output_t), intent(out) :: out
        type(config_t), intent(in)  :: cfg

        out%out_dir = cfg%out_dir
        allocate(out%prod(0), out%str(0), out%sched(0), out%group(0))
        if (len_trim(out%out_dir) > 0 .and. trim(out%out_dir) /= "." .and. trim(out%out_dir) /= "./") then
            call execute_command_line("mkdir -p '"//trim(out%out_dir)//"'")
        end if
        call nc_set_lossless(cfg%write_lossless)
    end subroutine output_init

    subroutine output_add_component(out, cfg, clk, comp, stat)
        ! Register the streams of instance `comp` (group = comp%name, template
        ! [instance] + [<kind>_<backend>]) and set comp%out_vars from `out_vars`
        ! (an unknown export stops). Call after the component's init, so its
        ! exports are declared, and before clock_start.
        type(output_t),     intent(inout)         :: out
        type(config_t),     intent(in)            :: cfg
        type(clock_t),      intent(inout)         :: clk
        class(component_t), intent(inout), target :: comp
        integer,            intent(out),   optional :: stat

        character(len=*), parameter :: proc = "kryos_output::output_add_component"
        character(len=NAME_LEN) :: group, extra, grp
        character(len=FIELD_NAME_LEN) :: vars(MAX_VARS)
        character(len=NAME_LEN), allocatable :: streams(:)
        integer :: i, ip

        if (present(stat)) stat = 0
        group = comp%name
        extra = trim(comp%kind)//"_"//trim(comp%backend)

        ! out_vars
        vars = ""
        call config_read(cfg, group, KEY_VARS, vars, "instance", extra)
        if (allocated(comp%out_vars)) deallocate(comp%out_vars)
        comp%out_vars = pack(vars, len_trim(vars) > 0)
        do i = 1, size(comp%out_vars)
            if (.not. fields_has(comp%exports, comp%out_vars(i))) then
                call fesm_error(proc, "out_vars of &"//trim(group)//" lists `"//trim(comp%out_vars(i))// &
                                "`, which is not an export of the instance", &
                                "exports: "//export_names(comp)//new_line("a")//"par file "//trim(cfg%par_file), stat)
                return
            end if
        end do

        ! streams: out_<stream> keys of the two templates
        streams = stream_names(cfg, extra)
        ip = 0
        do i = 1, size(streams)
            grp = ""
            call config_read(cfg, group, KEY_PREFIX//trim(streams(i)), grp, "instance", extra)
            if (is_none(grp)) cycle
            if (ip == 0) then
                call add_producer(out, group, comp=comp)
                ip = size(out%prod)
            end if
            call add_stream(out, cfg, clk, ip, trim(streams(i)), file_base(group, streams(i)), grp, stat)
            if (present(stat)) then
                if (stat /= 0) return
            end if
        end do
    end subroutine output_add_component

    subroutine output_add_writer(out, cfg, clk, name, stream, file_base, group, writer, ctx, stat)
        ! Register a stream of a producer that is not a component: `writer` is
        ! called with `ctx` and the file <out_dir>/<file_base> at the times of the
        ! schedule group `group` (an entry "none" registers nothing).
        type(output_t),   intent(inout)         :: out
        type(config_t),   intent(in)            :: cfg
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: name, stream, file_base, group
        procedure(output_writer_i)              :: writer
        class(*),         intent(inout), target :: ctx
        integer,          intent(out),   optional :: stat

        integer :: ip, k

        if (present(stat)) stat = 0
        if (is_none(group)) return

        ip = 0
        do k = 1, size(out%prod)         ! one producer entry per (name, ctx) pair
            if (trim(out%prod(k)%name) == trim(name) .and. associated(out%prod(k)%ctx)) then
                if (associated(out%prod(k)%ctx, ctx)) ip = k
            end if
        end do
        if (ip == 0) then
            call add_producer(out, name, writer=writer, ctx=ctx)
            ip = size(out%prod)
        end if
        call add_stream(out, cfg, clk, ip, stream, file_base, group, stat)
    end subroutine output_add_writer

    subroutine output_add_hub(out, cfg, clk, domain, writer, ctx, stat)
        ! The domain's hub: stream "hub", file <domain>_hub.nc, schedule group
        ! `out_2D` of the domain group (template [domain]).
        type(output_t),   intent(inout)         :: out
        type(config_t),   intent(in)            :: cfg
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: domain
        procedure(output_writer_i)              :: writer
        class(*),         intent(inout), target :: ctx
        integer,          intent(out),   optional :: stat

        character(len=NAME_LEN) :: grp

        grp = ""
        call config_read(cfg, domain, KEY_PREFIX//"2D", grp, "domain")
        call output_add_writer(out, cfg, clk, domain, "hub", trim(domain)//"_hub.nc", grp, writer, ctx, stat)
    end subroutine output_add_hub

    subroutine output_add_timing(out, cfg, clk, tmr, group, stat)
        ! The run's timing table (yelmox timing.txt): stream "timing", file
        ! timing.txt, written with timer_write_totals (seconds) when [kryos] timing
        ! is true. Schedule: `group`, by default that of the first registered 1D
        ! stream, so register it after the components; no 1D stream, no table.
        type(output_t),   intent(inout)         :: out
        type(config_t),   intent(in)            :: cfg
        type(clock_t),    intent(inout)         :: clk
        type(timer_class), intent(inout), target :: tmr
        character(len=*), intent(in),    optional :: group
        integer,          intent(out),   optional :: stat

        character(len=NAME_LEN) :: grp
        integer :: k

        if (present(stat)) stat = 0
        if (.not. cfg%timing) return
        grp = "none"
        if (present(group)) then
            grp = group
        else
            do k = 1, size(out%str)
                if (trim(out%str(k)%stream) == "1D") then
                    grp = out%str(k)%group
                    exit
                end if
            end do
        end if
        if (is_none(grp)) then
            call log_msg(LOG_INFO, "kryos_output", "timing table not written: no 1D stream to follow")
            return
        end if
        call output_add_writer(out, cfg, clk, "run", "timing", "timing.txt", grp, timing_writer, tmr, stat)
    end subroutine output_add_timing

    ! =========================================================================
    ! Dispatch
    ! =========================================================================

    subroutine output_write(out, clk)
        ! The output of the current clock time: every distinct schedule asked once,
        ! then the due streams written in registration order.
        type(output_t), intent(inout) :: out
        type(clock_t),  intent(in)    :: clk

        logical, allocatable :: due(:)
        real(dp) :: t
        integer :: i, k

        t = clock_time(clk)
        allocate(due(size(out%sched)))
        do i = 1, size(out%sched)
            due(i) = schedule_due(out%sched(i), clk)
        end do

        do k = 1, size(out%str)
            if (.not. due(out%str(k)%is)) cycle
            call write_stream(out, k, t)
        end do
    end subroutine output_write

    subroutine write_stream(out, k, t)
        type(output_t), intent(inout) :: out
        integer,        intent(in)    :: k
        real(dp),       intent(in)    :: t

        logical :: init
        integer :: ip

        init = .not. out%str(k)%created
        ip   = out%str(k)%ip
        if (associated(out%prod(ip)%comp)) then
            call component_output(out%prod(ip)%comp, trim(out%str(k)%file), t, trim(out%str(k)%stream), init)
        else
            call out%prod(ip)%write(out%prod(ip)%ctx, trim(out%str(k)%file), t, trim(out%str(k)%stream), init)
        end if
        out%str(k)%created = .true.
        out%str(k)%nwrite  = out%str(k)%nwrite + 1
    end subroutine write_stream

    ! =========================================================================
    ! Describe
    ! =========================================================================

    subroutine output_describe(out)
        ! One INFO line per stream: producer, stream, file and schedule group.
        type(output_t), intent(in) :: out
        integer :: k

        call log_msg(LOG_INFO, "kryos_output", "output: "//to_str(size(out%str))//" streams, "// &
                     to_str(size(out%sched))//" schedules, folder "//trim(out%out_dir))
        do k = 1, size(out%str)
            call log_msg(LOG_INFO, "kryos_output", "  "//trim(out%str(k)%producer)//" / "// &
                         trim(out%str(k)%stream)//" -> "//trim(out%str(k)%file)//" on &"//trim(out%str(k)%group))
        end do
    end subroutine output_describe

    function output_streams(out, producer) result(files)
        ! The files of the streams of `producer` (all when absent), in registration order.
        type(output_t),   intent(in)           :: out
        character(len=*), intent(in), optional :: producer
        character(len=PATH_LEN), allocatable :: files(:)
        integer :: k

        files = [character(len=PATH_LEN) ::]
        do k = 1, size(out%str)
            if (present(producer)) then
                if (trim(out%str(k)%producer) /= trim(producer)) cycle
            end if
            files = [files, out%str(k)%file]
        end do
    end function output_streams

    ! =========================================================================
    ! Internals
    ! =========================================================================

    subroutine add_producer(out, name, comp, writer, ctx)
        type(output_t),     intent(inout)                  :: out
        character(len=*),   intent(in)                     :: name
        class(component_t), intent(inout), target, optional :: comp
        procedure(output_writer_i),                optional :: writer
        class(*),           intent(inout), target, optional :: ctx
        type(producer_t), allocatable :: tmp(:)
        integer :: n

        n = size(out%prod)
        allocate(tmp(n + 1))
        tmp(1:n) = out%prod
        tmp(n + 1)%name = name
        if (present(comp)) tmp(n + 1)%comp => comp
        if (present(writer)) tmp(n + 1)%write => writer
        if (present(ctx)) tmp(n + 1)%ctx => ctx
        call move_alloc(tmp, out%prod)
    end subroutine add_producer

    subroutine add_stream(out, cfg, clk, ip, stream, base, group, stat)
        ! One more stream of producer ip on schedule group `group`, whose
        ! schedule is initialised the first time a stream names it.
        type(output_t),   intent(inout)         :: out
        type(config_t),   intent(in)            :: cfg
        type(clock_t),    intent(inout)         :: clk
        integer,          intent(in)            :: ip
        character(len=*), intent(in)            :: stream, base, group
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_output::add_stream"
        type(stream_t) :: s
        type(stream_t), allocatable :: tmp(:)
        integer :: k, is, n

        s%producer = out%prod(ip)%name
        s%stream   = stream
        s%group    = group
        s%file     = outfile_path(out%out_dir, base)
        s%ip       = ip

        do k = 1, size(out%str)
            if (trim(out%str(k)%file) == trim(s%file)) then
                call fesm_error(proc, "file "//trim(s%file)//" is the target of two streams", &
                                trim(out%str(k)%producer)//" / "//trim(out%str(k)%stream)//" and "// &
                                trim(s%producer)//" / "//trim(s%stream), stat)
                return
            end if
        end do

        is = schedule_index(out, group)
        if (is == 0) then
            call schedule_new(out, cfg, clk, group, stat)
            if (present(stat)) then
                if (stat /= 0) return
            end if
            is = size(out%sched)
        end if
        s%is = is

        n = size(out%str)
        allocate(tmp(n + 1))
        tmp(1:n) = out%str
        tmp(n + 1) = s
        call move_alloc(tmp, out%str)
    end subroutine add_stream

    integer function schedule_index(out, group) result(k)
        type(output_t),   intent(in) :: out
        character(len=*), intent(in) :: group
        integer :: i
        k = 0
        do i = 1, size(out%group)
            if (trim(out%group(i)) == trim(group)) then
                k = i
                return
            end if
        end do
    end function schedule_index

    subroutine schedule_new(out, cfg, clk, group, stat)
        ! The schedule of group `group` (in the par file, template [schedule]).
        type(output_t),   intent(inout)         :: out
        type(config_t),   intent(in)            :: cfg
        type(clock_t),    intent(inout)         :: clk
        character(len=*), intent(in)            :: group
        integer,          intent(out), optional :: stat

        type(schedule_t), allocatable :: tmp(:)
        character(len=NAME_LEN), allocatable :: gtmp(:)
        integer :: n

        call config_validate(cfg, group, "schedule", stat=stat)
        if (present(stat)) then
            if (stat /= 0) return
        end if

        n = size(out%sched)
        allocate(tmp(n + 1), gtmp(n + 1))
        tmp(1:n)  = out%sched
        gtmp(1:n) = out%group
        call schedule_init(tmp(n + 1), clk, group, trim(cfg%par_file), group, stat=stat)
        if (present(stat)) then
            if (stat /= 0) return
        end if
        gtmp(n + 1) = group
        call move_alloc(tmp, out%sched)
        call move_alloc(gtmp, out%group)
    end subroutine schedule_new

    function stream_names(cfg, extra) result(names)
        ! The stream names a component can have: `out_<stream>` keys of the
        ! [instance] template and of its backend template, 2D and 1D first.
        type(config_t),   intent(in) :: cfg
        character(len=*), intent(in) :: extra
        character(len=NAME_LEN), allocatable :: names(:)

        character(len=nml_key_len), allocatable :: k1(:), k2(:), keys(:)
        integer :: i, n

        k1 = nml_group_keys(trim(cfg%defaults_file), "instance")
        k2 = nml_group_keys(trim(cfg%defaults_file), trim(extra))
        keys = [k1, k2]
        names = [character(len=NAME_LEN) ::]
        do i = 1, size(keys)
            n = len(KEY_PREFIX)
            if (len_trim(keys(i)) <= n) cycle
            if (keys(i)(1:n) /= KEY_PREFIX) cycle
            if (trim(keys(i)) == KEY_VARS) cycle
            if (any(names == keys(i)(n+1:))) cycle
            names = [names, keys(i)(n+1:)]
        end do
    end function stream_names

    function file_base(instance, stream) result(base)
        ! <instance>.nc for 2D, <instance>_ts.nc for 1D, else <instance>_<stream>.nc.
        character(len=*), intent(in) :: instance, stream
        character(len=:), allocatable :: base
        select case (trim(stream))
        case ("2D")
            base = trim(instance)//".nc"
        case ("1D")
            base = trim(instance)//"_ts.nc"
        case default
            base = trim(instance)//"_"//trim(stream)//".nc"
        end select
    end function file_base

    logical function is_none(group)
        character(len=*), intent(in) :: group
        is_none = len_trim(group) == 0 .or. trim(group) == "none"
    end function is_none

    function export_names(c) result(s)
        class(component_t), intent(in) :: c
        character(len=:), allocatable :: s
        integer :: h
        s = ""
        do h = 1, c%exports%n
            s = s//" "//trim(c%exports%f(h)%name)
        end do
        s = trim(adjustl(s))
    end function export_names

    subroutine timing_writer(ctx, file, t, stream, init)
        ! timing.txt: the time and the wall time accumulated by every timer stage.
        class(*),         intent(inout) :: ctx
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        associate (unused => stream); end associate
        select type (ctx)
        type is (timer_class)
            call timer_write_totals(ctx, t, "s", file, init=init)
        class default
            call fesm_error("kryos_output::timing_writer", "context of stream `timing` is not a timer_class")
        end select
    end subroutine timing_writer

end module kryos_output
