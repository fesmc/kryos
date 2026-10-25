program test_run_clock
    ! kryos_run, milestone C3 (core-design section 12): the clock, the
    ! schedules, a loop over zero domains and the clock-only restart.
    !
    !  (a) output and restart times of both gate configs (yelmox_Greenland,
    !      yelmox_Antarctica; gate overrides time_end = 1000, tm_2D dt = 100,
    !      tm_2Dsm dt = 10, tm_1D dt = 10, tm_rst dt = 500) equal the lists of
    !      yelmox's fesm-utils `timeout` for the same groups, and a bundle is
    !      written at each restart time. The one difference is by design: a
    !      restart schedule excludes time_init (core-design section 4.2), yelmox's
    !      tm_rst list includes it, so the lists are compared above time_init.
    !  (b) continuation: 0 -> 500 writing a clock-only bundle, restart from it,
    !      500 -> 1000: the same tick sequence and the same output and restart
    !      times as one 0 -> 1000 run, and time_elapsed runs on.
    !  (c) a restart bundle at a time other than time_init is a logged branch.
    !
    ! The yelmox gate par files are read for their ctrl and tm_* keys when they
    ! are found ($KRYOS_YELMOX_PAR_DIR, else the topo-prep worktree of yelmox);
    ! otherwise the same values are built in. Par files and output go to
    ! tmp/test_run_clock/ (run from the repo root, as make check does). Abort
    ! paths run in a child process of this executable ("abort <case>").

    use, intrinsic :: iso_fortran_env, only: output_unit
    use precision,   only: dp
    use fesm_log,    only: log_set_unit
    use fesm_clock,  only: clock_time, clock_time_elapsed, clock_time_init, clock_time_end, clock_dtt
    use nml,         only: nml_read
    use timeout,     only: timeout_class, timeout_init
    use kryos_run,   only: run_t, run_init, run_step, run_done, run_end

    implicit none

    character(len=*), parameter :: dir = "tmp/test_run_clock/"
    integer :: nfail = 0

    type :: gate_t                          ! the timeline keys of a yelmox gate config
        character(len=16) :: name = ""
        character(len=16) :: method = "const"
        real(dp) :: time_init = 0.0_dp, dtt = 10.0_dp, time_const = 0.0_dp
        logical  :: from_file = .false.
    end type gate_t

    type :: trace_t                         ! what a run did, tick by tick
        real(dp), allocatable :: t(:)       ! time at init and after every step
        logical,  allocatable :: due(:,:)   ! (tick, schedule)
        character(len=64), allocatable :: group(:)
        integer  :: i_rst = 0
        logical  :: continued = .false., branched = .false.
        real(dp) :: elapsed = 0.0_dp, elapsed_init = 0.0_dp
    end type trace_t

    character(len=512) :: arg

    call get_command_argument(1, arg)
    if (trim(arg) == "abort") then
        call get_command_argument(2, arg)
        call abort_case(trim(arg))          ! must stop; falling through is a failure
        stop 0
    end if

    call execute_command_line("rm -rf "//dir//" && mkdir -p "//dir)

    call test_gate("Greenland")
    call test_gate("Antarctica")
    call test_continuation()
    call test_branch()
    call test_aborts()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_run_clock: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_run_clock: OK"

contains

    ! -------------------------------------------------------------------------
    ! Tests
    ! -------------------------------------------------------------------------

    subroutine test_gate(name)
        ! (a) the kryos schedules reproduce yelmox's timeout lists.
        character(len=*), intent(in) :: name

        character(len=*), parameter :: tm(4) = [character(len=8) :: "tm_2D", "tm_2Dsm", "tm_1D", "tm_rst"]
        character(len=*), parameter :: kr(4) = [character(len=8) :: "out_2D", "out_2Dsm", "out_1D", "out_rst"]
        real(dp), parameter :: time_end = 1000.0_dp
        real(dp), parameter :: dt(4) = [100.0_dp, 10.0_dp, 10.0_dp, 500.0_dp]
        type(gate_t)  :: g
        type(trace_t) :: tr
        type(timeout_class) :: tmo
        character(len=256) :: kpar, ypar
        real(dp), allocatable :: tk(:), ty(:)
        integer :: i, j
        logical :: found

        call gate_params(name, g)
        call check(trim(g%method) == "const" .and. g%dtt == 10.0_dp .and. g%time_init == 0.0_dp, &
                   name//": gate timeline (const, dtt = 10, time_init = 0)")

        kpar = dir//name//"_kryos.nml"
        ypar = dir//name//"_yelmox.nml"
        call write_kryos_par(kpar, dir//name//"/", "None", g%time_init, time_end, g%dtt, g%time_const, 500.0_dp, &
                             dt(1), dt(2), dt(3), "warn")
        call put(ypar, [character(len=100) :: &
            '&tm_2D',   '    method = "const"', '    dt = '//num(dt(1)), '/', &
            '&tm_2Dsm', '    method = "const"', '    dt = '//num(dt(2)), '/', &
            '&tm_1D',   '    method = "const"', '    dt = '//num(dt(3)), '/', &
            '&tm_rst',  '    method = "const"', '    dt = '//num(dt(4)), '/'])

        call run_trace(kpar, [character(len=8) :: "out_2D", "out_2Dsm", "out_1D"], tr)
        call check(tr%t(1) == g%time_init .and. tr%t(size(tr%t)) == time_end .and. size(tr%t) == 101, &
                   name//": ticks 0 .. 1000, 101 times")
        call check(all(tr%t(2:) - tr%t(:size(tr%t)-1) == g%dtt), name//": every step is dtt")
        call check(size(tr%group) == 4 .and. tr%i_rst == 4, name//": three output schedules and the restart")

        do i = 1, 4
            j = schedule_index(tr, trim(kr(i)))
            if (i == 4) j = tr%i_rst                         ! the restart schedule
            call check(j > 0, name//": schedule "//trim(kr(i))//" registered")
            if (j == 0) cycle
            tk = due_times(tr, j)

            call timeout_init(tmo, ypar, trim(tm(i)), "yelmox "//trim(tm(i)), g%time_init, time_end)
            ty = tmo%times
            if (i == 4) ty = pack(ty, ty > g%time_init)      ! kryos: no restart at time_init
            call check(size(tk) == size(ty), name//": "//trim(kr(i))//" vs "//trim(tm(i))//": count")
            if (size(tk) == size(ty)) call check(all(tk == ty), name//": "//trim(kr(i))//" vs "//trim(tm(i))//": times")
        end do

        ! a bundle at every restart time above time_init, none at time_init
        tk = due_times(tr, tr%i_rst)
        call check(size(tk) == 2 .and. tk(1) == 500.0_dp .and. tk(2) == 1000.0_dp, name//": restart times 500, 1000")
        do i = 1, size(tk)
            inquire(file=dir//name//"/restart-"//itoa(nint(tk(i)))//"/clock_restart.nc", exist=found)
            call check(found, name//": bundle written at "//itoa(nint(tk(i))))
        end do
        inquire(file=dir//name//"/restart-0/clock_restart.nc", exist=found)
        call check(.not. found, name//": no bundle at time_init")
        inquire(file=dir//name//"/restart-510/clock_restart.nc", exist=found)
        call check(.not. found, name//": no bundle off the schedule")
    end subroutine test_gate

    subroutine test_continuation()
        ! (b) 0 -> 500 + restart + 500 -> 1000 = 0 -> 1000.
        type(trace_t) :: full, c1, c2
        real(dp), allocatable :: t(:), tf(:)
        integer :: i, j, ng
        logical :: found, ok
        character(len=64) :: groups(3)

        groups = ["out_2D  ", "out_2Dsm", "out_1D  "]
        call write_kryos_par(dir//"full.nml", dir//"full/", "None", 0.0_dp, 1000.0_dp, 10.0_dp, 0.0_dp, &
                             500.0_dp, 100.0_dp, 10.0_dp, 10.0_dp, "warn")
        call write_kryos_par(dir//"c1.nml", dir//"c1/", "None", 0.0_dp, 500.0_dp, 10.0_dp, 0.0_dp, &
                             500.0_dp, 100.0_dp, 10.0_dp, 10.0_dp, "warn")
        call write_kryos_par(dir//"c2.nml", dir//"c2/", dir//"c1/restart-500", 500.0_dp, 1000.0_dp, 10.0_dp, &
                             0.0_dp, 500.0_dp, 100.0_dp, 10.0_dp, 10.0_dp, "info")

        call run_trace(dir//"full.nml", groups, full)
        call run_trace(dir//"c1.nml", groups, c1)
        inquire(file=dir//"c1/restart-500/clock_restart.nc", exist=found)
        call check(found, "continuation: chunk 1 wrote restart-500/clock_restart.nc")

        call capture_log(dir//"c2.log", .true.)
        call run_trace(dir//"c2.nml", groups, c2)
        call capture_log(dir//"c2.log", .false.)
        call check(log_has(dir//"c2.log", "continuation from"), "continuation: logged as a continuation")
        call check(.not. log_has(dir//"c2.log", "branch from"), "continuation: not logged as a branch")

        call check(c2%continued .and. .not. c2%branched, "continuation: clock restored as a continuation")
        call check(c2%elapsed_init == 500.0_dp .and. c2%elapsed == 1000.0_dp, &
                   "continuation: time_elapsed runs on (500 at restart, 1000 at the end)")

        ! the tick sequence: chunk 1, then chunk 2 without its first (restart) time
        t = [c1%t, c2%t(2:)]
        call check(size(t) == size(full%t), "continuation: same number of ticks")
        if (size(t) == size(full%t)) call check(all(t == full%t), "continuation: same tick sequence")
        call check(c1%t(size(c1%t)) == c2%t(1), "continuation: chunk 2 starts where chunk 1 ended")

        ! the output and restart times of the two chunks together
        ng = size(full%group)
        call check(size(c1%group) == ng .and. size(c2%group) == ng, "continuation: same schedules in every run")
        do i = 1, ng
            if (i == full%i_rst) cycle
            j = schedule_index(c1, trim(full%group(i)))
            tf = due_times(full, i)
            t = union(due_times(c1, j), due_times(c2, schedule_index(c2, trim(full%group(i)))))
            ok = size(t) == size(tf)
            if (ok) ok = all(t == tf)
            call check(ok, "continuation: times of "//trim(full%group(i))//" (chunk 1 + chunk 2 = full)")
        end do
        tf = due_times(full, full%i_rst)
        t = [due_times(c1, c1%i_rst), due_times(c2, c2%i_rst)]
        call check(size(t) == 2 .and. size(tf) == 2, "continuation: restart: one bundle per chunk")
        if (size(t) == 2 .and. size(tf) == 2) call check(all(t == tf), "continuation: restart times 500, 1000")
        inquire(file=dir//"c2/restart-1000/clock_restart.nc", exist=found)
        call check(found, "continuation: chunk 2 wrote restart-1000")
    end subroutine test_continuation

    subroutine test_branch()
        ! (c) the bundle time (500) differs from time_init (200): a branch.
        type(trace_t) :: br
        character(len=64) :: groups(1)
        real(dp), allocatable :: tk(:)

        groups = ["out_2D"]
        call write_kryos_par(dir//"br.nml", dir//"br/", dir//"c1/restart-500", 200.0_dp, 1000.0_dp, 10.0_dp, &
                             0.0_dp, 500.0_dp, 100.0_dp, 10.0_dp, 10.0_dp, "info")
        call capture_log(dir//"br.log", .true.)
        call run_trace(dir//"br.nml", groups, br)
        call capture_log(dir//"br.log", .false.)

        call check(log_has(dir//"br.log", "branch from"), "branch: logged as a branch")
        call check(.not. log_has(dir//"br.log", "continuation from"), "branch: not logged as a continuation")
        call check(br%branched .and. .not. br%continued, "branch: clock flags")
        call check(br%elapsed_init == 0.0_dp .and. br%elapsed == 800.0_dp, "branch: time_elapsed restarts at the branch")
        call check(size(br%t) == 81 .and. br%t(1) == 200.0_dp .and. br%t(81) == 1000.0_dp, "branch: ticks 200 .. 1000")
        tk = due_times(br, br%i_rst)
        call check(size(tk) == 2 .and. tk(1) == 500.0_dp .and. tk(2) == 1000.0_dp, "branch: restart times 500, 1000")
    end subroutine test_branch

    subroutine test_aborts()
        ! The stopping paths: the child process must exit non-zero.
        call write_kryos_par(dir//"ok.nml", dir//"ok/", "None", 0.0_dp, 50.0_dp, 10.0_dp, 0.0_dp, 50.0_dp, &
                             10.0_dp, 10.0_dp, 10.0_dp, "warn")
        call put(dir//"nokey.nml", [character(len=100) :: '&kryos', '    log_level = "warn"', '/', &
                                    '&clock', '    time_init = 0.0', '    dtt = 10.0', '/'])
        call put(dir//"domains.nml", [character(len=100) :: '&kryos', '    log_level = "warn"', &
                                      '    domains = "Greenland"', '/', &
                                      '&clock', '    time_init = 0.0', '    time_end = 50.0', '    dtt = 10.0', '/'])
        call write_kryos_par(dir//"nobundle.nml", dir//"nobundle/", dir//"no_such_bundle", 0.0_dp, 50.0_dp, &
                             10.0_dp, 0.0_dp, 50.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, "warn")
        call write_kryos_par(dir//"badend.nml", dir//"badend/", "None", 0.0_dp, 55.0_dp, 10.0_dp, 0.0_dp, 50.0_dp, &
                             10.0_dp, 10.0_dp, 10.0_dp, "warn")

        call check(.not. aborts("valid"),    "abort: a valid par file does not stop")
        call check(aborts("missing-key"),    "abort: [clock] without time_end stops")
        call check(aborts("domains"),        "abort: [kryos] domains stops until C6")
        call check(aborts("missing-bundle"), "abort: a restart folder without clock_restart.nc stops")
        call check(aborts("off-grid-end"),   "abort: time_end off the step grid stops")
    end subroutine test_aborts

    subroutine abort_case(name)
        ! Runs in the child process.
        character(len=*), intent(in) :: name
        type(run_t) :: run

        select case (name)
        case ("valid");          call run_init(run, dir//"ok.nml")
        case ("missing-key");    call run_init(run, dir//"nokey.nml")
        case ("domains");        call run_init(run, dir//"domains.nml")
        case ("missing-bundle"); call run_init(run, dir//"nobundle.nml")
        case ("off-grid-end");   call run_init(run, dir//"badend.nml")
        end select
    end subroutine abort_case

    logical function aborts(name)
        ! True if `<this program> abort <name>` exits non-zero.
        character(len=*), intent(in) :: name
        character(len=512) :: exe
        integer :: es

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" abort "//name//" > /dev/null 2>&1", exitstat=es)
        aborts = es /= 0
    end function aborts

    ! -------------------------------------------------------------------------
    ! The run under test
    ! -------------------------------------------------------------------------

    subroutine run_trace(par, groups, tr)
        ! run_init, run_step to the end, run_end; record every tick.
        character(len=*), intent(in)  :: par, groups(:)
        type(trace_t),    intent(out) :: tr

        type(run_t) :: run
        integer :: k, nt

        call run_init(run, par, groups)
        nt = nint((clock_time_end(run%clk) - clock_time_init(run%clk))/clock_dtt(run%clk)) + 1
        allocate(tr%t(nt), tr%due(nt, size(run%sched)))
        tr%group = run%sched_group
        tr%i_rst = run%i_rst
        tr%continued = run%clk%continued
        tr%branched  = run%clk%branched
        tr%elapsed_init = clock_time_elapsed(run%clk)

        k = 1
        tr%t(k)     = clock_time(run%clk)
        tr%due(k,:) = run%due
        do while (.not. run_done(run))
            call run_step(run)
            k = k + 1
            tr%t(k)     = clock_time(run%clk)
            tr%due(k,:) = run%due
        end do
        call check(k == nt, "run_trace: step count")
        tr%elapsed = clock_time_elapsed(run%clk)
        call run_end(run)
    end subroutine run_trace

    integer function schedule_index(tr, group)
        type(trace_t),    intent(in) :: tr
        character(len=*), intent(in) :: group
        integer :: i

        schedule_index = 0
        do i = 1, size(tr%group)
            if (i == tr%i_rst) cycle
            if (trim(tr%group(i)) == group) schedule_index = i
        end do
    end function schedule_index

    function due_times(tr, i) result(t)
        ! The times at which schedule i was due.
        type(trace_t), intent(in) :: tr
        integer,       intent(in) :: i
        real(dp), allocatable :: t(:)
        t = pack(tr%t, tr%due(:,i))
    end function due_times

    function union(a, b) result(u)
        ! Sorted union without duplicates.
        real(dp), intent(in) :: a(:), b(:)
        real(dp), allocatable :: u(:), w(:)
        real(dp) :: x
        integer :: i, j, n

        w = [a, b]
        n = size(w)
        do i = 2, n                                  ! insertion sort
            x = w(i)
            j = i - 1
            do while (j >= 1)
                if (w(j) <= x) exit
                w(j+1) = w(j)
                j = j - 1
            end do
            w(j+1) = x
        end do
        u = w(1:0)
        do i = 1, n
            if (i > 1) then
                if (w(i) == w(i-1)) cycle
            end if
            u = [u, w(i)]
        end do
    end function union

    ! -------------------------------------------------------------------------
    ! Gate parameters and par files
    ! -------------------------------------------------------------------------

    subroutine gate_params(name, g)
        ! The timeline of a yelmox gate config: from its par file when found,
        ! else the built-in values (identical in both files).
        character(len=*), intent(in)  :: name
        type(gate_t),     intent(out) :: g

        character(len=512) :: path, base, method
        integer :: i, st
        logical :: found, ok

        g%name = name
        base = ""
        call get_environment_variable("KRYOS_YELMOX_PAR_DIR", base, status=st)
        if (st /= 0 .or. len_trim(base) == 0) then
            call get_environment_variable("HOME", base)
            base = trim(base)//"/models/yelmox/.claude/worktrees/topo-prep/yelmox"
        end if
        path = trim(base)//"/yelmox_"//name//".nml"
        inquire(file=trim(path), exist=found)
        if (.not. found) then
            write(*,'(a)') "test_run_clock: "//trim(path)//" not found: built-in "//name//" gate values"
            return
        end if

        g%from_file = .true.
        call nml_read(trim(path), "ctrl", "time_init",    g%time_init)
        call nml_read(trim(path), "ctrl", "dtt",          g%dtt)
        call nml_read(trim(path), "ctrl", "tstep_method", g%method)
        call nml_read(trim(path), "ctrl", "tstep_const",  g%time_const)
        ok = .true.
        do i = 1, 4
            method = ""
            select case (i)
            case (1); call nml_read(trim(path), "tm_2D",   "method", method)
            case (2); call nml_read(trim(path), "tm_2Dsm", "method", method)
            case (3); call nml_read(trim(path), "tm_1D",   "method", method)
            case (4); call nml_read(trim(path), "tm_rst",  "method", method)
            end select
            if (trim(method) /= "const") ok = .false.
        end do
        call check(ok, name//": every tm_* group of "//trim(path)//" is method = const")
        write(*,'(a)') "test_run_clock: gate values of "//trim(path)
    end subroutine gate_params

    subroutine write_kryos_par(file, out_dir, restart, time_init, time_end, dtt, time_const, rst_dt, d2D, d2Dsm, d1D, level)
        ! [kryos] (no domain), [clock] (method const), three output schedules
        ! out_2D, out_2Dsm, out_1D and the restart schedule out_rst.
        character(len=*), intent(in) :: file, out_dir, restart, level
        real(dp),         intent(in) :: time_init, time_end, dtt, time_const, rst_dt, d2D, d2Dsm, d1D

        call put(file, [character(len=200) :: &
            '&kryos', &
            '    out_dir     = "'//trim(out_dir)//'"', &
            '    restart     = "'//trim(restart)//'"', &
            '    restart_out = "out_rst"', &
            '    log_level   = "'//trim(level)//'"', &
            '/', &
            '&clock', &
            '    method     = "const"', &
            '    time_const = '//num(time_const), &
            '    time_init  = '//num(time_init), &
            '    time_end   = '//num(time_end), &
            '    dtt        = '//num(dtt), &
            '/', &
            '&out_2D',   '    method = "const"', '    dt = '//num(d2D),   '/', &
            '&out_2Dsm', '    method = "const"', '    dt = '//num(d2Dsm), '/', &
            '&out_1D',   '    method = "const"', '    dt = '//num(d1D),   '/', &
            '&out_rst',  '    method = "const"', '    dt = '//num(rst_dt), '/'])
    end subroutine write_kryos_par

    subroutine put(file, lines)
        character(len=*), intent(in) :: file, lines(:)
        integer :: u, i

        open(newunit=u, file=file, status="replace", action="write")
        do i = 1, size(lines)
            write(u,'(a)') trim(lines(i))
        end do
        close(u)
    end subroutine put

    ! -------------------------------------------------------------------------
    ! Log capture and helpers
    ! -------------------------------------------------------------------------

    subroutine capture_log(file, on)
        ! Route the log (INFO) to `file` while on, back to stdout when off.
        character(len=*), intent(in) :: file
        logical,          intent(in) :: on
        integer, save :: u = 0          ! newunit values are negative; 0 = none open

        if (on) then
            open(newunit=u, file=file, status="replace", action="write")
            call log_set_unit(u)
        else
            call log_set_unit(output_unit)
            if (u /= 0) close(u)
            u = 0
        end if
    end subroutine capture_log

    logical function log_has(file, text)
        character(len=*), intent(in) :: file, text
        character(len=1024) :: line
        integer :: u, ios

        log_has = .false.
        open(newunit=u, file=file, status="old", action="read", iostat=ios)
        if (ios /= 0) return
        do
            read(u,'(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, text) > 0) then
                log_has = .true.
                exit
            end if
        end do
        close(u)
    end function log_has

    function num(x) result(s)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: s
        character(len=64) :: buf

        write(buf,'(f0.6)') x
        s = trim(adjustl(buf))
        if (s(1:1) == ".") s = "0"//s                ! gfortran f0.d drops the leading zero
        if (s(1:2) == "-.") s = "-0"//s(2:)
    end function num

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=32) :: buf

        write(buf,'(i0)') i
        s = trim(buf)
    end function itoa

    subroutine check(ok, msg)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            nfail = nfail + 1
            write(*,'(a)') "test_run_clock: FAIL: "//msg
        end if
    end subroutine check

end program test_run_clock
