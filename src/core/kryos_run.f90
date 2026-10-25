module kryos_run
    ! The run: clock, output and restart schedules, the loop over the domains,
    ! and the restart bundle (design decision 10, core-design sections 3, 4, 7).
    ! The driver only calls run_init, loops run_step until run_done, then
    ! run_end; tests and other hosts call the same sequence.
    !
    ! Milestone C3 (core-design section 12): the clock ([kryos] clock group),
    ! the schedules (one per distinct group, plus the restart schedule), a loop
    ! over zero domains, and the clock-only restart bundle. Places where the
    ! other milestones plug in are marked "C6", "C7", "C8", "C9", "C10":
    !
    !   C6 domains (kryos_domain)    run%dom, domain_init, domain loops
    !   C7 globe (kryos_globe)       run%globe, globe_init, global steps
    !   C8 output (kryos_output)     run%out, output dispatch at the due schedules
    !   C9 sequence                  the cold start and the step sequence
    !   C10 restart (kryos_restart)  manifest, the bundle contents of the instances
    !
    ! Style: run_t is a plain type with public components; every procedure is
    ! free. Errors stop through fesm_error (core-design section 8).

    use precision,    only: dp
    use fesm_log,     only: LOG_INFO, log_msg, log_set_level, fesm_error, to_str
    use ncio,         only: nc_set_lossless
    use fesm_clock,   only: clock_t, schedule_t, clock_init, clock_start, clock_advance, &
                            clock_finished, clock_time, clock_step_index, clock_restart_read, &
                            clock_restart_write, schedule_init, schedule_due
    use kryos_config, only: config_t, config_init, config_read, config_has, config_validate, &
                            CFG_NAME_LEN, CFG_PATH_LEN

    implicit none

    private
    public :: run_t, run_init, run_step, run_done, run_end

    character(len=*), parameter :: NO_RESTART = "None"      ! [kryos] restart: no bundle
    character(len=*), parameter :: NO_GROUP   = "none"      ! [kryos] restart_out: no restart schedule
    character(len=*), parameter :: CLOCK_FILE = "clock_restart.nc"

    type :: run_t
        type(config_t)     :: cfg                           ! par file, [kryos] keys
        character(len=CFG_PATH_LEN) :: out_dir = "./", restart = "None"
        type(clock_t)      :: clk
        ! Schedules by group name: the output schedules first, then (if
        ! [kryos] restart_out names a group) the restart schedule, index i_rst.
        ! due(i): schedule i is due at the current time; set once per time by
        ! run_init (time_init) and run_step, read by the output dispatch (C8)
        ! and by the restart writer.
        type(schedule_t),            allocatable :: sched(:)
        character(len=CFG_NAME_LEN), allocatable :: sched_group(:)
        logical,                     allocatable :: due(:)
        integer            :: i_rst = 0                     ! restart schedule in sched, 0 = none
        integer            :: n_dom = 0                     ! domains; 0 until C6
        ! C3 (step 1): type(phys_const_class) :: cnst   (phys_const_load, with C6: the domains' ctx)
        ! C6: type(domain_t), allocatable :: dom(:)      (n_dom = size(dom))
        ! C7: type(globe_t)  :: globe
        ! C8: type(output_t) :: out
        ! C9: type(timer_class) :: tmr                   (per-stage wall time)
    end type run_t

contains

    ! =========================================================================
    ! Lifecycle
    ! =========================================================================

    subroutine run_init(run, par_file, out_groups)
        ! core-design section 3.1. `out_groups`: schedule groups to register
        ! besides those the instances name (none until C6-C8) and
        ! [kryos] restart_out; for tests and hosts that dispatch their own output.
        type(run_t),      intent(out)          :: run
        character(len=*), intent(in)           :: par_file
        character(len=*), intent(in), optional :: out_groups(:)

        ! 1. config, logging, constants
        call config_init(run%cfg, par_file)
        call log_set_level(run%cfg%log_level)
        call nc_set_lossless(run%cfg%write_lossless)
        run%out_dir = with_slash(run%cfg%out_dir)
        run%restart = run%cfg%restart
        ! C3 (step 1): phys_const_load(run%cnst, cfg%phys_const_file, ...), phys_const_require

        ! 2. clock
        call run_clock_init(run)

        ! 3-4. globe, then the domains, each instance with its cadence (C7, C6, C9)
        if (size(run%cfg%domains) > 0) then
            call fesm_error("kryos_run::run_init", "[kryos] domains is not supported yet.", &
                "domains: "//join(run%cfg%domains)//"; the domain milestones (core-design section 12, C6-C9) "// &
                "are not in this build. par file "//trim(par_file))
        end if
        run%n_dom = 0
        ! C7: call globe_init(run%globe, ...)
        ! C6: allocate(run%dom(n)); call domain_init(run%dom(d), ...) for each name of cfg%domains

        ! 5. schedules
        call run_schedules_init(run, out_groups)

        ! 6. init_state of every instance (PASS_INIT), then restore or cold start
        ! C9: <comp>_init_state(inst, st0, PASS_INIT) for the globe, then the domains
        if (trim(run%restart) /= NO_RESTART) then
            call run_restore(run)
        else
            call run_cold_start(run)
        end if

        ! 7. start the clock (checks the restart times against every cadence, logs
        ! the table), describe, output at time_init
        call clock_start(run%clk)
        call run_describe(run)
        call run_schedules_due(run)
        ! C8: output_write_due(run%out, run%due)

    end subroutine run_init

    subroutine run_step(run)
        ! core-design section 3.4: advance, sequence, output + restart at the new time.
        type(run_t), intent(inout) :: run

        call clock_advance(run%clk)
        call run_sequence(run)

        ! 5. end of step
        ! C7: if (forcing_kill(run%globe%forcing)) call clock_finish(run%clk)
        call run_schedules_due(run)
        ! C8: output_write_due(run%out, run%due)
        if (run%i_rst > 0) then
            if (run%due(run%i_rst)) call run_restart_write(run)
        end if

        call log_msg(LOG_INFO, "kryos_run", "step "//to_str(clock_step_index(run%clk))// &
                     ": time = "//year_label(clock_time(run%clk))//" yr")

    end subroutine run_step

    logical function run_done(run)
        type(run_t), intent(in) :: run
        run_done = clock_finished(run%clk)
    end function run_done

    subroutine run_end(run)
        ! core-design section 3.5. No output: it was written by the last step.
        type(run_t), intent(inout) :: run

        ! C9: <comp>_end for every instance (Yelmo last), timer_print_totals(run%tmr, ...)
        call log_msg(LOG_INFO, "kryos_run", "end: time = "//year_label(clock_time(run%clk))//" yr after "// &
                     to_str(clock_step_index(run%clk))//" steps")

    end subroutine run_end

    ! =========================================================================
    ! Init pieces
    ! =========================================================================

    subroutine run_clock_init(run)
        ! The [kryos] clock group (template [clock]) through config_read, then
        ! clock_init in its argument form, so the group may have any name. The
        ! template holds placeholders for time_init, time_end and dtt: they are
        ! required in the group. time_const / time_const_cal are passed only when
        ! the group sets them (fesm_clock needs exactly one for method = const).
        type(run_t), intent(inout) :: run

        character(len=CFG_NAME_LEN) :: group
        character(len=32)           :: method, calendar, year
        real(dp) :: time_init, time_end, dtt, time_ref, time_const, time_const_cal
        integer  :: ticks_per_year

        group = run%cfg%clock
        call config_validate(run%cfg, group, "clock", required=[character(len=9) :: "time_init", "time_end", "dtt"])

        call config_read(run%cfg, group, "time_init",      time_init,      "clock")
        call config_read(run%cfg, group, "time_end",       time_end,       "clock")
        call config_read(run%cfg, group, "dtt",            dtt,            "clock")
        call config_read(run%cfg, group, "method",         method,         "clock")
        call config_read(run%cfg, group, "calendar",       calendar,       "clock")
        call config_read(run%cfg, group, "year",           year,           "clock")
        call config_read(run%cfg, group, "time_ref",       time_ref,       "clock")
        call config_read(run%cfg, group, "ticks_per_year", ticks_per_year, "clock")

        if (config_has(run%cfg, group, "time_const")) then
            call config_read(run%cfg, group, "time_const", time_const, "clock")
            if (config_has(run%cfg, group, "time_const_cal")) then
                call config_read(run%cfg, group, "time_const_cal", time_const_cal, "clock")
                call clock_init(run%clk, time_init, time_end, dtt, calendar, year, method=method, &
                                time_const=time_const, time_const_cal=time_const_cal, &
                                time_ref=time_ref, ticks_per_year=ticks_per_year)
            else
                call clock_init(run%clk, time_init, time_end, dtt, calendar, year, method=method, &
                                time_const=time_const, time_ref=time_ref, ticks_per_year=ticks_per_year)
            end if
        else if (config_has(run%cfg, group, "time_const_cal")) then
            call config_read(run%cfg, group, "time_const_cal", time_const_cal, "clock")
            call clock_init(run%clk, time_init, time_end, dtt, calendar, year, method=method, &
                            time_const_cal=time_const_cal, time_ref=time_ref, ticks_per_year=ticks_per_year)
        else
            call clock_init(run%clk, time_init, time_end, dtt, calendar, year, method=method, &
                            time_ref=time_ref, ticks_per_year=ticks_per_year)
        end if

    end subroutine run_clock_init

    subroutine run_schedules_init(run, out_groups)
        ! One schedule per distinct group, in order: `out_groups`, then the
        ! groups the instances name (C8); then the restart schedule from
        ! [kryos] restart_out (restart = .true.: no time_init, always time_end).
        ! A group is a par file group of template [schedule]; method is required.
        type(run_t),      intent(inout)        :: run
        character(len=*), intent(in), optional :: out_groups(:)

        character(len=CFG_NAME_LEN), allocatable :: names(:)
        character(len=CFG_NAME_LEN) :: group
        integer :: i, n

        allocate(names(0))
        if (present(out_groups)) then
            do i = 1, size(out_groups)
                call add_group(names, out_groups(i))
            end do
        end if
        ! C8: call add_group(names, <instance out_* key>) for every instance, and [domain] out_2D

        n = size(names)
        group = run%cfg%restart_out
        if (len_trim(group) > 0 .and. trim(group) /= NO_GROUP) n = n + 1
        allocate(run%sched(n), run%sched_group(n), run%due(n))
        run%due = .false.

        do i = 1, size(names)
            call config_validate(run%cfg, names(i), "schedule", required=["method"])
            call schedule_init(run%sched(i), run%clk, trim(names(i)), run%cfg%par_file, trim(names(i)))
            run%sched_group(i) = names(i)
        end do

        if (n > size(names)) then
            run%i_rst = n
            call config_validate(run%cfg, group, "schedule", required=["method"])
            call schedule_init(run%sched(n), run%clk, "restart", run%cfg%par_file, trim(group), restart=.true.)
            run%sched_group(n) = group
        end if

    end subroutine run_schedules_init

    subroutine add_group(names, group)
        ! Append `group` to `names` unless it is blank, "none" or already there.
        character(len=CFG_NAME_LEN), allocatable, intent(inout) :: names(:)
        character(len=*),                          intent(in)    :: group

        if (len_trim(group) == 0 .or. trim(group) == NO_GROUP) return
        if (any(names == trim(group))) return
        names = [character(len=CFG_NAME_LEN) :: names, trim(group)]
    end subroutine add_group

    subroutine run_schedules_due(run)
        ! Evaluate every schedule once at the current time (schedule_due moves
        ! its cursor, so it must not be asked twice).
        type(run_t), intent(inout) :: run
        integer :: i

        do i = 1, size(run%sched)
            run%due(i) = schedule_due(run%sched(i), run%clk)
        end do
    end subroutine run_schedules_due

    subroutine run_describe(run)
        ! core-design section 8: the run line (the clock table is logged by clock_start).
        type(run_t), intent(in) :: run
        character(len=1), parameter :: nl = new_line("a")

        call log_msg(LOG_INFO, "kryos_run", "run: "//to_str(run%n_dom)//" domains", &
                     "par file "//trim(run%cfg%par_file)//nl// &
                     "out dir  "//trim(run%out_dir)//nl// &
                     "restart  "//trim(run%restart))
    end subroutine run_describe

    ! =========================================================================
    ! Cold start, restore, step sequence
    ! =========================================================================

    subroutine run_cold_start(run)
        ! core-design section 3.2, in the order of yelmox domain_init_state.
        type(run_t), intent(inout) :: run
        integer :: d

        do d = 1, run%n_dom
            ! C9: a) bed and sea level   b) boundary fields from the hub
            !     c) ice state           d) the surface on the final geometry
        end do
    end subroutine run_cold_start

    subroutine run_restore(run)
        ! core-design section 3.3: restore, never step. C3 restores the clock only.
        type(run_t), intent(inout) :: run

        character(len=CFG_PATH_LEN) :: file
        logical :: found
        integer :: d

        file = trim(strip_slash(run%restart))//"/"//CLOCK_FILE
        inquire(file=trim(file), exist=found)
        if (.not. found) then
            call fesm_error("kryos_run::run_restore", "restart bundle has no "//CLOCK_FILE, &
                            "[kryos] restart = "//trim(run%restart)//"; looked for "//trim(file))
        end if

        ! continuation if the bundle time = time_init, else a logged branch (fesm_clock)
        call clock_restart_read(run%clk, trim(file))
        ! C6-C9: cadence_restart_read(inst%cadence, run%clk, file) for every instance
        ! C10: manifest check against the par file; globe, then per domain the
        !      instances ice, geo, climate, surface, shelf, geothermal, sediments

        do d = 1, run%n_dom
            ! C10: domain_check_reference, domain_hub_update (the hub is derived)
        end do
    end subroutine run_restore

    subroutine run_sequence(run)
        ! core-design section 3.4, steps 1-4 (yelmox.f90 order, N domains).
        type(run_t), intent(inout) :: run
        integer :: d

        ! 1. globe, start of step (C7, C9): geo_global_begin, gather_into V_ice, forcing_update
        ! 2. hooks: none in phase 1
        ! 3. geo, every domain, then the global reduction (C9)
        do d = 1, run%n_dom
            ! C9: exchange ice -> geo, globe -> geo; geo_update when its cadence is due
        end do
        ! C9: geo_global_reduce(run%globe%geo, run%dom, t)
        ! 4. per domain: ice, hub, climate, surface, shelf (C9)
        do d = 1, run%n_dom
            ! C9: exchange_into ice; ice_update; domain_hub_update; climate, surface, shelf
            !     each when its cadence is due, followed by cadence_mark
        end do
    end subroutine run_sequence

    ! =========================================================================
    ! Restart bundle
    ! =========================================================================

    subroutine run_restart_write(run)
        ! core-design section 7: <out_dir>/restart-<time>/. C3 writes
        ! clock_restart.nc (the clock; the cadences of the instances go in the
        ! same file); C10 adds the manifest and the instances' files.
        type(run_t), intent(in) :: run

        character(len=CFG_PATH_LEN) :: dir

        dir = trim(run%out_dir)//"restart-"//year_label(clock_time(run%clk))
        call execute_command_line("mkdir -p '"//trim(dir)//"'")
        call clock_restart_write(run%clk, trim(dir)//"/"//CLOCK_FILE)
        ! C6-C9: cadence_restart_write(inst%cadence, file) for every instance
        ! C10: manifest.nml; globe, then per domain the instances (write order = restore order)

        call log_msg(LOG_INFO, "kryos_run", "restart bundle "//trim(dir))
    end subroutine run_restart_write

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    function year_label(t) result(s)
        ! `t` [yr] as the shortest decimal that reads back to the same value:
        ! 15000, 0.5, -21000 (restart-<time> folder names, log lines).
        real(dp), intent(in) :: t
        character(len=:), allocatable :: s

        character(len=64) :: buf, fmt
        real(dp) :: x
        integer  :: nd

        buf = ""
        do nd = 0, 17
            write(fmt, '(a,i0,a)') "(f0.", nd, ")"
            write(buf, fmt) t
            read(buf, *) x
            if (x == t) exit
        end do
        s = trim(adjustl(buf))
        if (s(1:1) == ".") s = "0"//s                      ! gfortran f0.d drops the leading zero
        if (len(s) > 1) then
            if (s(1:2) == "-.") s = "-0"//s(2:)
        end if
        if (s(len(s):len(s)) == ".") s = s(1:len(s)-1)
        if (s == "-0") s = "0"
    end function year_label

    function with_slash(dir) result(s)
        character(len=*), intent(in) :: dir
        character(len=:), allocatable :: s
        s = trim(dir)
        if (len(s) == 0) then
            s = "./"
        else if (s(len(s):len(s)) /= "/") then
            s = s//"/"
        end if
    end function with_slash

    function strip_slash(dir) result(s)
        character(len=*), intent(in) :: dir
        character(len=:), allocatable :: s
        s = trim(dir)
        do while (len(s) > 1)
            if (s(len(s):len(s)) /= "/") exit
            s = s(1:len(s)-1)
        end do
    end function strip_slash

    function join(names) result(s)
        character(len=*), intent(in) :: names(:)
        character(len=:), allocatable :: s
        integer :: i
        s = ""
        do i = 1, size(names)
            s = s//" "//trim(names(i))
        end do
        s = trim(adjustl(s))
    end function join

end module kryos_run
