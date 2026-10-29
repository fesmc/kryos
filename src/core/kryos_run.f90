module kryos_run
    ! The run: physical constants, clock, the globe and the domains with their
    ! instances, output streams, the restart schedule and bundle (design
    ! decision 10, core-design sections 3, 4, 6, 7, 9). The driver only calls
    ! run_init, loops run_step until run_done, then run_end; tests and other
    ! hosts call the same sequence.
    !
    ! This module holds the lifecycle and its setup; the sequence itself (the
    ! cold start of core-design 3.2, the restore of 3.3 and the step of 3.4,
    ! in yelmox order) is the submodule kryos_run_sequence. The restart
    ! bundles (their files and manifest) are kryos_restart's: the run writes
    ! one when its restart schedule is due, and restores from the bundle that
    ! [kryos] restart names (a folder relative to the run directory, or an
    ! absolute path).
    !
    ! Style: run_t is a plain type with public components; every procedure is
    ! free. Errors stop through fesm_error (core-design section 8). The output
    ! streams and the gather sources point into run%globe and run%dom: the run
    ! must be a target that stays in place (the driver's variable).

    use precision,        only: dp
    use fesm_log,         only: LOG_INFO, log_msg, log_set_level, fesm_error, to_str
    use ncio,             only: nc_set_lossless
    use nml,              only: nml_has_group
    use timer,            only: timer_class, timer_init, timer_step, timer_print_totals
    use phys_constants,   only: phys_const_class, phys_const_load, phys_const_require, phys_const_log
    use fesm_clock,       only: clock_t, schedule_t, clock_init, clock_start, clock_advance, clock_finish, &
                                clock_finished, clock_time, clock_step_index, clock_restart_allowed, &
                                schedule_init, schedule_due
    use kryos_config,     only: config_t, config_init, config_read, config_has, config_validate, &
                                CFG_NAME_LEN, CFG_PATH_LEN
    use kryos_component,  only: component_t, census_t, PASS_INIT, NAME_LEN
    use kryos_wiring,     only: wiring_describe
    use kryos_forcing,    only: forcing_present, forcing_kill
    use kryos_geo,        only: geo_present
    use kryos_globe,      only: globe_t, globe_init, globe_init_state, globe_end, globe_nodes
    use kryos_domain,     only: domain_t, domain_init, domain_create, domain_wire, domain_component, &
                                domain_pass, domain_end, domain_hub_write, exchange_static, &
                                NCOMP, INIT_ORDER
    use kryos_experiment, only: experiment_list_t, experiment_proto_t, experiments_create, experiments_describe, &
                                experiments_end
    use kryos_output,     only: output_t, output_init, output_add_component, output_add_hub, &
                                output_add_timing, output_write, output_describe
    use kryos_restart,    only: year_label, restart_dir, restart_write

    implicit none

    private
    public :: run_t, run_init, run_step, run_done, run_end

    character(len=*), parameter :: NO_RESTART = "None"      ! [kryos] restart: no bundle
    character(len=*), parameter :: NO_GROUP   = "none"      ! [kryos] restart_out: no restart schedule

    ! Stages of the timer (yelmox timer table; timing.txt).
    integer, parameter :: TT_OUTPUT = 1, TT_FORCING = 2, TT_GEO = 3, TT_COUPLE_ICE = 4, TT_ICE = 5, &
                          TT_HUB = 6, TT_CLIMATE = 7, TT_SURFACE = 8, TT_SHELF = 9, TT_EXPER = 10
    character(len=16), parameter :: TT_LABELS(10) = [character(len=16) :: "output", "forcing", "geo", &
        "couple_ice", "ice", "hub", "climate", "surface", "shelf", "experiments"]

    type :: run_t
        type(config_t)     :: cfg                           ! par file, [kryos] keys
        character(len=CFG_PATH_LEN) :: out_dir = "./", restart = "None"
        type(phys_const_class) :: cnst                      ! the run's one record (core-design 9)
        type(clock_t)      :: clk
        ! Schedules of the run itself (the output streams have theirs in out):
        ! the groups a host asks for (run_init out_groups), then the restart
        ! schedule ([kryos] restart_out), index i_rst. due(i): schedule i is
        ! due at the current time; set once per time by run_init and run_step.
        type(schedule_t),            allocatable :: sched(:)
        character(len=CFG_NAME_LEN), allocatable :: sched_group(:)
        logical,                     allocatable :: due(:)
        integer            :: i_rst = 0                     ! restart schedule in sched, 0 = none
        integer            :: n_dom = 0                     ! size(dom)
        type(census_t)     :: census                        ! admitted instances (backend constraints)
        type(globe_t)      :: globe
        type(domain_t), allocatable :: dom(:)               ! in [kryos] domains order
        type(experiment_list_t), allocatable :: exps(:)     ! the experiments of each domain, indexed like dom (core-design 3.9)
        integer, allocatable :: dom_order(:)                ! domain indices in name order (gathers)
        type(output_t)     :: out
        type(timer_class)  :: tmr                           ! wall time per stage
    end type run_t

    interface                                               ! submodule kryos_run_sequence
        module subroutine run_cold_start(run)
            ! core-design 3.2 (with 3.6 items 1-3): every domain's initial state.
            type(run_t), target, intent(inout) :: run
        end subroutine run_cold_start

        module subroutine run_restore(run)
            ! core-design 3.3 (with 3.6 item 4): every instance from the bundle
            ! [kryos] restart names; restore, never step.
            type(run_t), target, intent(inout) :: run
        end subroutine run_restore

        module subroutine run_sequence(run)
            ! core-design 3.4 steps 1-4 at the clock's (advanced) time.
            type(run_t), target, intent(inout) :: run
        end subroutine run_sequence
    end interface

contains

    ! =========================================================================
    ! Lifecycle
    ! =========================================================================

    subroutine run_init(run, par_file, out_groups, protos)
        ! core-design section 3.1. `out_groups`: schedule groups the run
        ! evaluates besides [kryos] restart_out (run%due), for tests and hosts
        ! that dispatch their own output; the instances' streams are output's.
        ! `protos`: the experiments of kind "host" the host provides, selected
        ! by the name of their par-file group (core-design 3.9).
        type(run_t),      intent(out), target  :: run
        character(len=*), intent(in)           :: par_file
        character(len=*), intent(in), optional :: out_groups(:)
        type(experiment_proto_t), intent(in), optional :: protos(:)

        integer :: d

        ! 1. config, logging, constants
        call config_init(run%cfg, par_file)
        call log_set_level(run%cfg%log_level)
        call nc_set_lossless(run%cfg%write_lossless)
        run%out_dir = with_slash(run%cfg%out_dir)
        run%restart = run%cfg%restart
        call run_const_load(run)

        ! 2. clock
        call run_clock_init(run)

        ! 3. globe (forcing, geo global), before the domains: the climate checks
        !    the coefficients it names against the forcing
        if (nml_has_group(trim(run%cfg%par_file), trim(run%cfg%globe))) then
            call globe_init(run%globe, run%cfg, run%clk, run%cnst, run%census)
        else
            call log_msg(LOG_INFO, "kryos_run", "no group &"//trim(run%cfg%globe)// &
                         " in the par file: the run has no global instances")
        end if

        ! 4. domains: groups, grids, hub, refs; instances; rows (with the globe's nodes)
        run%n_dom = size(run%cfg%domains)
        allocate(run%dom(run%n_dom))
        do d = 1, run%n_dom
            call domain_init(run%dom(d), run%cfg, trim(run%cfg%domains(d)), run%cnst, clock_time(run%clk))
            call domain_create(run%dom(d), run%cfg, run%clk, run%cnst, run%census)
        end do
        do d = 1, run%n_dom
            call domain_wire(run%dom(d), run%clk, extra=globe_nodes(run%globe, run%clk))
        end do
        allocate(run%exps(run%n_dom))
        do d = 1, run%n_dom
            if (present(protos)) then
                call experiments_create(run%exps(d), run%cfg, run%dom(d), run%clk, run%cnst, protos)
            else
                call experiments_create(run%exps(d), run%cfg, run%dom(d), run%clk, run%cnst)
            end if
        end do
        call run_check_names(run)
        run%dom_order = name_order(run%dom)

        ! 5. schedules: the run's (restart), the output streams (instances, hubs, timing)
        call run_schedules_init(run, out_groups)
        call run_output_init(run)

        ! 6. the exports of every instance as initialised at t0 (PASS_INIT), the
        !    static rows; then restore or cold start
        call run_pass_init(run)
        if (trim(run%restart) /= NO_RESTART) then
            call run_restore(run)
        else
            call run_cold_start(run)
        end if

        ! 7. describe, start the clock (checks the restart times against every
        !    cadence, logs the table), output at time_init
        call run_describe(run)
        call clock_start(run%clk)
        call run_schedules_due(run)
        call timer_init(run%tmr, TT_LABELS)
        call output_write(run%out, run%clk)
        call timer_step(run%tmr, TT_OUTPUT)

    end subroutine run_init

    subroutine run_step(run)
        ! core-design section 3.4: advance, sequence, output + restart at the new time.
        type(run_t), intent(inout), target :: run

        call clock_advance(run%clk)
        call run_sequence(run)

        ! 5. end of step: the kill switch (yelmox coefs_kill), output, restart
        if (forcing_present(run%globe%forcing)) then
            if (forcing_kill(run%globe%forcing)) then
                call log_msg(LOG_INFO, "kryos_run", "forcing kill switch tripped at time = "// &
                             year_label(clock_time(run%clk))//" yr: the run ends after this step")
                call clock_finish(run%clk)
            end if
        end if
        call run_schedules_due(run)
        call output_write(run%out, run%clk)
        if (run%i_rst > 0) then
            if (run%due(run%i_rst)) call run_restart_write(run)
        end if
        call timer_step(run%tmr, TT_OUTPUT)

        call log_msg(LOG_INFO, "kryos_run", "step "//to_str(clock_step_index(run%clk))// &
                     ": time = "//year_label(clock_time(run%clk))//" yr")

    end subroutine run_step

    logical function run_done(run)
        type(run_t), intent(in) :: run
        run_done = clock_finished(run%clk)
    end function run_done

    subroutine run_end(run)
        ! core-design section 3.5: every instance's end (the globe, then each
        ! domain with its ice last), the timer totals. No output: it was
        ! written by the last step.
        type(run_t), intent(inout) :: run
        integer :: d

        call globe_end(run%globe)
        do d = 1, run%n_dom
            call experiments_end(run%exps(d))
            call domain_end(run%dom(d))
        end do
        call timer_print_totals(run%tmr, "s")
        call log_msg(LOG_INFO, "kryos_run", "end: time = "//year_label(clock_time(run%clk))//" yr after "// &
                     to_str(clock_step_index(run%clk))//" steps")

    end subroutine run_end

    ! =========================================================================
    ! Init pieces
    ! =========================================================================

    subroutine run_const_load(run)
        ! core-design 9: the run's physical constants, group [kryos] phys_const
        ! of [kryos] phys_const_file (yelmox: phys_const_load of the domain's
        ! group of input/yelmo_phys_const.nml); handed to every instance and to
        ! the hubs.
        type(run_t), intent(inout) :: run
        logical :: found

        inquire(file=trim(run%cfg%phys_const_file), exist=found)
        if (.not. found) then
            call fesm_error("kryos_run::run_init", "physical constants file not found: "// &
                            trim(run%cfg%phys_const_file), "[kryos] phys_const_file of "//trim(run%cfg%par_file))
        end if
        call phys_const_load(run%cnst, trim(run%cfg%phys_const_file), group=trim(run%cfg%phys_const))
        call phys_const_require(run%cnst, "kryos_run::run_init")
    end subroutine run_const_load

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

    subroutine run_check_names(run)
        ! Domain names, instance names and experiment names are unique in a run
        ! (core-design 5, 6): output files and restart bundles are named by them.
        type(run_t), intent(in), target :: run

        character(len=NAME_LEN), allocatable :: names(:), owners(:)
        class(component_t), pointer :: c
        integer :: d, j, i

        allocate(names(0), owners(0))
        if (forcing_present(run%globe%forcing)) call add(run%globe%forcing%c%name, "the globe")
        if (geo_present(run%globe%geo))         call add(run%globe%geo%c%name, "the globe")
        do d = 1, run%n_dom
            do j = 1, d - 1
                if (trim(run%dom(j)%name) == trim(run%dom(d)%name)) then
                    call fesm_error("kryos_run::run_init", "domain `"//trim(run%dom(d)%name)// &
                                    "` is listed twice in [kryos] domains", "par file "//trim(run%cfg%par_file))
                end if
            end do
            do j = 1, NCOMP
                c => domain_component(run%dom(d), j)
                if (associated(c)) call add(c%name, "domain "//trim(run%dom(d)%name))
            end do
            if (allocated(run%exps(d)%s)) then
                do j = 1, size(run%exps(d)%s)
                    call add(run%exps(d)%s(j)%c%name, "an experiment of domain "//trim(run%dom(d)%name))
                end do
            end if
        end do

    contains

        subroutine add(name, owner)
            character(len=*), intent(in) :: name, owner
            do i = 1, size(names)
                if (trim(names(i)) /= trim(name)) cycle
                call fesm_error("kryos_run::run_init", "instance `"//trim(name)//"` is used twice in the run", &
                                "by "//trim(owners(i))//" and by "//trim(owner)//": instance names name the "// &
                                "output files and restart bundles, one group per instance")
            end do
            names  = [character(len=NAME_LEN) :: names, name]
            owners = [character(len=NAME_LEN) :: owners, owner]
        end subroutine add

    end subroutine run_check_names

    subroutine run_schedules_init(run, out_groups)
        ! One schedule per distinct group of `out_groups`, then the restart
        ! schedule from [kryos] restart_out (restart = .true.: no time_init,
        ! always time_end). A group is a par file group of template [schedule];
        ! method is required.
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

    subroutine run_output_init(run)
        ! core-design 6: the streams of every instance (the globe's, then each
        ! domain's in INIT_ORDER, then its hub, <domain>_hub.nc on [domain]
        ! out_2D) and the timing table (on the schedule of the first 1D stream).
        ! Before clock_start: the schedules register with the clock.
        type(run_t), intent(inout), target :: run

        class(component_t), pointer :: c
        integer :: d, j

        call output_init(run%out, run%cfg)
        if (forcing_present(run%globe%forcing)) call output_add_component(run%out, run%cfg, run%clk, run%globe%forcing%c)
        if (geo_present(run%globe%geo))         call output_add_component(run%out, run%cfg, run%clk, run%globe%geo%c)
        do d = 1, run%n_dom
            do j = 1, NCOMP
                c => domain_component(run%dom(d), INIT_ORDER(j))
                if (associated(c)) call output_add_component(run%out, run%cfg, run%clk, c)
            end do
            call output_add_hub(run%out, run%cfg, run%clk, run%dom(d)%name, hub_writer, run%dom(d))
        end do
        call output_add_timing(run%out, run%cfg, run%clk, run%tmr)
        call output_describe(run%out)
    end subroutine run_output_init

    subroutine hub_writer(ctx, file, t, stream, init)
        ! The writer of a domain's hub stream (output_writer_i): domain_hub_write.
        class(*),         intent(inout) :: ctx
        character(len=*), intent(in)    :: file, stream
        real(dp),         intent(in)    :: t
        logical,          intent(in)    :: init
        associate (unused => stream); end associate
        select type (ctx)
        type is (domain_t)
            call domain_hub_write(ctx, file, t, init)
        class default
            call fesm_error("kryos_run::hub_writer", "context of stream `hub` is not a domain_t")
        end select
    end subroutine hub_writer

    subroutine run_pass_init(run)
        ! core-design 3.1 step 6 and 3.6 item 1: PASS_INIT of every instance
        ! (the globe, then each domain in INIT_ORDER: each puts the exports of
        ! its state as initialised at t0), then each domain's static rows (ice
        ! reference -> geo, geothermal and sediments -> ice, hub z_srf_ref ->
        ! climate).
        type(run_t), intent(inout), target :: run
        integer :: d, j

        call globe_init_state(run%globe, run%clk, PASS_INIT)
        do d = 1, run%n_dom
            do j = 1, NCOMP
                call domain_pass(run%dom(d), INIT_ORDER(j), run%clk, PASS_INIT)
            end do
            call exchange_static(run%dom(d), clock_time(run%clk), run%clk)
        end do
    end subroutine run_pass_init

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
        ! core-design section 8: the run line, the constants, the rows of each
        ! domain (after the field states, logged by domain_create); the clock
        ! table follows from clock_start.
        type(run_t), intent(in) :: run
        character(len=1), parameter :: nl = new_line("a")
        integer :: d

        call log_msg(LOG_INFO, "kryos_run", "run: "//to_str(run%n_dom)//" domains", &
                     "par file "//trim(run%cfg%par_file)//nl// &
                     "out dir  "//trim(run%out_dir)//nl// &
                     "restart  "//trim(run%restart))
        if (run%cnst%initialized) call phys_const_log(run%cnst)
        do d = 1, run%n_dom
            call log_msg(LOG_INFO, "kryos_run", "domain "//trim(run%dom(d)%name)//":")
            call wiring_describe(run%dom(d)%wires)
            call experiments_describe(run%exps(d))
        end do
    end subroutine run_describe

    ! =========================================================================
    ! Restart bundle
    ! =========================================================================

    subroutine run_restart_write(run)
        ! core-design section 7: the bundle <out_dir>/restart-<time>/ (kryos_restart).
        ! A run ended by the kill switch off the cadence grids writes none
        ! (fesm_clock clock_restart_allowed) and says why.
        type(run_t), intent(inout), target :: run

        if (.not. clock_restart_allowed(run%clk)) then
            call log_msg(LOG_INFO, "kryos_run", "no restart bundle at time = "//year_label(clock_time(run%clk))// &
                         " yr: the time is not on the grid of every cadence (the run was killed between them)")
            return
        end if
        call restart_write(restart_dir(run%out_dir, clock_time(run%clk)), run%clk, run%globe, run%dom)
    end subroutine run_restart_write

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    function name_order(dom) result(idx)
        ! Indices of the domains in name order (gathers sum in this order,
        ! contracts-design 4.2).
        type(domain_t), intent(in) :: dom(:)
        integer, allocatable :: idx(:)
        integer :: i, j, k

        idx = [(i, i = 1, size(dom))]
        do i = 2, size(idx)
            k = idx(i)
            j = i - 1
            do while (j >= 1)
                if (llt(trim(dom(k)%name), trim(dom(idx(j))%name))) then
                    idx(j + 1) = idx(j)
                    j = j - 1
                else
                    exit
                end if
            end do
            idx(j + 1) = k
        end do
    end function name_order

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

end module kryos_run
