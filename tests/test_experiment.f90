module test_experiment_toys
    ! Experiments for the mechanism tests (kryos_experiment): a recorder
    ! (logs every hook call), one that needs a toy ice (the select type stop),
    ! and a stateful counter (the restart bindings), and the toy ice types
    ! they look for in the domain.
    use precision,         only: dp
    use fesm_log,          only: fesm_error
    use kryos_component,   only: init_ctx_t, step_time_t, NAME_LEN, PATH_LEN
    use kryos_ice,         only: ice_t
    use kryos_domain,      only: domain_t
    use kryos_experiment
    implicit none
    private

    type, extends(ice_t), public :: toy_ice_t
    contains
        procedure :: init => toy_init
        procedure :: update => toy_update
        procedure :: restart_write => toy_rw
        procedure :: restart_read => toy_rw
    end type toy_ice_t

    type, extends(toy_ice_t), public :: other_ice_t       ! a different dynamic type
    end type other_ice_t

    ! recorder: what each hook call saw, in call order
    integer, parameter, public :: MAXLOG = 64
    integer,                 public :: n_log = 0
    character(len=NAME_LEN), public :: log_name(MAXLOG) = ""
    integer,                 public :: log_id(MAXLOG) = 0
    type(hook_time_t),       public :: log_st(MAXLOG)

    type, extends(experiment_t), public :: recorder_t     ! rec_a wants every hook, rec_b only ICE_BND, others none
    contains
        procedure :: init => rec_init
        procedure :: hook => rec_hook
    end type recorder_t

    type, extends(experiment_t), public :: need_toy_t     ! stops at init unless the ice is a toy_ice_t
    contains
        procedure :: init => need_init
        procedure :: hook => need_hook
    end type need_toy_t

    type, extends(experiment_t), public :: need_step_t    ! stops at init unless the ice runs every step
        real(dp) :: dt_seen = -1.0_dp                    ! ctx%dt_ice as init saw it
    contains
        procedure :: init => step_init
        procedure :: hook => step_hook
    end type need_step_t

    type, extends(experiment_t), public :: counter_t      ! stateful: n, in <name>_toy.dat
        integer :: n = 0
    contains
        procedure :: init => cnt_init
        procedure :: hook => cnt_hook
        procedure :: restart_write => cnt_write
        procedure :: restart_read  => cnt_read
        procedure :: restart_files => cnt_files
    end type counter_t

contains

    subroutine toy_init(self, ctx)
        class(toy_ice_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        associate (u => self%nx + ctx%time); end associate
    end subroutine toy_init
    subroutine toy_update(self, st)
        class(toy_ice_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => self%nx + st%t1); end associate
    end subroutine toy_update
    subroutine toy_rw(self, dir, t)
        class(toy_ice_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_rw

    subroutine rec_init(self, ctx, dom)
        class(recorder_t), intent(inout) :: self
        type(exp_ctx_t),   intent(in)    :: ctx
        type(domain_t),    intent(inout) :: dom
        associate (u => ctx%time + len(dom%name)); end associate
        select case (trim(self%name))
        case ("rec_a");  self%wants = .true.
        case ("rec_b");  self%wants(HOOK_ICE_BND) = .true.
        end select
    end subroutine rec_init
    subroutine rec_hook(self, dom, id, st)
        class(recorder_t), intent(inout) :: self
        type(domain_t),    intent(inout) :: dom
        integer,           intent(in)    :: id
        type(hook_time_t), intent(in)    :: st
        associate (u => len(dom%name)); end associate
        n_log = n_log + 1
        log_name(n_log) = self%name
        log_id(n_log)   = id
        log_st(n_log)   = st
    end subroutine rec_hook

    subroutine need_init(self, ctx, dom)
        class(need_toy_t), intent(inout) :: self
        type(exp_ctx_t),   intent(in)    :: ctx
        type(domain_t),    intent(inout) :: dom
        associate (u => ctx%time); end associate
        if (.not. allocated(dom%ice%c)) then
            call fesm_error("test_experiment::need_toy", "experiment "//trim(self%name)//" needs an ice instance", &
                            "domain "//trim(dom%name)//" has none")
        end if
        select type (c => dom%ice%c)
        type is (toy_ice_t)
            self%wants(HOOK_STEP_BEGIN) = .true.
        class default
            call fesm_error("test_experiment::need_toy", "experiment "//trim(self%name)//" needs a toy_ice_t", &
                            "domain "//trim(dom%name)//": the ice instance has another type")
        end select
    end subroutine need_init
    subroutine need_hook(self, dom, id, st)
        class(need_toy_t), intent(inout) :: self
        type(domain_t),    intent(inout) :: dom
        integer,           intent(in)    :: id
        type(hook_time_t), intent(in)    :: st
        associate (u => len(dom%name) + id + st%t1 + len_trim(self%name)); end associate
    end subroutine need_hook

    subroutine step_init(self, ctx, dom)
        class(need_step_t), intent(inout) :: self
        type(exp_ctx_t),    intent(in)    :: ctx
        type(domain_t),     intent(inout) :: dom
        associate (u => len(dom%name)); end associate
        self%dt_seen = ctx%dt_ice
        call experiment_require_ice_every_step(ctx, "the toy")
        self%wants(HOOK_STEP_BEGIN) = .true.
    end subroutine step_init
    subroutine step_hook(self, dom, id, st)
        class(need_step_t), intent(inout) :: self
        type(domain_t),     intent(inout) :: dom
        integer,            intent(in)    :: id
        type(hook_time_t),  intent(in)    :: st
        associate (u => len(dom%name) + id + st%t1 + len_trim(self%name)); end associate
    end subroutine step_hook

    subroutine cnt_init(self, ctx, dom)
        class(counter_t), intent(inout) :: self
        type(exp_ctx_t),  intent(in)    :: ctx
        type(domain_t),   intent(inout) :: dom
        associate (u => ctx%time + len(dom%name)); end associate
        self%restart = .true.
        self%wants(HOOK_STEP_BEGIN) = .true.
    end subroutine cnt_init
    subroutine cnt_hook(self, dom, id, st)
        class(counter_t),  intent(inout) :: self
        type(domain_t),    intent(inout) :: dom
        integer,           intent(in)    :: id
        type(hook_time_t), intent(in)    :: st
        associate (u => len(dom%name) + id + st%t1); end associate
        self%n = self%n + 1
    end subroutine cnt_hook
    function cnt_files(self) result(files)
        class(counter_t), intent(in) :: self
        character(len=PATH_LEN), allocatable :: files(:)
        files = [character(len=PATH_LEN) :: trim(self%name)//"_toy.dat"]
    end function cnt_files
    subroutine cnt_write(self, dir, t)
        class(counter_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        integer :: u
        associate (unused => t); end associate
        open(newunit=u, file=trim(dir)//"/"//trim(self%name)//"_toy.dat", status="replace", action="write")
        write(u,*) self%n
        close(u)
    end subroutine cnt_write
    subroutine cnt_read(self, dir, t)
        class(counter_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        integer :: u
        associate (unused => t); end associate
        open(newunit=u, file=trim(dir)//"/"//trim(self%name)//"_toy.dat", status="old", action="read")
        read(u,*) self%n
        close(u)
    end subroutine cnt_read

end module test_experiment_toys


program test_experiment
    ! kryos_experiment (core-design 3.9): the mechanism without a run, on a
    ! domain_t built by hand (toy ice instances) and a par file written to
    ! tmp/test_experiment/ (run from the repo root, as make check does).
    !
    !   * creation: the groups named by the domain's `experiments`, in order,
    !     the host prototypes selected by name (one recorder adopted three
    !     times), kind, domain, wants;
    !   * hooks: order per id, an experiment that does not want a hook is skipped,
    !     the hook_time_t arrives unchanged; an empty list does nothing;
    !   * the stops, each in a child process of this executable ("abort <case>"),
    !     with the message: an experiment on an ice of another type (select
    !     type), unknown kind, unknown key, a name listed twice, a host
    !     experiment without a prototype, a restart file missing from the bundle;
    !   * restart: a stateful toy writes the files restart_files declares and
    !     reads its state back; a stateless one writes nothing.

    use precision,       only: dp
    use fesm_clock,      only: clock_t, clock_init, clock_start, cadence_init
    use fesm_log,        only: log_set_level, LOG_WARN
    use phys_constants,  only: phys_const_class, phys_const_set
    use kryos_component, only: PATH_LEN
    use kryos_config,    only: config_t, config_init
    use kryos_domain,    only: domain_t
    use kryos_experiment
    use test_experiment_toys

    implicit none

    character(len=*), parameter :: dir = "tmp/test_experiment/"
    character(len=*), parameter :: par = dir//"par.nml"
    integer :: nfail = 0
    character(len=512) :: arg
    type(phys_const_class) :: cnst
    type(config_t) :: cfg
    type(clock_t), target :: clk

    call log_set_level(LOG_WARN)
    call phys_const_set(cnst, label="Earth", rho_ice=910.0_dp, rho_sw=1028.0_dp, rho_w=1000.0_dp, g=9.81_dp)
    call clock_init(clk, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")
    call clock_start(clk)

    call get_command_argument(1, arg)
    if (trim(arg) == "abort") then
        call get_command_argument(2, arg)
        call config_init(cfg, par)
        call abort_case(trim(arg))          ! must stop; falling through is a failure
        stop 0
    end if

    call execute_command_line("mkdir -p "//dir)
    call write_par()
    call config_init(cfg, par)

    call test_create_and_hooks()
    call test_restart()
    call test_stops()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_experiment: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_experiment: OK"

contains

    subroutine test_create_and_hooks()
        type(domain_t) :: d
        type(experiment_list_t) :: list
        type(experiment_proto_t), allocatable :: protos(:)
        type(hook_time_t) :: st
        integer :: id

        write(*,"(a)") "creation and hooks"
        d%name = "Greenland"
        protos = host_protos()
        call experiments_create(list, cfg, d, clk, cnst, protos)
        call check("three experiments", experiments_count(list) == 3)
        call check("in the order of the key", trim(list%s(1)%c%name) == "rec_a" .and. &
                   trim(list%s(2)%c%name) == "rec_b" .and. trim(list%s(3)%c%name) == "rec_c")
        call check("kind host, domain set by the factory", trim(list%s(1)%c%kind) == "host" .and. &
                   trim(list%s(1)%c%domain) == "Greenland")
        call check("wants set in init", all(list%s(1)%c%wants) .and. count(list%s(2)%c%wants) == 1 .and. &
                   list%s(2)%c%wants(HOOK_ICE_BND) .and. .not. any(list%s(3)%c%wants))

        st%t0 = 10.0_dp;  st%t1 = 20.0_dp;  st%t_eval = 1970.0_dp;  st%t_elapsed = 20.0_dp;  st%dtt = 10.0_dp
        st%clk => clk
        n_log = 0
        do id = 1, NHOOK
            call experiments_hook(list, d, id, st)
        end do
        call check("calls: rec_a at every hook, rec_b at ICE_BND only, rec_c never", n_log == 5)
        call check("order: STEP_BEGIN a; ICE_BND a, b; COLD_PRE a; COLD_POST a", &
                   trim(log_name(1)) == "rec_a" .and. log_id(1) == HOOK_STEP_BEGIN .and. &
                   trim(log_name(2)) == "rec_a" .and. log_id(2) == HOOK_ICE_BND .and. &
                   trim(log_name(3)) == "rec_b" .and. log_id(3) == HOOK_ICE_BND .and. &
                   trim(log_name(4)) == "rec_a" .and. log_id(4) == HOOK_COLD_PRE_STATE .and. &
                   trim(log_name(5)) == "rec_a" .and. log_id(5) == HOOK_COLD_POST_STATE)
        call check("hook_time_t arrives unchanged", log_st(3)%t0 == 10.0_dp .and. log_st(3)%t1 == 20.0_dp .and. &
                   log_st(3)%t_eval == 1970.0_dp .and. log_st(3)%t_elapsed == 20.0_dp .and. log_st(3)%dtt == 10.0_dp .and. &
                   associated(log_st(3)%clk))
        call experiments_describe(list)
        call experiments_end(list)

        d%name = "Quiet"                          ! no `experiments` key: an empty list
        call experiments_create(list, cfg, d, clk, cnst)
        call check("no experiments: empty list", experiments_count(list) == 0)
        n_log = 0
        call experiments_hook(list, d, HOOK_STEP_BEGIN, st)
        call check("no experiments: a hook does nothing", n_log == 0)

        d%name = "Toyland"                        ! a host experiment that finds its ice
        allocate(toy_ice_t :: d%ice%c)
        call experiments_create(list, cfg, d, clk, cnst, [proto("need_toy", "needs")])
        call check("select type on the right ice type: created, wants STEP_BEGIN", &
                   experiments_count(list) == 1 .and. list%s(1)%c%wants(HOOK_STEP_BEGIN))

        d%name = "StepEvery"                      ! the ice due every step: ctx%dt_ice = dtt, the guard passes
        call toy_ice_domain(d, "every")
        call experiments_create(list, cfg, d, clk, cnst, [proto("need_step", "needstep")])
        call check("ice at every step: the guard passes, ctx%dt_ice = dtt", experiments_count(list) == 1 .and. &
                   list%s(1)%c%wants(HOOK_STEP_BEGIN))
        select type (c => list%s(1)%c)
        type is (need_step_t)
            call check("ctx%dt_ice is the ice cadence step [yr]", c%dt_seen == 10.0_dp)
        end select
    end subroutine test_create_and_hooks

    subroutine test_restart()
        type(domain_t) :: d
        type(experiment_list_t) :: list, list2
        character(len=PATH_LEN), allocatable :: files(:)
        character(len=*), parameter :: bundle = dir//"bundle"
        logical :: found
        integer :: n

        write(*,"(a)") "restart"
        call execute_command_line("rm -rf "//bundle//" && mkdir -p "//bundle)
        d%name = "Antarctica"
        call experiments_create(list, cfg, d, clk, cnst, host_protos())
        call check("a stateful and a stateless experiment", experiments_count(list) == 2)
        select type (c => list%s(1)%c)
        type is (counter_t)
            c%n = 7
        end select
        call experiments_restart_write(list, bundle, 100.0_dp)

        files = experiments_restart_files(list, 1)
        inquire(file=bundle//"/cnt_toy.dat", exist=found)
        call check("stateful: restart_files names the file written", size(files) == 1 .and. &
                   trim(files(1)) == "cnt_toy.dat" .and. found)
        files = experiments_restart_files(list, 2)
        inquire(file=bundle//"/rec_a_restart.nc", exist=found)
        call check("stateless: no file declared, none written", size(files) == 0 .and. .not. found)

        call experiments_create(list2, cfg, d, clk, cnst, host_protos())
        n = -1
        select type (c => list2%s(1)%c)
        type is (counter_t)
            call experiments_restart_read(list2, bundle, 100.0_dp)
            n = c%n
        end select
        call check("read back into a new experiment", n == 7)
    end subroutine test_restart

    subroutine test_stops()
        write(*,"(a)") "stops"
        call stops("select_mismatch", "needs a toy_ice_t", "select type on another ice type")
        call stops("unknown_kind",    "has no template", "unknown kind (no template)", expect_text="nonexistent")
        call stops("unknown_key",     "bogus", "unknown key in the experiment group")
        call stops("duplicate",       "listed twice", "a name listed twice")
        call stops("no_proto",        "no host experiment named", "a host experiment without a prototype")
        call stops("missing_file",    "is not in the bundle", "a restart file missing from the bundle")
        call stops("dt_slow",         "needs the ice to run at every step", "the ice on a 20 yr cadence, dtt = 10")
        call stops("dt_off",          "needs the ice to run at every step", "the ice off")
        call stops("dt_none",         "needs the ice to run at every step", "no ice instance")
        call stops("ok",              "", "a valid case does not stop", ok_expected=.true.)
    end subroutine test_stops

    subroutine abort_case(name)
        ! Runs in the child process.
        character(len=*), intent(in) :: name
        type(domain_t) :: d
        type(experiment_list_t) :: list

        select case (name)
        case ("select_mismatch")
            d%name = "Toyland"
            allocate(other_ice_t :: d%ice%c)
            call experiments_create(list, cfg, d, clk, cnst, [proto("need_toy", "needs")])
        case ("ok")
            d%name = "Toyland"
            allocate(toy_ice_t :: d%ice%c)
            call experiments_create(list, cfg, d, clk, cnst, [proto("need_toy", "needs")])
        case ("dt_slow")
            d%name = "StepSlow"
            call toy_ice_domain(d, "slow")
            call experiments_create(list, cfg, d, clk, cnst, [proto("need_step", "needstep")])
        case ("dt_off")
            d%name = "StepOff"
            call toy_ice_domain(d, "off")
            call experiments_create(list, cfg, d, clk, cnst, [proto("need_step", "needstep")])
        case ("dt_none")
            d%name = "StepOff"
            call experiments_create(list, cfg, d, clk, cnst, [proto("need_step", "needstep")])
        case ("unknown_kind")
            d%name = "BadKind"
            call experiments_create(list, cfg, d, clk, cnst, host_protos())
        case ("unknown_key")
            d%name = "BadKey"
            call experiments_create(list, cfg, d, clk, cnst, host_protos())
        case ("duplicate")
            d%name = "Twice"
            call experiments_create(list, cfg, d, clk, cnst, host_protos())
        case ("no_proto")
            d%name = "Greenland"
            call experiments_create(list, cfg, d, clk, cnst, [proto("rec_a", "recorder")])
        case ("missing_file")
            d%name = "Antarctica"
            call experiments_create(list, cfg, d, clk, cnst, host_protos())
            call experiments_restart_read(list, dir//"empty", 100.0_dp)
        end select
    end subroutine abort_case

    subroutine toy_ice_domain(d, mode)
        ! A toy ice with a cadence of the given kind, registered on a clock of the
        ! same step that is not started (the run's clock is).
        type(domain_t),   intent(inout) :: d
        character(len=*), intent(in)    :: mode
        type(clock_t) :: c2

        if (allocated(d%ice%c)) deallocate(d%ice%c)
        allocate(toy_ice_t :: d%ice%c)
        call clock_init(c2, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")
        select case (mode)
        case ("every"); call cadence_init(d%ice%c%cadence, c2, "toy_ice", .true.)
        case ("slow");  call cadence_init(d%ice%c%cadence, c2, "toy_ice", 20.0_dp)
        case ("off");   call cadence_init(d%ice%c%cadence, c2, "toy_ice", .false.)
        end select
    end subroutine toy_ice_domain

    subroutine stops(name, text, label, expect_text, ok_expected)
        ! The child `abort <name>` exits non-zero and its output has `text`
        ! (and `expect_text`); ok_expected: exits zero.
        character(len=*), intent(in)           :: name, text, label
        character(len=*), intent(in), optional :: expect_text
        logical,          intent(in), optional :: ok_expected

        character(len=512) :: exe
        integer :: es, es2
        logical :: want_ok, good

        want_ok = .false.
        if (present(ok_expected)) want_ok = ok_expected
        call execute_command_line("mkdir -p "//dir//"empty")
        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" abort "//name//" > "//dir//name//".log 2>&1", exitstat=es)
        if (want_ok) then
            call check(label, es == 0)
            return
        end if
        good = es /= 0
        if (len_trim(text) > 0) then
            call execute_command_line("grep -q '"//text//"' "//dir//name//".log", exitstat=es2)
            good = good .and. es2 == 0
        end if
        if (present(expect_text)) then
            call execute_command_line("grep -q '"//expect_text//"' "//dir//name//".log", exitstat=es2)
            good = good .and. es2 == 0
        end if
        call check(label//" ("//dir//name//".log)", good)
    end subroutine stops

    ! -------------------------------------------------------------------------
    ! Prototypes and par file
    ! -------------------------------------------------------------------------

    function proto(name, what) result(p)
        ! A prototype named `name`: a recorder, the ice-needing experiment or the counter.
        character(len=*), intent(in) :: name, what
        type(experiment_proto_t) :: p
        p%name = name
        select case (what)
        case ("recorder"); allocate(recorder_t :: p%c)
        case ("needs");    allocate(need_toy_t :: p%c)
        case ("needstep"); allocate(need_step_t :: p%c)
        case ("counter");  allocate(counter_t :: p%c)
        end select
    end function proto

    function host_protos() result(protos)
        ! rec_a, rec_b, rec_c: one recorder type adopted three times; cnt: the counter.
        type(experiment_proto_t), allocatable :: protos(:)
        allocate(protos(4))
        protos(1) = proto("rec_a", "recorder")
        protos(2) = proto("rec_b", "recorder")
        protos(3) = proto("rec_c", "recorder")
        protos(4) = proto("cnt",   "counter")
    end function host_protos

    subroutine write_par()
        integer :: u

        open(newunit=u, file=par, status="replace", action="write")
        write(u,'(a)') '&kryos', '    domains = "Greenland"', '/'
        write(u,'(a)') '&Greenland', '    experiments = "rec_a" "rec_b" "rec_c"', '/'
        write(u,'(a)') '&Antarctica', '    experiments = "cnt" "rec_a"', '/'
        write(u,'(a)') '&Quiet', '/'
        write(u,'(a)') '&Toyland', '    experiments = "need_toy"', '/'
        write(u,'(a)') '&BadKind', '    experiments = "exp_unknown"', '/'
        write(u,'(a)') '&BadKey', '    experiments = "exp_badkey"', '/'
        write(u,'(a)') '&Twice', '    experiments = "rec_a" "rec_a"', '/'
        write(u,'(a)') '&rec_a', '    kind = "host"', '/'
        write(u,'(a)') '&rec_b', '    kind = "host"', '/'
        write(u,'(a)') '&rec_c', '    kind = "host"', '/'
        write(u,'(a)') '&cnt', '    kind = "host"', '/'
        write(u,'(a)') '&need_toy', '    kind = "host"', '/'
        write(u,'(a)') '&StepEvery', '    experiments = "need_step"', '/'
        write(u,'(a)') '&StepSlow', '    experiments = "need_step"', '/'
        write(u,'(a)') '&StepOff', '    experiments = "need_step"', '/'
        write(u,'(a)') '&need_step', '    kind = "host"', '/'
        write(u,'(a)') '&exp_unknown', '    kind = "nonexistent"', '/'
        write(u,'(a)') '&exp_badkey', '    kind  = "host"', '    bogus = 1.0', '/'
        close(u)
    end subroutine write_par

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

end program test_experiment
