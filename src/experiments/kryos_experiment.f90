module kryos_experiment
    ! The experiment plugin mechanism (phase 2 WP1, core-design section 3.9,
    ! design decision 1): an experiment sits between the domain and the run, takes the
    ! domain_t and calls a hook at one of four fixed points of the run sequence
    ! (HOOK_*; the sequence is kryos_run_sequence). It reaches backend internals
    ! only by `select type` on the domain's instances and stops on a mismatch;
    ! the core never names an experiment type.
    !
    ! Same pattern as the components: the abstract experiment_t with deferred
    ! bindings called from one place each, the free-procedure dispatchers below
    ! (experiments_*), a factory by `kind` (the submodule kryos_experiment_factory,
    ! the only place that names the concrete kinds) and an adopt for experiments
    ! a host program provides (kind = "host": a prototype of the group's name,
    ! passed to run_init). An experiment is configured by one par-file group,
    ! named by the domain's `experiments` key (one list per domain; the order is
    ! the call order), template [experiment] + [experiment_<kind>].
    !
    ! The experiments of a domain are held in run_t, one experiment_list_t per
    ! domain (not in domain_t: hook takes a domain_t). Stateless by default; a
    ! stateful one sets `restart` and overrides the restart bindings.

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, log_msg, log_set_prefix, fesm_error, to_str
    use fesm_clock,      only: clock_t, clock_time, clock_time_rel, clock_dtt
    use phys_constants,  only: phys_const_class
    use kryos_component, only: NAME_LEN, PATH_LEN, NLIB
    use kryos_config,    only: config_t, config_read, config_validate_experiment
    use kryos_domain,    only: domain_t

    implicit none

    private

    integer, parameter, public :: HOOK_STEP_BEGIN = 1, HOOK_ICE_BND = 2, &
                                  HOOK_COLD_PRE_STATE = 3, HOOK_COLD_POST_STATE = 4, NHOOK = 4
    character(len=24), parameter, public :: HOOK_NAME(NHOOK) = [character(len=24) :: &
        "STEP_BEGIN", "ICE_BND", "COLD_PRE_STATE", "COLD_POST_STATE"]

    integer, parameter :: MAX_EXPERIMENTS = 16      ! per domain

    type, public :: exp_ctx_t                       ! like init_ctx_t; built by experiments_create
        character(len=PATH_LEN) :: par_file = "", library(NLIB) = ""
        character(len=NAME_LEN) :: instance = "", kind = "", domain = ""   ! instance = the experiment's group
        type(phys_const_class)  :: cnst
        real(dp) :: time = 0.0_dp, time_rel = 0.0_dp, dtt = 0.0_dp, sec_year = 0.0_dp   ! at t0
    end type exp_ctx_t

    type, public :: hook_time_t                     ! the `st` of a hook (the RUN's step)
        real(dp) :: t0 = 0.0_dp, t1 = 0.0_dp        ! model-time interval [yr]; cold start: t0 = t1 = time_init
        real(dp) :: t_eval = 0.0_dp                 ! clock_time_rel at t1
        real(dp) :: t_elapsed = 0.0_dp              ! clock_time_elapsed at t1 (window time); cold start 0
        real(dp) :: dtt = 0.0_dp                    ! run step [yr]
        type(clock_t), pointer :: clk => null()     ! read-only
    end type hook_time_t

    type, abstract, public :: experiment_t
        character(len=NAME_LEN) :: name = "", kind = "", domain = ""
        logical :: wants(NHOOK) = .false.           ! set in init; the dispatcher skips the other hooks
        logical :: restart = .false.                ! stateful: kryos_restart keeps its files
    contains
        procedure(x_init), deferred :: init         ! (ctx, dom): read the group, check the domain's instances, set wants
        procedure(x_hook), deferred :: hook         ! (dom, id, st): select case (id) inside
        procedure :: restart_write => experiment_restart_write    ! (dir, t): default nothing
        procedure :: restart_read  => experiment_restart_read
        procedure :: restart_files => experiment_restart_files    ! () -> base names restart_write creates
        procedure :: end           => experiment_end
    end type experiment_t

    abstract interface
        subroutine x_init(self, ctx, dom)
            import :: experiment_t, exp_ctx_t, domain_t
            class(experiment_t), intent(inout) :: self
            type(exp_ctx_t),     intent(in)    :: ctx
            type(domain_t),      intent(in)    :: dom    ! keep no pointer into it
        end subroutine x_init

        subroutine x_hook(self, dom, id, st)
            import :: experiment_t, domain_t, hook_time_t
            class(experiment_t), intent(inout) :: self
            type(domain_t),      intent(inout) :: dom
            integer,             intent(in)    :: id
            type(hook_time_t),   intent(in)    :: st
        end subroutine x_hook
    end interface

    type, public :: experiment_slot_t
        class(experiment_t), allocatable :: c
    end type experiment_slot_t

    type, public :: experiment_list_t               ! the experiments of one domain, in call order
        type(experiment_slot_t), allocatable :: s(:)
    end type experiment_list_t

    type, public :: experiment_proto_t              ! a host's experiment (kind = "host"), selected by name
        character(len=NAME_LEN) :: name = ""
        class(experiment_t), allocatable :: c
    end type experiment_proto_t

    interface                                       ! factory: submodule kryos_experiment_factory
        module subroutine experiment_create(slot, ctx, protos)
            type(experiment_slot_t),  intent(out) :: slot
            type(exp_ctx_t),          intent(in)  :: ctx      ! instance, kind, domain
            type(experiment_proto_t), intent(in)  :: protos(:)
        end subroutine experiment_create
    end interface

    public :: experiment_create, experiment_adopt
    public :: experiments_create, experiments_hook, experiments_count, experiments_describe, experiments_end
    public :: experiments_restart_write, experiments_restart_read, experiments_restart_files

contains

    ! =========================================================================
    ! Default bindings
    ! =========================================================================

    subroutine experiment_restart_write(self, dir, t)
        ! Default: nothing (stateless).
        class(experiment_t), intent(inout) :: self
        character(len=*),    intent(in)    :: dir
        real(dp),            intent(in)    :: t
        associate (unused => t + len(dir)); end associate
        associate (unused => self%restart); end associate
    end subroutine experiment_restart_write

    subroutine experiment_restart_read(self, dir, t)
        class(experiment_t), intent(inout) :: self
        character(len=*),    intent(in)    :: dir
        real(dp),            intent(in)    :: t
        associate (unused => t + len(dir)); end associate
        associate (unused => self%restart); end associate
    end subroutine experiment_restart_read

    function experiment_restart_files(self) result(files)
        ! The base names (no folder) of the files restart_write creates:
        ! <name>_restart.nc for a stateful experiment, none otherwise. An
        ! experiment whose files differ overrides it, in step with its
        ! restart_write.
        class(experiment_t), intent(in) :: self
        character(len=PATH_LEN), allocatable :: files(:)

        if (self%restart) then
            files = [character(len=PATH_LEN) :: trim(self%name)//"_restart.nc"]
        else
            allocate(files(0))
        end if
    end function experiment_restart_files

    subroutine experiment_end(self)
        ! Default: nothing to release.
        class(experiment_t), intent(inout) :: self
        associate (unused => self%restart); end associate
    end subroutine experiment_end

    ! =========================================================================
    ! Factory support
    ! =========================================================================

    subroutine experiment_adopt(slot, proto, ctx)
        ! A host-provided extension of experiment_t, identified as the factory
        ! would (same shape as ice_adopt).
        type(experiment_slot_t),     intent(out) :: slot
        class(experiment_t), allocatable, intent(in) :: proto
        type(exp_ctx_t),             intent(in)  :: ctx
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(slot%c, mold=proto)
        slot%c%name   = ctx%instance
        slot%c%kind   = ctx%kind
        slot%c%domain = ctx%domain
    end subroutine experiment_adopt

    ! =========================================================================
    ! Dispatchers: the only call sites of the bindings
    ! =========================================================================

    subroutine experiments_create(list, cfg, dom, clk, cnst, protos)
        ! The experiments of domain `dom`: the groups named by its `experiments`
        ! key, in order, each validated (template [experiment] + [experiment_<kind>]),
        ! made by the factory and initialised. Names are unique in the domain
        ! (kryos_run checks them against the instances). `protos`: the host's
        ! experiments of kind "host", selected by the group name.
        type(experiment_list_t),  intent(out)          :: list
        type(config_t),           intent(in)           :: cfg
        type(domain_t),           intent(in)           :: dom
        type(clock_t),            intent(in)           :: clk
        type(phys_const_class),   intent(in)           :: cnst
        type(experiment_proto_t), intent(in), optional :: protos(:)

        character(len=*), parameter :: proc = "kryos_experiment::experiments_create"
        character(len=NAME_LEN) :: names(MAX_EXPERIMENTS), kind
        type(experiment_proto_t), allocatable :: none(:)
        type(exp_ctx_t) :: ctx
        integer :: i, j, n

        names = ""
        call config_read(cfg, dom%name, "experiments", names, "domain")
        n = count(len_trim(names) > 0)
        allocate(list%s(n))
        if (n == 0) return

        if (.not. present(protos)) allocate(none(0))
        do i = 1, n
            do j = 1, i - 1
                if (trim(names(j)) == trim(names(i))) then
                    call fesm_error(proc, "domain &"//trim(dom%name)//": experiment `"//trim(names(i))// &
                                    "` is listed twice in `experiments`", "par file "//trim(cfg%par_file))
                end if
            end do
            call config_validate_experiment(cfg, trim(names(i)), kind)

            ctx%par_file = cfg%par_file
            ctx%library  = cfg%library
            ctx%instance = names(i)
            ctx%kind     = kind
            ctx%domain   = dom%name
            ctx%cnst     = cnst
            ctx%time     = clock_time(clk)
            ctx%time_rel = clock_time_rel(clk)
            ctx%dtt      = clock_dtt(clk)
            ctx%sec_year = clk%sec_year
            if (present(protos)) then
                call experiment_create(list%s(i), ctx, protos)
            else
                call experiment_create(list%s(i), ctx, none)
            end if

            call log_set_prefix(trim(names(i)))
            call list%s(i)%c%init(ctx, dom)
            call log_set_prefix("")
        end do
    end subroutine experiments_create

    integer function experiments_count(list) result(n)
        type(experiment_list_t), intent(in) :: list
        n = 0
        if (allocated(list%s)) n = size(list%s)
    end function experiments_count

    subroutine experiments_hook(list, dom, id, st)
        ! Hook `id` of every experiment of the list that wants it, in order.
        type(experiment_list_t), intent(inout) :: list
        type(domain_t),          intent(inout) :: dom
        integer,                 intent(in)    :: id
        type(hook_time_t),       intent(in)    :: st

        integer :: i

        if (.not. allocated(list%s)) return
        do i = 1, size(list%s)
            if (.not. list%s(i)%c%wants(id)) cycle
            call log_set_prefix(trim(list%s(i)%c%name))
            call list%s(i)%c%hook(dom, id, st)
            call log_set_prefix("")
        end do
    end subroutine experiments_hook

    subroutine experiments_restart_write(list, dir, t)
        ! restart_write of every stateful experiment into the bundle folder;
        ! each file it declares (restart_files) must exist afterwards.
        type(experiment_list_t), intent(inout) :: list
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t

        character(len=PATH_LEN), allocatable :: files(:)
        logical :: found
        integer :: i, j

        if (.not. allocated(list%s)) return
        do i = 1, size(list%s)
            if (.not. list%s(i)%c%restart) cycle
            call log_set_prefix(trim(list%s(i)%c%name))
            call list%s(i)%c%restart_write(dir, t)
            call log_set_prefix("")
            files = experiments_restart_files(list, i)
            do j = 1, size(files)
                inquire(file=trim(dir)//"/"//trim(files(j)), exist=found)
                if (found) cycle
                call fesm_error("kryos_experiment::experiments_restart_write", "experiment `"// &
                                trim(list%s(i)%c%name)//"` declared the restart file "//trim(files(j))// &
                                " but did not write it", "bundle "//trim(dir))
            end do
        end do
    end subroutine experiments_restart_write

    subroutine experiments_restart_read(list, dir, t)
        ! restart_read of every stateful experiment; every file it declares
        ! must be in the bundle folder.
        type(experiment_list_t), intent(inout) :: list
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t

        character(len=PATH_LEN), allocatable :: files(:)
        logical :: found
        integer :: i, j

        if (.not. allocated(list%s)) return
        do i = 1, size(list%s)
            if (.not. list%s(i)%c%restart) cycle
            files = experiments_restart_files(list, i)
            do j = 1, size(files)
                inquire(file=trim(dir)//"/"//trim(files(j)), exist=found)
                if (found) cycle
                call fesm_error("kryos_experiment::experiments_restart_read", "restart file "//trim(files(j))// &
                                " of experiment `"//trim(list%s(i)%c%name)//"` is not in the bundle", &
                                "bundle "//trim(dir))
            end do
            call log_set_prefix(trim(list%s(i)%c%name))
            call list%s(i)%c%restart_read(dir, t)
            call log_set_prefix("")
        end do
    end subroutine experiments_restart_read

    function experiments_restart_files(list, i) result(files)
        ! The base names the restart_write of experiment i creates.
        type(experiment_list_t), intent(in) :: list
        integer,                 intent(in) :: i
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => list%s(i)%c)                    ! (as ice_restart_files: no call through the slot itself)
            files = c%restart_files()
        end associate
    end function experiments_restart_files

    subroutine experiments_describe(list)
        ! core-design 8: one line per experiment, its kind and the hooks it wants.
        type(experiment_list_t), intent(in) :: list

        character(len=:), allocatable :: hooks
        integer :: i, h

        if (.not. allocated(list%s)) return
        do i = 1, size(list%s)
            hooks = ""
            do h = 1, NHOOK
                if (list%s(i)%c%wants(h)) hooks = hooks//" "//trim(HOOK_NAME(h))
            end do
            if (len(hooks) == 0) hooks = " (none)"
            call log_msg(LOG_INFO, "kryos_experiment", "experiment "//trim(list%s(i)%c%name)//" ("// &
                         trim(list%s(i)%c%kind)//", domain "//trim(list%s(i)%c%domain)//"), hooks:"//hooks)
        end do
    end subroutine experiments_describe

    subroutine experiments_end(list)
        type(experiment_list_t), intent(inout) :: list

        integer :: i

        if (.not. allocated(list%s)) return
        do i = 1, size(list%s)
            call log_set_prefix(trim(list%s(i)%c%name))
            call list%s(i)%c%end()
            call log_set_prefix("")
        end do
    end subroutine experiments_end

end module kryos_experiment
