module kryos_globe
    ! The globe (core-design sections 1, 3.1 step 3, 3.4; contracts-design 4.2,
    ! 6, 7): the instances that belong to the run, not to a domain: the
    ! forcing (tsgen: coefficients on the forcing time) and the global part of
    ! geo (fastisostasy bsl: sea level and ocean area, the reduction of the
    ! domains' volume changes). A global geo that is resolved on the host's loads
    ! (opt%host, vilma2) has its own grid: the globe owns the coupler of the global grids
    ! (cpl: that grid, the background grid, the domains' ice grids) and the global ice
    ! field gi (kryos_globe_ice); the domains' couplers get the grid by
    ! globe_register_grid.
    !
    ! Plain type, free procedures. The globe never sees a domain: the
    ! exchanges are helpers on the rows a domain resolved (wire_resolve with
    ! the globe nodes, globe_nodes) and on the component field states, so the
    ! run sequence calls them once per domain:
    !
    !     exchange_globe    globe export -> a domain instance's import (bsl, A_ocn,
    !                       coefficients); rows kind globe -> kind
    !     gather_into       domain instances' 0-D exports summed into a globe
    !                       import (dV_bsl, M_ice), in the order of the sources
    !     exchange_to_globe a domain instance's 0-D export -> a globe import
    !                       (a derived coefficient, one producer)
    !
    ! plus the globe's own step: globe_update (bsl_update, coefs_update on the
    ! instance's cadence and clock), globe_reduce (geo_reduce at the time of
    ! the last update) and globe_init / globe_init_state / globe_end; for the
    ! restart bundles (kryos_restart) one instance's init_state (globe_pass),
    ! its restart files (globe_restart_write/read) and its component
    ! (globe_component). The
    ! wrappers of the components (geo_reduce, geo_scatter, forcing_set,
    ! forcing_kill) stay the call sites of the bindings; geo_scatter acts on a
    ! domain's geo slot and is called by the domain sequence.
    !
    ! Kinds: "forcing" and "geo" name the globe side of a row; "geo" is the
    ! table kind "geo_g" (kryos_wiring).

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t
    use fesm_clock,      only: clock_t, cadence_due, cadence_mark, cadence_t_last, cadence_dt_years, &
                               clock_time, clock_time_rel, clock_time_init, clock_time_end, clock_dtt, &
                               clock_to_years
    use coupler,         only: coupler_class, coupler_init, coupler_add_grid
    use phys_constants,  only: phys_const_class
    use kryos_config,    only: config_t, config_read, config_validate, config_validate_instance, config_cadence
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, step_time_of, NAME_LEN, PATH_LEN
    use kryos_wiring,    only: wire_t, wire_node_t, wire_node, wire_exchange, wire_gather, wire_runs
    use kryos_globe_ice, only: globe_ice_t, globe_ice_init, globe_ice_end, GI_NAME_LEN
    use kryos_forcing,   only: forcing_slot_t, forcing_create, forcing_init, forcing_init_state, &
                               forcing_update, forcing_end, forcing_present, forcing_restart_write, &
                               forcing_restart_read, forcing_restart_files
    use kryos_geo,       only: geo_slot_t, geo_create, geo_init, geo_init_state, geo_update, &
                               geo_reduce, geo_end, geo_present, geo_restart_write, geo_restart_read, &
                               geo_restart_files

    implicit none

    private

    type, public :: globe_t
        type(forcing_slot_t) :: forcing                 ! tsgen: coefficients on the forcing time
        type(geo_slot_t)     :: geo                     ! global part: fastisostasy bsl (scope "global")
        character(len=8)     :: clock_forcing = "model" ! the instances' `clock` key: model | forcing
        character(len=8)     :: clock_geo     = "model"
        real(dp)             :: geo_dt = 0.0_dp, geo_phase = 0.0_dp   ! the geo's dt and phase keys (the regional ports take them)
        type(coupler_class)  :: cpl                     ! global grids by name, map cache (a host-resolved geo)
        type(globe_ice_t)    :: gi                      ! the global ice field (a host-resolved geo)
        logical              :: has_gi = .false.        ! the global geo is host-resolved: cpl and gi are in use
    end type globe_t

    type, public :: gather_src_t                        ! one domain's contribution to a gather
        logical           :: active = .false.           ! the domain has a resolved, active gather row
        type(wire_t)      :: w                          ! that row (rank 0: no workspace)
        type(wire_node_t) :: node                       ! the producing instance
    end type gather_src_t

    public :: globe_init, globe_init_state, globe_update, globe_time, globe_reduce, globe_end
    public :: globe_pass, globe_component, globe_restart_write, globe_restart_read, globe_restart_files
    public :: globe_nodes, globe_register_grid
    public :: exchange_globe, exchange_to_globe, gather_src, gather_into

contains

    ! =========================================================================
    ! Init
    ! =========================================================================

    subroutine globe_init(globe, cfg, clk, cnst, census)
        ! The [kryos] globe group (template [globe]) names the instances; each
        ! present one, forcing then geo global (yelmox: bsl_init + bsl_update,
        ! then coefs_init): validate its group, factory (admission), init with
        ! its init_ctx_t, cadence_init on the clock from its keys dt, phase
        ! (dt = 0: every step; dt < 0: off) and clock (model | forcing). A key
        ! left "" in [globe] means no such instance. Call before clock_start.
        type(globe_t),         intent(out)   :: globe
        type(config_t),        intent(in)    :: cfg
        type(clock_t),         intent(inout) :: clk
        type(phys_const_class), intent(in)   :: cnst
        type(census_t),        intent(inout) :: census

        character(len=NAME_LEN) :: name, backend
        type(init_ctx_t) :: ctx

        call config_validate(cfg, cfg%globe, "globe")

        name = ""
        call config_read(cfg, cfg%globe, "forcing", name, "globe")
        if (len_trim(name) > 0) then
            call config_validate_instance(cfg, name, "forcing", backend)
            call globe_ctx(ctx, cfg, name, "forcing", backend, clk, cnst)
            call forcing_create(globe%forcing, ctx, census)
            call forcing_init(globe%forcing, ctx)
            call config_cadence(cfg, name, "forcing", backend, clk, globe%forcing%c%cadence, globe%clock_forcing)
            call log_msg(LOG_INFO, "kryos_globe", "forcing `"//trim(name)//"` ("//trim(backend)//"), clock "// &
                         trim(globe%clock_forcing)//", "//cadence_word(cadence_dt_years(globe%forcing%c%cadence, clk)))
        end if

        name = ""
        call config_read(cfg, cfg%globe, "geo", name, "globe")
        if (len_trim(name) > 0) then
            call config_validate_instance(cfg, name, "geo", backend)
            call globe_ctx(ctx, cfg, name, "geo", backend, clk, cnst)
            call geo_create(globe%geo, ctx, census)
            call geo_init(globe%geo, ctx)
            call config_cadence(cfg, name, "geo", backend, clk, globe%geo%c%cadence, globe%clock_geo, &
                                dt_key=globe%geo_dt, phase_key=globe%geo_phase)
            call log_msg(LOG_INFO, "kryos_globe", "geo `"//trim(name)//"` ("//trim(backend)//", global), clock "// &
                         trim(globe%clock_geo)//", "//cadence_word(cadence_dt_years(globe%geo%c%cadence, clk)))
            if (globe%geo%c%opt%host) call globe_gia_init(globe, cfg, clk, cnst, name, backend)
        end if
    end subroutine globe_init

    subroutine globe_gia_init(globe, cfg, clk, cnst, name, backend)
        ! The global ice field of a host-resolved geo (phase4-plan 3.4, 8.2): the coupler of the global grids
        ! (its map folder is [kryos] map_fldr) and the group the instance's `ice_bg` names ([globe_ice]):
        ! the datasrc of the background and of the external history, on the clock of the run (calendar year of
        ! model time zero, the run's span). The domains' ice grids join at gia_reference.
        type(globe_t),          intent(inout) :: globe
        type(config_t),         intent(in)    :: cfg
        type(clock_t),          intent(in)    :: clk
        type(phys_const_class), intent(in)    :: cnst
        character(len=*),       intent(in)    :: name, backend

        character(len=GI_NAME_LEN) :: group

        if (.not. allocated(globe%geo%c%grid_obj)) then
            call fesm_error("kryos_globe::globe_init", "host-resolved geo `"//trim(name)//"` has no grid object", &
                            "its init must set grid_obj with grid, nx and ny (opt%host)")
        end if
        group = "globe_ice"
        call config_read(cfg, name, "ice_bg", group, "instance", "geo_"//trim(backend))
        call coupler_init(globe%cpl, trim(cfg%map_fldr))
        call globe_ice_init(globe%gi, cfg, globe%cpl, cnst, globe%geo%c%grid_obj, group=trim(group), &
                            time_ref=clock_to_years(clk, clk%n_ref), time_span=[clock_time_init(clk), clock_time_end(clk)])
        globe%has_gi = .true.
    end subroutine globe_gia_init

    subroutine globe_register_grid(globe, cpl)
        ! The global geo's own grid (a host-resolved geo: GAUSS-L<lmax>) into the coupler of a domain, which
        ! resolves the rows geo_g -> geo on it; no-op without one.
        type(globe_t),       intent(in)    :: globe
        type(coupler_class), intent(inout) :: cpl
        if (.not. geo_present(globe%geo)) return
        if (.not. allocated(globe%geo%c%grid_obj)) return
        call coupler_add_grid(cpl, trim(globe%geo%c%grid), globe%geo%c%grid_obj)
    end subroutine globe_register_grid

    subroutine globe_ctx(ctx, cfg, name, comp, backend, clk, cnst)
        ! The init_ctx_t of a global instance (no domain, no grid), at t0.
        type(init_ctx_t),       intent(out) :: ctx
        type(config_t),         intent(in)  :: cfg
        character(len=*),       intent(in)  :: name, comp, backend
        type(clock_t),          intent(in)  :: clk
        type(phys_const_class), intent(in)  :: cnst

        ctx%par_file  = cfg%par_file
        ctx%library   = cfg%library
        ctx%instance  = name
        ctx%component = comp
        ctx%backend   = backend
        ctx%domain    = ""
        ctx%grid_name = ""
        ctx%scope     = "global"
        ctx%map_fldr  = cfg%map_fldr
        ctx%cnst      = cnst
        ctx%time      = clock_time(clk)
        ctx%time_rel  = clock_time_rel(clk)
        ctx%dtt       = clock_dtt(clk)
        ctx%sec_year  = clk%sec_year
    end subroutine globe_ctx

    function cadence_word(dt) result(s)
        real(dp), intent(in) :: dt
        character(len=:), allocatable :: s
        if (dt > 0.0_dp) then
            s = "every "//to_str(dt)//" yr"
        else
            s = "off"
        end if
    end function cadence_word

    subroutine globe_init_state(globe, clk, pass)
        ! init_state of every present instance, forcing then geo (core-design
        ! 3.1 step 6): the exports as initialised at t0 (PASS_INIT).
        type(globe_t), intent(inout) :: globe
        type(clock_t), intent(in)    :: clk
        integer,       intent(in)    :: pass

        if (forcing_present(globe%forcing)) call globe_pass(globe, "forcing", clk, pass)
        if (geo_present(globe%geo))         call globe_pass(globe, "geo", clk, pass)
    end subroutine globe_init_state

    subroutine globe_pass(globe, which, clk, pass)
        ! init_state of one instance at the clock's time (t0 = t1 = now):
        ! PASS_INIT, or PASS_RESTART after its restore (kryos_restart).
        type(globe_t),    intent(inout) :: globe
        character(len=*), intent(in)    :: which           ! "forcing" | "geo"
        type(clock_t),    intent(in)    :: clk
        integer,          intent(in)    :: pass

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            call forcing_init_state(globe%forcing, globe_time(globe, which, clk, init=.true.), pass)
        case ("geo")
            call require(geo_present(globe%geo), which)
            call geo_init_state(globe%geo, globe_time(globe, which, clk, init=.true.), pass)
        case default
            call unknown(which)
        end select
    end subroutine globe_pass

    function globe_component(globe, which) result(c)
        ! The instance `which` ("forcing" | "geo") as a component, or a
        ! disassociated pointer when the globe has none.
        type(globe_t), target, intent(in) :: globe
        character(len=*),      intent(in) :: which
        class(component_t), pointer :: c

        c => null()
        select case (trim(which))
        case ("forcing");  if (forcing_present(globe%forcing)) c => globe%forcing%c
        case ("geo");      if (geo_present(globe%geo))         c => globe%geo%c
        case default;      call unknown(which)
        end select
    end function globe_component

    subroutine globe_restart_write(globe, which, dir, t)
        ! <comp>_restart_write of one instance into the folder `dir` at t (the
        ! backend's state; the field states are the core's, kryos_restart).
        type(globe_t),    intent(inout) :: globe
        character(len=*), intent(in)    :: which           ! "forcing" | "geo"
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            call forcing_restart_write(globe%forcing, dir, t)
        case ("geo")
            call require(geo_present(globe%geo), which)
            call geo_restart_write(globe%geo, dir, t)
        case default
            call unknown(which)
        end select
    end subroutine globe_restart_write

    subroutine globe_restart_read(globe, which, dir, t)
        ! <comp>_restart_read of one instance from the folder `dir`.
        type(globe_t),    intent(inout) :: globe
        character(len=*), intent(in)    :: which           ! "forcing" | "geo"
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            call forcing_restart_read(globe%forcing, dir, t)
        case ("geo")
            call require(geo_present(globe%geo), which)
            call geo_restart_read(globe%geo, dir, t)
        case default
            call unknown(which)
        end select
    end subroutine globe_restart_read

    function globe_restart_files(globe, which) result(files)
        ! <comp>_restart_files of one instance: the base names its restart_write
        ! creates in the folder it is given.
        type(globe_t),    intent(in) :: globe
        character(len=*), intent(in) :: which              ! "forcing" | "geo"
        character(len=PATH_LEN), allocatable :: files(:)

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            files = forcing_restart_files(globe%forcing)
        case ("geo")
            call require(geo_present(globe%geo), which)
            files = geo_restart_files(globe%geo)
        case default
            call unknown(which)
            allocate(files(0))
        end select
    end function globe_restart_files

    subroutine globe_end(globe)
        ! forcing, then geo (run_end).
        type(globe_t), intent(inout) :: globe
        if (forcing_present(globe%forcing)) call forcing_end(globe%forcing)
        if (geo_present(globe%geo)) call geo_end(globe%geo)
        if (globe%has_gi) call globe_ice_end(globe%gi)
    end subroutine globe_end

    function globe_nodes(globe, clk) result(nodes)
        ! The wire nodes of the present instances (forcing, geo_g), unchecked
        ! (the core-filled and gathered imports of a global node are not any one
        ! domain's to wire), except the host-resolved geo: its imports all have a
        ! producer, a row or the core (wire_resolve's rule). A domain's
        ! wire_resolve gets them next to its own nodes. The nodes point into the
        ! globe's field states: `globe` must be a target at the caller, and live
        ! as long as the rows.
        type(globe_t), target, intent(in) :: globe
        type(clock_t),         intent(in) :: clk
        type(wire_node_t), allocatable :: nodes(:)

        allocate(nodes(0))
        if (forcing_present(globe%forcing)) &
            nodes = [nodes, wire_node(globe%forcing%c, cadence_dt_years(globe%forcing%c%cadence, clk), check=.false.)]
        if (geo_present(globe%geo)) &
            nodes = [nodes, wire_node(globe%geo%c, cadence_dt_years(globe%geo%c%cadence, clk), &
                                      check=globe%geo%c%opt%host)]
    end function globe_nodes

    ! =========================================================================
    ! The globe's own step
    ! =========================================================================

    function globe_time(globe, which, clk, init) result(st)
        ! The step_time_t of an instance: its cadence interval [t0, t1] (init:
        ! t0 = t1 = now), t_eval = t1, or the forcing time of the clock for an
        ! instance with clock = "forcing". After the mark of the update the
        ! interval is [t, t] (the time globe_reduce works at).
        type(globe_t),    intent(in) :: globe
        character(len=*), intent(in) :: which           ! "forcing" | "geo"
        type(clock_t),    intent(in) :: clk
        logical,          intent(in), optional :: init
        type(step_time_t) :: st

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            st = step_time_of(globe%forcing%c%cadence, clk, globe%clock_forcing, init)
        case ("geo")
            call require(geo_present(globe%geo), which)
            st = step_time_of(globe%geo%c%cadence, clk, globe%clock_geo, init)
        case default
            call unknown(which)
        end select

    end function globe_time

    subroutine globe_update(globe, which, clk)
        ! One instance's update when its cadence is due, then its mark:
        ! "geo" = bsl_update (A_ocean from bsl) at the start of the step,
        ! "forcing" = coefs_update (after gather_into of M_ice and the derived
        ! rows, exchange_to_globe, of the previous step).
        type(globe_t),    intent(inout) :: globe
        character(len=*), intent(in)    :: which           ! "forcing" | "geo"
        type(clock_t),    intent(in)    :: clk

        type(step_time_t) :: st

        select case (trim(which))
        case ("forcing")
            call require(forcing_present(globe%forcing), which)
            if (.not. cadence_due(globe%forcing%c%cadence, clk)) return
            st = globe_time(globe, which, clk)
            call forcing_update(globe%forcing, st)
            call cadence_mark(globe%forcing%c%cadence, clk)
        case ("geo")
            call require(geo_present(globe%geo), which)
            if (.not. cadence_due(globe%geo%c%cadence, clk)) return
            st = globe_time(globe, which, clk)
            call geo_update(globe%geo, st)
            call cadence_mark(globe%geo%c%cadence, clk)
        case default
            call unknown(which)
        end select
    end subroutine globe_update

    subroutine globe_reduce(globe, clk)
        ! geo_reduce of the global part after gather_into(dV_bsl): bsl -=
        ! sum(dV)/A_ocn, bsl put again at the time of the last update (the
        ! instance's cadence must be due every step, or have run at this
        ! time). At the cold start: after every domain's geo PASS_GEOMETRY.
        type(globe_t), intent(inout) :: globe
        type(clock_t), intent(in)    :: clk

        call require(geo_present(globe%geo), "geo")
        call geo_reduce(globe%geo, globe_time(globe, "geo", clk))
    end subroutine globe_reduce

    ! =========================================================================
    ! Exchanges
    ! =========================================================================

    subroutine exchange_globe(globe, src_kind, wires, dst_kind, dst, clk, init, cpl)
        ! The domain's resolved rows globe `src_kind` -> `dst_kind`, in table
        ! order: the published globe exports into the destination instance's
        ! imports `dst` (its field state). The lag check reads the globe
        ! instance's cadence (t_last, due) and the clock; init: cold start.
        ! Rows are 0-D (every export of the forcing and of a barystatic geo)
        ! or rank 2 (the anomalies of a host-resolved geo into a port: bilinear
        ! onto the domain's grid): those remap through `cpl`, the domain's
        ! coupler, in which wire_resolve primed the map.
        type(globe_t),       intent(in)    :: globe
        character(len=*),    intent(in)    :: src_kind     ! "forcing" | "geo"
        type(wire_t),        intent(inout) :: wires(:)     ! the domain's rows (domain%wires)
        character(len=*),    intent(in)    :: dst_kind     ! kind of the destination instance
        type(field_state_t), intent(inout) :: dst
        type(clock_t),       intent(in)    :: clk
        logical,             intent(in), optional :: init
        type(coupler_class), intent(inout), optional :: cpl

        type(coupler_class) :: cpl0                        ! 0-D rows do not remap
        character(len=16)   :: kind
        real(dp) :: t, t_last
        logical  :: due, rank2
        integer  :: i

        kind = globe_kind(src_kind)
        t = clock_time(clk)
        do i = 1, size(wires)
            if (.not. wire_runs(wires(i), kind, dst_kind)) cycle
            rank2 = dst%f(wires(i)%hd)%rank > 0
            if (rank2 .and. .not. present(cpl)) then
                call fesm_error("kryos_globe", "row "//trim(wires(i)%src_kind)//" -> "//trim(wires(i)%dst_kind)// &
                                ": import `"//trim(dst%f(wires(i)%hd)%name)//"` is not 0-D", &
                                "a rank-2 row remaps: pass the domain's coupler (cpl)")
            end if
            select case (trim(kind))
            case ("forcing")
                call require_0d(dst, wires(i)%hd, wires(i)%src_kind, wires(i)%dst_kind)
                t_last = cadence_t_last(globe%forcing%c%cadence, clk)
                due    = cadence_due(globe%forcing%c%cadence, clk)
                call wire_exchange(wires(i), globe%forcing%c%exports, dst, cpl0, t, t_last, due=due, init=init)
            case ("geo_g")
                t_last = cadence_t_last(globe%geo%c%cadence, clk)
                due    = cadence_due(globe%geo%c%cadence, clk)
                if (rank2) then
                    call wire_exchange(wires(i), globe%geo%c%exports, dst, cpl, t, t_last, due=due, init=init)
                else
                    call wire_exchange(wires(i), globe%geo%c%exports, dst, cpl0, t, t_last, due=due, init=init)
                end if
            end select
        end do
    end subroutine exchange_globe

    subroutine exchange_to_globe(globe, wires, src_kind, src, t_last, dst_kind, clk, due, init)
        ! The domain's resolved rows `src_kind` -> globe `dst_kind` that are
        ! not gathers (a derived coefficient, one producer: climate dTo ->
        ! forcing), in table order: the published exports `src` of the
        ! producing instance into the globe instance's imports. t_last: the
        ! producer's cadence_t_last, due: it is due at the clock time. A
        ! lag 1 row runs before the globe update of the next step (it reads
        ! the producer's previous step), at the cold start with init.
        type(globe_t),       intent(inout) :: globe
        type(wire_t),        intent(inout) :: wires(:)
        character(len=*),    intent(in)    :: src_kind
        type(field_state_t), intent(in)    :: src
        real(dp),            intent(in)    :: t_last
        character(len=*),    intent(in)    :: dst_kind     ! "forcing" | "geo"
        type(clock_t),       intent(in)    :: clk
        logical,             intent(in), optional :: due, init

        type(coupler_class) :: cpl
        character(len=16)   :: kind
        real(dp) :: t
        integer  :: i

        kind = globe_kind(dst_kind)
        t = clock_time(clk)
        do i = 1, size(wires)
            if (.not. wire_runs(wires(i), src_kind, kind)) cycle
            select case (trim(kind))
            case ("forcing")
                call require_0d(globe%forcing%c%imports, wires(i)%hd, wires(i)%src_kind, wires(i)%dst_kind)
                call wire_exchange(wires(i), src, globe%forcing%c%imports, cpl, t, t_last, due=due, init=init)
            case ("geo_g")
                call require_0d(globe%geo%c%imports, wires(i)%hd, wires(i)%src_kind, wires(i)%dst_kind)
                call wire_exchange(wires(i), src, globe%geo%c%imports, cpl, t, t_last, due=due, init=init)
            end select
        end do
    end subroutine exchange_to_globe

    function gather_src(wires, node, dst_kind) result(g)
        ! The contribution of one domain to a gather: its resolved gather row
        ! from the instance `node` (wire_node of a domain instance) into the
        ! globe's `dst_kind`; inactive when the domain has no such active row
        ! (no instance, or the globe has no such import).
        type(wire_t),      intent(in) :: wires(:)           ! the domain's rows
        type(wire_node_t), intent(in) :: node
        character(len=*),  intent(in) :: dst_kind           ! "forcing" | "geo"
        type(gather_src_t) :: g

        character(len=16) :: kind
        integer :: i

        kind = globe_kind(dst_kind)
        g%node = node
        do i = 1, size(wires)
            if (.not. (wires(i)%gather .and. trim(wires(i)%src_kind) == trim(node%kind) &
                       .and. trim(wires(i)%dst_kind) == trim(kind))) cycle
            g%w = wires(i)
            g%active = wires(i)%active
            return
        end do
    end function gather_src

    subroutine gather_into(globe, dst_kind, srcs, t)
        ! The sum of the 0-D exports of the domains' instances into one
        ! import of the globe instance `dst_kind`, in the order of srcs (the
        ! caller passes domain-name order): the values published at exactly t
        ! count, stale ones count 0 (a geo that did not run, an ice not yet
        ! published); the sum is formed in the destination's precision and
        ! put at t. t is the time of the producers' last update: the clock
        ! time for dV_bsl, the previous step for the lag 1 M_ice. Nothing is
        ! put when no source is active (no domain feeds this import).
        type(globe_t),      intent(inout) :: globe
        character(len=*),   intent(in)    :: dst_kind       ! "forcing" | "geo"
        type(gather_src_t), intent(in)    :: srcs(:)
        real(dp),           intent(in)    :: t

        type(wire_t),      allocatable :: ws(:)
        type(wire_node_t), allocatable :: nodes(:)
        integer :: i, n, k

        n = count(srcs%active)
        if (n == 0) return
        allocate(ws(n), nodes(n))
        k = 0
        do i = 1, size(srcs)
            if (.not. srcs(i)%active) cycle
            k = k + 1
            ws(k)    = srcs(i)%w
            nodes(k) = srcs(i)%node
        end do

        select case (trim(globe_kind(dst_kind)))
        case ("forcing")
            call require(forcing_present(globe%forcing), dst_kind)
            call wire_gather(ws, nodes, globe%forcing%c%imports, t)
        case ("geo_g")
            call require(geo_present(globe%geo), dst_kind)
            call wire_gather(ws, nodes, globe%geo%c%imports, t)
        end select
    end subroutine gather_into

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    function globe_kind(which) result(kind)
        ! The table kind of a globe instance.
        character(len=*), intent(in) :: which
        character(len=16) :: kind
        select case (trim(which))
        case ("forcing");        kind = "forcing"
        case ("geo", "geo_g");   kind = "geo_g"
        case default
            kind = ""
            call unknown(which)
        end select
    end function globe_kind

    subroutine require(present_, which)
        logical,          intent(in) :: present_
        character(len=*), intent(in) :: which
        if (present_) return
        call fesm_error("kryos_globe", "the globe has no "//trim(which)//" instance", &
                        "key `"//trim(which)//"` of the [globe] group is empty")
    end subroutine require

    subroutine unknown(which)
        character(len=*), intent(in) :: which
        call fesm_error("kryos_globe", "unknown globe instance '"//trim(which)//"'", "use forcing or geo")
    end subroutine unknown

    subroutine require_0d(dst, h, src_kind, dst_kind)
        type(field_state_t), intent(in) :: dst
        integer,             intent(in) :: h
        character(len=*),    intent(in) :: src_kind, dst_kind
        if (dst%f(h)%rank == 0) return
        call fesm_error("kryos_globe", "row "//trim(src_kind)//" -> "//trim(dst_kind)//": import `"// &
                        trim(dst%f(h)%name)//"` is not 0-D", &
                        "the globe exchanges scalars; fields cross through the domain exchange")
    end subroutine require_0d

end module kryos_globe
