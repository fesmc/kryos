module kryos_domain
    ! One domain (core-design sections 1, 2, 5; contracts-design section 5): its
    ! groups, its grids and the default grid chain, the hub (fesm-utils htopo:
    ! the prepared reference topography on the finest grid and the current
    ! geometry), the coupler, the reference topography on every grid in use
    ! (refs), the slots of the component instances, the resolved wiring rows,
    ! and the three procedures the run sequence calls on it: exchange,
    ! domain_hub_update and domain_check_reference.
    !
    ! domain_init does everything that needs no instance; domain_create then
    ! makes the instances (the init_ctx_t of each: its grid, the domain's
    ! regions on it, the reference on it, the run's constants and calendar
    ! year, the axes of the producers made before it; factory, init, cadence),
    ! and domain_wire resolves the rows. The run sequence (kryos_run) drives
    ! the instances through domain_pass, domain_update, domain_time and
    ! domain_due, and the rows through exchange, exchange_into and
    ! exchange_static; this module is the one place that names the seven
    ! slots (domain_component and the dispatchers below).
    !
    ! Plain types and free procedures; no OOP. Hub precision: the hub is htopo's
    ! (kind wp, sp in phase 1); the field states hold it widened to dp.
    !
    ! References (core R4, contracts-design 5.1): the hub prepares the reference
    ! once (htopo [topo] keys z_bed_f_sd, remove_lakes, H_ice_min, grad_lim_zb;
    ! decision 2026-10-24), and every instance grid's reference is made from that
    ! prepared hub reference, once, by the coupler's `con` map through its
    ! hub-precision path (the remap Yelmo's topo_init/topo_pd came from in
    ! yelmox), and stored widened to dp. The hub's own reference is the prepared
    ! one widened. Only grids in use (an instance on them) get one.
    !
    ! Domain, instance and run groups are read from the par file; the hub's
    ! groups ([topo] and [regions], named by the domain's `topo` and `regions`
    ! keys) must be in the par file too, because htopo_init takes one file.
    !
    ! Errors stop through fesm_error naming the domain and key; the procedures
    ! that can fail on a par-file or consistency error take an optional `stat`.

    use precision,       only: dp, wp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t, fields_init, fields_add, fields_lock, fields_put, &
                               fields_has, fields_handle, fields_describe, FIELD_INSTANT, FIELD_INTENSIVE, &
                               FIELD_EXTENSIVE
    use fesm_clock,      only: clock_t, clock_dtt, clock_time, clock_time_rel, cadence_t_last, cadence_due, &
                               cadence_dt_years, cadence_mark
    use coords,          only: grid_cdo_read_desc
    use regions,         only: regions_init_nml
    use nml,             only: nml_has_group
    use ncio,            only: nc_open, nc_close, nc_write, nc_time_index
    use coupler,         only: coupler_class, coupler_init, coupler_prime, remap
    use htopo,           only: htopo_class, htopo_init, htopo_update, htopo_write_init, htopo_write_step
    use phys_constants,  only: phys_const_class
    use kryos_component, only: component_t, grid_ref_t, init_ctx_t, step_time_t, census_t, axis_t, &
                               step_time_of, axis_find, kryos_mask, component_prepare_call, NAME_LEN, PATH_LEN
    use kryos_config,    only: config_t, config_read, config_validate, config_validate_instance, config_cadence
    use kryos_wiring,    only: wire_t, wire_node_t, wire_node, wire_node_states, wire_resolve, &
                               wire_prime, wire_exchange, wire_runs, wire_label, LAG_INIT, T_GIA
    use kryos_ice,       only: ice_slot_t, ice_create, ice_init, ice_init_state, ice_update, ice_end, &
                               ice_restart_write, ice_restart_read, ice_restart_files
    use kryos_geo,       only: geo_slot_t, geo_create, geo_init, geo_init_state, geo_update, geo_end, &
                               geo_restart_write, geo_restart_read, geo_restart_files
    use kryos_climate,   only: climate_slot_t, climate_create, climate_init, climate_init_state, &
                               climate_update, climate_end, climate_restart_write, climate_restart_read, climate_restart_files
    use kryos_surface,   only: surface_slot_t, surface_create, surface_init, surface_init_state, &
                               surface_update, surface_end, surface_restart_write, surface_restart_read, surface_restart_files
    use kryos_shelf,     only: shelf_slot_t, shelf_create, shelf_init, shelf_init_state, shelf_update, shelf_end, &
                               shelf_restart_write, shelf_restart_read, shelf_restart_files
    use kryos_geothermal, only: geothermal_slot_t, geothermal_create, geothermal_init, geothermal_init_state, &
                               geothermal_update, geothermal_end, geothermal_restart_write, geothermal_restart_read, &
                               geothermal_restart_files
    use kryos_sediments, only: sediments_slot_t, sediments_create, sediments_init, sediments_init_state, &
                               sediments_update, sediments_end, sediments_restart_write, sediments_restart_read, &
                               sediments_restart_files

    implicit none

    private

    ! The components of a domain, in the order of the keys of the [domain]
    ! template (core-design section 2).
    integer, parameter, public :: C_ICE = 1, C_GEO = 2, C_CLIMATE = 3, C_SURFACE = 4, C_SHELF = 5, &
                                  C_GEOTHERMAL = 6, C_SEDIMENTS = 7, NCOMP = 7

    ! Default grid chain (yelmox kryos.f90:574-578): the component whose grid a
    ! component follows when its instance sets no `grid`; 0 = the hub.
    integer, parameter, public :: GRID_DEFAULT(NCOMP) = [0, C_ICE, C_ICE, C_CLIMATE, 0, C_ICE, C_ICE]

    ! Component names: the [domain] keys, the wiring kinds and component_t%kind.
    character(len=16), parameter, public :: COMP_NAME(NCOMP) = &
        [character(len=16) :: "ice", "geo", "climate", "surface", "shelf", "geothermal", "sediments"]

    ! The order the instances are made in (yelmox kryos.f90 domain_init;
    ! core-design 3.1 step 4): a producer of an axis (the climate's depth)
    ! before its consumers (surface, shelf).
    integer, parameter, public :: INIT_ORDER(NCOMP) = [C_ICE, C_GEO, C_GEOTHERMAL, C_SEDIMENTS, &
                                                       C_CLIMATE, C_SURFACE, C_SHELF]

    ! The order the instances are written to a restart bundle and restored from
    ! it (core-design 3.3, 7): the ice before geo (FastIsostasy's restore takes
    ! the restored ice geometry), the climate before surface and shelf.
    integer, parameter, public :: RESTORE_ORDER(NCOMP) = [C_ICE, C_GEO, C_CLIMATE, C_SURFACE, C_SHELF, &
                                                          C_GEOTHERMAL, C_SEDIMENTS]

    ! exchange_into: the sources of the rows into one consumer, in this order
    ! (yelmox couple_to_yelmo: isostasy, surface, shelf, climate; then the
    ! static producers, whose rows have run once at init).
    character(len=16), parameter :: INTO_ORDER(8) = [character(len=16) :: "geo", "surface", "shelf", &
        "climate", "geothermal", "sediments", "ice", "hub"]

    real(dp), parameter :: REF_TOL = 1.0_dp     ! [m] domain_check_reference: round-off only

    type :: hub_handles_t                       ! handles into hub_imports / hub_exports
        integer :: hi_z_bed = 0, hi_H_ice = 0, hi_z_srf = 0, hi_f_ice = 0, hi_f_grnd = 0, hi_z_sl = 0
        integer :: hi_dz_bed = 0, hi_dH_ice = 0
        integer :: h_z_bed = 0, h_H_ice = 0, h_z_srf = 0, h_f_ice = 0, h_f_grnd = 0, h_z_sl = 0
        integer :: h_z_srf_ref = 0
    end type hub_handles_t

    type, public :: domain_t
        character(len=NAME_LEN) :: name = ""            ! its par-file group, e.g. "Greenland" (= the FesmData domain)
        character(len=NAME_LEN) :: grid_hub = ""
        character(len=NAME_LEN) :: topo = "", regions = "", basins = "", out_2D = "none"   ! hub groups and schedule
        character(len=NAME_LEN) :: inst(NCOMP) = ""     ! instance (group) name per component; "" = absent
        character(len=NAME_LEN) :: backend(NCOMP) = ""  ! its backend, from its group
        character(len=NAME_LEN) :: grid(NCOMP) = ""     ! resolved native grid; "" for an absent component
        character(len=8)        :: clock(NCOMP) = "model"   ! the instance's `clock` key: model | forcing
        logical                 :: grid_default(NCOMP) = .true.   ! follows the chain (no `grid` in its group)
        type(htopo_class)       :: hub                  ! sp in phase 1
        logical                 :: hub_mirror = .false. ! grid_hub == grid(ice): the hub mirrors the ice geometry
        type(field_state_t)     :: hub_imports          ! from the ice: mirror set, or anomaly set (finer hub)
        type(field_state_t)     :: hub_exports          ! z_bed, H_ice, z_srf, f_ice, f_grnd, z_sl, z_srf_ref
        type(hub_handles_t)     :: hh
        real(dp)                :: hub_t = 0.0_dp       ! time of the last hub update (lag checks)
        type(coupler_class)     :: cpl                  ! grids by name, map cache (maps/)
        type(grid_ref_t), allocatable :: refs(:)        ! one per distinct grid in use, hub first
        type(ice_slot_t)        :: ice                  ! instance slots (domain_create fills them)
        type(geo_slot_t)        :: geo
        type(climate_slot_t)    :: climate
        type(surface_slot_t)    :: surface
        type(shelf_slot_t)      :: shelf
        type(geothermal_slot_t) :: geothermal
        type(sediments_slot_t)  :: sediments
        type(wire_t), allocatable :: wires(:)           ! rows resolved for this domain (domain_wire)
    end type domain_t

    public :: domain_init, domain_create, domain_describe
    public :: domain_ref_index, domain_comp_index, domain_component, domain_present
    public :: domain_time, domain_due, domain_pass, domain_prepare, domain_update, domain_end
    public :: domain_restart_write, domain_restart_read, domain_restart_files
    public :: domain_nodes, domain_wire
    public :: exchange, exchange_into, exchange_static, domain_hub_update, domain_check_reference
    public :: domain_hub_write

contains

    ! =========================================================================
    ! Init
    ! =========================================================================

    subroutine domain_init(dom, cfg, name, cnst, t0, stat)
        ! The domain `name` (a [domain]-template group of the par file) without
        ! its instances: groups, grids and the chain, the hub loaded and
        ! prepared, the coupler, refs, the hub exports published at t0.
        !   1. the group and the instance groups validated (templates);
        !   2. the grid of every present component: its group's `grid`, else the
        !      chain of GRID_DEFAULT;
        !   3. htopo_init (prepared reference), coupler_init;
        !   4. refs: hub first, then each other grid in use (con, hub path).
        type(domain_t),         intent(out)   :: dom
        type(config_t),         intent(in)    :: cfg
        character(len=*),       intent(in)    :: name
        type(phys_const_class), intent(in)    :: cnst        ! of the run (one record)
        real(dp),               intent(in)    :: t0          ! time_init [yr]
        integer,                intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_domain::domain_init"
        character(len=NAME_LEN) :: explicit(NCOMP)
        integer :: k, j, s

        if (present(stat)) stat = 0
        dom%name = name

        ! --- the domain group ---
        call config_validate(cfg, name, "domain", required=[character(len=16) :: "grid_hub"], stat=stat)
        if (failed(stat)) return
        call config_read(cfg, name, "grid_hub", dom%grid_hub, "domain")
        call config_read(cfg, name, "topo",     dom%topo,     "domain")
        call config_read(cfg, name, "regions",  dom%regions,  "domain")
        call config_read(cfg, name, "basins",   dom%basins,   "domain")
        call config_read(cfg, name, "out_2D",   dom%out_2D,   "domain")
        do k = 1, NCOMP
            call config_read(cfg, name, trim(COMP_NAME(k)), dom%inst(k), "domain")
        end do
        if (len_trim(dom%grid_hub) == 0) then
            call fesm_error(proc, "domain &"//trim(name)//": grid_hub is empty", &
                            "par file "//trim(cfg%par_file), stat)
            return
        end if

        ! the hub's groups: htopo_init reads them all from the par file
        if (.not. nml_has_group(trim(cfg%par_file), trim(dom%topo))) then
            call fesm_error(proc, "domain &"//trim(name)//": group &"//trim(dom%topo)// &
                            " (topo) is missing from the par file", &
                            "the hub reads its topography group from the par file "//trim(cfg%par_file), stat)
            return
        end if
        if (.not. nml_has_group(trim(cfg%par_file), trim(dom%regions))) then
            call fesm_error(proc, "domain &"//trim(name)//": group &"//trim(dom%regions)// &
                            " (regions) is missing from the par file", &
                            "the hub reads its regions group from the par file "//trim(cfg%par_file), stat)
            return
        end if

        ! --- instances: unique names, groups valid, grids read ---
        explicit = ""
        do k = 1, NCOMP
            if (len_trim(dom%inst(k)) == 0) cycle
            do j = 1, k - 1
                if (trim(dom%inst(j)) /= trim(dom%inst(k))) cycle
                call fesm_error(proc, "domain &"//trim(name)//": instance `"//trim(dom%inst(k))// &
                                "` is named by two components", trim(COMP_NAME(j))//" and "// &
                                trim(COMP_NAME(k))//": one group per component", stat)
                return
            end do
            call config_validate_instance(cfg, dom%inst(k), trim(COMP_NAME(k)), dom%backend(k), stat=s)
            if (s /= 0) then
                call fesm_error(proc, "domain &"//trim(name)//": "//trim(COMP_NAME(k))//" instance `"// &
                                trim(dom%inst(k))//"` is not a valid instance group", &
                                "par file "//trim(cfg%par_file), stat)
                return
            end if
            call config_read(cfg, dom%inst(k), "grid", explicit(k), "instance")
        end do

        ! --- the grids ---
        call resolve_grids(dom, explicit)
        dom%hub_mirror = (len_trim(dom%inst(C_ICE)) > 0) .and. (trim(dom%grid(C_ICE)) == trim(dom%grid_hub))

        ! --- the hub and the coupler ---
        call htopo_init(dom%hub, trim(cfg%par_file), trim(name), trim(dom%topo), trim(dom%regions), &
                        trim(name), trim(dom%grid_hub), cnst, map_fldr=trim(cfg%map_fldr))
        call coupler_init(dom%cpl, trim(cfg%map_fldr))

        ! --- the reference on every grid in use ---
        call make_refs(dom)

        ! --- hub field states, published at t0 (the reference is the initial geometry) ---
        call hub_declare(dom)
        call hub_put(dom, t0, hub_f_ice(dom))
        call fields_put(dom%hub_exports, dom%hh%h_z_srf_ref, dom%hub%ref%z_srf, t1=t0)   ! static: once
        dom%hub_t = t0

        call domain_describe(dom)

        ! The instances: domain_create (needs the clock for their cadences),
        ! then domain_wire. Until then every slot is unallocated.
    end subroutine domain_init

    subroutine resolve_grids(dom, explicit)
        ! The grid of each present component: its own `grid`, else the grid of
        ! the component it follows (GRID_DEFAULT; the hub for 0). The chain runs
        ! through absent components too (an absent ice still stands for the hub
        ! grid), so a geo without an ice follows the hub.
        type(domain_t),          intent(inout) :: dom
        character(len=NAME_LEN), intent(in)    :: explicit(NCOMP)

        character(len=NAME_LEN) :: g(NCOMP)
        integer :: k

        do k = 1, NCOMP
            if (len_trim(explicit(k)) > 0) then
                g(k) = explicit(k)
            else if (GRID_DEFAULT(k) == 0) then
                g(k) = dom%grid_hub
            else
                g(k) = g(GRID_DEFAULT(k))           ! parents have lower indices
            end if
            dom%grid_default(k) = (len_trim(explicit(k)) == 0)
            dom%grid(k) = ""
            if (len_trim(dom%inst(k)) > 0) dom%grid(k) = g(k)
        end do
    end subroutine resolve_grids

    subroutine make_refs(dom)
        ! refs(1): the prepared hub reference widened. refs(i > 1): each other
        ! grid in use, from the hub reference by the coupler's `con` map in the
        ! hub's precision, widened (core R4). Stored in order of first use.
        type(domain_t), intent(inout) :: dom

        character(len=NAME_LEN) :: list(NCOMP + 1)
        integer :: n, k, i

        n = 1
        list(1) = dom%grid_hub
        do k = 1, NCOMP
            if (len_trim(dom%grid(k)) == 0) cycle
            if (any(list(1:n) == dom%grid(k))) cycle
            n = n + 1
            list(n) = dom%grid(k)
        end do

        allocate(dom%refs(n))
        associate (ref => dom%hub%ref)
            dom%refs(1)%grid    = dom%grid_hub
            dom%refs(1)%z_bed   = real(ref%z_bed,    dp)
            dom%refs(1)%H_ice   = real(ref%H_ice,    dp)
            dom%refs(1)%z_srf   = real(ref%z_srf,    dp)
            dom%refs(1)%z_bed_sd = real(ref%z_bed_sd, dp)
            do i = 2, n
                dom%refs(i)%grid = list(i)
                call coupler_prime(dom%cpl, trim(dom%grid_hub), trim(list(i)), method="con")
                call ref_remap(dom, ref%z_bed,    list(i), dom%refs(i)%z_bed)
                call ref_remap(dom, ref%H_ice,    list(i), dom%refs(i)%H_ice)
                call ref_remap(dom, ref%z_srf,    list(i), dom%refs(i)%z_srf)
                call ref_remap(dom, ref%z_bed_sd, list(i), dom%refs(i)%z_bed_sd)
            end do
        end associate
    end subroutine make_refs

    subroutine ref_remap(dom, a, grid, b)
        ! b = widened con remap of the hub-precision field a, hub -> grid.
        type(domain_t),        intent(inout) :: dom
        real(wp),              intent(in)    :: a(:,:)
        character(len=*),      intent(in)    :: grid
        real(dp), allocatable, intent(out)   :: b(:,:)
        real(wp), allocatable :: tmp(:,:)

        call remap(dom%cpl, a, trim(dom%grid_hub), tmp, trim(grid), method="con")
        b = real(tmp, dp)
    end subroutine ref_remap

    ! =========================================================================
    ! Instances
    ! =========================================================================

    subroutine domain_create(dom, cfg, clk, cnst, census)
        ! The instances of the domain (core-design 3.1 step 4), in INIT_ORDER:
        ! each present component's init_ctx_t (instance_ctx), its factory
        ! (admission against the run's census), <comp>_init, its cadence on
        ! the clock (config_cadence: keys dt, phase, clock), one describe line
        ! and its field states. The axes an instance declares (the climate's
        ! month and depth) reach the instances made after it in ctx%axes.
        ! Call after domain_init and before domain_wire and clock_start.
        type(domain_t),  target, intent(inout) :: dom
        type(config_t),          intent(in)    :: cfg
        type(clock_t),           intent(inout) :: clk
        type(phys_const_class),  intent(in)    :: cnst
        type(census_t),          intent(inout) :: census

        type(init_ctx_t) :: ctx
        type(axis_t), allocatable :: axes(:)
        class(component_t), pointer :: c
        integer :: j, k, i

        allocate(axes(0))
        do j = 1, NCOMP
            k = INIT_ORDER(j)
            if (len_trim(dom%inst(k)) == 0) cycle

            call instance_ctx(dom, k, cfg, clk, cnst, axes, ctx)
            call create_init(dom, k, ctx, census)

            c => domain_component(dom, k)
            call config_cadence(cfg, dom%inst(k), trim(COMP_NAME(k)), dom%backend(k), clk, c%cadence, dom%clock(k))
            if (allocated(c%axes)) then
                do i = 1, size(c%axes)
                    if (axis_find(axes, c%axes(i)%name) == 0) axes = [axes, c%axes(i)]
                end do
            end if
            call instance_describe(dom, k, clk)
        end do
    end subroutine domain_create

    subroutine instance_ctx(dom, k, cfg, clk, cnst, axes, ctx)
        ! The init_ctx_t of component k (core-design 2, 3.6 item 5): the
        ! instance, its grid (maps/grid_<name>.txt), the domain's regions on it
        ! (the [regions] group of the hub grid, remapped by nearest neighbour,
        ! as yelmox domain_regions_load), the reference on it, the run's
        ! constants, library and calendar year, and the time of the clock (t0).
        type(domain_t),         intent(in)  :: dom
        integer,                intent(in)  :: k
        type(config_t),         intent(in)  :: cfg
        type(clock_t),          intent(in)  :: clk
        type(phys_const_class), intent(in)  :: cnst
        type(axis_t),           intent(in)  :: axes(:)
        type(init_ctx_t),       intent(out) :: ctx

        integer :: ir

        ctx%par_file  = cfg%par_file
        ctx%library   = cfg%library
        ctx%instance  = dom%inst(k)
        ctx%component = COMP_NAME(k)
        ctx%backend   = dom%backend(k)
        ctx%domain    = dom%name
        ctx%grid_name = dom%grid(k)
        ctx%scope     = "regional"
        call grid_cdo_read_desc(ctx%grid, trim(dom%grid(k)), trim(cfg%map_fldr))
        call regions_init_nml(ctx%reg, trim(cfg%par_file), trim(dom%regions), domain=trim(dom%name), &
                              grid_name=trim(dom%grid_hub), grid=ctx%grid, verbose=.false.)
        ir = domain_ref_index(dom, dom%grid(k))
        if (ir == 0) then
            call fesm_error("kryos_domain::domain_create", "domain "//trim(dom%name)//": no reference on grid "// &
                            trim(dom%grid(k)), "instance "//trim(dom%inst(k)))
        end if
        ctx%ref      = dom%refs(ir)
        ctx%cnst     = cnst
        ctx%time     = clock_time(clk)
        ctx%time_rel = clock_time_rel(clk)
        ctx%dtt      = clock_dtt(clk)
        ctx%sec_year = clk%sec_year
        ctx%axes     = axes
    end subroutine instance_ctx

    subroutine create_init(dom, k, ctx, census)
        ! Factory and init of component k.
        type(domain_t),   intent(inout) :: dom
        integer,          intent(in)    :: k
        type(init_ctx_t), intent(in)    :: ctx
        type(census_t),   intent(inout) :: census

        select case (k)
        case (C_ICE)
            call ice_create(dom%ice, ctx, census);                call ice_init(dom%ice, ctx)
        case (C_GEO)
            call geo_create(dom%geo, ctx, census);                call geo_init(dom%geo, ctx)
        case (C_CLIMATE)
            call climate_create(dom%climate, ctx, census);        call climate_init(dom%climate, ctx)
        case (C_SURFACE)
            call surface_create(dom%surface, ctx, census);        call surface_init(dom%surface, ctx)
        case (C_SHELF)
            call shelf_create(dom%shelf, ctx, census);            call shelf_init(dom%shelf, ctx)
        case (C_GEOTHERMAL)
            call geothermal_create(dom%geothermal, ctx, census);  call geothermal_init(dom%geothermal, ctx)
        case (C_SEDIMENTS)
            call sediments_create(dom%sediments, ctx, census);    call sediments_init(dom%sediments, ctx)
        end select
    end subroutine create_init

    subroutine instance_describe(dom, k, clk)
        ! core-design 8: backend, capabilities, cadence and clock, then the
        ! field states (fields_describe of imports and exports).
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        type(clock_t),          intent(in) :: clk

        class(component_t), pointer :: c
        real(dp) :: dt
        character(len=:), allocatable :: cad, caps

        c => domain_component(dom, k)
        dt = cadence_dt_years(c%cadence, clk)
        if (dt > 0.0_dp) then
            cad = "every "//to_str(dt)//" yr"
        else
            cad = "off (static)"
        end if
        caps = ""
        if (c%caps%has_ocn)       caps = caps//" has_ocn"
        if (c%caps%has_ocn_shelf) caps = caps//" has_ocn_shelf"
        if (c%caps%has_smb)       caps = caps//" has_smb"
        if (c%caps%has_Qd)        caps = caps//" has_Qd"
        if (c%caps%has_flotation) caps = caps//" has_flotation"
        if (len(caps) == 0) caps = " none"
        call log_msg(LOG_INFO, "kryos_domain", "domain "//trim(dom%name)//": "//trim(COMP_NAME(k))//" `"// &
                     trim(c%name)//"` ("//trim(c%backend)//") on "//trim(c%grid)//", "//cad//", clock "// &
                     trim(dom%clock(k))//", caps:"//caps)
        call fields_describe(c%imports)
        call fields_describe(c%exports)
    end subroutine instance_describe

    function domain_time(dom, k, clk, init) result(st)
        ! The step_time_t of component k (step_time_of with its cadence and
        ! clock key); init: t0 = t1 = now (the init passes).
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        type(clock_t),          intent(in) :: clk
        logical,                intent(in), optional :: init
        type(step_time_t) :: st

        class(component_t), pointer :: c

        c => require_component(dom, k, "domain_time")
        st = step_time_of(c%cadence, clk, dom%clock(k), init)
    end function domain_time

    logical function domain_present(dom, k)
        ! Is component k's instance created? (presence, never a backend name)
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        domain_present = associated(domain_component(dom, k))
    end function domain_present

    logical function domain_due(dom, k, clk)
        ! Component k is present and its cadence is due at the clock's time.
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        type(clock_t),          intent(in) :: clk

        class(component_t), pointer :: c

        domain_due = .false.
        c => domain_component(dom, k)
        if (associated(c)) domain_due = cadence_due(c%cadence, clk)
    end function domain_due

    subroutine domain_pass(dom, k, clk, pass)
        ! <comp>_init_state of component k at the clock's time (t0 = t1 = now),
        ! pass PASS_INIT .. PASS_RESTART (core-design 3.2). Absent: nothing.
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(in)    :: k
        type(clock_t),          intent(in)    :: clk
        integer,                intent(in)    :: pass

        type(step_time_t) :: st

        if (.not. domain_present(dom, k)) return
        st = domain_time(dom, k, clk, init=.true.)
        select case (k)
        case (C_ICE);        call ice_init_state(dom%ice, st, pass)
        case (C_GEO);        call geo_init_state(dom%geo, st, pass)
        case (C_CLIMATE);    call climate_init_state(dom%climate, st, pass)
        case (C_SURFACE);    call surface_init_state(dom%surface, st, pass)
        case (C_SHELF);      call shelf_init_state(dom%shelf, st, pass)
        case (C_GEOTHERMAL); call geothermal_init_state(dom%geothermal, st, pass)
        case (C_SEDIMENTS);  call sediments_init_state(dom%sediments, st, pass)
        end select
    end subroutine domain_pass

    subroutine domain_prepare(dom, k, clk, init)
        ! <comp>_prepare of component k: its imports into the backend's
        ! boundary, over its cadence interval (init: t0 = t1 = now, the cold
        ! start). The run sequence calls it when the component is due, before
        ! the experiment hook and domain_update. Absent: nothing.
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(in)    :: k
        type(clock_t),          intent(in)    :: clk
        logical,                intent(in), optional :: init

        class(component_t), pointer :: c

        c => domain_component(dom, k)
        if (.not. associated(c)) return
        call component_prepare_call(c, domain_time(dom, k, clk, init))
    end subroutine domain_prepare

    subroutine domain_update(dom, k, clk)
        ! <comp>_update of component k over its cadence interval, then the
        ! cadence mark (core-design 3.4 `mark`). The caller tests domain_due.
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(in)    :: k
        type(clock_t),          intent(in)    :: clk

        type(step_time_t) :: st

        st = domain_time(dom, k, clk)
        select case (k)
        case (C_ICE);        call ice_update(dom%ice, st);               call cadence_mark(dom%ice%c%cadence, clk)
        case (C_GEO);        call geo_update(dom%geo, st);               call cadence_mark(dom%geo%c%cadence, clk)
        case (C_CLIMATE);    call climate_update(dom%climate, st);       call cadence_mark(dom%climate%c%cadence, clk)
        case (C_SURFACE);    call surface_update(dom%surface, st);       call cadence_mark(dom%surface%c%cadence, clk)
        case (C_SHELF);      call shelf_update(dom%shelf, st);           call cadence_mark(dom%shelf%c%cadence, clk)
        case (C_GEOTHERMAL); call geothermal_update(dom%geothermal, st); call cadence_mark(dom%geothermal%c%cadence, clk)
        case (C_SEDIMENTS);  call sediments_update(dom%sediments, st);   call cadence_mark(dom%sediments%c%cadence, clk)
        end select
    end subroutine domain_update

    subroutine domain_restart_write(dom, k, dir, t)
        ! <comp>_restart_write of component k: the backend's own state into the
        ! folder `dir` at time t (nothing for a backend without a restart file,
        ! constraints restart = .false.). Absent: nothing. The field states are
        ! the core's (kryos_restart).
        type(domain_t),   intent(inout) :: dom
        integer,          intent(in)    :: k
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t

        if (.not. domain_present(dom, k)) return
        select case (k)
        case (C_ICE);        call ice_restart_write(dom%ice, dir, t)
        case (C_GEO);        call geo_restart_write(dom%geo, dir, t)
        case (C_CLIMATE);    call climate_restart_write(dom%climate, dir, t)
        case (C_SURFACE);    call surface_restart_write(dom%surface, dir, t)
        case (C_SHELF);      call shelf_restart_write(dom%shelf, dir, t)
        case (C_GEOTHERMAL); call geothermal_restart_write(dom%geothermal, dir, t)
        case (C_SEDIMENTS);  call sediments_restart_write(dom%sediments, dir, t)
        end select
    end subroutine domain_restart_write

    subroutine domain_restart_read(dom, k, dir, t)
        ! <comp>_restart_read of component k from the folder `dir` (the
        ! counterpart of domain_restart_write). Absent: nothing.
        type(domain_t),   intent(inout) :: dom
        integer,          intent(in)    :: k
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t

        if (.not. domain_present(dom, k)) return
        select case (k)
        case (C_ICE);        call ice_restart_read(dom%ice, dir, t)
        case (C_GEO);        call geo_restart_read(dom%geo, dir, t)
        case (C_CLIMATE);    call climate_restart_read(dom%climate, dir, t)
        case (C_SURFACE);    call surface_restart_read(dom%surface, dir, t)
        case (C_SHELF);      call shelf_restart_read(dom%shelf, dir, t)
        case (C_GEOTHERMAL); call geothermal_restart_read(dom%geothermal, dir, t)
        case (C_SEDIMENTS);  call sediments_restart_read(dom%sediments, dir, t)
        end select
    end subroutine domain_restart_read

    function domain_restart_files(dom, k) result(files)
        ! <comp>_restart_files of component k: the base names its restart_write
        ! creates in the folder it is given. Absent: none.
        type(domain_t), intent(in) :: dom
        integer,        intent(in) :: k
        character(len=PATH_LEN), allocatable :: files(:)

        allocate(files(0))
        if (.not. domain_present(dom, k)) return
        select case (k)
        case (C_ICE);        files = ice_restart_files(dom%ice)
        case (C_GEO);        files = geo_restart_files(dom%geo)
        case (C_CLIMATE);    files = climate_restart_files(dom%climate)
        case (C_SURFACE);    files = surface_restart_files(dom%surface)
        case (C_SHELF);      files = shelf_restart_files(dom%shelf)
        case (C_GEOTHERMAL); files = geothermal_restart_files(dom%geothermal)
        case (C_SEDIMENTS);  files = sediments_restart_files(dom%sediments)
        end select
    end function domain_restart_files

    subroutine domain_end(dom)
        ! <comp>_end of every instance, the ice last (Yelmo deallocates its
        ! state; core-design 3.5).
        type(domain_t), intent(inout) :: dom

        if (allocated(dom%geo%c))        call geo_end(dom%geo)
        if (allocated(dom%climate%c))    call climate_end(dom%climate)
        if (allocated(dom%surface%c))    call surface_end(dom%surface)
        if (allocated(dom%shelf%c))      call shelf_end(dom%shelf)
        if (allocated(dom%geothermal%c)) call geothermal_end(dom%geothermal)
        if (allocated(dom%sediments%c))  call sediments_end(dom%sediments)
        if (allocated(dom%ice%c))        call ice_end(dom%ice)
    end subroutine domain_end

    function require_component(dom, k, who) result(c)
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        character(len=*),       intent(in) :: who
        class(component_t), pointer :: c
        c => domain_component(dom, k)
        if (.not. associated(c)) then
            call fesm_error("kryos_domain::"//who, "domain "//trim(dom%name)//": no "//trim(COMP_NAME(k))// &
                            " instance", "test domain_present first")
        end if
    end function require_component

    ! =========================================================================
    ! Lookup
    ! =========================================================================

    integer function domain_ref_index(dom, grid) result(i)
        ! Index of the reference on `grid` in dom%refs; 0 if no instance of the
        ! domain is on it (and it is not the hub grid).
        type(domain_t),   intent(in) :: dom
        character(len=*), intent(in) :: grid
        integer :: j

        i = 0
        if (.not. allocated(dom%refs)) return
        do j = 1, size(dom%refs)
            if (trim(dom%refs(j)%grid) == trim(grid)) then
                i = j
                return
            end if
        end do
    end function domain_ref_index

    integer function domain_comp_index(kind) result(k)
        ! Component index of a kind name ("ice", ...); 0 if unknown (e.g. "hub").
        character(len=*), intent(in) :: kind
        integer :: j

        k = 0
        do j = 1, NCOMP
            if (trim(COMP_NAME(j)) == trim(kind)) then
                k = j
                return
            end if
        end do
    end function domain_comp_index

    function domain_component(dom, k) result(c)
        ! The instance of component k, or a disassociated pointer if its slot
        ! is not allocated. The one place that names the seven slots.
        type(domain_t), target, intent(in) :: dom
        integer,                intent(in) :: k
        class(component_t), pointer :: c

        c => null()
        select case (k)
        case (C_ICE);        if (allocated(dom%ice%c))        c => dom%ice%c
        case (C_GEO);        if (allocated(dom%geo%c))        c => dom%geo%c
        case (C_CLIMATE);    if (allocated(dom%climate%c))    c => dom%climate%c
        case (C_SURFACE);    if (allocated(dom%surface%c))    c => dom%surface%c
        case (C_SHELF);      if (allocated(dom%shelf%c))      c => dom%shelf%c
        case (C_GEOTHERMAL); if (allocated(dom%geothermal%c)) c => dom%geothermal%c
        case (C_SEDIMENTS);  if (allocated(dom%sediments%c))  c => dom%sediments%c
        end select
    end function domain_component

    ! =========================================================================
    ! Describe
    ! =========================================================================

    subroutine domain_describe(dom)
        ! Grids per component ("(default)" where the chain decided), the hub
        ! and the references, logged once (core-design section 8).
        type(domain_t), intent(in) :: dom

        character(len=*), parameter :: who = "kryos_domain"
        integer :: k, i

        call log_msg(LOG_INFO, who, "domain "//trim(dom%name)//": hub "//trim(dom%grid_hub)// &
                     " ("//to_str(dom%hub%nx)//" x "//to_str(dom%hub%ny)//"), basins '"//trim(dom%basins)//"'")
        do k = 1, NCOMP
            if (len_trim(dom%inst(k)) == 0) then
                call log_msg(LOG_INFO, who, "  "//pad(COMP_NAME(k))//" absent")
                cycle
            end if
            call log_msg(LOG_INFO, who, "  "//pad(COMP_NAME(k))//" "//pad(dom%inst(k))//" ("// &
                         trim(dom%backend(k))//") on "//trim(dom%grid(k))// &
                         merge(" (default)", "          ", dom%grid_default(k)))
        end do
        do i = 1, size(dom%refs)
            call log_msg(LOG_INFO, who, "  reference on "//trim(dom%refs(i)%grid)// &
                         merge(" (hub)", "      ", i == 1)//": z_bed, H_ice, z_srf, z_bed_sd, widened")
        end do
        if (dom%hub_mirror) then
            call log_msg(LOG_INFO, who, "  hub mirrors the ice geometry (hub grid = ice grid)")
        else
            call log_msg(LOG_INFO, who, "  hub is the reference plus the ice anomalies")
        end if
    end subroutine domain_describe

    function pad(s) result(p)
        ! s left-justified in at least 12 characters (log alignment).
        character(len=*), intent(in) :: s
        character(len=:), allocatable :: p
        p = trim(s)//repeat(" ", max(0, 12 - len_trim(s)))
    end function pad

    ! =========================================================================
    ! The hub: field states and refresh
    ! =========================================================================

    subroutine hub_declare(dom)
        ! hub_exports: the geometry every consumer reads (contracts-design 2.7).
        ! hub_imports: what the ice puts in, the mirror set when the hub is on
        ! the ice grid, else the anomaly set (empty without an ice instance).
        type(domain_t), intent(inout) :: dom

        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [dom%hub%nx, dom%hub%ny]
        g  = trim(dom%grid_hub)

        associate (h => dom%hh)
            call fields_init(dom%hub_exports, trim(dom%name)//":hub_exports")
            call fields_add(dom%hub_exports, "z_bed",     "m", FIELD_INSTANT, FIELD_INTENSIVE, h%h_z_bed,     grid=g, dims=d2)
            call fields_add(dom%hub_exports, "H_ice",     "m", FIELD_INSTANT, FIELD_EXTENSIVE, h%h_H_ice,     grid=g, dims=d2)
            call fields_add(dom%hub_exports, "z_srf",     "m", FIELD_INSTANT, FIELD_INTENSIVE, h%h_z_srf,     grid=g, dims=d2)
            call fields_add(dom%hub_exports, "f_ice",     "1", FIELD_INSTANT, FIELD_INTENSIVE, h%h_f_ice,     grid=g, dims=d2)
            call fields_add(dom%hub_exports, "f_grnd",    "1", FIELD_INSTANT, FIELD_INTENSIVE, h%h_f_grnd,    grid=g, dims=d2)
            call fields_add(dom%hub_exports, "z_sl",      "m", FIELD_INSTANT, FIELD_INTENSIVE, h%h_z_sl,      grid=g, dims=d2)
            call fields_add(dom%hub_exports, "z_srf_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, h%h_z_srf_ref, grid=g, dims=d2)
            call fields_lock(dom%hub_exports)

            call fields_init(dom%hub_imports, trim(dom%name)//":hub_imports")
            if (len_trim(dom%inst(C_ICE)) > 0) then
                if (dom%hub_mirror) then
                    call fields_add(dom%hub_imports, "z_bed",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_z_bed,  grid=g, dims=d2)
                    call fields_add(dom%hub_imports, "H_ice",  "m", FIELD_INSTANT, FIELD_EXTENSIVE, h%hi_H_ice,  grid=g, dims=d2)
                    call fields_add(dom%hub_imports, "z_srf",  "m", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_z_srf,  grid=g, dims=d2)
                    call fields_add(dom%hub_imports, "f_ice",  "1", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_f_ice,  grid=g, dims=d2)
                    call fields_add(dom%hub_imports, "f_grnd", "1", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_f_grnd, grid=g, dims=d2)
                else
                    call fields_add(dom%hub_imports, "dz_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_dz_bed, grid=g, dims=d2)
                    call fields_add(dom%hub_imports, "dH_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, h%hi_dH_ice, grid=g, dims=d2)
                end if
                call fields_add(dom%hub_imports, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, h%hi_z_sl, grid=g, dims=d2)
            end if
            call fields_lock(dom%hub_imports)
        end associate
    end subroutine hub_declare

    function hub_f_ice(dom) result(f)
        ! Ice-covered fraction of the hub's current geometry: each hub cell is
        ! ice-covered or not (the hub of an ice sheet finer than the ice is
        ! never fractional; a mirror hub takes the ice's own f_ice instead).
        type(domain_t), intent(in) :: dom
        real(dp) :: f(dom%hub%nx, dom%hub%ny)
        f = merge(1.0_dp, 0.0_dp, dom%hub%H_ice > 0.0_wp)
    end function hub_f_ice

    subroutine hub_put(dom, t, f_ice)
        ! Publish the hub's current geometry (kind wp, widened exactly) at t.
        type(domain_t), intent(inout) :: dom
        real(dp),       intent(in)    :: t
        real(dp),       intent(in)    :: f_ice(:,:)

        associate (h => dom%hh, e => dom%hub_exports, p => dom%hub)
            call fields_put(e, h%h_z_bed,  p%z_bed,  t1=t)
            call fields_put(e, h%h_H_ice,  p%H_ice,  t1=t)
            call fields_put(e, h%h_z_srf,  p%z_srf,  t1=t)
            call fields_put(e, h%h_f_ice,  f_ice,    t1=t)
            call fields_put(e, h%h_f_grnd, p%f_grnd, t1=t)
            call fields_put(e, h%h_z_sl,   p%z_sl,   t1=t)
        end associate
    end subroutine hub_put

    subroutine domain_hub_update(dom, t, clk, init, stat)
        ! Refresh the hub after the ice step (couple_yelmo_to_htopo, contracts-
        ! design 5.2): the ice rows into hub_imports, then
        !   mirror (grid_hub == grid(ice)): the hub is the ice geometry, its
        !     fractional f_grnd and f_ice kept;
        !   else: the hub keeps its reference and takes the ice's anomalies
        !     (bilinear; the bed displacement against refs(ice)%z_bed, the
        !     thickness change against refs(ice)%H_ice), and htopo_update
        !     recomputes the grounding and the surface on the hub.
        ! Publishes hub_exports at t and sets hub_t. Without an ice instance
        ! the hub stays the reference and is republished.
        type(domain_t), target, intent(inout) :: dom
        real(dp),               intent(in)    :: t
        type(clock_t),          intent(in)    :: clk
        logical,                intent(in),  optional :: init    ! cold start: every row counts as lag 0
        integer,                intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_domain::domain_hub_update"
        real(dp), allocatable :: f_ice(:,:)
        integer :: s

        if (present(stat)) stat = 0

        if (len_trim(dom%inst(C_ICE)) > 0) then
            if (.not. allocated(dom%ice%c)) then
                call fesm_error(proc, "domain "//trim(dom%name)//": ice instance `"//trim(dom%inst(C_ICE))// &
                                "` is configured but not created", "create the instances first", stat)
                return
            end if
            call exchange(dom, "ice", "hub", t, clk, init=init, stat=s)
            if (s /= 0) then
                call fesm_error(proc, "domain "//trim(dom%name)//": the ice rows into the hub failed", "", stat)
                return
            end if

            associate (h => dom%hh, im => dom%hub_imports, p => dom%hub)
                if (dom%hub_mirror) then
                    p%z_bed  = real(im%f(h%hi_z_bed)%v(:,:,1),  wp)
                    p%H_ice  = real(im%f(h%hi_H_ice)%v(:,:,1),  wp)
                    p%z_srf  = real(im%f(h%hi_z_srf)%v(:,:,1),  wp)
                    p%f_grnd = real(im%f(h%hi_f_grnd)%v(:,:,1), wp)
                    p%z_sl   = real(im%f(h%hi_z_sl)%v(:,:,1),   wp)
                    f_ice    = im%f(h%hi_f_ice)%v(:,:,1)
                else
                    call htopo_update(p, real(im%f(h%hi_dz_bed)%v(:,:,1), wp), &
                                         real(im%f(h%hi_dH_ice)%v(:,:,1), wp), &
                                         real(im%f(h%hi_z_sl)%v(:,:,1),   wp))
                    f_ice = hub_f_ice(dom)
                end if
            end associate
        else
            f_ice = hub_f_ice(dom)
        end if

        call hub_put(dom, t, f_ice)
        dom%hub_t = t
    end subroutine domain_hub_update

    ! =========================================================================
    ! Reference check (decision 4, contracts-design 5.3)
    ! =========================================================================

    subroutine domain_check_reference(dom, stat)
        ! The reference bed geo and the ice work from is the domain's: geo's
        ! z_bed_ref export against refs(geo grid)%z_bed, and the ice's
        ! z_bed_ref (if declared) against refs(ice grid)%z_bed. Logs the mean
        ! and max |difference|; stops above 1 m (yelmox check_isostasy_
        ! reference): on a restart the geo or ice restart came from another
        ! topography. Exact at the cold start, round-off after steps.
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_domain::domain_check_reference"
        integer, parameter :: kinds(2) = [C_GEO, C_ICE]
        class(component_t), pointer :: c
        integer :: j, h, i
        real(dp) :: dmean, dmax

        if (present(stat)) stat = 0

        do j = 1, size(kinds)
            c => domain_component(dom, kinds(j))
            if (.not. associated(c)) cycle
            if (.not. fields_has(c%exports, "z_bed_ref")) cycle
            h = fields_handle(c%exports, "z_bed_ref")
            associate (f => c%exports%f(h))
                if (.not. f%published) then
                    call fesm_error(proc, trim(c%name)//": z_bed_ref was never published", &
                                    "the check runs after the instance's PASS_GEOMETRY", stat)
                    return
                end if
                i = domain_ref_index(dom, f%grid)
                if (i == 0) then
                    call fesm_error(proc, trim(c%name)//": no reference topography on grid "//trim(f%grid), &
                                    "domain "//trim(dom%name), stat)
                    return
                end if
                dmean = sum(f%v(:,:,1) - dom%refs(i)%z_bed) / real(size(dom%refs(i)%z_bed), dp)
                dmax  = maxval(abs(f%v(:,:,1) - dom%refs(i)%z_bed))
            end associate

            call log_msg(LOG_INFO, proc, trim(c%name)//" z_bed_ref against the domain reference on "// &
                         trim(c%grid)//" [m]: mean diff "//to_str(dmean)//", max |diff| "//to_str(dmax)// &
                         " (tolerance "//to_str(REF_TOL)//")")

            if (dmax > REF_TOL) then
                call fesm_error(proc, trim(c%name)//": its reference bedrock does not match the domain's", &
                    "max |diff| = "//to_str(dmax)//" m > "//to_str(REF_TOL)//" m on grid "//trim(c%grid)// &
                    new_line("a")//"z_bed is rebuilt as z_bed_ref + displacement from the two independently, "// &
                    "so both must describe the same bedrock."//new_line("a")//"On a restart this usually means the "// &
                    "restart of "//trim(c%name)//" came from a run with another topography ("//trim(dom%topo)//").", &
                    stat)
                return
            end if
        end do
    end subroutine domain_check_reference

    ! =========================================================================
    ! Wiring: nodes, resolution, exchange
    ! =========================================================================

    function domain_nodes(dom, clk) result(nodes)
        ! The wire nodes of the domain: the hub, then every created instance
        ! (cadence step from the clock). The nodes point into dom (target).
        type(domain_t), target, intent(in) :: dom
        type(clock_t),          intent(in) :: clk
        type(wire_node_t), allocatable :: nodes(:)

        class(component_t), pointer :: c
        integer :: k

        allocate(nodes(0))
        nodes = [nodes, wire_node_states("hub", dom%name, dom%hub_imports, dom%hub_exports, wp, clock_dtt(clk))]
        do k = 1, NCOMP
            c => domain_component(dom, k)
            if (.not. associated(c)) cycle
            nodes = [nodes, wire_node(c, cadence_dt_years(c%cadence, clk))]
        end do
    end function domain_nodes

    subroutine domain_wire(dom, clk, extra, table, check, stat)
        ! Resolve the wiring rows of this domain (wire_resolve: gaps, links,
        ! contracts; stops naming the row) and prime every remapped pair
        ! (wire_prime). extra: the globe's nodes the domain's rows reach
        ! (forcing, global geo; with check = .false.). table: the rows to resolve
        ! (default: the phase 1 table). check = .false.: the created instances
        ! may leave required imports unwired (tests of a partial domain; a run
        ! never sets it). dom%wires holds all rows, active or not; exchange runs
        ! the active ones.
        type(domain_t), target, intent(inout) :: dom
        type(clock_t),          intent(in)    :: clk
        type(wire_node_t),      intent(in),  optional :: extra(:)
        type(wire_t),           intent(in),  optional :: table(:)
        logical,                intent(in),  optional :: check
        integer,                intent(out), optional :: stat

        type(wire_node_t), allocatable :: nodes(:)

        nodes = domain_nodes(dom, clk)
        if (present(check)) then
            if (.not. check) nodes%check = .false.
        end if
        if (present(extra)) nodes = [nodes, extra]
        call wire_resolve(dom%wires, nodes, dom%cpl, table, stat)
        if (present(stat)) then
            if (stat /= 0) then
                if (allocated(dom%wires)) deallocate(dom%wires)    ! half-resolved rows must not run
                return
            end if
        end if
        call wire_prime(dom%wires, dom%cpl)
    end subroutine domain_wire

    subroutine exchange(dom, src, dst, t, clk, init, stat)
        ! The domain's active rows src -> dst, in table order (core-design 3.4):
        ! `src` and `dst` are component kinds ("ice", "geo", ...) or "hub" (a
        ! source, or "hub" as the destination of the ice rows). Each row picks the
        ! two field states, the producer's last update (hub: hub_t; instance: its
        ! cadence) and whether it is due, and the references on the two grids,
        ! and runs wire_exchange (which checks the lag and remaps).
        ! init: cold start, every row counts as lag 0.
        type(domain_t), target, intent(inout) :: dom
        character(len=*),       intent(in)    :: src, dst
        real(dp),               intent(in)    :: t
        type(clock_t),          intent(in)    :: clk
        logical,                intent(in),  optional :: init
        integer,                intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_domain::exchange"
        logical  :: ini
        integer  :: i, s

        if (present(stat)) stat = 0
        ini = .false.
        if (present(init)) ini = init
        if (.not. wired(dom, src//" -> "//dst, stat)) return

        do i = 1, size(dom%wires)
            if (.not. wire_runs(dom%wires(i), src, dst)) cycle
            call run_row(dom, i, t, clk, ini, s)
            if (s /= 0) then
                call fesm_error(proc, "domain "//trim(dom%name)//": row "//wire_label(dom%wires(i))//" failed", &
                                "t = "//to_str(t), stat)
                return
            end if
        end do
    end subroutine exchange

    subroutine exchange_into(dom, dst, t, clk, init, stat)
        ! Every row into `dst`, by source in INTO_ORDER (yelmox couple_to_yelmo:
        ! geo, surface, shelf, climate), each source's rows in table order.
        ! Rows of static producers (LAG_INIT) have run once, at init, and are
        ! skipped by wire_exchange.
        type(domain_t), target, intent(inout) :: dom
        character(len=*),       intent(in)    :: dst
        real(dp),               intent(in)    :: t
        type(clock_t),          intent(in)    :: clk
        logical,                intent(in),  optional :: init
        integer,                intent(out), optional :: stat

        integer :: j, s

        if (present(stat)) stat = 0
        do j = 1, size(INTO_ORDER)
            if (trim(INTO_ORDER(j)) == trim(dst)) cycle
            call exchange(dom, trim(INTO_ORDER(j)), dst, t, clk, init=init, stat=s)
            if (s /= 0) then
                call fesm_error("kryos_domain::exchange_into", "domain "//trim(dom%name)//": rows into "// &
                                trim(dst)//" failed", "", stat)
                return
            end if
        end do
    end subroutine exchange_into

    subroutine exchange_static(dom, t, clk, stat)
        ! Every active LAG_INIT row not yet exchanged, in table order: the
        ! static fields after PASS_INIT (core-design 3.6 item 1: ice z_bed_ref,
        ! H_ice_ref -> geo; geothermal Q_geo, sediments H_sed -> ice; hub
        ! z_srf_ref -> climate).
        type(domain_t), target, intent(inout) :: dom
        real(dp),               intent(in)    :: t
        type(clock_t),          intent(in)    :: clk
        integer,                intent(out), optional :: stat

        integer :: i, s

        if (present(stat)) stat = 0
        if (.not. wired(dom, "static rows", stat)) return

        do i = 1, size(dom%wires)
            associate (w => dom%wires(i))
                if (.not. (w%active .and. w%lag == LAG_INIT .and. .not. w%done .and. .not. w%gather &
                           .and. w%transform /= T_GIA)) cycle
            end associate
            call run_row(dom, i, t, clk, .true., s)
            if (s /= 0) then
                call fesm_error("kryos_domain::exchange_static", "domain "//trim(dom%name)//": row "// &
                                wire_label(dom%wires(i))//" failed", "t = "//to_str(t), stat)
                return
            end if
        end do
    end subroutine exchange_static

    subroutine run_row(dom, i, t, clk, init, stat)
        ! Row i of the domain: the two field states, the producer's last update
        ! and whether it is due, the references on the two grids; wire_exchange
        ! checks the lag and remaps.
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(in)    :: i
        real(dp),               intent(in)    :: t
        type(clock_t),          intent(in)    :: clk
        logical,                intent(in)    :: init
        integer,                intent(out)   :: stat

        type(field_state_t), pointer :: s_exp, d_imp
        type(grid_ref_t),    pointer :: rs, rd
        real(dp) :: t_last
        logical  :: due
        integer  :: ir

        associate (w => dom%wires(i))
            call kind_states(dom, w%src_kind, exports=s_exp)
            call kind_states(dom, w%dst_kind, imports=d_imp)
            call producer_time(dom, w%src_kind, clk, t_last, due)

            rs => null();  rd => null()
            ir = domain_ref_index(dom, w%src_grid)
            if (ir > 0) rs => dom%refs(ir)
            ir = domain_ref_index(dom, w%dst_grid)
            if (ir > 0) rd => dom%refs(ir)

            call wire_exchange(w, s_exp, d_imp, dom%cpl, t, t_last, ref_src=rs, ref_dst=rd, &
                               due=due, init=init, stat=stat)
        end associate
    end subroutine run_row

    logical function wired(dom, what, stat)
        ! The rows are resolved (domain_wire), else an error naming `what`.
        type(domain_t),   intent(in) :: dom
        character(len=*), intent(in) :: what
        integer,          intent(inout), optional :: stat
        wired = allocated(dom%wires)
        if (wired) return
        call fesm_error("kryos_domain::exchange", "domain "//trim(dom%name)//": the wiring is not resolved (domain_wire)", &
                        "row group "//trim(what), stat)
    end function wired

    subroutine kind_states(dom, kind, imports, exports)
        ! The field states of a kind: the hub's, or the created instance's.
        type(domain_t), target, intent(in) :: dom
        character(len=*),       intent(in) :: kind
        type(field_state_t), pointer, intent(out), optional :: imports, exports

        class(component_t), pointer :: c

        if (present(imports)) imports => null()
        if (present(exports)) exports => null()
        if (trim(kind) == "hub") then
            if (present(imports)) imports => dom%hub_imports
            if (present(exports)) exports => dom%hub_exports
            return
        end if
        c => domain_component(dom, domain_comp_index(kind))
        if (.not. associated(c)) &
            call fesm_error("kryos_domain::exchange", "domain "//trim(dom%name)//": no "//trim(kind)//" instance", &
                            "a resolved row needs it")
        if (present(imports)) imports => c%imports
        if (present(exports)) exports => c%exports
    end subroutine kind_states

    subroutine producer_time(dom, kind, clk, t_last, due)
        ! The producer's last update and whether it is due at the clock's time.
        ! The hub is refreshed every step (hub_t); an instance by its cadence.
        type(domain_t), target, intent(in) :: dom
        character(len=*),       intent(in) :: kind
        type(clock_t),          intent(in) :: clk
        real(dp),               intent(out) :: t_last
        logical,                intent(out) :: due

        class(component_t), pointer :: c

        if (trim(kind) == "hub") then
            t_last = dom%hub_t
            due    = .true.
            return
        end if
        c => domain_component(dom, domain_comp_index(kind))
        t_last = cadence_t_last(c%cadence, clk)
        due    = cadence_due(c%cadence, clk)
    end subroutine producer_time

    ! =========================================================================
    ! Hub output (core-design section 6): <domain>_hub.nc
    ! =========================================================================

    subroutine domain_hub_write(dom, file, t, init)
        ! The hub's current geometry at t: z_bed, H_ice, z_srf, f_grnd, z_sl
        ! (htopo_write_step) plus f_ice and the four-value mask (kryos_mask on
        ! the hub exports). init: first write, creates the file with the regions
        ! and basins (htopo_write_init). The caller (output, core C8) owns the
        ! schedule `dom%out_2D` and the file name.
        type(domain_t),   intent(in) :: dom
        character(len=*), intent(in) :: file
        real(dp),         intent(in) :: t
        logical,          intent(in) :: init

        integer :: ncid, n, nx, ny
        integer, allocatable :: mask(:,:)

        nx = dom%hub%nx
        ny = dom%hub%ny
        if (init) call htopo_write_init(dom%hub, trim(file), t)
        call htopo_write_step(dom%hub, trim(file), t)

        mask = kryos_mask(dom%hub_exports%f(dom%hh%h_f_ice)%v(:,:,1), dom%hub_exports%f(dom%hh%h_f_grnd)%v(:,:,1))
        call nc_open(trim(file), ncid, writable=.true.)
        n = nc_time_index(trim(file), "time", t, ncid)
        call nc_write(trim(file), "f_ice", real(dom%hub_exports%f(dom%hh%h_f_ice)%v(:,:,1), wp), &
                      dim1="xc", dim2="yc", dim3="time", start=[1,1,n], count=[nx,ny,1], ncid=ncid, &
                      units="1", long_name="Ice-covered fraction", prec_abs=1e-4_wp)
        call nc_write(trim(file), "mask", mask, dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], count=[nx,ny,1], ncid=ncid, units="1", &
                      long_name="Surface type (0 ocean, 1 land, 2 grounded ice, 3 floating ice)")
        call nc_close(ncid)
    end subroutine domain_hub_write

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    logical function failed(stat)
        integer, intent(in), optional :: stat
        failed = .false.
        if (present(stat)) failed = stat /= 0
    end function failed

end module kryos_domain
