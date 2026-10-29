submodule (kryos_run) kryos_run_sequence
    ! The coupling sequence of the run (core-design 3.2-3.4, 3.6): the cold
    ! start (yelmox kryos_startup.f90 domain_init_state, 131-195), the restore
    ! from a bundle (domain_startup / domain_restart_read, 95-103, 470-530)
    ! and one step (yelmox.f90 main loop, 176-208), for N domains. The order
    ! is the bit-identity contract with yelmox (and, for the restore, with the
    ! continuous run); each block names the yelmox call it stands for.
    ! Presence is tested (domain_present, *_present), never a backend name;
    ! the wiring rows check every lag (kryos_wiring).
    !
    ! Globe <-> domain traffic in one place each:
    !   bsl:     globe geo -> every geo (exchange_globe), the dV_bsl of the
    !            domains gathered, reduced, scattered back (reduce_bsl);
    !   forcing: M_ice of the domains gathered (lag 1), the climate's derived
    !            coefficient (climate_to_forcing, lag 1; at t0 after the
    !            climate's PASS_BOUNDARY), coefficients into the climate
    !            (forcing_to_climate, lag 0).

    use fesm_clock,      only: clock_to_years, clock_time_rel, clock_time_elapsed, clock_dtt, cadence_t_last, &
                               cadence_due, cadence_dt_years
    use kryos_component, only: PASS_GEOMETRY, PASS_BOUNDARY, PASS_STATE, PASS_SPINUP, PASS_FINAL
    use kryos_wiring,    only: wire_node
    use kryos_geo,       only: geo_scatter
    use kryos_globe,     only: gather_src_t, gather_src, gather_into, globe_update, globe_reduce, &
                               exchange_globe, exchange_to_globe
    use kryos_domain,    only: domain_present, domain_due, domain_prepare, domain_update, domain_time, exchange, exchange_into, &
                               domain_hub_update, domain_check_reference, &
                               C_ICE, C_GEO, C_CLIMATE, C_SURFACE, C_SHELF, RESTORE_ORDER
    use kryos_experiment, only: hook_time_t, experiments_hook, HOOK_STEP_BEGIN, HOOK_ICE_BND, &
                                HOOK_COLD_PRE_STATE, HOOK_COLD_POST_STATE
    use kryos_restart,   only: GLOBE_ORDER, restart_read_clock, restart_check_manifest, restart_read_global, &
                               restart_read_instance, restart_read_experiments

    implicit none

contains

    ! =========================================================================
    ! Cold start (core-design 3.2)
    ! =========================================================================

    module subroutine run_cold_start(run)
        ! After PASS_INIT and the static rows (run_pass_init). Per domain, in
        ! the order of yelmox domain_init_state; the bed-and-sea-level block a)
        ! runs for every domain first, so the global reduction of the initial
        ! volume changes (3.6 item 2) sits between geo's PASS_GEOMETRY and the
        ! geo -> ice rows.
        type(run_t), target, intent(inout) :: run

        real(dp) :: t0
        logical, allocatable :: ran(:)
        integer :: d

        t0 = clock_time(run%clk)

        ! a) bed and sea level (domain_init_isostasy): the geo state from the
        !    ice's reference and geometry and the global bsl
        allocate(ran(run%n_dom))
        do d = 1, run%n_dom
            associate (dom => run%dom(d))
                ran(d) = domain_present(dom, C_GEO)
                if (.not. ran(d)) cycle
                call exchange(dom, "ice", "geo", t0, run%clk, init=.true.)        ! H_ice (reference rows: static)
                call bsl_to_geo(run, d, init=.true.)                              ! bsl after bsl_update(time_rel0)
                call domain_pass(dom, C_GEO, run%clk, PASS_GEOMETRY)              ! isos_init_ref + isos_init_state
                call domain_check_reference(dom)                                  ! check_isostasy_reference, 1 m
            end associate
        end do
        call reduce_bsl(run, t0, ran, init=.true.)                                ! isos_init_state's bsl change

        do d = 1, run%n_dom
            associate (dom => run%dom(d))
                call exchange(dom, "geo", "ice", t0, run%clk, init=.true.)        ! couple_isostasy_to_yelmo
                call domain_pass(dom, C_ICE, run%clk, PASS_GEOMETRY)              ! calc_ytopo_diagnostic

                ! b) boundary fields from the hub (couple_yelmo_to_htopo, steps with init = .true.)
                call domain_hub_update(dom, t0, run%clk, init=.true.)
                if (domain_present(dom, C_CLIMATE)) then                          ! step_climate(init)
                    call exchange(dom, "hub", "climate", t0, run%clk, init=.true.)
                    call forcing_to_climate(run, d, init=.true.)
                    call domain_pass(dom, C_CLIMATE, run%clk, PASS_BOUNDARY)
                    call climate_to_forcing(run, d, init=.true.)                  ! 3.6 item 3: dTo at t0
                end if
                if (domain_present(dom, C_SURFACE)) then                          ! step_surface(init)
                    call exchange(dom, "hub", "surface", t0, run%clk, init=.true.)
                    call exchange(dom, "climate", "surface", t0, run%clk, init=.true.)
                    call domain_pass(dom, C_SURFACE, run%clk, PASS_BOUNDARY)      ! spin-up + one cycle
                end if
                if (domain_present(dom, C_SHELF)) then                            ! step_shelf
                    call exchange(dom, "hub", "shelf", t0, run%clk, init=.true.)
                    call exchange(dom, "climate", "shelf", t0, run%clk, init=.true.)
                    call domain_pass(dom, C_SHELF, run%clk, PASS_BOUNDARY)
                end if
                call exchange_into(dom, "ice", t0, run%clk, init=.true.)          ! couple_to_yelmo
                call domain_prepare(dom, C_ICE, run%clk, init=.true.)             ! imports -> Yelmo's bnd

                ! c) the ice state: yelmo_init_state (robin-cold), then domain_init_ice;
                !    the experiments' hooks around the state (the cb_ref guess, recon, ...)
                call run_hook(run, d, HOOK_COLD_PRE_STATE)
                call domain_pass(dom, C_ICE, run%clk, PASS_STATE)
                call run_hook(run, d, HOOK_COLD_POST_STATE)
                call domain_pass(dom, C_ICE, run%clk, PASS_SPINUP)

                ! d) the surface on the final initial geometry: one more cycle, no spin-up
                call domain_hub_update(dom, t0, run%clk, init=.true.)            ! couple_yelmo_to_htopo
                if (domain_present(dom, C_SURFACE)) then
                    call exchange(dom, "hub", "surface", t0, run%clk, init=.true.)
                    call domain_pass(dom, C_SURFACE, run%clk, PASS_FINAL)         ! step_surface
                    call exchange(dom, "surface", "ice", t0, run%clk, init=.true.) ! couple_surface_to_yelmo
                end if
            end associate
        end do
    end subroutine run_cold_start

    ! =========================================================================
    ! Restore (core-design 3.3, 3.6 item 4)
    ! =========================================================================

    module subroutine run_restore(run)
        ! After PASS_INIT and the static rows (run_pass_init), instead of the
        ! cold start; restore, never step. The bundle is [kryos] restart
        ! (relative to the run directory, or absolute).
        !
        !   1. the clock, then every instance's cadence (continuation: the
        !      bundle time is time_init; a branch, another time, restores the
        !      clock only: a run with instances stops);
        !   2. the manifest against the run's instances;
        !   3. the globe: forcing, then geo global (bsl_now, A_ocean_now);
        !   4. per domain, in RESTORE_ORDER: ice; then the ice -> geo rows and
        !      the global bsl into geo (FastIsostasy rebuilds its state from the
        !      restored H_ice and bsl, yelmox domain_restart_read); geo, climate,
        !      surface, shelf, geothermal, sediments. Each instance: its
        !      restart_read, its imports and exports as the bundle holds them,
        !      PASS_RESTART (kryos_restart). Then the domain's stateful
        !      experiments (restart_read; no hook runs in a restore). The rows into the other instances
        !      do not run: their imports are restored as received, and the
        !      static (LAG_INIT) rows stay done from run_pass_init;
        !   5. per domain, domain_check_reference (the restored references
        !      against the domain's) and the hub (derived, never stored) from
        !      the restored ice.
        type(run_t), target, intent(inout) :: run

        character(len=:), allocatable :: dir
        real(dp) :: t
        integer  :: d, j, k

        dir = strip_slash(run%restart)

        ! 1. clock and cadences
        call restart_read_clock(dir, run%clk, run%globe, run%dom)
        if (run%clk%branched .and. run_has_instances(run)) then
            call fesm_error("kryos_run::run_restore", "a branch from a restart bundle restores the clock only", &
                            "[kryos] restart = "//trim(run%restart)//": the bundle time is not [clock] time_init = "// &
                            year_label(clock_time(run%clk))//" yr. Set time_init to the bundle time to continue "// &
                            "the run; a branch of the instances (their field times moved to time_init) is not "// &
                            "implemented")
        end if
        t = clock_time(run%clk)

        ! 2. the manifest
        call restart_check_manifest(dir, run%globe, run%dom, run%exps)

        ! 3. the globe (yelmox bsl_startup, coefs_restart_read)
        do j = 1, size(GLOBE_ORDER)
            call restart_read_global(dir, run%globe, trim(GLOBE_ORDER(j)), run%clk)
        end do

        ! 4. the domains (domain_restart_read)
        do d = 1, run%n_dom
            associate (dom => run%dom(d))
                do j = 1, size(RESTORE_ORDER)
                    k = RESTORE_ORDER(j)
                    if (.not. domain_present(dom, k)) cycle
                    if (k == C_GEO) then                                      ! geo rebuilds on the restored ice, bsl
                        call exchange(dom, "ice", "geo", t, run%clk, init=.true.)
                        call bsl_to_geo(run, d, init=.true.)
                    end if
                    call restart_read_instance(dir, dom, k, run%clk)
                end do
                call restart_read_experiments(dir, dom, run%exps(d), run%clk)         ! stateful ones; no hook runs

                ! 5. reference check (check_isostasy_reference), hub (couple_yelmo_to_htopo)
                call domain_check_reference(dom)
                call domain_hub_update(dom, t, run%clk)
            end associate
        end do

        call log_msg(LOG_INFO, "kryos_run", "restored from "//dir//" at time = "//year_label(t)//" yr")
    end subroutine run_restore

    logical function run_has_instances(run)
        ! Does the run have any instance (globe or domain)?
        type(run_t), target, intent(in) :: run
        integer :: d, k

        run_has_instances = forcing_present(run%globe%forcing) .or. geo_present(run%globe%geo)
        do d = 1, run%n_dom
            do k = 1, NCOMP
                if (associated(domain_component(run%dom(d), k))) run_has_instances = .true.
            end do
        end do
    end function run_has_instances

    function strip_slash(dir) result(s)
        ! A folder name without trailing slashes ("restart-500/" -> "restart-500").
        character(len=*), intent(in) :: dir
        character(len=:), allocatable :: s
        s = trim(dir)
        do while (len(s) > 1)
            if (s(len(s):len(s)) /= "/") exit
            s = s(1:len(s)-1)
        end do
    end function strip_slash

    ! =========================================================================
    ! One step (core-design 3.4; yelmox.f90:176-208)
    ! =========================================================================

    module subroutine run_sequence(run)
        ! Steps 1-4 at the clock's time t (advanced by run_step); step 5 (kill
        ! switch, output, restart) is run_step's. `mark` is in domain_update.
        type(run_t), target, intent(inout) :: run

        real(dp) :: t
        logical, allocatable :: ran(:)
        logical :: ice_due
        integer :: d

        t = clock_time(run%clk)

        ! 1. globe, start of the step: bsl_update(time_rel) (yelmox.f90:180),
        !    then coefs_update(time_rel, V_ice) (:185-186) on the previous step's
        !    ice mass and derived coefficients
        if (geo_present(run%globe%geo)) call globe_update(run%globe, "geo", run%clk)
        if (forcing_present(run%globe%forcing)) then
            call gather_forcing(run)
            do d = 1, run%n_dom
                call climate_to_forcing(run, d, init=.false.)
            end do
            call globe_update(run%globe, "forcing", run%clk)
        end if
        call timer_step(run%tmr, TT_FORCING)

        ! 2. experiment hooks, every domain in order (yelmox step_relax, step_optimize)
        do d = 1, run%n_dom
            call run_hook(run, d, HOOK_STEP_BEGIN)
        end do
        call timer_step(run%tmr, TT_EXPER)

        ! 3. geo: every domain on its cadence (step_isostasy, :194), then the
        !    global reduction (gather -> reduce -> scatter, decision 5)
        allocate(ran(run%n_dom))
        do d = 1, run%n_dom
            associate (dom => run%dom(d))
                ran(d) = domain_due(dom, C_GEO, run%clk)
                if (.not. ran(d)) cycle
                call exchange(dom, "ice", "geo", t, run%clk)                       ! H_ice of the previous step
                call bsl_to_geo(run, d, init=.false.)                              ! bsl at the start of the step
                call domain_update(dom, C_GEO, run%clk)                            ! isos_update
            end associate
        end do
        call reduce_bsl(run, t, ran, init=.false.)
        call timer_step(run%tmr, TT_GEO)

        ! 4. per domain: ice, hub, climate, surface, shelf (:196-208)
        do d = 1, run%n_dom
            associate (dom => run%dom(d))
                call exchange_into(dom, "ice", t, run%clk)                        ! couple_to_yelmo
                ice_due = domain_due(dom, C_ICE, run%clk)
                if (ice_due) then
                    call domain_prepare(dom, C_ICE, run%clk)                      ! imports -> Yelmo's bnd
                    call run_hook(run, d, HOOK_ICE_BND)                           ! experiments on the prepared boundary
                end if
                call timer_step(run%tmr, TT_COUPLE_ICE)
                if (ice_due) call domain_update(dom, C_ICE, run%clk)              ! step_icesheet
                call timer_step(run%tmr, TT_ICE)
                call domain_hub_update(dom, t, run%clk)                            ! couple_yelmo_to_htopo
                call timer_step(run%tmr, TT_HUB)
                if (domain_due(dom, C_CLIMATE, run%clk)) then                      ! step_climate (dt_clim)
                    call exchange(dom, "hub", "climate", t, run%clk)
                    call forcing_to_climate(run, d, init=.false.)
                    call domain_update(dom, C_CLIMATE, run%clk)
                end if
                call timer_step(run%tmr, TT_CLIMATE)
                if (domain_due(dom, C_SURFACE, run%clk)) then                      ! step_surface
                    call exchange(dom, "hub", "surface", t, run%clk)
                    call exchange(dom, "climate", "surface", t, run%clk)
                    call domain_update(dom, C_SURFACE, run%clk)
                end if
                call timer_step(run%tmr, TT_SURFACE)
                if (domain_due(dom, C_SHELF, run%clk)) then                        ! step_shelf
                    call exchange(dom, "hub", "shelf", t, run%clk)
                    call exchange(dom, "climate", "shelf", t, run%clk)
                    call domain_update(dom, C_SHELF, run%clk)
                end if
                call timer_step(run%tmr, TT_SHELF)
            end associate
        end do
    end subroutine run_sequence

    subroutine run_hook(run, d, id)
        ! Hook `id` of the experiments of domain d, at the clock's time: the
        ! run's step [t1 - dtt, t1] (the cold start: t0 = t1), the forcing
        ! time, and the elapsed time of the clock (core-design 3.9).
        type(run_t), target, intent(inout) :: run
        integer,             intent(in)    :: d, id

        type(hook_time_t) :: st
        logical :: cold

        if (.not. allocated(run%exps(d)%s)) return
        cold = id == HOOK_COLD_PRE_STATE .or. id == HOOK_COLD_POST_STATE
        st%t1        = clock_time(run%clk)
        st%t0        = st%t1
        if (.not. cold) st%t0 = clock_to_years(run%clk, run%clk%n - run%clk%dt)
        st%t_eval    = clock_time_rel(run%clk)
        st%t_elapsed = clock_time_elapsed(run%clk)
        st%dtt       = clock_dtt(run%clk)
        st%clk       => run%clk
        call experiments_hook(run%exps(d), run%dom(d), id, st)
    end subroutine run_hook

    ! =========================================================================
    ! Globe <-> domain
    ! =========================================================================

    subroutine bsl_to_geo(run, d, init)
        ! The global geo's exports (bsl, A_ocn) into domain d's geo imports.
        type(run_t), target, intent(inout) :: run
        integer,             intent(in)    :: d
        logical,             intent(in)    :: init

        if (.not. geo_present(run%globe%geo)) return
        if (.not. domain_present(run%dom(d), C_GEO)) return
        call exchange_globe(run%globe, "geo", run%dom(d)%wires, "geo", run%dom(d)%geo%c%imports, run%clk, init=init)
    end subroutine bsl_to_geo

    subroutine reduce_bsl(run, t, ran, init)
        ! contracts-design 6: the dV_bsl of the domains whose geo ran at t,
        ! gathered in name order (stale ones count 0), the global reduction
        ! bsl -= sum(dV)/A_ocn (geo_reduce), and the reduced bsl back into those
        ! domains' geo, each exporting it (geo_scatter): what z_sl at the ice
        ! is formed with. Nothing without a global geo or a gather source.
        type(run_t), target, intent(inout) :: run
        real(dp),            intent(in)    :: t
        logical,             intent(in)    :: ran(:)
        logical,             intent(in)    :: init

        type(gather_src_t), allocatable :: srcs(:)
        integer :: i, d

        if (.not. geo_present(run%globe%geo)) return
        allocate(srcs(run%n_dom))
        do i = 1, run%n_dom
            d = run%dom_order(i)
            if (.not. domain_present(run%dom(d), C_GEO)) cycle
            srcs(i) = gather_src(run%dom(d)%wires, wire_node(run%dom(d)%geo%c, &
                                 cadence_dt_years(run%dom(d)%geo%c%cadence, run%clk), check=.false.), "geo")
        end do
        if (.not. any(srcs%active)) return

        call gather_into(run%globe, "geo", srcs, t)
        call globe_reduce(run%globe, run%clk)
        do d = 1, run%n_dom
            if (.not. ran(d)) cycle
            call bsl_to_geo(run, d, init)
            call geo_scatter(run%dom(d)%geo, domain_time(run%dom(d), C_GEO, run%clk, init=init))
        end do
    end subroutine reduce_bsl

    subroutine gather_forcing(run)
        ! M_ice of every domain's ice at the previous step (lag 1), summed in
        ! name order into the forcing's import (yelmox.f90:185: the feedback
        ! variable of coefs_update).
        type(run_t), target, intent(inout) :: run

        type(gather_src_t), allocatable :: srcs(:)
        integer :: i, d

        allocate(srcs(run%n_dom))
        do i = 1, run%n_dom
            d = run%dom_order(i)
            if (.not. domain_present(run%dom(d), C_ICE)) cycle
            srcs(i) = gather_src(run%dom(d)%wires, wire_node(run%dom(d)%ice%c, &
                                 cadence_dt_years(run%dom(d)%ice%c%cadence, run%clk), check=.false.), "forcing")
        end do
        call gather_into(run%globe, "forcing", srcs, clock_to_years(run%clk, run%clk%n - run%clk%dt))
    end subroutine gather_forcing

    subroutine climate_to_forcing(run, d, init)
        ! Domain d's climate -> forcing rows (derived coefficients, decision 1
        ! of C9: the wiring row, lag 1): at t0 after the climate's
        ! PASS_BOUNDARY (init), then at the start of every step.
        type(run_t), target, intent(inout) :: run
        integer,             intent(in)    :: d
        logical,             intent(in)    :: init

        if (.not. forcing_present(run%globe%forcing)) return
        if (.not. domain_present(run%dom(d), C_CLIMATE)) return
        associate (c => run%dom(d)%climate%c)
            call exchange_to_globe(run%globe, run%dom(d)%wires, "climate", c%exports, &
                                   cadence_t_last(c%cadence, run%clk), "forcing", run%clk, &
                                   due=cadence_due(c%cadence, run%clk), init=init)
        end associate
    end subroutine climate_to_forcing

    subroutine forcing_to_climate(run, d, init)
        ! The forcing's coefficients into domain d's climate imports coef_<name>
        ! (lag 0: coefs_update of this step).
        type(run_t), target, intent(inout) :: run
        integer,             intent(in)    :: d
        logical,             intent(in)    :: init

        if (.not. forcing_present(run%globe%forcing)) return
        if (.not. domain_present(run%dom(d), C_CLIMATE)) return
        call exchange_globe(run%globe, "forcing", run%dom(d)%wires, "climate", run%dom(d)%climate%c%imports, &
                            run%clk, init=init)
    end subroutine forcing_to_climate

end submodule kryos_run_sequence
