submodule (kryos_run) kryos_run_sequence
    ! The coupling sequence of the run (core-design 3.2, 3.4, 3.6): the cold
    ! start (yelmox kryos_startup.f90 domain_init_state, 131-195) and one step
    ! (yelmox.f90 main loop, 176-208), for N domains. The order is the
    ! bit-identity contract with yelmox; each block names the yelmox call it
    ! stands for. Presence is tested (domain_present, *_present), never a
    ! backend name; the wiring rows check every lag (kryos_wiring).
    !
    ! Globe <-> domain traffic in one place each:
    !   bsl:     globe geo -> every geo (exchange_globe), the dV_bsl of the
    !            domains gathered, reduced, scattered back (reduce_bsl);
    !   forcing: M_ice of the domains gathered (lag 1), the climate's derived
    !            coefficient (climate_to_forcing, lag 1; at t0 after the
    !            climate's PASS_BOUNDARY), coefficients into the climate
    !            (forcing_to_climate, lag 0).

    use fesm_clock,      only: clock_to_years, cadence_t_last, cadence_due, cadence_dt_years
    use kryos_component, only: PASS_GEOMETRY, PASS_BOUNDARY, PASS_STATE, PASS_FINAL
    use kryos_wiring,    only: wire_node
    use kryos_geo,       only: geo_scatter
    use kryos_globe,     only: gather_src_t, gather_src, gather_into, globe_update, globe_reduce, &
                               exchange_globe, exchange_to_globe
    use kryos_domain,    only: domain_present, domain_due, domain_update, domain_time, exchange, exchange_into, &
                               domain_hub_update, domain_check_reference, &
                               C_ICE, C_GEO, C_CLIMATE, C_SURFACE, C_SHELF

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

                ! c) the ice state: yelmo_init_state (robin-cold), domain_init_ice
                call domain_pass(dom, C_ICE, run%clk, PASS_STATE)

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
    ! One step (core-design 3.4; yelmox.f90:176-208)
    ! =========================================================================

    module subroutine run_sequence(run)
        ! Steps 1-4 at the clock's time t (advanced by run_step); step 5 (kill
        ! switch, output, restart) is run_step's. `mark` is in domain_update.
        type(run_t), target, intent(inout) :: run

        real(dp) :: t
        logical, allocatable :: ran(:)
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

        ! 2. hooks: none in phase 1 (yelmox step_relax, step_optimize)

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
                call timer_step(run%tmr, TT_COUPLE_ICE)
                if (domain_due(dom, C_ICE, run%clk)) call domain_update(dom, C_ICE, run%clk)   ! step_icesheet
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
