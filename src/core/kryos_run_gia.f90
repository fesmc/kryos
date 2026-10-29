submodule (kryos_run) kryos_run_gia
    ! The global ice field's side of the run sequence (phase4-plan 4, 5.1, 7.1, 8.1; core-design 3.2, 3.4):
    ! for a host-resolved global geo (vilma2) the core, not the wiring rows, turns the domains' ice into the
    ! Earth's load and builds the reference it is measured against.
    !
    !   gia_reference   after PASS_INIT: the domains' static exports (footprint, H_ice_pd, z_bed_ref), the
    !                   background X and the external history are merged into the reference of the Earth
    !                   (kryos_globe_ice); the global geo's imports z_bed_ref, H_ice_ref, host_w,
    !                   host_Ig_ref, host_C_ref at t0
    !   gia_gather      the host terms of the ice of the previous step (lag 1; at the cold start: the final
    !                   ice at t0) into the imports H_ice_ext, host_w, host_Ig, host_C, host_F
    !   gia_seed        cold start, between the ice's PASS_SPINUP and the reconcile: gather at t0, the Earth's
    !                   PASS_STATE (its own spin-up under the t0 load, then the dt = 0 seed), the anomalies
    !                   into the ports
    !
    ! The scatter at a step is the generic globe -> geo rows of bsl_to_geo and the port's update (block 3 of
    ! the sequence): the rows geo_g -> geo are rank 2 for a port. The T_GIA rows of the wiring table exist for
    ! resolution, logging and the lag check: gia_gather checks each domain's row against the ice's last
    ! update (wire_lag_status) before it reads the exports.
    !
    ! Domains are visited in name order (run%dom_order); the merge itself is ordered by the ice grid names and
    ! does not depend on the order (kryos_globe_ice).

    use fesm_clock,      only: cadence_t_last, cadence_due
    use fesm_fields,     only: fields_get, fields_put
    use coordinates,     only: grid_class
    use coords,          only: grid_cdo_read_desc
    use kryos_component, only: PASS_STATE, NAME_LEN
    use kryos_geo,       only: geo_present
    use kryos_globe,     only: globe_pass
    use kryos_globe_ice, only: globe_ice_add_domain, globe_ice_reference_begin, globe_ice_reference_add, &
                               globe_ice_reference_end, globe_ice_gather_begin, globe_ice_gather_add, &
                               globe_ice_gather_end
    use kryos_wiring,    only: wire_lag_status, T_GIA
    use kryos_domain,    only: domain_t, domain_pass, C_ICE, C_GEO

    implicit none

contains

    ! =========================================================================
    ! Reference (plan 5.1)
    ! =========================================================================

    module subroutine gia_reference(run)
        type(run_t), target, intent(inout) :: run

        character(len=*), parameter :: proc = "kryos_run::gia_reference"
        real(dp), allocatable :: fp(:,:), hpd(:,:), zbr(:,:)
        type(grid_class) :: g
        real(dp) :: t0
        integer :: i, d

        if (.not. run%globe%has_gi) return
        t0 = clock_time(run%clk)

        associate (gi => run%globe%gi, cpl => run%globe%cpl)
            do i = 1, run%n_dom
                d = run%dom_order(i)
                associate (ice => run%dom(d)%ice)
                    allocate(fp(ice%c%nx, ice%c%ny), hpd(ice%c%nx, ice%c%ny), zbr(ice%c%nx, ice%c%ny))
                    call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                    call fields_get(ice%c%exports, ice%c%h_H_ice_pd,  hpd)
                    call fields_get(ice%c%exports, ice%c%h_z_bed_ref, zbr)
                    call grid_cdo_read_desc(g, trim(run%dom(d)%grid(C_ICE)), trim(run%cfg%map_fldr))
                    call globe_ice_add_domain(gi, cpl, trim(run%dom(d)%name), g, fp, hpd, zbr)
                    deallocate(fp, hpd, zbr)
                end associate
            end do

            call globe_ice_reference_begin(gi, cpl)
            do i = 1, run%n_dom
                d = run%dom_order(i)
                call globe_ice_reference_add(gi, cpl, trim(run%dom(d)%name), run%dom(d)%ice)
            end do
            call globe_ice_reference_end(gi, cpl)

            associate (c => run%globe%geo%c)
                call fields_put(c%imports, c%hi_z_bed_ref,   gi%z_bed_ref,   t1=t0)
                call fields_put(c%imports, c%hi_H_ice_ref,   gi%H_ice_ref,   t1=t0)
                call fields_put(c%imports, c%hi_host_w,      gi%host_w,      t1=t0)
                call fields_put(c%imports, c%hi_host_Ig_ref, gi%host_Ig_ref, t1=t0)
                call fields_put(c%imports, c%hi_host_C_ref,  gi%host_C_ref,  t1=t0)
            end associate
        end associate
        associate (unused => proc); end associate
    end subroutine gia_reference

    ! =========================================================================
    ! Gather (plan 4.1)
    ! =========================================================================

    module subroutine gia_gather(run, init)
        type(run_t), target, intent(inout) :: run
        logical,             intent(in)    :: init

        real(dp), allocatable :: H(:,:), fg(:,:), zb(:,:), zs(:,:), fp(:,:)
        real(dp) :: t, t_load, t_dom
        integer  :: i, d

        if (.not. run%globe%has_gi) return
        t = clock_time(run%clk)

        ! the ice rows: published at the producer's last update (lag 1; the cold start: now)
        t_load = t
        do i = 1, run%n_dom
            d = run%dom_order(i)
            t_dom = gia_row_check(run, d, t, init)
            if (i == 1) then
                t_load = t_dom
            else if (t_dom /= t_load) then
                call fesm_error("kryos_run::gia_gather", "the ice instances were last updated at different times", &
                                "domain "//trim(run%dom(run%dom_order(1))%name)//" at t = "//to_str(t_load)//", domain "// &
                                trim(run%dom(d)%name)//" at t = "//to_str(t_dom)//": the load of the Earth is one "// &
                                "time; give the ice instances one cadence")
            end if
        end do

        associate (gi => run%globe%gi, cpl => run%globe%cpl)
            call globe_ice_gather_begin(gi, cpl, t_load)
            do i = 1, run%n_dom
                d = run%dom_order(i)
                associate (ice => run%dom(d)%ice)
                    allocate(H(ice%c%nx, ice%c%ny), fg(ice%c%nx, ice%c%ny), zb(ice%c%nx, ice%c%ny), &
                             zs(ice%c%nx, ice%c%ny), fp(ice%c%nx, ice%c%ny))
                    call fields_get(ice%c%exports, ice%c%h_H_ice,     H)
                    call fields_get(ice%c%exports, ice%c%h_f_grnd,    fg)
                    call fields_get(ice%c%exports, ice%c%h_z_bed,     zb)
                    call fields_get(ice%c%exports, ice%c%h_z_sl,      zs)
                    call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                    call globe_ice_gather_add(gi, cpl, trim(run%dom(d)%name), ice, H, fg, zb, zs, fp)
                    deallocate(H, fg, zb, zs, fp)
                end associate
            end do
            call globe_ice_gather_end(gi, cpl)

            associate (c => run%globe%geo%c)
                call fields_put(c%imports, c%hi_H_ice_ext, gi%H_ice_ext, t1=t_load)
                call fields_put(c%imports, c%hi_host_w,    gi%host_w,    t1=t_load)
                call fields_put(c%imports, c%hi_host_Ig,   gi%host_Ig,   t1=t_load)
                call fields_put(c%imports, c%hi_host_C,    gi%host_C,    t1=t_load)
                call fields_put(c%imports, c%hi_host_F,    gi%host_F,    t1=t_load)
            end associate
        end associate
    end subroutine gia_gather

    function gia_row_check(run, d, t, init) result(t_last)
        ! The T_GIA rows of domain d into the global geo (host_w stands for them: they have one source and
        ! one lag): active, and the ice's exports (H_ice, f_grnd, z_bed, z_sl) published at its last update,
        ! before t (lag 1) or at t (the cold start). Returns that time. Stops naming the row.
        type(run_t), target, intent(in) :: run
        integer,             intent(in) :: d
        real(dp),            intent(in) :: t
        logical,             intent(in) :: init
        real(dp) :: t_last

        character(len=:), allocatable :: msg
        integer :: i

        t_last = cadence_t_last(run%dom(d)%ice%c%cadence, run%clk)
        do i = 1, size(run%dom(d)%wires)
            associate (w => run%dom(d)%wires(i))
                if (w%transform /= T_GIA .or. trim(w%dst_field) /= "host_w") cycle
                if (.not. w%active) then
                    call fesm_error("kryos_run::gia_gather", "domain "//trim(run%dom(d)%name)// &
                                    ": the gather row into `host_w` is inactive", trim(w%reason))
                end if
                msg = wire_lag_status(w, run%dom(d)%ice%c%exports, run%globe%geo%c%imports, t, t_last, &
                                      due=cadence_due(run%dom(d)%ice%c%cadence, run%clk), init=init)
                if (len(msg) > 0) then
                    call fesm_error("kryos_run::gia_gather", "domain "//trim(run%dom(d)%name)//": the load of the ice: "// &
                                    msg, "t = "//to_str(t))
                end if
                return
            end associate
        end do
        call fesm_error("kryos_run::gia_gather", "domain "//trim(run%dom(d)%name)//" has no gather row into `host_w`", &
                        "the wiring table's T_GIA rows (ice -> geo_g)")
    end function gia_row_check

    ! =========================================================================
    ! Seed (plan 8.1, 7.4)
    ! =========================================================================

    module subroutine gia_seed(run)
        type(run_t), target, intent(inout) :: run

        integer :: d

        if (.not. run%globe%has_gi) return
        call gia_gather(run, init=.true.)                                   ! the t0 load: final ice
        call globe_pass(run%globe, "geo", run%clk, PASS_STATE)              ! spin-up under it, the dt = 0 seed
        do d = 1, run%n_dom
            if (.not. domain_present(run%dom(d), C_GEO)) cycle
            call bsl_to_geo(run, d, init=.true.)                            ! rows geo_g -> geo: anomalies, dphi, bsl
            call domain_pass(run%dom(d), C_GEO, run%clk, PASS_STATE)        ! the port holds them
        end do
    end subroutine gia_seed

end submodule kryos_run_gia
