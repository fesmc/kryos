module test_globe_toy
    ! A regional geo that only declares its fields (the test puts its exports).
    use precision,       only: sp, dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_geo,       only: geo_t, geo_opt_t

    implicit none
    private
    public :: toy_geo_t

    type, extends(geo_t) :: toy_geo_t
    contains
        procedure :: init          => toy_init
        procedure :: update        => toy_update
        procedure :: restart_write => toy_rst_w
        procedure :: restart_read  => toy_rst_r
    end type toy_geo_t

contains

    subroutine toy_init(self, ctx)
        class(toy_geo_t), intent(inout) :: self
        type(init_ctx_t), intent(in)    :: ctx
        associate (unused => ctx%time); end associate
        self%prec = sp
        self%nx = 2;  self%ny = 2
        self%opt = geo_opt_t(A_ocn=.true., dV_bsl=.true.)
    end subroutine toy_init

    subroutine toy_update(self, st)
        class(toy_geo_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (unused => st%t1 + self%nx); end associate
    end subroutine toy_update

    subroutine toy_rst_w(self, dir, t)
        class(toy_geo_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (unused => t + self%nx + len(dir)); end associate
    end subroutine toy_rst_w

    subroutine toy_rst_r(self, dir, t)
        class(toy_geo_t), intent(inout) :: self
        character(len=*), intent(in)    :: dir
        real(dp),         intent(in)    :: t
        associate (unused => t + self%nx + len(dir)); end associate
    end subroutine toy_rst_r

end module test_globe_toy

program test_globe
    ! kryos_globe (core C7, make check): a globe from a par file (tsgen
    ! forcing with a ramp coefficient, a feedback coefficient on M_ice and a
    ! derived one; the fastisostasy global part) driven through two toy
    ! domains the way the run sequence will: M_ice and dV_bsl gathered from
    ! the domains' exports, the global bsl reduced, bsl scattered to the
    ! domains' geo through geo_scatter, coefficients and bsl exchanged into the
    ! domains' imports, a derived coefficient sent to the forcing. Every value
    ! is checked bit for bit against the libraries driven directly (tsgen,
    ! bsl_class), as yelmox drives them. With fastisostasy=0 the test is not
    ! built (test_skip in config/common.mk).
    !
    ! Run from the repo root (input/kryos_defaults.nml); files in tmp/test_globe/.

    use precision,        only: sp, dp, wp
    use fesm_clock,       only: clock_t, clock_init, clock_start, clock_advance, clock_time, clock_time_rel, &
                                CADENCE_STEP
    use fesm_fields,      only: field_state_t, fields_init, fields_add, fields_lock, fields_put, fields_get, &
                                FIELD_INSTANT, FIELD_INTENSIVE, FIELD_EXTENSIVE
    use phys_constants,   only: phys_const_class
    use coupler,          only: coupler_class
    use tsgen,            only: tsgen_class, tsgen_init, tsgen_update
    use fastisostasy,     only: bsl_class, bsl_init, bsl_update
    use kryos_config,     only: config_t, config_init
    use kryos_component,  only: init_ctx_t, step_time_t, census_t, constraints_t, PASS_INIT
    use kryos_wiring,     only: wire_t, wire_node_t, wire_node, wire_node_states, wire_resolve, wire_label
    use kryos_geo,        only: geo_slot_t, geo_adopt, geo_init, geo_scatter
    use test_globe_toy,   only: toy_geo_t
    use kryos_globe,      only: globe_t, gather_src_t, globe_init, globe_init_state, globe_update, globe_reduce, &
                                globe_end, globe_nodes, exchange_globe, exchange_to_globe, &
                                gather_src, gather_into

    implicit none

    integer,  parameter :: ND = 2, NSTEP = 8
    real(dp), parameter :: DTT = 10.0_dp, TIME_REF = -100.0_dp
    character(len=*), parameter :: DIR = "tmp/test_globe/", PAR = DIR//"par.nml", &
                                   DEFAULTS = "input/forcing_defaults.nml"

    type :: toy_dom_t                                ! the instances of a domain the globe talks to
        character(len=16)    :: name = ""
        type(geo_slot_t)     :: geo
        type(field_state_t)  :: ice_i, ice_e, clim_i, clim_e
        integer              :: h_M = 0, hc_alpha = 0, hc_beta = 0, h_dTo = 0
        type(wire_t), allocatable :: w(:)
    end type toy_dom_t

    type(config_t)         :: cfg
    type(clock_t)          :: clk
    type(phys_const_class) :: cnst
    type(census_t)         :: census
    type(globe_t), target  :: globe
    type(toy_dom_t), target :: dom(ND)
    type(gather_src_t)     :: srcs(ND)
    type(coupler_class)    :: cpl
    type(tsgen_class)      :: ts_a, ts_b
    type(bsl_class)        :: bsl
    real(dp) :: rel0, t, tprev
    real(sp) :: x, msum, dvsum, f_dTo
    integer  :: n, d, h
    logical  :: stale

    call execute_command_line("mkdir -p "//DIR)
    call par_write()

    ! ---- globe_init: forcing then geo, cadences on the clock, before clock_start ----
    call config_init(cfg, PAR)
    call clock_init(clk, 0.0_dp, 200.0_dp, DTT, "365_day", "tropical", method="cal", time_ref=TIME_REF)
    call globe_init(globe, cfg, clk, cnst, census)
    call clock_start(clk)
    rel0 = clock_time_rel(clk)
    call check(rel0 == 100.0_dp, "forcing time at t0 = t - time_ref")
    call check(globe%clock_forcing == "forcing" .and. globe%clock_geo == "forcing", "clock keys")
    call check(globe%forcing%c%cadence%mode == CADENCE_STEP .and. globe%geo%c%cadence%mode == CADENCE_STEP, "every-step cadences")
    call check(trim(globe%forcing%c%name) == "forcing" .and. trim(globe%geo%c%name) == "geo_glob" .and. &
               trim(globe%geo%c%scope) == "global", "instances")

    ! ---- the libraries, directly, as yelmox ----
    call bsl_init(bsl, PAR, rel0)
    call bsl_update(bsl, rel0)
    call tsgen_init(ts_a, PAR, rel0, group="coef_alpha", defaults_file=DEFAULTS, defaults_group="coef")
    call tsgen_init(ts_b, PAR, rel0, group="coef_beta",  defaults_file=DEFAULTS, defaults_group="coef")

    ! ---- two toy domains, rows resolved with the globe nodes ----
    call build_domains()
    do d = 1, ND
        call put_M(d, 0, 0.0_dp)
    end do
    call globe_init_state(globe, clk, PASS_INIT)

    ! cold start: the init values reach the climate and the geo imports
    call exchange_globe(globe, "forcing", dom(1)%w, "climate", dom(1)%clim_i, clk, init=.true.)
    call get(dom(1)%clim_i, dom(1)%hc_alpha, x);  call check(x == ts_a%f_now, "init: coef_alpha")
    call get(dom(1)%clim_i, dom(1)%hc_beta,  x);  call check(x == ts_b%f_now, "init: coef_beta")
    do d = 1, ND
        call exchange_globe(globe, "geo", dom(d)%w, "geo", dom(d)%geo%c%imports, clk, init=.true.)
        call get(dom(d)%geo%c%imports, dom(d)%geo%c%hi_bsl, x)
        call check(x == bsl%bsl_now, "init: bsl")
        call get(dom(d)%geo%c%imports, dom(d)%geo%c%hi_A_ocn, x)
        call check(x == bsl%A_ocean_now, "init: A_ocn")
    end do

    ! cold start reduction (core 3.6 item 2): the geo PASS_GEOMETRY changed bsl
    call put_dV(1, 3.0e11_dp, 0.0_dp);  call put_dV(2, -1.0e11_dp, 0.0_dp)
    call build_srcs("geo")
    call gather_into(globe, "geo", srcs, 0.0_dp)
    call globe_reduce(globe, clk)
    bsl%bsl_now = bsl%bsl_now - (real(3.0e11_dp, sp) + real(-1.0e11_dp, sp)) / bsl%A_ocean_now
    do d = 1, ND
        call exchange_globe(globe, "geo", dom(d)%w, "geo", dom(d)%geo%c%imports, clk, init=.true.)
        call geo_scatter(dom(d)%geo, step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=DTT, dtt=DTT))
        call get(dom(d)%geo%c%exports, dom(d)%geo%c%h_bsl, x)
        call check(x == bsl%bsl_now, "cold start: reduced bsl scattered")
    end do
    call check(bsl%bsl_now /= 0.0_sp, "cold start: the reduction moved bsl")

    ! the climate's derived coefficient at t0 reaches the forcing (core 3.6 item 3)
    call put_dTo(0, 0.0_dp)
    call exchange_to_globe(globe, dom(1)%w, "climate", dom(1)%clim_e, 0.0_dp, "forcing", clk, init=.true.)
    f_dTo = real(0.5_dp, sp)

    ! ---- steps ----
    do n = 1, NSTEP
        call clock_advance(clk)
        t     = clock_time(clk)
        tprev = t - DTT

        ! start of the step: lag 1 rows (M_ice of the domains, the derived coefficient), then the globe
        call build_srcs("forcing")
        call gather_into(globe, "forcing", srcs, tprev)
        call exchange_to_globe(globe, dom(1)%w, "climate", dom(1)%clim_e, tprev, "forcing", clk)
        call globe_update(globe, "geo", clk)
        call globe_update(globe, "forcing", clk)

        call bsl_update(bsl, clock_time_rel(clk))
        msum = real(real(M_val(1, n - 1), sp) + real(M_val(2, n - 1), sp), sp)
        call tsgen_update(ts_a, clock_time_rel(clk), var=0.0_wp)
        call tsgen_update(ts_b, clock_time_rel(clk), var=msum)
        f_dTo = real(dTo_val(n - 1), sp)

        ! coefficients and bsl into the domains
        call exchange_globe(globe, "forcing", dom(1)%w, "climate", dom(1)%clim_i, clk)
        call get(dom(1)%clim_i, dom(1)%hc_alpha, x);  call check(x == ts_a%f_now, "step: coef_alpha")
        call get(dom(1)%clim_i, dom(1)%hc_beta,  x);  call check(x == ts_b%f_now, "step: coef_beta")
        call get(globe%forcing%c%exports, globe%forcing%c%h_coef(3), x)
        call check(x == f_dTo, "step: derived coefficient copied from its import")
        do d = 1, ND
            call exchange_globe(globe, "geo", dom(d)%w, "geo", dom(d)%geo%c%imports, clk)
            call get(dom(d)%geo%c%imports, dom(d)%geo%c%hi_bsl, x)
            call check(x == bsl%bsl_now, "step: bsl before the reduction")
            call get(dom(d)%geo%c%imports, dom(d)%geo%c%hi_A_ocn, x)
            call check(x == bsl%A_ocean_now, "step: A_ocn")
        end do

        ! geo of the domains ran: dV_bsl; in step 3 the second domain did not (stale: counts 0)
        stale = n == 3
        call put_dV(1, dV_val(1, n), t)
        if (.not. stale) call put_dV(2, dV_val(2, n), t)
        dvsum = real(dV_val(1, n), sp)
        if (.not. stale) dvsum = dvsum + real(dV_val(2, n), sp)
        call build_srcs("geo")
        call gather_into(globe, "geo", srcs, t)
        call globe_reduce(globe, clk)
        bsl%bsl_now = bsl%bsl_now - dvsum / bsl%A_ocean_now

        do d = 1, ND
            call exchange_globe(globe, "geo", dom(d)%w, "geo", dom(d)%geo%c%imports, clk)
            call geo_scatter(dom(d)%geo, step_time_t(t0=tprev, t1=t, t_eval=t, dt=DTT, dtt=DTT))
            call get(dom(d)%geo%c%exports, dom(d)%geo%c%h_bsl, x)
            call check(x == bsl%bsl_now, "step: reduced bsl scattered")
        end do

        ! end of the step: the domains publish M_ice and the climate its derived coefficient
        do d = 1, ND
            call put_M(d, n, t)
        end do
        call put_dTo(n, t)
        write(*,'(a,i2,a,3es14.6)') "test_globe: step ", n, ": alpha, beta, bsl = ", &
            real(ts_a%f_now, dp), real(ts_b%f_now, dp), real(bsl%bsl_now, dp)
    end do
    call check(ts_b%f_now /= ts_b%par%f_min, "the feedback coefficient moved (M_ice reached the forcing)")

    call globe_end(globe)
    call execute_command_line("rm -rf "//DIR)
    write(*,'(a)') "test_globe: OK"

contains

    ! ---- values the toy domains publish ----
    function M_val(d, n) result(v)
        integer, intent(in) :: d, n
        real(dp) :: v
        v = merge(1000.1_dp + 9.37_dp*real(n, dp), 400.3_dp - 3.21_dp*real(n, dp), d == 1)
    end function M_val

    function dV_val(d, n) result(v)
        integer, intent(in) :: d, n
        real(dp) :: v
        v = merge(1.0e11_dp*real(n, dp) + 7.0e9_dp, -4.0e10_dp*real(n, dp) + 1.3e9_dp, d == 1)
    end function dV_val

    function dTo_val(n) result(v)
        integer, intent(in) :: n
        real(dp) :: v
        v = 0.5_dp + 0.05_dp*real(n, dp)
    end function dTo_val

    subroutine put_M(d, n, t)
        integer,  intent(in) :: d, n
        real(dp), intent(in) :: t
        call fields_put(dom(d)%ice_e, dom(d)%h_M, M_val(d, n), t1=t)
    end subroutine put_M

    subroutine put_dV(d, v, t)
        integer,  intent(in) :: d
        real(dp), intent(in) :: v, t
        call fields_put(dom(d)%geo%c%exports, dom(d)%geo%c%h_dV_bsl, v, t1=t)
    end subroutine put_dV

    subroutine put_dTo(n, t)
        integer,  intent(in) :: n
        real(dp), intent(in) :: t
        call fields_put(dom(1)%clim_e, dom(1)%h_dTo, dTo_val(n), t1=t)
    end subroutine put_dTo

    subroutine get(st, h, x)
        type(field_state_t), intent(in)  :: st
        integer,             intent(in)  :: h
        real(sp),            intent(out) :: x
        call fields_get(st, h, x)
    end subroutine get

    subroutine build_srcs(dst_kind)
        ! Per-domain gather contributions, domain order.
        character(len=*), intent(in) :: dst_kind
        integer :: k
        do k = 1, ND
            if (dst_kind == "forcing") then
                srcs(k) = gather_src(dom(k)%w, wire_node_states("ice", "ice_"//trim(dom(k)%name), dom(k)%ice_i, &
                                     dom(k)%ice_e, sp, DTT, check=.false.), dst_kind)
            else
                srcs(k) = gather_src(dom(k)%w, wire_node(dom(k)%geo%c, DTT, check=.false.), dst_kind)
            end if
        end do
    end subroutine build_srcs

    subroutine build_domains()
        ! The instances of a domain the globe reaches, as toy states, and their
        ! rows from the real table: wire_resolve with the globe nodes.
        type(toy_geo_t)     :: proto
        type(init_ctx_t)    :: ctx
        type(wire_node_t), allocatable :: nodes(:)
        integer :: k

        do k = 1, ND
            dom(k)%name = merge("d1", "d2", k == 1)

            ctx%par_file  = PAR
            ctx%instance  = "geo_"//trim(dom(k)%name)
            ctx%component = "geo"
            ctx%backend   = "toy"
            ctx%domain    = dom(k)%name
            ctx%grid_name = "TOY"
            ctx%scope     = "regional"
            call geo_adopt(dom(k)%geo, proto, constraints_t(), ctx, census)
            call geo_init(dom(k)%geo, ctx)

            call fields_init(dom(k)%ice_i, "ice_"//trim(dom(k)%name)//":imports");  call fields_lock(dom(k)%ice_i)
            call fields_init(dom(k)%ice_e, "ice_"//trim(dom(k)%name)//":exports")
            call fields_add(dom(k)%ice_e, "M_ice", "Gt", FIELD_INSTANT, FIELD_EXTENSIVE, dom(k)%h_M)
            call fields_lock(dom(k)%ice_e)

            call fields_init(dom(k)%clim_i, "climate_"//trim(dom(k)%name)//":imports")
            call fields_init(dom(k)%clim_e, "climate_"//trim(dom(k)%name)//":exports")
            if (k == 1) then
                call fields_add(dom(k)%clim_i, "coef_alpha", "1", FIELD_INSTANT, FIELD_INTENSIVE, dom(k)%hc_alpha)
                call fields_add(dom(k)%clim_i, "coef_beta",  "1", FIELD_INSTANT, FIELD_INTENSIVE, dom(k)%hc_beta)
                call fields_add(dom(k)%clim_e, "dTo", "K", FIELD_INSTANT, FIELD_INTENSIVE, dom(k)%h_dTo)
            end if
            call fields_lock(dom(k)%clim_i);  call fields_lock(dom(k)%clim_e)
        end do

        do k = 1, ND
            nodes = [wire_node_states("ice", "ice_"//trim(dom(k)%name), dom(k)%ice_i, dom(k)%ice_e, sp, DTT, &
                                      check=.false.), &
                     wire_node(dom(k)%geo%c, DTT, check=.false.)]
            if (k == 1) nodes = [nodes, wire_node_states("climate", "climate_d1", dom(k)%clim_i, dom(k)%clim_e, sp, DTT)]
            nodes = [nodes, globe_nodes(globe, clk)]
            call wire_resolve(dom(k)%w, nodes, cpl)
        end do

        ! the rows the globe uses are active, the others are not
        call check(row_active(dom(1)%w, "geo_g.bsl -> geo.bsl") .and. row_active(dom(1)%w, "geo_g.A_ocn -> geo.A_ocn") &
                   .and. row_active(dom(1)%w, "geo.dV_bsl -> geo_g.dV_bsl") .and. &
                   row_active(dom(1)%w, "ice.M_ice -> forcing.M_ice") .and. &
                   row_active(dom(1)%w, "forcing.alpha -> climate.coef_alpha") .and. &
                   row_active(dom(1)%w, "forcing.beta -> climate.coef_beta") .and. &
                   row_active(dom(1)%w, "climate.dTo -> forcing.dTo"), "rows of domain 1 resolve")
        call check(.not. row_active(dom(2)%w, "climate.dTo -> forcing.dTo") .and. &
                   row_active(dom(2)%w, "geo.dV_bsl -> geo_g.dV_bsl"), "rows of domain 2 resolve")
    end subroutine build_domains

    logical function row_active(w, label)
        type(wire_t),     intent(in) :: w(:)
        character(len=*), intent(in) :: label
        integer :: i
        row_active = .false.
        do i = 1, size(w)
            if (wire_label(w(i)) == label) row_active = w(i)%active
        end do
    end function row_active

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            write(*,'(a)') "test_globe: FAILED ("//what//")"
            error stop 1
        end if
    end subroutine check

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&kryos", ' log_level = "warn"', "/", &
                       "&globe", ' forcing = "forcing"', ' geo = "geo_glob"', "/", &
                       "&forcing", ' backend = "tsgen"', ' clock = "forcing"', &
                       " coefs = 'alpha' 'beta' 'dTo'", "/", &
                       "&geo_glob", ' backend = "fastisostasy"', ' clock = "forcing"', ' bsl = "barysealevel"', "/", &
                       "&coef_alpha", ' method = "ramp-slope"', " f_min = 0.1", " f_max = 0.9", &
                       " df_dt_max = 0.003", "/", &
                       "&coef_beta", ' method = "exp"', " f_min = 0.0", " f_max = 1.0", " df_dt_max = 0.01", &
                       " dt_ave = 30.0", " tol = 0.5", "/", &
                       "&coef_dTo", ' method = "derived"', ' from = "air_anomaly"', " factor = 1.0", &
                       ' units = "K"', "/", &
                       "&barysealevel", ' method = "fastiso"', " bsl_init = 0.0", ' sl_path = "None"', &
                       ' sl_name = "none"', " A_ocean_pd = 3.625e14", ' A_ocean_path = "None"', &
                       ' restart = "None"', "/"
        close(u)
    end subroutine par_write

end program test_globe
