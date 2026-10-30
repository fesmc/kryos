module g2b_ice
    ! The mock ice of gate G2b: an ice instance (ice_adopt, as the tests do) that only declares its fields and has Yelmo's
    ! flotation rule (calc_H_grnd and the gl_sep = 1 rule of the science case, in sp, as yelmo_flotation_). Its thickness is
    ! replayed by the harness (g2b_replay); its bed and sea level are what the Earth's scatter gives it.
    use precision,       only: sp, dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_ice,       only: ice_t, ice_opt_t
    use topography,      only: calc_H_grnd, calc_f_grnd_subgrid_linear
    implicit none
    private

    type, extends(ice_t), public :: replay_ice_t
        real(sp) :: rho_ice = 0.0_sp, rho_sw = 0.0_sp
    contains
        procedure :: init => replay_init
        procedure :: update => replay_update
        procedure :: restart_write => replay_rw
        procedure :: restart_read => replay_rr
        procedure :: flotation => replay_flotation
    end type replay_ice_t

contains

    subroutine replay_init(self, ctx)
        class(replay_ice_t), intent(inout) :: self
        type(init_ctx_t),    intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt = ice_opt_t(z_bed_ref=.true., footprint=.true., H_ice_pd=.true.)
        self%caps%has_flotation = .true.
        self%rho_ice = real(ctx%cnst%rho_ice, sp)
        self%rho_sw  = real(ctx%cnst%rho_sw, sp)
    end subroutine replay_init

    subroutine replay_flotation(self, H, z_bed, z_sl, f_grnd)
        ! Yelmo's grounded fraction of ice H over z_bed and z_sl (kryos_ice_yelmo: yelmo_flotation_ with gl_sep = 1 and
        ! boundaries = "zeros"): calc_H_grnd (no f_ice) with the run's densities, then the binary rule of gl_sep = 1.
        class(replay_ice_t), intent(in)  :: self
        real(dp),            intent(in)  :: H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp),            intent(out) :: f_grnd(:,:)
        real(sp), allocatable :: H_s(:,:), f_ice(:,:), H_grnd(:,:), f_g(:,:), f_x(:,:), f_y(:,:)

        allocate(H_s(self%nx, self%ny), f_ice(self%nx, self%ny), H_grnd(self%nx, self%ny), f_g(self%nx, self%ny), &
                 f_x(self%nx, self%ny), f_y(self%nx, self%ny))
        H_s   = real(H, sp)
        f_ice = 1.0_sp
        H_grnd = 0.0_sp
        call calc_H_grnd(H_grnd, H_s, f_ice, real(z_bed, sp), real(z_sl, sp), self%rho_ice, self%rho_sw, use_f_ice=.false.)
        call calc_f_grnd_subgrid_linear(f_g, f_x, f_y, H_grnd, "zeros")
        f_grnd = real(f_g, dp)
    end subroutine replay_flotation

    subroutine replay_update(self, st)
        class(replay_ice_t), intent(inout) :: self
        type(step_time_t),   intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine replay_update

    subroutine replay_rw(self, dir, t)
        class(replay_ice_t), intent(inout) :: self
        character(len=*),    intent(in)    :: dir
        real(dp),            intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine replay_rw

    subroutine replay_rr(self, dir, t)
        class(replay_ice_t), intent(inout) :: self
        character(len=*),    intent(in)    :: dir
        real(dp),            intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine replay_rr

end module g2b_ice

module g2b_mod
    ! Gate G2b (phase4-plan 13.2): a mock ice replaying a PaleoMIST deglaciation on the NH-32KM and ANT-32KM grids, through the
    ! globe of the run (kryos_globe: the global VILMA2 geo, the coupler of the global grids, the global ice field kryos_globe_ice)
    ! and the domain side of the run sequence written out here as kryos_run_gia and bsl_to_geo do it (they are private to the
    ! run): the reference after PASS_INIT, the cold start (gather at t0, VILMA2's spin-up and seed, the anomalies into the ports,
    ! the ice's reconcile), then per step the gather of the previous step's ice (lag 1), VILMA2's update, the bilinear scatter
    ! into the ports and the ice with its new bed and sea level. The standalone vilma.x gets the same history from the global
    ! 0.5 degree file (tools/levante/g2b_prep.jl), lagged by one step, and the same &vilma.
    ! Not part of `make check` (a Levante job: lmax 85 to 170, the 3D viscosity and the spin-up of the science case).
    !
    !   g2b_replay.x <lmax> <nstep> [t0]     in a run folder with input/ (kryos), fesm-utils/, fesmdata/ and maps/ (writable)
    !
    ! Environment: G2B_RATE_TOL (equil_rate_tol, 1e-5), G2B_TIME_MAX (equil_time_max, 2e4), G2B_SPINUP=0 (no spin-up: for the maps
    ! and smoke runs), G2B_NOUT (record the Gauss grid and the ice grids every NOUT steps, 10), G2B_DATA (folder of the standalone
    ! inputs, ../data).
    !
    ! Writes kryos_out.nc (Gauss grid: dz_bed, dz_ss [VILMA2's own: the export plus dphi], rsl, C_ocn at the records; host_w,
    ! area; the series per step: bsl, dphi, bsl_ive, mass_resid and the sea-level equivalent of each domain's gathered ice and of
    ! the external part), g2b_<domain>.nc (ice grid, at the records: H_ice, z_bed, z_sl, f_grnd, dz_bed of the port, lon, lat) and
    ! vilma_g2b.nml (the &ctl and &vilma of the standalone run, written at the start; G2B_NMLONLY=1 writes it and the par file only).
    ! Timings go to stdout ("g2b:" lines).

    use precision,       only: sp, dp
    use fesm_log,        only: LOG_INFO, log_msg, log_set_level, to_str
    use coords,          only: grid_cdo_read_desc
    use coordinates,     only: grid_class
    use coupler,         only: coupler_class, coupler_init, coupler_add_grid, remap
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_read, nc_size, nc_set_lossless
    use phys_constants,  only: phys_const_class, phys_const_load, phys_const_require
    use fesm_clock,      only: clock_t, clock_init, clock_start, clock_advance, clock_time, clock_time_init, clock_time_end
    use fesm_fields,     only: fields_get, fields_put
    use kryos_config,    only: config_t, config_init
    use kryos_component, only: census_t, constraints_t, init_ctx_t, PASS_INIT, PASS_GEOMETRY, PASS_STATE
    use kryos_ice,       only: ice_t, ice_slot_t, ice_adopt, ice_init, ice_flotation
    use kryos_globe,     only: globe_t, globe_init, globe_pass, globe_update, globe_register_grid, globe_end
    use kryos_globe_ice, only: globe_ice_add_domain, globe_ice_reference_begin, globe_ice_reference_add, &
                               globe_ice_reference_end, globe_ice_gather_begin, globe_ice_gather_add, globe_ice_gather_end
    use kryos_geo_vilma2, only: vilma_earth_t
    use g2b_ice,         only: replay_ice_t

    implicit none

    integer, parameter :: ND = 2                        ! domains, in the name order of the run (Antarctica, North)
    character(len=*), parameter :: PAR = "g2b_par.nml"
    character(len=16),  parameter :: DNAME(ND) = [character(len=16) :: "Antarctica", "North"]
    character(len=16),  parameter :: DGRID(ND) = [character(len=16) :: "ANT-32KM", "NH-32KM"]
    character(len=128), parameter :: DFILE(ND) = [character(len=128) :: &
        "fesmdata/Antarctica/ANT-32KM/ANT-32KM_IceHistory-PaleoMIST.nc", "fesmdata/North/NH-32KM/NH-32KM_IceHistory-PaleoMIST.nc"]
    real(dp), parameter :: DTT = 10.0_dp

    type :: dom_t                                       ! one mock domain
        character(len=16) :: name = "", gname = ""
        type(grid_class) :: g
        type(coupler_class) :: cpl                      ! the domain's coupler: its ice grid and the Earth's (the scatter's)
        type(ice_slot_t) :: ice
        integer :: nx = 0, ny = 0, nt = 0
        real(dp), allocatable :: ages(:), hist(:,:,:)   ! the history on the ice grid (years ago), H_ice
        real(dp), allocatable :: fp(:,:), hpd(:,:), zbr(:,:)       ! footprint, present-day ice, reference bed (static)
        real(dp), allocatable :: H(:,:), zb(:,:), zs(:,:), fg(:,:) ! the ice's state now
        real(dp), allocatable :: dzb(:,:)               ! the port's dz_bed
        real(sp), allocatable :: o_H(:,:,:), o_zb(:,:,:), o_zs(:,:,:), o_fg(:,:,:), o_dzb(:,:,:)
    end type dom_t

    type(config_t),         save :: cfg
    type(phys_const_class), save :: cnst
    type(clock_t),          save :: clk
    type(census_t),         save :: census
    type(globe_t),          save :: globe
    type(dom_t),            save :: dom(ND)
    integer, save :: lmax = 170, nstep = 2100, nout = 10, nrec = 0, nphi = 0, nlat = 0
    real(dp), save :: t0 = -21000.0_dp
    logical, save :: spinup = .true.
    character(len=16), save :: rate_tol = "1.0e-5", time_max = "2.0e4"
    character(len=256), save :: datadir = "../data"
    real(dp), allocatable, save :: lon(:), lat(:), area(:,:), tt(:), ser(:,:), esl(:,:)
    real(sp), allocatable, save :: o_bed(:,:,:), o_ss(:,:,:), o_rsl(:,:,:), o_C(:,:,:)
    real(dp), allocatable, save :: tr(:)
    integer(8), save :: c0, c1, rate
    real(dp), save :: tm_gather = 0.0_dp, tm_geo = 0.0_dp, tm_scatter = 0.0_dp, tm_ice = 0.0_dp, tm_out = 0.0_dp

contains

    subroutine run()
        integer :: n, d
        real(dp) :: tgeo
        character(len=64) :: arg
        character(len=8)  :: v

        call get_command_argument(1, arg); read(arg, *) lmax
        call get_command_argument(2, arg); read(arg, *) nstep
        if (command_argument_count() > 2) then
            call get_command_argument(3, arg); read(arg, *) t0
        end if
        call get_environment_variable("G2B_RATE_TOL", v); if (len_trim(v) > 0) rate_tol = v
        call get_environment_variable("G2B_TIME_MAX", v); if (len_trim(v) > 0) time_max = v
        call get_environment_variable("G2B_SPINUP", v);   if (len_trim(v) > 0) spinup = trim(v) /= "0"
        call get_environment_variable("G2B_NOUT", v);     if (len_trim(v) > 0) read(v, *) nout
        call get_environment_variable("G2B_DATA", datadir)
        if (len_trim(datadir) == 0) datadir = "../data"
        nphi = 4*lmax; nlat = 2*lmax + 2
        nrec = nstep/nout + 1
        if (mod(nstep, nout) /= 0) nrec = nrec + 1                 ! the last step is recorded too
        call system_clock(count_rate=rate)

        call par_write()
        call config_init(cfg, PAR)
        call log_set_level(cfg%log_level)
        call nc_set_lossless(.true.)
        call phys_const_load(cnst, trim(cfg%phys_const_file), group=trim(cfg%phys_const))
        call phys_const_require(cnst, "g2b_replay")
        call write_standalone_inputs()
        call get_environment_variable("G2B_NMLONLY", v)
        if (len_trim(v) > 0) then
            write(*,'(a)') "g2b: G2B_NMLONLY: wrote "//PAR//" and vilma_g2b.nml"
            return
        end if
        call clock_init(clk, t0, 0.0_dp, DTT, "365_day", "tropical", method="rel", time_ref=1950.0_dp, ticks_per_year=36000)

        ! the globe: VILMA2 (geo_glob), the coupler of the global grids, the global ice field
        call system_clock(c0)
        call globe_init(globe, cfg, clk, cnst, census)
        do d = 1, ND
            call dom_setup(dom(d), DNAME(d), DGRID(d), DFILE(d))
        end do
        call allocate_records()
        call system_clock(c1); write(*,'(a,f9.2,a)') "g2b: init (grids, VILMA2 setup incl. 3D viscosity, mock ice) ", real(c1-c0, dp)/rate, " s"

        ! PASS_INIT, then the cold start (kryos_run_sequence: gia_reference, PASS_GEOMETRY, the t0 gather, gia_seed)
        call system_clock(c0)
        call globe_pass(globe, "geo", clk, PASS_INIT)
        do d = 1, ND
            call dom_publish(dom(d), clock_time(clk), reference=.true.)
        end do
        call gia_reference()
        call globe_pass(globe, "geo", clk, PASS_GEOMETRY)
        call system_clock(c1); write(*,'(a,f9.2,a)') "g2b: reference ", real(c1-c0, dp)/rate, " s"
        call system_clock(c0)
        call gia_gather(clock_time(clk))
        if (.not. spinup) then
            call log_msg(LOG_INFO, "g2b_replay", "G2B_SPINUP = 0: the seed without VILMA2's spin-up")
        end if
        call globe_pass(globe, "geo", clk, PASS_STATE)
        do d = 1, ND
            call scatter(dom(d))
            call dom_publish(dom(d), clock_time(clk), reference=.false.)     ! the ice on the Earth's solution (reconcile)
        end do
        call system_clock(c1); write(*,'(a,f9.2,a)') "g2b: t0 gather, spin-up and seed ", real(c1-c0, dp)/rate, " s"
        call clock_start(clk)
        call coords_of_earth()
        tt(0) = clock_time(clk)
        call record(0)

        do n = 1, nstep
            call clock_advance(clk)
            tt(n) = clock_time(clk)
            call system_clock(c0)
            call gia_gather(tt(n) - DTT)                              ! the ice of the previous step (lag 1)
            call system_clock(c1); tm_gather = tm_gather + real(c1-c0, dp)/rate
            call system_clock(c0)
            call globe_update(globe, "geo", clk)                      ! VILMA2 over [t - dt, t] under that load
            call system_clock(c1); tm_geo = tm_geo + real(c1-c0, dp)/rate
            tgeo = real(c1-c0, dp)/rate
            do d = 1, ND
                call system_clock(c0)
                call scatter(dom(d))
                call system_clock(c1); tm_scatter = tm_scatter + real(c1-c0, dp)/rate
                call system_clock(c0)
                call dom_publish(dom(d), tt(n), reference=.false.)
                call system_clock(c1); tm_ice = tm_ice + real(c1-c0, dp)/rate
            end do
            call record(n)
            if (mod(n, 10) == 0 .or. n == 1) then
                write(*,'(a,i0,a,f9.1,a,f8.3,a,f9.3,a)') "g2b: step ", n, " t = ", tt(n), " yr, geo ", tgeo, &
                    " s, bsl = ", ser(2, n), " m"
            end if
        end do
        call system_clock(c0)
        call write_outputs()
        call system_clock(c1); tm_out = real(c1-c0, dp)/rate
        write(*,'(a,5f10.2)') "g2b: wall per phase [s]: gather, geo (update), scatter, ice (publish), output ", &
            tm_gather, tm_geo, tm_scatter, tm_ice, tm_out
        write(*,'(a,f9.3,a)') "g2b: mean per step: gather ", tm_gather/max(nstep, 1), " s"
        write(*,'(a,f9.3,a)') "g2b: mean per step: geo update ", tm_geo/max(nstep, 1), " s"
        call globe_end(globe)
    end subroutine run

    ! =========================================================================
    ! The mock domains
    ! =========================================================================

    subroutine dom_setup(dm, name, gname, file)
        type(dom_t),      intent(inout) :: dm
        character(len=*), intent(in)    :: name, gname, file

        class(ice_t), allocatable :: proto
        type(init_ctx_t) :: ctx
        real(dp), allocatable :: z2(:,:)

        dm%name = name; dm%gname = gname
        call grid_cdo_read_desc(dm%g, trim(gname), trim(cfg%map_fldr))
        dm%nx = dm%g%G%nx; dm%ny = dm%g%G%ny
        dm%nt = nc_size(file, "time")
        allocate(dm%ages(dm%nt), dm%hist(dm%nx, dm%ny, dm%nt), z2(dm%nx, dm%ny))
        call nc_read(file, "time", dm%ages)
        call nc_read(file, "H_ice", dm%hist)
        call nc_read(file, "z_bed", z2, start=[1,1,dm%nt], count=[dm%nx, dm%ny, 1])      ! the record at 0 ka: the reference bed
        if (dm%ages(dm%nt) /= 0.0_dp .or. dm%ages(1) <= dm%ages(dm%nt)) then
            write(*,'(a)') "g2b: the history of "//trim(name)//" does not end at 0 years ago"
            error stop 1
        end if
        allocate(dm%fp(dm%nx, dm%ny), dm%hpd(dm%nx, dm%ny), dm%zbr(dm%nx, dm%ny))
        dm%fp  = merge(1.0_dp, 0.0_dp, maxval(dm%hist, dim=3) > 0.0_dp)   ! where the ice model can hold ice: the history's own extent
        dm%hpd = dm%hist(:,:,dm%nt)
        dm%zbr = z2
        allocate(dm%H(dm%nx, dm%ny), dm%zb(dm%nx, dm%ny), dm%zs(dm%nx, dm%ny), dm%fg(dm%nx, dm%ny), dm%dzb(dm%nx, dm%ny))
        dm%zb = dm%zbr; dm%zs = 0.0_dp; dm%dzb = 0.0_dp

        ctx%instance  = "ice_"//trim(name)
        ctx%component = "ice"
        ctx%backend   = "replay"
        ctx%domain    = name
        ctx%grid_name = gname
        ctx%grid      = dm%g
        ctx%cnst      = cnst
        allocate(replay_ice_t :: proto)
        call ice_adopt(dm%ice, proto, constraints_t(), ctx, census)
        call ice_init(dm%ice, ctx)

        call coupler_init(dm%cpl, trim(cfg%map_fldr))
        call coupler_add_grid(dm%cpl, trim(gname), dm%g)
        call globe_register_grid(globe, dm%cpl)                           ! the Earth's grid, for the scatter
        write(*,'(a,a,a,i0,a,i0,a,i0,a)') "g2b: domain ", trim(name), " on ", dm%nx, " x ", dm%ny, ", ", dm%nt, " records"
        write(*,'(a,i0,a,i0,a)') "g2b:   footprint ", nint(sum(dm%fp)), " cells, ice at 0 ka ", count(dm%hpd > 0.0_dp), " cells"
    end subroutine dom_setup

    subroutine hist_at(dm, t, H)
        ! The history at model time t [yr] (age = -t): linear between the records, held outside.
        type(dom_t), intent(in)  :: dm
        real(dp),    intent(in)  :: t
        real(dp),    intent(out) :: H(:,:)
        real(dp) :: age, w
        integer :: k

        age = -t
        if (age >= dm%ages(1)) then
            H = dm%hist(:,:,1)
        else if (age <= dm%ages(dm%nt)) then
            H = dm%hist(:,:,dm%nt)
        else
            do k = 1, dm%nt - 1
                if (age <= dm%ages(k) .and. age > dm%ages(k+1)) exit
            end do
            w = (dm%ages(k) - age)/(dm%ages(k) - dm%ages(k+1))
            H = (1.0_dp - w)*dm%hist(:,:,k) + w*dm%hist(:,:,k+1)
        end if
    end subroutine hist_at

    subroutine dom_publish(dm, t, reference)
        ! The ice's exports at t (what its update and exports_geometry put): the replayed thickness, the bed and sea level it
        ! holds, and the grounded fraction of its own rule. reference: the first publication (PASS_INIT, the reference bed and
        ! sea level 0), with the static exports (footprint, H_ice_pd, z_bed_ref).
        type(dom_t), intent(inout) :: dm
        real(dp),    intent(in)    :: t
        logical,     intent(in)    :: reference

        call hist_at(dm, t, dm%H)
        call ice_flotation(dm%ice, dm%H, dm%zb, dm%zs, dm%fg)
        associate (c => dm%ice%c)
            if (reference) then
                call fields_put(c%exports, c%h_footprint, dm%fp,  t1=t)
                call fields_put(c%exports, c%h_H_ice_pd,  dm%hpd, t1=t)
                call fields_put(c%exports, c%h_z_bed_ref, dm%zbr, t1=t)
            end if
            call fields_put(c%exports, c%h_H_ice,  dm%H,  t1=t)
            call fields_put(c%exports, c%h_z_bed,  dm%zb, t1=t)
            call fields_put(c%exports, c%h_z_sl,   dm%zs, t1=t)
            call fields_put(c%exports, c%h_f_grnd, dm%fg, t1=t)
        end associate
    end subroutine dom_publish

    subroutine scatter(dm)
        ! The Earth's anomalies into the port (rows geo_g -> geo: dz_bed and dz_ss bilinear onto the ice grid, dphi), then into
        ! the ice (rows geo -> ice: z_bed = z_bed_ref + dz_bed, z_sl = dphi + dz_ss).
        type(dom_t), intent(inout) :: dm
        real(dp), allocatable :: g(:,:), r(:,:)
        real(dp) :: dphi

        associate (c => globe%geo%c)
            allocate(g(c%nx, c%ny))
            call fields_get(c%exports, c%h_dz_bed, g)
            call remap(dm%cpl, g, trim(c%grid), r, trim(dm%gname), method="bilin")
            dm%dzb = r
            dm%zb  = dm%zbr + r
            call fields_get(c%exports, c%h_dz_ss, g)
            call remap(dm%cpl, g, trim(c%grid), r, trim(dm%gname), method="bilin")
            call fields_get(c%exports, c%h_dphi, dphi)
            dm%zs = r + dphi
        end associate
    end subroutine scatter

    ! =========================================================================
    ! The run sequence's global ice field (kryos_run_gia)
    ! =========================================================================

    subroutine gia_reference()
        ! After PASS_INIT: the domains' static exports, the background X and the external history into the reference of the
        ! Earth (kryos_run_gia gia_reference).
        real(dp), allocatable :: fp(:,:), hpd(:,:), zbr(:,:)
        integer :: d

        associate (gi => globe%gi, cpl => globe%cpl)
            do d = 1, ND
                associate (ice => dom(d)%ice)
                    allocate(fp(ice%c%nx, ice%c%ny), hpd(ice%c%nx, ice%c%ny), zbr(ice%c%nx, ice%c%ny))
                    call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                    call fields_get(ice%c%exports, ice%c%h_H_ice_pd,  hpd)
                    call fields_get(ice%c%exports, ice%c%h_z_bed_ref, zbr)
                    call globe_ice_add_domain(gi, cpl, trim(dom(d)%name), dom(d)%g, fp, hpd, zbr)
                    deallocate(fp, hpd, zbr)
                end associate
            end do
            call globe_ice_reference_begin(gi, cpl)
            do d = 1, ND
                call globe_ice_reference_add(gi, cpl, trim(dom(d)%name), dom(d)%ice)
            end do
            call globe_ice_reference_end(gi, cpl)
            associate (c => globe%geo%c)
                call fields_put(c%imports, c%hi_z_bed_ref,   gi%z_bed_ref,   t1=clock_time(clk))
                call fields_put(c%imports, c%hi_H_ice_ref,   gi%H_ice_ref,   t1=clock_time(clk))
                call fields_put(c%imports, c%hi_host_w,      gi%host_w,      t1=clock_time(clk))
                call fields_put(c%imports, c%hi_host_Ig_ref, gi%host_Ig_ref, t1=clock_time(clk))
                call fields_put(c%imports, c%hi_host_C_ref,  gi%host_C_ref,  t1=clock_time(clk))
            end associate
        end associate
    end subroutine gia_reference

    subroutine gia_gather(t_load)
        ! The host terms of the ice published at t_load (kryos_run_gia gia_gather) into the Earth's imports.
        real(dp), intent(in) :: t_load
        real(dp), allocatable :: H(:,:), fg(:,:), zb(:,:), zs(:,:), fp(:,:)
        real(dp) :: tp
        integer :: d

        associate (gi => globe%gi, cpl => globe%cpl)
            call globe_ice_gather_begin(gi, cpl, t_load)
            do d = 1, ND
                associate (ice => dom(d)%ice)
                    allocate(H(ice%c%nx, ice%c%ny), fg(ice%c%nx, ice%c%ny), zb(ice%c%nx, ice%c%ny), &
                             zs(ice%c%nx, ice%c%ny), fp(ice%c%nx, ice%c%ny))
                    call fields_get(ice%c%exports, ice%c%h_H_ice, H, t1=tp)
                    if (tp /= t_load) then
                        write(*,'(a,f10.2,a,f10.2)') "g2b: the ice was published at ", tp, ", the load is read at ", t_load
                        error stop 1
                    end if
                    call fields_get(ice%c%exports, ice%c%h_f_grnd,    fg)
                    call fields_get(ice%c%exports, ice%c%h_z_bed,     zb)
                    call fields_get(ice%c%exports, ice%c%h_z_sl,      zs)
                    call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                    call globe_ice_gather_add(gi, cpl, trim(dom(d)%name), ice, H, fg, zb, zs, fp)
                    deallocate(H, fg, zb, zs, fp)
                end associate
            end do
            call globe_ice_gather_end(gi, cpl)
            associate (c => globe%geo%c)
                call fields_put(c%imports, c%hi_H_ice_ext, gi%H_ice_ext, t1=t_load)
                call fields_put(c%imports, c%hi_host_w,    gi%host_w,    t1=t_load)
                call fields_put(c%imports, c%hi_host_Ig,   gi%host_Ig,   t1=t_load)
                call fields_put(c%imports, c%hi_host_C,    gi%host_C,    t1=t_load)
                call fields_put(c%imports, c%hi_host_F,    gi%host_F,    t1=t_load)
            end associate
        end associate
    end subroutine gia_gather

    ! =========================================================================
    ! Records and output
    ! =========================================================================

    subroutine allocate_records()
        integer :: d
        allocate(tt(0:nstep), ser(7, 0:nstep), esl(3, 0:nstep))
        allocate(lon(nphi), lat(nlat), area(nphi, nlat))
        allocate(o_bed(nphi, nlat, nrec), o_ss(nphi, nlat, nrec), o_rsl(nphi, nlat, nrec), o_C(nphi, nlat, nrec), tr(nrec))
        ser = 0.0_dp; esl = 0.0_dp
        do d = 1, ND
            allocate(dom(d)%o_H(dom(d)%nx, dom(d)%ny, nrec), dom(d)%o_zb(dom(d)%nx, dom(d)%ny, nrec), &
                     dom(d)%o_zs(dom(d)%nx, dom(d)%ny, nrec), dom(d)%o_fg(dom(d)%nx, dom(d)%ny, nrec), &
                     dom(d)%o_dzb(dom(d)%nx, dom(d)%ny, nrec))
        end do
    end subroutine allocate_records

    subroutine coords_of_earth()
        integer :: j
        associate (gg => globe%geo%c%grid_obj)
            if (gg%G%nx /= nphi .or. gg%G%ny /= nlat) then
                write(*,'(a)') "g2b: the Earth's grid is not the Gauss grid of this lmax"
                error stop 1
            end if
            lon = gg%G%x
            lat = gg%G%y
            do j = 1, nlat
                area(:,j) = gg%area(:,j)*gg%cs%xy_conv**2
            end do
        end associate
    end subroutine coords_of_earth

    subroutine record(k)
        integer, intent(in) :: k
        real(dp), allocatable :: a(:,:)
        real(dp) :: s, dphi
        integer :: ir, d, i

        allocate(a(nphi, nlat))
        associate (c => globe%geo%c)
            call fields_get(c%exports, c%h_dphi, dphi);       ser(1, k) = dphi
            call fields_get(c%exports, c%h_bsl, s);           ser(2, k) = s
            call fields_get(c%exports, c%h_bsl_ive, s);       ser(3, k) = s
            call fields_get(c%exports, c%h_mass_resid, s);    ser(4, k) = s
            call fields_get(c%exports, c%h_A_ocn, s);         ser(5, k) = s
            ser(6, k) = globe%gi%diag%clip_w
            ser(7, k) = globe%gi%diag%cons_I
            if (globe%gi%diag%nd >= 1) then
                do i = 1, min(globe%gi%diag%nd, 2)
                    esl(i, k) = globe%gi%diag%esl(i)
                end do
                esl(3, k) = globe%gi%diag%ext_esl
            end if
            if (mod(k, nout) /= 0 .and. k /= nstep .and. k /= 0) return
            ir = k/nout + 1
            if (k == nstep .and. mod(k, nout) /= 0) ir = nrec          ! nstep not a multiple of nout: the last record
            tr(ir) = tt(k)
            call fields_get(c%exports, c%h_dz_bed, a);  o_bed(:,:,ir) = real(a, sp)
            call fields_get(c%exports, c%h_dz_ss,  a);  o_ss(:,:,ir)  = real(a + dphi, sp)       ! VILMA2's own dz_ss
            call fields_get(c%exports, c%h_C_ocn,  a);  o_C(:,:,ir)   = real(a, sp)
            select type (e => c)
            type is (vilma_earth_t)
                o_rsl(:,:,ir) = real(e%se%rsl, sp)
            end select
        end associate
        do d = 1, ND
            dom(d)%o_H(:,:,ir)   = real(dom(d)%H, sp)
            dom(d)%o_zb(:,:,ir)  = real(dom(d)%zb, sp)
            dom(d)%o_zs(:,:,ir)  = real(dom(d)%zs, sp)
            dom(d)%o_fg(:,:,ir)  = real(dom(d)%fg, sp)
            dom(d)%o_dzb(:,:,ir) = real(dom(d)%dzb, sp)
        end do
    end subroutine record

    subroutine write_outputs()
        character(len=*), parameter :: f = "kryos_out.nc"
        character(len=64) :: fd
        real(dp), allocatable :: x(:), y(:)
        integer :: d, nr

        nr = nrec
        call nc_create(f)
        call nc_write_dim(f, "lon", x=lon, units="degrees_east")
        call nc_write_dim(f, "lat", x=lat, units="degrees_north")
        call nc_write_dim(f, "rec", x=tr(1:nr), units="years", unlimited=.false.)
        call nc_write_dim(f, "step", x=tt, units="years", unlimited=.false.)
        call nc_write(f, "dz_bed", o_bed(:,:,1:nr), dim1="lon", dim2="lat", dim3="rec", units="m")
        call nc_write(f, "dz_ss",  o_ss(:,:,1:nr),  dim1="lon", dim2="lat", dim3="rec", units="m")
        call nc_write(f, "rsl",    o_rsl(:,:,1:nr), dim1="lon", dim2="lat", dim3="rec", units="m")
        call nc_write(f, "C_ocn",  o_C(:,:,1:nr),   dim1="lon", dim2="lat", dim3="rec", units="1")
        call nc_write(f, "host_w", globe%gi%host_w, dim1="lon", dim2="lat", units="1")
        call nc_write(f, "area",   area, dim1="lon", dim2="lat", units="m2")
        call nc_write(f, "dphi",       ser(1,:), dim1="step", units="m")
        call nc_write(f, "bsl",        ser(2,:), dim1="step", units="m")
        call nc_write(f, "bsl_ive",    ser(3,:), dim1="step", units="m")
        call nc_write(f, "mass_resid", ser(4,:), dim1="step", units="1")
        call nc_write(f, "A_ocn",      ser(5,:), dim1="step", units="m2")
        call nc_write(f, "clip_w",     ser(6,:), dim1="step", units="1")
        call nc_write(f, "cons_I",     ser(7,:), dim1="step", units="1")
        call nc_write(f, "esl_domain1", esl(1,:), dim1="step", units="m")
        call nc_write(f, "esl_domain2", esl(2,:), dim1="step", units="m")
        call nc_write(f, "esl_ext",     esl(3,:), dim1="step", units="m")
        do d = 1, ND
            fd = "g2b_"//trim(dom(d)%name)//".nc"
            call nc_create(fd)
            allocate(x(dom(d)%nx), y(dom(d)%ny))
            x = dom(d)%g%G%x; y = dom(d)%g%G%y
            call nc_write_dim(fd, "x", x=x, units="km")
            call nc_write_dim(fd, "y", x=y, units="km")
            call nc_write_dim(fd, "rec", x=tr(1:nr), units="years", unlimited=.false.)
            call nc_write(fd, "lon", dom(d)%g%lon, dim1="x", dim2="y", units="degrees_east")
            call nc_write(fd, "lat", dom(d)%g%lat, dim1="x", dim2="y", units="degrees_north")
            call nc_write(fd, "footprint", dom(d)%fp, dim1="x", dim2="y", units="1")
            call nc_write(fd, "z_bed_ref", dom(d)%zbr, dim1="x", dim2="y", units="m")
            call nc_write(fd, "H_ice",   dom(d)%o_H(:,:,1:nr),   dim1="x", dim2="y", dim3="rec", units="m")
            call nc_write(fd, "z_bed",   dom(d)%o_zb(:,:,1:nr),  dim1="x", dim2="y", dim3="rec", units="m")
            call nc_write(fd, "z_sl",    dom(d)%o_zs(:,:,1:nr),  dim1="x", dim2="y", dim3="rec", units="m")
            call nc_write(fd, "f_grnd",  dom(d)%o_fg(:,:,1:nr),  dim1="x", dim2="y", dim3="rec", units="1")
            call nc_write(fd, "dz_bed",  dom(d)%o_dzb(:,:,1:nr), dim1="x", dim2="y", dim3="rec", units="m")
            deallocate(x, y)
        end do
        write(*,'(a)') "g2b: domain order of the series esl_domain1/2: "//trim(globe%gi%diag%name(1))//", "// &
            trim(globe%gi%diag%name(min(2, globe%gi%diag%nd)))
    end subroutine write_outputs

    ! =========================================================================
    ! The parameters, the standalone run's
    ! =========================================================================

    subroutine vilma_group(u)
        ! &vilma of the science case (par/kryos_NH-ANT_vilma2.nml) at this lmax, with the spin-up of G2c's choice; the same
        ! lines go to the kryos par file and to the standalone's.
        integer, intent(in) :: u
        write(u,'(a)') "&vilma", " lmax = "//trim(itoa(lmax)), ' earth = "PREM"', " n_layer = 11", &
            " r_bot = 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0", &
            " r_top = 6371.0e3, 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3", &
            " eta = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 0.0", &
            " rheology = 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2", " l_visc_3d = .true.", &
            ' visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"', " visc_3d_structure = 6", &
            " pre_spinup_1d = "//merge(".true. ", ".false.", spinup)
        if (spinup) write(u,'(a)') " equil_rate_tol = "//trim(rate_tol), " equil_time_max = "//trim(time_max)
    end subroutine vilma_group

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: b
        write(b, '(i0)') i
        s = trim(b)
    end function itoa

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&kryos", ' globe = "globe"', ' clock = "clock"', ' map_fldr = "maps"', ' phys_const = "Earth"', &
            ' phys_const_file = "input/phys_const.nml"', ' log_level = "info"', "/"
        write(u,'(a)') "&clock", " time_init = "//ftoa(t0), " time_end = 0.0", " dtt = 10.0", ' method = "rel"', &
            " time_ref = 1950.0", ' calendar = "365_day"', ' year = "tropical"', " ticks_per_year = 36000", "/"
        write(u,'(a)') "&globe", ' forcing = ""', ' geo = "geo_glob"', "/"
        write(u,'(a)') "&geo_glob", ' backend = "vilma2"', ' grid = ""', " dt = 10.0", " phase = 0.0", ' clock = "model"', &
            ' out_2D = "none"', ' out_1D = "none"', ' out_vars = ""', ' vilma = "vilma"', ' ice_bg = "globe_ice"', &
            " allow_no_spinup = "//merge(".false.", ".true. ", spinup), ' out_gia = "none"', ' out_vilma = "none"', "/"
        call vilma_group(u)
        write(u,'(a)') "/"
        ! the global ice field of the science case: X = PaleoMIST at 0 ka on GLOBAL-0.5DEG, the history as an anomaly
        write(u,'(a)') "&globe_ice", ' grid = "GLOBAL-0.5DEG"', ' z_bed_ref = "X_z_bed"', ' H_ice_ref = "X_H_ice"', &
            ' ext = "file"', ' ext_field = "X_dH_ice"', " ext_tref = 0.0", ' out_2D = "none"', "/"
        write(u,'(a)') "&X_z_bed", ' filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"', &
            ' name = "z_bed"', ' mode = "transient"', "/"
        write(u,'(a)') "&X_H_ice", ' filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"', &
            ' name = "H_ice"', ' mode = "transient"', "/"
        write(u,'(a)') "&X_dH_ice", ' filename = "fesmdata/Global/{grid_name}/{grid_name}_IceHistory-PaleoMIST.nc"', &
            ' name = "H_ice"', ' mode = "transient"', ' anomaly_ref = "slice"', " ref_time = 0.0", &
            ' out_of_range = "error"', "/"
        close(u)
    end subroutine par_write

    function ftoa(x) result(s)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: s
        character(len=32) :: b
        write(b, '(es24.16)') x
        s = trim(adjustl(b))
    end function ftoa

    subroutine write_standalone_inputs()
        ! The standalone run: the global 0.5 degree history at the steps of this run, lagged by one (g2b_prep.jl), the present-day
        ! reference of the same product, i_eq = 1 (reference = present day), the same &vilma and the run's densities.
        integer :: u
        open(newunit=u, file="vilma_g2b.nml", status="replace", action="write")
        write(u,'(a)') "&ctl", ' file_forcing = "'//trim(datadir)//'/g2b_force_'//trim(itoa(nstep))//'.nc"', ' name_ice = "h_ice"', &
            ' name_time = "time"', ' file_ref = ""', ' name_zbed_eq = "z_bed_eq"', ' name_hice_ref = "h_ice_ref"', &
            ' file_out = "vilma_out.nc"', ' file_hor = ""', " write_lossless = .true.", " time_init = -1.0e30", &
            " time_end = 1.0e30", " remap_input = .true.", ' name_lon = "lon"', ' name_lat = "lat"', " i_eq = 1", &
            ' z_bed_ref_file = "'//trim(datadir)//'/g2b_ref.nc"', ' h_ice_ref_file = "'//trim(datadir)//'/g2b_ref.nc"', &
            ' z_bed_eq_file = ""', ' h_ice_eq_file = ""', ' rsl_restart_file = ""', &
            ' name_z_bed_ref = "z_bed"', ' name_h_ice_ref = "H_ice"', ' name_rsl = "rsl"', ' restart_in_file = ""', "/"
        call vilma_group(u)
        write(u,'(a,es24.16)') " rho_ice = ", cnst%rho_ice
        write(u,'(a,es24.16)') " rho_water = ", cnst%rho_sw
        write(u,'(a)') "/"
        close(u)
    end subroutine write_standalone_inputs

end module g2b_mod

program g2b_replay
    use g2b_mod, only: run
    implicit none
    call run()
end program g2b_replay
