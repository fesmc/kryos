program test_gia_sequence
    ! Phase 4 WP4: VILMA2 in the run sequence (phase4-plan 4-8, 13.1), Yelmo on the local test grids (GRL-16KM,
    ! ANT-32KM from ice_data), VILMA2 at lmax 16 on GAUSS-L16, the background X and the external history synthetic
    ! (a 4 degree lon-lat grid, written here). The runs are children of this executable (one thread, in
    ! tmp/test_gia_sequence/<run>/); the inspecting child drives run_init / run_step itself and checks the run
    ! from inside:
    !
    !   cold start  the two-pass cold start gives a consistent reference: the t0 gather is the reference (host terms
    !               equal, no F), VILMA2's seed is at rest (no bed change, no sea-level shift), the ice's z_bed and
    !               z_sl are the reference's; the global instance has its grid (GAUSS-L16) and the ports step on its
    !               cadence;
    !   gather      at every step the host terms VILMA2 got are those of a globe_ice_t built here from the same files
    !               and fed with the ice's exports of the previous step (bitwise), and the load time is that step;
    !   scatter     the port holds the Gauss-grid anomalies remapped bilinearly (bitwise, an independent remap),
    !               the ice's z_sl is dphi + dz_ss and its z_bed the reference bed + dz_bed on the ice grid, and
    !               the anomalies move once the load changes (an external history that thins);
    !   cadence     the Earth stepping every second step (geo dt = 20 yr, dtt = 10 yr): the gather only when it is due, the
    !               ports and the ice hold the last solution in between;
    !   spin-up     with VILMA2's own spin-up (pre_spinup_1d) the seed holds the whole t0 load, the external ice included;
    !   order       two domains (Greenland, Antarctica): swapping them in [kryos] domains gives identical files
    !               (KryosTools.Gate, `identical`).
    !
    ! Inputs: the hub topographies and basins of GRL-16KM and ANT-32KM from ice_data ($KRYOS_ICE_DATA, default
    ! ./ice_data); without them the test reports and passes. Julia with the project of tools/ ($KRYOS_TOOLS) is
    ! needed for the merge of the two single-domain par files and the comparison; without it the test reports OPEN
    ! and runs the Greenland run only. Built only with vilma=1 yelmo=1 fastisostasy=1 chion=1 (config/common.mk).
    !     test_gia_sequence.x inspect <par file> <steps> <dt_geo> <spin|nospin>   the run, checked from inside
    !     test_gia_sequence.x plain <par file>               the run as the driver does it

    use precision,       only: sp, dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_write_attr, nc_read
    use coords,          only: grid_class, grid_init, grid_cdo_read_desc, grid_cdo_write_desc_short
    use coupler,         only: coupler_class, coupler_init, coupler_add_grid, remap
    use fesm_fields,     only: fields_get
    use fesm_clock,      only: clock_time, clock_to_years, clock_time_init, clock_time_end, clock_dtt, cadence_dt_years
    use kryos_run,       only: run_t, run_init, run_step, run_done, run_end
    use kryos_component, only: NAME_LEN
    use kryos_domain,    only: C_ICE, C_GEO
    use kryos_globe_ice, only: globe_ice_t, globe_ice_init, globe_ice_add_domain, globe_ice_reference_begin, &
                               globe_ice_reference_add, globe_ice_reference_end, globe_ice_gather_begin, &
                               globe_ice_gather_add, globe_ice_gather_end, globe_ice_end

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_gia_sequence/"
    character(len=*), parameter :: BGN = "TBG-4DEG"
    integer, parameter :: NZ = 5, NDOM = 2, NXB = 90, NYB = 45
    logical  :: allow_nospin = .true.              ! allow_no_spinup in the par file being written (the runs without VILMA2's spin-up)
    real(dp) :: restart_dt = 0.0_dp                ! > 0: the restart schedule of the par file being written, every restart_dt yr
    real(dp) :: dt_geo = 10.0_dp                    ! [yr] cadence of the Earth in the par file being written
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    ! the two domains: name, instance suffix, hub grid, topography file
    character(len=16), parameter :: DNAME(NDOM) = [character(len=16) :: "Greenland", "Antarctica"]
    character(len=4),  parameter :: DSFX(NDOM)  = [character(len=4)  :: "grl", "ant"]
    character(len=16), parameter :: DGRID(NDOM) = [character(len=16) :: "GRL-16KM", "ANT-32KM"]
    character(len=32), parameter :: DTOPO(NDOM) = [character(len=32) :: "TOPO-M17-v5", "TOPO-BedMachine"]
    integer, parameter :: DNX(NDOM) = [106, 191], DNY(NDOM) = [181, 191]
    character(len=*), parameter :: PCT = 'skip = ["*:pc_tau_max"]'

    character(len=512) :: arg, par_arg, steps_arg, dt_arg, spin_arg, repo, root, tools
    character(len=:), allocatable :: d
    integer :: nfail = 0, es, nstep
    logical :: ok, have_tools

    call get_command_argument(1, arg)
    if (trim(arg) == "plain") then
        call get_command_argument(2, par_arg)
        call plain_child(trim(par_arg))
        stop 0
    else if (trim(arg) == "inspect") then
        call get_command_argument(2, par_arg)
        call get_command_argument(3, steps_arg)
        read(steps_arg, *) nstep
        call get_command_argument(4, dt_arg)
        read(dt_arg, *) dt_geo
        call get_command_argument(5, spin_arg)
        call inspect_child(trim(par_arg), nstep, trim(spin_arg) == "spin")
        stop 0
    end if

    call get_environment_variable("PWD", repo)
    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = trim(repo)//"/ice_data"
    if (root(1:1) /= "/") root = trim(repo)//"/"//root
    if (.not. (exists(topo_file(1)) .and. exists(topo_file(2)))) then
        write(*,'(a)') "test_gia_sequence: ice_data topographies of GRL-16KM and ANT-32KM not found under "// &
                       trim(root)//": skipped"
        stop 0
    end if
    call get_environment_variable("KRYOS_TOOLS", tools)
    if (len_trim(tools) == 0) tools = trim(repo)//"/tools"
    call execute_command_line("julia --project='"//trim(tools)//"' -e 'using KryosTools' > /dev/null 2>&1", exitstat=es)
    have_tools = es == 0

    d = trim(repo)//"/"//DIR
    call sh("rm -rf '"//d//"'")
    call sh("mkdir -p '"//d//"maps'")
    call sh("cp '"//trim(repo)//"'/maps/grid_*.txt '"//d//"maps/'")
    call write_data(1, d//"clim_grl.nc")
    call write_data(2, d//"clim_ant.nc")
    call write_background()

    call greenland_ref()                                 ! one domain, ice outside the footprint = X (ext = "ref")
    call greenland_slow()                                ! the Earth every second step
    call no_spinup_stop()                                ! no VILMA2 spin-up, no opt-out: the run stops at init
    if (have_tools) then
        call two_domains()                               ! both, an external history that thins; spin-up; the order swap
    else
        write(*,'(a)') "test_gia_sequence: OPEN: Julia or the project "//trim(tools)//" is not available (the merge "// &
                       "step and the comparison need it; set KRYOS_TOOLS): the two-domain runs are skipped"
    end if

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_gia_sequence: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_gia_sequence: OK"

contains

    ! =========================================================================
    ! The run sets
    ! =========================================================================

    subroutine greenland_ref()
        ! Greenland alone, X outside the footprint, 3 steps: the checks of the inspecting child.
        write(*,'(a)') "--- Greenland, ext = ref"
        call setup_run_dir("g1")
        call write_domain_par(d//"g1/grl.nml", 1, 0.0_dp, 30.0_dp, "None")
        call sh("cp '"//d//"g1/grl.nml' '"//d//"g1/par.nml'")
        call append_globe(d//"g1/par.nml", "ref", 0.0_dp)
        ok = child("g1", "inspect par.nml 3 10 nospin", "run.log")
        call check(ok, "inspecting run g1 (checks in "//DIR//"g1/run.log)")
    end subroutine greenland_ref

    subroutine no_spinup_stop()
        ! gia_check: without pre_spinup_1d / equil_time_max the run stops unless allow_no_spinup (plan 7.4).
        write(*,'(a)') "--- no VILMA2 spin-up"
        call setup_run_dir("ns")
        allow_nospin = .false.
        call write_domain_par(d//"ns/grl.nml", 1, 0.0_dp, 10.0_dp, "None")
        allow_nospin = .true.
        call sh("cp '"//d//"ns/grl.nml' '"//d//"ns/par.nml'")
        call append_globe(d//"ns/par.nml", "ref", 0.0_dp)
        ok = child("ns", "plain par.nml", "run.log", quiet=.true.)
        call check(.not. ok .and. in_file(d//"ns/run.log", "has no VILMA2 spin-up before its seed") .and. &
                   in_file(d//"ns/run.log", "allow_no_spinup"), &
                   "a run without VILMA2's spin-up stops at init, naming the key that accepts it ("//DIR//"ns/run.log)")
    end subroutine no_spinup_stop

    subroutine greenland_slow()
        ! The Earth steps every second ice step: dt = 20 yr, dtt = 10 yr, 4 steps.
        write(*,'(a)') "--- Greenland, ext = ref, the Earth every 20 yr"
        dt_geo = 20.0_dp
        call setup_run_dir("g2")
        call write_domain_par(d//"g2/grl.nml", 1, 0.0_dp, 40.0_dp, "None")
        dt_geo = 10.0_dp
        call sh("cp '"//d//"g2/grl.nml' '"//d//"g2/par.nml'")
        call append_globe(d//"g2/par.nml", "ref", 0.0_dp)
        ok = child("g2", "inspect par.nml 4 20 nospin", "run.log")
        call check(ok, "inspecting run g2 (checks in "//DIR//"g2/run.log)")
    end subroutine greenland_slow

    subroutine two_domains()
        ! Greenland + Antarctica, ext = file: the external ice thins over [-100, -70] yr, so the Earth is loaded
        ! from the start; the same run with the domains swapped.
        write(*,'(a)') "--- Greenland + Antarctica, ext = file; the order swap"
        call two_run("s2", .false., "inspect par.nml 3 10 nospin", spin=.false.)
        allow_nospin = .false.                           ! with the spin-up the opt-out is not needed
        call two_run("s3", .false., "inspect par.nml 2 10 spin", spin=.true.)
        allow_nospin = .true.
        call two_run("s2sw", .true., "plain par.nml", spin=.false.)
        call write_gate("gate_swap.toml", checks= &
            [character(len=400) :: '[[check]]', 'type = "identical"', 'ref = "s2"', 'run = "s2sw"'], &
            chunks="[-100.0, -70.0]", rtimes="[]")
        call gate_check("gate_swap.toml", "swap", "the domain order swapped gives identical files (KryosTools.Gate)")
        call chunks_set()
    end subroutine two_domains

    subroutine chunks_set()
        ! Restarts: the continuous run -100 -> -60 against the chunks -100 -> -80 -> -60 (the bundle at -80: VILMA2's
        ! restart file, the Earth's imports and exports on the Gauss grid, the ports, the domains), with VILMA2's
        ! spin-up. The reference is rebuilt before the restore (plan 8.2); every file and the bundles are compared.
        write(*,'(a)') "--- restarts: chunks == continuous run"
        allow_nospin = .false.
        restart_dt = 20.0_dp
        call two_run_t("c",       -100.0_dp, -60.0_dp, "None",           "plain par.nml", spin=.true.)
        call two_run_t("c2",      -100.0_dp, -60.0_dp, "None",           "plain par.nml", spin=.true.)
        call two_run_t("chunk/c00", -100.0_dp, -80.0_dp, "None",         "plain par.nml", spin=.true.)
        call two_run_t("chunk/c01", -80.0_dp, -60.0_dp, "../c00/restart--80", "plain par.nml", spin=.true.)
        restart_dt = 0.0_dp
        allow_nospin = .true.
        call write_gate("gate_chunk.toml", checks= &
            [character(len=400) :: &
             '[[check]]', 'type = "identical"', 'ref = "c"', 'run = "c2"', &
             '[[check]]', 'type = "chunks"', 'ref = "c"', 'run = "chunk"'], &
            chunks="[-100.0, -80.0, -60.0]", rtimes="[-80.0, -60.0]")
        call gate_check("gate_chunk.toml", "chunk", "cont2 == cont (determinism), chunks == cont: every file and the bundles "// &
                        "at -80 and -60 (KryosTools.Gate)")
    end subroutine chunks_set

    subroutine two_run_t(sub, t0, t1, restart, args, spin)
        ! Greenland + Antarctica (ext = file) over [t0, t1] from `restart`, run as `args`.
        character(len=*), intent(in) :: sub, restart, args
        real(dp),         intent(in) :: t0, t1
        logical,          intent(in) :: spin
        integer :: es
        call setup_run_dir(sub)
        call write_domain_par(d//sub//"/grl.nml", 1, t0, t1, restart)
        call write_domain_par(d//sub//"/ant.nml", 2, t0, t1, restart)
        call execute_command_line("cd '"//d//sub//"' && julia --project='"//trim(tools)//"' '"//trim(repo)// &
                                  "/tools/bin/kryos-merge.jl' grl.nml ant.nml --out par.nml --report merge.md > merge.log 2>&1", &
                                  exitstat=es)
        if (es /= 0) call sh("tail -n 20 '"//d//sub//"/merge.log'")
        call check(es == 0, "merge of the two single-domain par files ("//DIR//sub//"/merge.log)")
        call append_globe(d//sub//"/par.nml", "file", t0, spin=spin)
        ok = child(sub, args, "run.log")
        call check(ok, "two-domain run "//sub//" ("//DIR//sub//"/run.log)")
    end subroutine two_run_t

    subroutine two_run(sub, swap, args, spin)
        character(len=*), intent(in) :: sub, args
        logical,          intent(in) :: swap, spin
        integer :: es
        call setup_run_dir(sub)
        call write_domain_par(d//sub//"/grl.nml", 1, -100.0_dp, -70.0_dp, "None")
        call write_domain_par(d//sub//"/ant.nml", 2, -100.0_dp, -70.0_dp, "None")
        call execute_command_line("cd '"//d//sub//"' && julia --project='"//trim(tools)//"' '"//trim(repo)// &
                                  "/tools/bin/kryos-merge.jl' grl.nml ant.nml --out par.nml --report merge.md"// &
                                  trim(merge(" --swap", "       ", swap))//" > merge.log 2>&1", exitstat=es)
        if (es /= 0) call sh("tail -n 20 '"//d//sub//"/merge.log'")
        call check(es == 0, "merge of the two single-domain par files ("//DIR//sub//"/merge.log)")
        call append_globe(d//sub//"/par.nml", "file", -100.0_dp, spin=spin)
        ok = child(sub, args, "run.log")
        call check(ok, "two-domain run "//sub//trim(merge(" (domains swapped)", "                  ", swap))// &
                   " ("//DIR//sub//"/run.log)")
    end subroutine two_run

    subroutine write_gate(file, checks, chunks, rtimes)
        ! A kryos-against-kryos gate file for the runs of this test.
        character(len=*),  intent(in) :: file, chunks, rtimes
        character(len=*),  intent(in) :: checks(:)
        integer :: u, i
        open(newunit=u, file=d//file, status="replace", action="write")
        write(u,'(a)') 'gate = "test_gia_sequence"', 'kind = "kryos"', "chunks = "//chunks, &
                       "restart_times = "//rtimes, 'domain_order = ["Greenland", "Antarctica"]', PCT, &
                       '[domains.Greenland.instances]', 'domain = "Greenland"', 'ice = "ice_grl"', 'geo = "geo_grl"', &
                       'climate = "climate_grl"', 'surface = "surface_grl"', 'shelf = "shelf_grl"', &
                       'geothermal = "ghf_grl"', 'sediments = "sed_grl"', 'forcing = "forcing"', 'geo_glob = "geo_glob"', &
                       '[domains.Antarctica.instances]', 'domain = "Antarctica"', 'ice = "ice_ant"', 'geo = "geo_ant"', &
                       'climate = "climate_ant"', 'surface = "surface_ant"', 'shelf = "shelf_ant"', &
                       'geothermal = "ghf_ant"', 'sediments = "sed_ant"', 'forcing = "forcing"', 'geo_glob = "geo_glob"'
        do i = 1, size(checks)
            write(u,'(a)') trim(checks(i))
        end do
        close(u)
    end subroutine write_gate

    subroutine gate_check(file, tag, what)
        ! tools/bin/kryos-gate-kk.jl check: the checks of the gate file over the run folders of this test.
        character(len=*), intent(in) :: file, tag, what
        integer :: es
        call execute_command_line("cd '"//d//"' && julia --project='"//trim(tools)//"' '"//trim(repo)// &
                                  "/tools/bin/kryos-gate-kk.jl' check "//file//" . --md gate_"//tag//".md > gate_"//tag// &
                                  ".log 2>&1", exitstat=es)
        if (es /= 0) call sh("grep -E '\| fail \|' '"//d//"gate_"//tag//".md' | head -n 30; tail -n 10 '"//d//"gate_"//tag//".log'")
        call check(es == 0, what//" (report "//DIR//"gate_"//tag//".md)")
    end subroutine gate_check

    ! =========================================================================
    ! The runs (child processes, in the run directory)
    ! =========================================================================

    subroutine plain_child(par_file)
        character(len=*), intent(in) :: par_file
        type(run_t), target :: run
        call run_init(run, par_file)
        do while (.not. run_done(run))
            call run_step(run)
        end do
        call run_end(run)
        write(*,'(a)') "test_gia_sequence (plain): OK"
    end subroutine plain_child

    subroutine inspect_child(par_file, nsteps, spin)
        ! The run, from inside: run_init, the checks of the cold start, then every step with the standalone
        ! gather beforehand and the checks of the gather and the scatter afterwards.
        character(len=*), intent(in) :: par_file
        integer,          intent(in) :: nsteps
        logical,          intent(in) :: spin         ! the run has VILMA2's own spin-up before the seed

        type(run_t), target :: run
        type(globe_ice_t) :: gi2
        type(coupler_class) :: cpl2, cpl3
        real(dp) :: t_geo, dtt
        logical :: due
        integer :: k

        call run_init(run, par_file)
        call coupler_init(cpl3, trim(run%cfg%map_fldr))
        call coupler_add_grid(cpl3, trim(run%globe%geo%c%grid), run%globe%geo%c%grid_obj)
        call standalone_init(run, gi2, cpl2)

        call check_cold_start(run, cpl3, spin)
        dtt = clock_dtt(run%clk)
        t_geo = clock_time(run%clk)                                          ! the Earth's last solution
        do k = 1, nsteps
            due = abs((clock_time(run%clk) + dtt)/dt_geo - real(nint((clock_time(run%clk) + dtt)/dt_geo), dp)) < 1.0e-9_dp
            if (due) call standalone_gather(run, gi2, cpl2, clock_time(run%clk))      ! the ice of the previous step
            call run_step(run)
            if (due) then
                t_geo = clock_time(run%clk)
                call check_gather(run, gi2, k)
            end if
            call check_scatter(run, cpl3, k, t_geo)
        end do
        call check_moved(run)
        call run_end(run)
        call globe_ice_end(gi2)

        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_gia_sequence (inspect): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_gia_sequence (inspect): OK"
    end subroutine inspect_child

    ! =========================================================================
    ! The checks
    ! =========================================================================

    subroutine check_cold_start(run, cpl, spin)
        ! The two-pass cold start: after run_init the Earth stands on the reference, loaded with the t0 ice, seeded, and
        ! the domains are on its solution. Yelmo's initial ice is not exactly its present-day ice (its own init moves
        ! some cells), so the seed is not at rest; what is checked is that the t0 load is the gathered final ice, the
        ! reference is the one the domains have, and the anomalies are consistent with that load.
        type(run_t), target, intent(in)    :: run
        type(coupler_class), intent(inout) :: cpl
        logical,             intent(in)    :: spin
        real(dp), allocatable :: a(:,:), h(:,:), fg(:,:)
        type(grid_class) :: g
        real(dp) :: dI, dC, dF, dphi, t0, t1, bsl_ive, esl, v_ice, fac, a_now
        integer :: i, k

        write(*,'(a)') "cold start"
        t0 = clock_time(run%clk)
        associate (c => run%globe%geo%c, gi => run%globe%gi)
            call check(run%globe%has_gi .and. trim(c%grid) == "GAUSS-L16" .and. c%nx == 64 .and. c%ny == 34, &
                       "the global geo is on GAUSS-L16 (64 x 34)")
            call check(allocated(c%grid_obj) .and. trim(gi%tgt) == "GAUSS-L16", "its grid object is the target of globe_ice")
            do i = 1, run%n_dom
                call check(abs(cadence_dt_years(run%dom(i)%geo%c%cadence, run%clk) - cadence_dt_years(c%cadence, run%clk)) &
                           == 0.0_dp, "the port of "//trim(run%dom(i)%name)//" steps on the global cadence")
            end do

            ! the t0 load: the final ice at t0, gathered once; the host terms differ from the reference where Yelmo's
            ! initial ice differs from its present day
            call check(gi%t_load == t0, "the t0 load is the gather at t0")
            allocate(a(c%nx, c%ny))
            call fields_get(c%imports, c%hi_host_Ig, a, t1=t1)
            call check(all(a == gi%host_Ig) .and. t1 == t0, "the Earth's host_Ig is that gather, published at t0")
            call fields_get(c%imports, c%hi_host_Ig_ref, a, t1=t1)
            call check(all(a == gi%host_Ig_ref) .and. t1 == t0, "the Earth's reference host_Ig_ref is the reference of globe_ice (at t0)")
            dI = maxabs_diff(gi%host_Ig, gi%host_Ig_ref)
            dC = maxabs_diff(gi%host_C, gi%host_C_ref)
            dF = maxval(abs(gi%host_F))
            write(*,'(a,3es11.3)') "    t0 gather minus the reference: max |I_g| [m], |C|, |F| [m] = ", dI, dC, dF
            call check(dI < 0.25_dp*maxval(gi%host_Ig_ref) .and. dC < 0.5_dp, "the t0 load is near the reference")
            call check(maxval(gi%host_w) > 0.0_dp .and. maxval(gi%host_Ig_ref) > 0.0_dp, "the host region holds ice")

            ! the grounded ice of the load is the ice model's own H_ice f_grnd
            esl = 0.0_dp
            do k = 1, gi%diag%nd
                do i = 1, run%n_dom
                    if (trim(run%dom(i)%name) /= trim(gi%diag%name(k))) cycle
                    associate (ice => run%dom(i)%ice%c)
                        allocate(h(ice%nx, ice%ny), fg(ice%nx, ice%ny))
                        call fields_get(ice%exports, ice%h_f_grnd, fg)
                        call grid_cdo_read_desc(g, trim(run%dom(i)%grid(C_ICE)), trim(run%cfg%map_fldr))
                        call fields_get(ice%exports, ice%h_H_ice, h)     ! after the reconcile: same ice (kept)
                        v_ice = sum(h*fg*g%area)*g%cs%xy_conv**2
                        ! (after a large seed the reconcile moves the grounding line: the comparison holds for the small seeds)
                        if (.not. spin) call check(abs(gi%diag%V_g(k) - v_ice) <= 1.0e-3_dp*v_ice, trim(gi%diag%name(k))// &
                                   ": the gathered grounded volume is the ice model's H_ice f_grnd (true areas, 0.1 %)")
                        deallocate(h, fg)
                    end associate
                end do
                esl = esl + gi%diag%esl(k)
            end do

            ! the seed's ice-volume equivalent bsl_ive is the t0 load's sea-level equivalent over VILMA2's ocean area
            ! (the synthetic bed has less ocean than the Earth): the domains' ice, and the external ice change too
            ! when VILMA2's seed applies it (its dt = 0 branch solves with the previous load, known issue F7: without
            ! a spin-up the external ice enters with the first call)
            call fields_get(c%exports, c%h_bsl_ive, bsl_ive)
            fac = gi%area_ocn/gi%diag%a_ocn_ref
            write(*,'(a,3es12.4)') "    t0 sea-level equivalent [m] over VILMA2's ocean area: domains, + external, VILMA2 bsl_ive = ", &
                  esl*fac, (esl + gi%diag%ext_esl)*fac, bsl_ive
            if (spin) then                                   ! the spin-up holds the whole t0 load: the seed has it
                call fields_get(c%exports, c%h_A_ocn, a_now)     ! tens of metres of sea level: the coastline has moved
                esl = (esl + gi%diag%ext_esl)*gi%area_ocn/a_now
                write(*,'(a,es12.4)') "    the whole t0 load over A_ocn of the seed: ", esl
                call check(abs(bsl_ive - esl) <= 2.0e-3_dp*abs(esl) + 1.0e-9_dp, &
                           "with VILMA2's spin-up the seed has the whole t0 load, external ice included (0.2 %)")
            else
                call check(abs(bsl_ive - esl*fac) <= 0.01_dp*abs(esl*fac) + 1.0e-9_dp .or. &
                           abs(bsl_ive - (esl + gi%diag%ext_esl)*fac) <= 0.01_dp*abs((esl + gi%diag%ext_esl)*fac) + 1.0e-9_dp, &
                           "VILMA2's seed has the gathered ice-volume equivalent (1 %: its ocean area, partial cells)")
            end if

            call fields_get(c%exports, c%h_dz_bed, a, t1=t1)
            write(*,'(a,es11.3)') "    seed: max |dz_bed| [m] = ", maxval(abs(a))
            call fields_get(c%exports, c%h_dphi, dphi, t0=t0)
            call check(t1 == clock_time(run%clk) .and. t0 == t1 .and. maxval(abs(a)) < merge(500.0_dp, 1.0_dp, spin), &
                       "the seed: the exports are at t0, a bed change of mm to cm (the whole load with a spin-up: metres)")
            write(*,'(a,es11.3)') "    seed: dphi [m] = ", dphi
        end associate

        call check_scatter(run, cpl, 0, clock_time(run%clk))   ! the ports and the ice stand on the seed
    end subroutine check_cold_start

    subroutine check_gather(run, gi2, k)
        ! The host terms VILMA2 got at step k are the standalone gather's of the previous step's ice.
        type(run_t), target, intent(in) :: run
        type(globe_ice_t),   intent(in) :: gi2
        integer,             intent(in) :: k

        real(dp), allocatable :: a(:,:)
        real(dp) :: t1, bsl_ive, esl, a_now
        character(len=16) :: tag

        write(tag, '(a,i0)') "step ", k
        write(*,'(a)') trim(tag)//": gather"
        associate (c => run%globe%geo%c)
            allocate(a(c%nx, c%ny))
            call fields_get(c%imports, c%hi_host_Ig, a, t1=t1)
            call check(all(a == gi2%host_Ig), trim(tag)//": host_Ig is the standalone gather's, bitwise")
            call fields_get(c%imports, c%hi_host_C, a)
            call check(all(a == gi2%host_C), trim(tag)//": host_C bitwise")
            call fields_get(c%imports, c%hi_host_F, a)
            call check(all(a == gi2%host_F), trim(tag)//": host_F bitwise")
            call fields_get(c%imports, c%hi_host_w, a)
            call check(all(a == gi2%host_w), trim(tag)//": host_w bitwise")
            call fields_get(c%imports, c%hi_H_ice_ext, a)
            call check(all(a == gi2%H_ice_ext), trim(tag)//": H_ice_ext bitwise")
            call check(t1 == gi2%t_load .and. t1 == clock_time(run%clk) - clock_dtt(run%clk), &
                       trim(tag)//": the load is the previous step's (lag 1)")
            write(*,'(a,es11.3,a,es11.3)') "    ext. ice: max |H_ext - X| [m] = ", maxval(abs(a - run%globe%gi%H_ice_ref)), &
                  ", max |I_g - ref| [m] = ", maxabs_diff(run%globe%gi%host_Ig, run%globe%gi%host_Ig_ref)
            ! the call ramped to this load: VILMA2's ice-volume equivalent is its sea-level equivalent over its ocean
            ! area at the end of the call (the A_ocn export: with tens of metres of sea-level change the coastline of the
            ! synthetic bed moves, so the reference's area would not do)
            call fields_get(c%exports, c%h_bsl_ive, bsl_ive)
            call fields_get(c%exports, c%h_A_ocn, a_now)
            esl = (sum(gi2%diag%esl) + gi2%diag%ext_esl)*gi2%area_ocn/a_now
            write(*,'(a,2es12.4)') "    sea-level equivalent of the load [m]: gathered over A_ocn, VILMA2 bsl_ive = ", esl, bsl_ive
            call check(abs(bsl_ive - esl) <= 2.0e-3_dp*abs(esl) + 1.0e-9_dp, &
                       trim(tag)//": VILMA2's ice-volume equivalent is the gathered load's (0.2 %)")
        end associate
    end subroutine check_gather

    subroutine check_scatter(run, cpl, k, t_geo)
        ! The port holds the Gauss-grid anomalies bilinearly remapped; the ice's z_sl = dphi + dz_ss, z_bed = reference + dz_bed.
        type(run_t), target, intent(in)    :: run
        type(coupler_class), intent(inout) :: cpl
        integer,             intent(in)    :: k
        real(dp),            intent(in)    :: t_geo        ! the Earth's last solution (the step's time when it is due)

        real(dp), allocatable :: g(:,:), r(:,:), p(:,:), v(:,:)
        real(dp) :: dphi, bsl, t1, dphi2
        character(len=NAME_LEN) :: nm
        character(len=16) :: tag
        integer :: i

        write(tag, '(a,i0)') "step ", k
        write(*,'(a)') trim(tag)//": scatter"
        associate (c => run%globe%geo%c)
            allocate(g(c%nx, c%ny))
            call fields_get(c%exports, c%h_dphi, dphi, t1=t1)
            call fields_get(c%exports, c%h_bsl, bsl)
            call check(t1 == t_geo, trim(tag)//": the Earth's exports are those of its last solution (t = "// &
                       trim(adjustl(ftoa(t_geo)))//")")
            write(*,'(a,es11.3,a,es11.3)') "    dphi [m] = ", dphi, ", bsl [m] = ", bsl
            do i = 1, run%n_dom
                nm = trim(run%dom(i)%name)
                associate (geo => run%dom(i)%geo%c, ice => run%dom(i)%ice%c)
                    allocate(p(geo%nx, geo%ny), v(geo%nx, geo%ny))
                    call fields_get(c%exports, c%h_dz_bed, g)
                    call remap(cpl, g, trim(c%grid), r, trim(geo%grid), method="bilin")
                    call fields_get(geo%exports, geo%h_dz_bed, p)
                    call check(all(p == r), trim(nm)//": the port's dz_bed is the bilinear remap of the Earth's, bitwise")
                    write(*,'(a,es11.3)') "    "//trim(nm)//": max |dz_bed| on the ice grid [m] = ", maxval(abs(p))
                    call fields_get(ice%imports, ice%hi_z_bed, v)
                    call fields_get(geo%exports, geo%h_z_bed_ref, r)
                    call check(maxval(abs(v - (r + p))) <= 2.0e-6_dp*maxval(abs(v)) + 1.0e-6_dp, &
                               trim(nm)//": the ice's z_bed is the reference bed + dz_bed")

                    call fields_get(c%exports, c%h_dz_ss, g)
                    call remap(cpl, g, trim(c%grid), r, trim(geo%grid), method="bilin")
                    call fields_get(geo%exports, geo%h_dz_ss, p)
                    call check(all(p == r), trim(nm)//": the port's dz_ss is the bilinear remap of the Earth's, bitwise")
                    call fields_get(ice%imports, ice%hi_z_sl, v)
                    call check(maxval(abs(v - (p + dphi))) <= 2.0e-6_dp*maxval(abs(v)) + 1.0e-6_dp, &
                               trim(nm)//": the ice's z_sl = dphi + dz_ss on the ice grid")
                    call fields_get(geo%exports, geo%h_dphi, dphi2, t1=t1)
                    call check(dphi2 == dphi .and. t1 == t_geo, trim(nm)//": the port's dphi is the Earth's, at its last solution")
                    deallocate(p, v, r)

                end associate
            end do
        end associate
    end subroutine check_scatter

    subroutine check_moved(run)
        ! The Earth responded to the load of the run: a bed change on the ice grids and a sea-level shift.
        type(run_t), target, intent(in) :: run
        real(dp), allocatable :: p(:,:)
        real(dp) :: dphi
        integer :: i

        write(*,'(a)') "response"
        associate (c => run%globe%geo%c)
            call fields_get(c%exports, c%h_dphi, dphi)
            do i = 1, run%n_dom
                associate (geo => run%dom(i)%geo%c)
                    allocate(p(geo%nx, geo%ny))
                    call fields_get(geo%exports, geo%h_dz_bed, p)
                    call check(maxval(abs(p)) > 0.0_dp .or. abs(dphi) > 0.0_dp, &
                               trim(run%dom(i)%name)//": the bed or the sea level has moved (the load changed)")
                    deallocate(p)
                end associate
            end do
        end associate
    end subroutine check_moved

    ! =========================================================================
    ! The standalone gather: a globe_ice_t of the same files, fed by the ice's exports
    ! =========================================================================

    subroutine standalone_init(run, gi2, cpl2)
        type(run_t), target, intent(in)    :: run
        type(globe_ice_t),   intent(inout) :: gi2
        type(coupler_class), intent(inout) :: cpl2

        real(dp), allocatable :: fp(:,:), hpd(:,:), zbr(:,:)
        type(grid_class) :: g
        integer :: i

        call coupler_init(cpl2, trim(run%cfg%map_fldr))
        call globe_ice_init(gi2, run%cfg, cpl2, run%cnst, run%globe%geo%c%grid_obj, group="globe_ice", &
                            time_ref=clock_to_years(run%clk, run%clk%n_ref), &
                            time_span=[clock_time_init(run%clk), clock_time_end(run%clk)])
        do i = 1, run%n_dom
            associate (ice => run%dom(i)%ice)
                allocate(fp(ice%c%nx, ice%c%ny), hpd(ice%c%nx, ice%c%ny), zbr(ice%c%nx, ice%c%ny))
                call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                call fields_get(ice%c%exports, ice%c%h_H_ice_pd,  hpd)
                call fields_get(ice%c%exports, ice%c%h_z_bed_ref, zbr)
                call grid_cdo_read_desc(g, trim(run%dom(i)%grid(C_ICE)), trim(run%cfg%map_fldr))
                call globe_ice_add_domain(gi2, cpl2, trim(run%dom(i)%name), g, fp, hpd, zbr)
                deallocate(fp, hpd, zbr)
            end associate
        end do
        call globe_ice_reference_begin(gi2, cpl2)
        do i = 1, run%n_dom
            call globe_ice_reference_add(gi2, cpl2, trim(run%dom(i)%name), run%dom(i)%ice)
        end do
        call globe_ice_reference_end(gi2, cpl2)
    end subroutine standalone_init

    subroutine standalone_gather(run, gi2, cpl2, t_load)
        type(run_t), target, intent(in)    :: run
        type(globe_ice_t),   intent(inout) :: gi2
        type(coupler_class), intent(inout) :: cpl2
        real(dp),            intent(in)    :: t_load

        real(dp), allocatable :: H(:,:), fg(:,:), zb(:,:), zs(:,:), fp(:,:)
        integer :: i

        call globe_ice_gather_begin(gi2, cpl2, t_load)
        do i = 1, run%n_dom
            associate (ice => run%dom(i)%ice)
                allocate(H(ice%c%nx, ice%c%ny), fg(ice%c%nx, ice%c%ny), zb(ice%c%nx, ice%c%ny), zs(ice%c%nx, ice%c%ny), &
                         fp(ice%c%nx, ice%c%ny))
                call fields_get(ice%c%exports, ice%c%h_H_ice,     H)
                call fields_get(ice%c%exports, ice%c%h_f_grnd,    fg)
                call fields_get(ice%c%exports, ice%c%h_z_bed,     zb)
                call fields_get(ice%c%exports, ice%c%h_z_sl,      zs)
                call fields_get(ice%c%exports, ice%c%h_footprint, fp)
                call globe_ice_gather_add(gi2, cpl2, trim(run%dom(i)%name), ice, H, fg, zb, zs, fp)
                deallocate(H, fg, zb, zs, fp)
            end associate
        end do
        call globe_ice_gather_end(gi2, cpl2)
    end subroutine standalone_gather

    real(dp) function maxabs_diff(a, b)
        real(dp), intent(in) :: a(:,:), b(:,:)
        maxabs_diff = maxval(abs(a - b))
    end function maxabs_diff

    ! =========================================================================
    ! Inputs
    ! =========================================================================

    subroutine write_background()
        ! The background grid (description in the shared map folder), and the history file TBG-4DEG_IceHistory-TEST.nc:
        ! the bed (constant) and the ice of the external sheet only (a cap over northern Europe, none in the domains),
        ! at 150, 100, 50 and 0 yr ago; the sheet is 1 + age/100 times its present size, so it thins towards 0 ka.
        character(len=*), parameter :: ages_name = "time"
        type(grid_class) :: g
        real(dp) :: ages(4), zb(NXB,NYB,4), h(NXB,NYB,4), lon, lat, dist
        character(len=:), allocatable :: f
        integer :: i, j, k

        call grid_init(g, name=BGN, mtype="latlon", units="degrees", planet="Spherical Earth", &
                       x0=2.0_dp, dx=4.0_dp, nx=NXB, y0=-88.0_dp, dy=4.0_dp, ny=NYB)
        call grid_cdo_write_desc_short(g, d//"maps")
        ages = [150.0_dp, 100.0_dp, 50.0_dp, 0.0_dp]
        do k = 1, 4
            do j = 1, NYB
                do i = 1, NXB
                    lon = 2.0_dp + 4.0_dp*(i - 1); lat = -88.0_dp + 4.0_dp*(j - 1)
                    zb(i,j,k) = -1500.0_dp*sin(0.07_dp*i)*cos(0.11_dp*j) + 200.0_dp
                    dist = gc_deg(lon, lat, 20.0_dp, 62.0_dp)
                    h(i,j,k) = max(0.0_dp, 1800.0_dp*(1.0_dp - dist/14.0_dp))*(1.0_dp + ages(k)/100.0_dp)
                end do
            end do
        end do
        f = d//BGN//"_IceHistory-TEST.nc"
        call nc_create(f)
        call nc_write_dim(f, "lon", x=2.0_dp, dx=4.0_dp, nx=NXB, units="degrees_east")
        call nc_write_dim(f, "lat", x=-88.0_dp, dx=4.0_dp, nx=NYB, units="degrees_north")
        call nc_write_dim(f, ages_name, x=ages, units="years ago", unlimited=.false.)
        call nc_write(f, "z_bed", zb, dim1="lon", dim2="lat", dim3=ages_name, units="m")
        call nc_write(f, "H_ice", h, dim1="lon", dim2="lat", dim3=ages_name, units="m")
        call nc_write_attr(f, "grid_name", BGN)
    end subroutine write_background

    real(dp) function gc_deg(lon1, lat1, lon2, lat2)
        ! great-circle distance [deg]
        real(dp), intent(in) :: lon1, lat1, lon2, lat2
        real(dp) :: a, b, dl
        a = lat1*PI/180.0_dp; b = lat2*PI/180.0_dp; dl = (lon1 - lon2)*PI/180.0_dp
        gc_deg = acos(max(-1.0_dp, min(1.0_dp, sin(a)*sin(b) + cos(a)*cos(b)*cos(dl))))*180.0_dp/PI
    end function gc_deg

    subroutine append_globe(file, ext, tref, spin)
        ! The groups of the global part (VILMA2 at lmax 16, the global ice field on the synthetic background) after the
        ! domains' groups. ext = "ref": the ice outside the footprints is X; "file": X plus the history's anomaly.
        character(len=*), intent(in) :: file, ext
        real(dp),         intent(in) :: tref       ! [yr] unused for ref; the window starts here
        logical,          intent(in), optional :: spin   ! VILMA2's own spin-up before the seed (pre_spinup_1d)
        integer :: u
        associate (unused => tref); end associate
        open(newunit=u, file=file, status="old", position="append", action="write")
        write(u,'(a)') "&vilma", "  lmax = 16"
        if (present(spin)) then
            if (spin) write(u,'(a)') "  pre_spinup_1d = .true."
        end if
        write(u,'(a)') "/"
        call w(u, [character(len=400) :: &
            "&globe_ice", "  grid = '"//BGN//"'", "  z_bed_ref = 'gi_X_zb'", "  H_ice_ref = 'gi_X_H'", &
            "  ext = '"//trim(ext)//"'", "  ext_field = 'gi_hist'", "  ext_tref = 0.0", "/", &
            "&gi_X_zb", "  filename = '"//d//"{grid_name}_IceHistory-TEST.nc'", "  name = 'z_bed'", &
            "  mode = 'transient'", "/", &
            "&gi_X_H", "  filename = '"//d//"{grid_name}_IceHistory-TEST.nc'", "  name = 'H_ice'", &
            "  mode = 'transient'", "/", &
            "&gi_hist", "  filename = '"//d//"{grid_name}_IceHistory-TEST.nc'", "  name = 'H_ice'", &
            "  mode = 'transient'", "  anomaly_ref = 'slice'", "  ref_time = 0.0", "  out_of_range = 'error'", "/"])
        close(u)
    end subroutine append_globe

    function topo_file(id) result(f)
        integer, intent(in) :: id
        character(len=:), allocatable :: f
        f = trim(root)//"/"//trim(DNAME(id))//"/"//trim(DGRID(id))//"/"//trim(DGRID(id))//"_"//trim(DTOPO(id))//".nc"
    end function topo_file
    subroutine setup_run_dir(sub)
        ! tmp/test_run_domains/<sub>: input/ (links), maps, the climate data (linked).
        character(len=*), intent(in) :: sub
        call sh("mkdir -p '"//d//sub//"/input'")
        call sh("ln -s '"//d//"maps' '"//d//sub//"/maps'")                 ! shared: grids, the weights made once
        call sh("for f in '"//trim(repo)//"'/input/*; do ln -s ""$f"" '"//d//sub//"/input/'; done")
        call sh("for f in '"//trim(repo)//"'/yelmo/input/*.nml '"//trim(repo)//"'/yelmo/input/*.md; do "// &
                "ln -sf ""$f"" '"//d//sub//"/input/'; done")
        call sh("ln -sf '"//trim(repo)//"/chion/input/chion_defaults.nml' '"//d//sub//"/input/'")
        call sh("ln -sf '"//trim(repo)//"/fesm-utils/par/datasrc_defaults.nml' '"//d//sub//"/input/'")
        call sh("for f in '"//trim(repo)//"'/insol/input/LA2004.*.txt; do ln -sf ""$f"" '"//d//sub//"/input/'; done")
        call sh("ln -s '"//d//"clim_grl.nc' '"//d//sub//"/clim_grl.nc'; ln -s '"//d//"clim_ant.nc' '"//d//sub//"/clim_ant.nc'")
    end subroutine setup_run_dir
    subroutine write_data(id, file)
        ! Synthetic snapshots on the grid of domain id: the reference pd (tas from a lapse rate, the latitude and a
        ! seasonal cycle, pr varying, zs = the topography's surface, ocean profiles), ctl = pd and alt = 10 K
        ! colder, half the pr.
        integer,          intent(in) :: id
        character(len=*), intent(in) :: file

        type(grid_class) :: grid
        real(sp), allocatable :: zs(:,:), t(:,:,:), p(:,:,:), o(:,:,:), s(:,:,:)
        real(sp), parameter :: prof(NZ) = [273.5_sp, 273.0_sp, 274.0_sp, 274.5_sp, 272.5_sp]
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: t0, dlat
        integer :: i, j, m, nx, ny

        nx = DNX(id); ny = DNY(id)
        allocate(zs(nx,ny), t(nx,ny,12), p(nx,ny,12), o(nx,ny,NZ), s(nx,ny,NZ))
        call grid_cdo_read_desc(grid, trim(DGRID(id)), trim(repo)//"/maps")
        call nc_read(topo_file(id), "z_srf", zs)
        zs = max(zs, 0.0_sp)
        t0 = 281.15_dp; dlat = 60.0_dp
        if (id == 2) then
            t0 = 262.15_dp; dlat = 70.0_dp
        end if
        do j = 1, ny
            do i = 1, nx
                do m = 1, 12
                    t(i,j,m) = real(t0 - 0.8_dp*(abs(grid%lat(i,j)) - dlat) - 0.0065_dp*zs(i,j) &
                                    + 10.0_dp*cos(2.0_dp*pi*real(m - 7, dp)/12.0_dp), sp)
                    p(i,j,m) = real(0.3_dp + 0.013_dp*mod(i, 40) + 0.007_dp*mod(j, 40) + 0.11_dp*m, sp)
                end do
                o(i,j,:) = prof
                s(i,j,:) = 34.5_sp
            end do
        end do
        call nc_create(file)
        call nc_write_dim(file, "xc",    x=grid%G%x)
        call nc_write_dim(file, "yc",    x=grid%G%y)
        call nc_write_dim(file, "month", x=1.0_dp, dx=1.0_dp, nx=12)
        call nc_write_dim(file, "depth", x=[0.0_dp, 200.0_dp, 500.0_dp, 1000.0_dp, 3000.0_dp])
        call nc_write(file, "t_pd",  t,                dim1="xc", dim2="yc", dim3="month")
        call nc_write(file, "p_pd",  p,                dim1="xc", dim2="yc", dim3="month")
        call nc_write(file, "t_alt", t - 10.0_sp,      dim1="xc", dim2="yc", dim3="month")
        call nc_write(file, "p_alt", 0.5_sp*p,         dim1="xc", dim2="yc", dim3="month")
        call nc_write(file, "zs",    zs,               dim1="xc", dim2="yc")
        call nc_write(file, "to",    o,                dim1="xc", dim2="yc", dim3="depth")
        call nc_write(file, "so",    s,                dim1="xc", dim2="yc", dim3="depth")
    end subroutine write_data

    subroutine write_domain_par(file, id, time_init, time_end, restart, exps, dto)
        ! The single-domain par file of domain id as the converter writes one: the domain's instance groups carry
        ! its suffix, the physics and snapshot groups do not (they collide in a file with both domains: the merge
        ! step renames them). Library defaults except where the two domains differ (values of yelmox_Greenland.nml and
        ! yelmox_Antarctica.nml where they matter). `exps`: the experiment of the domain; `dto`: a term on the
        ! derived coefficient dTo (air anomaly) in the climate.
        character(len=*), intent(in) :: file, restart
        integer,          intent(in) :: id
        real(dp),         intent(in) :: time_init, time_end
        character(len=*), intent(in), optional :: exps
        logical,          intent(in), optional :: dto

        character(len=:), allocatable :: nm, sx, grid, topo, basins, exp_key, terms, coefs
        logical :: with_dto, has_terms, has_ocean
        integer :: u, nz_aa
        real(dp) :: f_p

        nm = trim(DNAME(id)); sx = trim(DSFX(id)); grid = trim(DGRID(id))
        topo   = trim(root)//"/{domain}/{grid_name}/{grid_name}_"//trim(DTOPO(id))//".nc"
        basins = trim(root)//"/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"
        exp_key = ""
        if (present(exps)) exp_key = "  experiments = '"//exps//"'"
        with_dto = .false.
        if (present(dto)) with_dto = dto
        has_terms = .true.; has_ocean = .false.
        terms = "  terms = 'lgm'"
        coefs = "  coefs = 'alpha'"
        if (id == 1 .or. with_dto) then         ! as yelmox_North: alpha const + the derived dTo (term ocean)
            has_ocean = .true.
            terms = "  terms = 'lgm' 'ocean'"
            coefs = "  coefs = 'alpha' 'dTo'"
        else                                    ! the Antarctica base: no forcing coefficients, no terms
            has_terms = .false.
            terms = ""
            coefs = "  coefs = ''"
        end if
        nz_aa = 6; f_p = 0.06_dp
        if (id == 2) then
            nz_aa = 5; f_p = 0.04_dp
        end if

        open(newunit=u, file=file, status="replace", action="write")
        call w(u, [character(len=400) :: &
            "&kryos", "  domains = '"//nm//"'", "  restart = '"//restart//"'", "  restart_out = '"//trim(merge("out_rst", "none   ", restart_dt > 0.0_dp))//"'", &
            "  write_lossless = True", "  log_level = 'info'", "/", &
            "&clock", "  method = 'const'", "  time_const = 0.0", "  time_init = "//ftoa(time_init), &
            "  time_end = "//ftoa(time_end), "  dtt = 10.0", "/", &
            "&out_2D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_1D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_rst", "  method = 'const'", "  dt = "//ftoa(max(restart_dt, 1.0_dp)), "/", &
            "&globe", "  forcing = 'forcing'", "  geo = 'geo_glob'", "/", &
            "&forcing", "  backend = 'tsgen'", "  clock = 'forcing'", coefs, "  out_1D = 'out_1D'", "/", &
            "&geo_glob", "  backend = 'vilma2'", "  dt = "//ftoa(dt_geo), "  allow_no_spinup = "//trim(merge("True ", "False", allow_nospin)), &
            "  out_1D = 'out_1D'", "/", &
            "&"//nm, "  grid_hub = '"//grid//"'", "  topo = 'topo'", "  regions = 'regions'", &
            "  basins = 'basin'", "  out_2D = 'out_2D'", "  ice = 'ice_"//sx//"'", "  geo = 'geo_"//sx//"'", &
            "  climate = 'climate_"//sx//"'", "  surface = 'surface_"//sx//"'", "  shelf = 'shelf_"//sx//"'", &
            "  geothermal = 'ghf_"//sx//"'", "  sediments = 'sed_"//sx//"'", exp_key, "/", &
            "&topo", "  path = '"//topo//"'", "  remove_lakes = True", "  H_ice_min = 1.0", &
            "  grad_lim_zb = 0.5", "/", &
            "&regions", "  layers = 'basin'", "  path_basin = '"//basins//"'", "/", &
            "&ice_"//sx, "  backend = 'yelmo'", "  yelmo = 'yelmo'", "  init_equil_time = 0.0", &
            "  out_2D = 'out_2D'", "  out_1D = 'out_1D'", "/", &
            "&geo_"//sx, "  backend = 'vilma2'", "  out_2D = 'out_2D'", "  out_1D = 'out_1D'", "/", &
            "&climate_"//sx, "  backend = 'snapshots'", "  clock = 'forcing'", "  dt = 20.0", &
            "  snapshots = 'snapshots'", "  reference = 'pd'", terms, "  out_2D = 'out_2D'", "/", &
            "&surface_"//sx, "  backend = 'chion'", "  clock = 'forcing'", "  surface_chion = 'surface_chion'", &
            "  chion = 'chion'", "  orbit = 'orbit'", "  out_2D = 'out_2D'", "/", &
            "&shelf_"//sx, "  backend = 'marine_shelf'", "  marine_shelf = 'marine_shelf'", "  corr = 'shelf_corr'", &
            "  out_2D = 'out_2D'", "/", &
            "&ghf_"//sx, "  backend = 'data'", "  dt = -1.0", "  Q_geo = 'ghf_const'", "  out_2D = 'out_2D'", "/", &
            "&ghf_const", "  mode = 'const'", "  value = 0.05", "/", &
            "&sed_"//sx, "  backend = 'data'", "  dt = -1.0", "  H_sed = 'sed_none'", "/", &
            "&sed_none", "  mode = 'const'", "  value = 0.0", "/"])
        if (present(exps)) then
            select case (exps)
            case ("exp_recon")
                call w(u, [character(len=400) :: "&exp_recon", "  kind = 'recon'", &
                    "  recon_path = '"//trim(root)//"/{domain}/{grid_name}/{grid_name}_TOPO-ICE-6G_C.nc'", &
                    "  recon_var = 'dz'", "  recon_time = -20000.0", "  recon_regions = 'all'", "/"])
            case ("exp_relax")
                call w(u, [character(len=400) :: "&exp_relax", "  kind = 'relax'", "  topo_rel = 3", "  tau1 = 100.0", &
                    "  tau2 = 2000.0", "  time1 = 20.0", "  time2 = 60.0", "  m = 2.0", "/"])
            end select
        end if
        ! forcing coefficients (only what the domain uses: the converted Antarctica base has none)
        if (has_terms) call w(u, [character(len=400) :: "&coef_alpha", "  method = 'const'", "  f_min = 1.0", "/"])
        if (has_ocean) call w(u, [character(len=400) :: &
            "&coef_dTo", "  method = 'derived'", "  from = 'air_anomaly'", "  factor = 0.25", "  region = 'all'", &
            "  units = 'K'", "/"])
        ! climate: snapshots, terms, field groups
        call w(u, [character(len=400) :: &
            "&snapshots", "  f_p = "//ftoa(f_p), "/", &
            "&snapshot_pd", "  tas = 'v_t_pd'", "  pr = 'v_p_pd'", "  zs = 'v_zs'", "  to = 'v_to'", "  so = 'v_so'", "/"])
        if (has_terms) call w(u, [character(len=400) :: &
            "&snapshot_ctl", "  tas = 'v_t_pd'", "  pr = 'v_p_pd'", "  zs = 'v_zs'", "/", &
            "&snapshot_alt", "  tas = 'v_t_alt'", "  pr = 'v_p_alt'", "  zs = 'v_zs'", "/", &
            "&term_lgm", "  coef = 'alpha'", "  snap0 = 'ctl'", "  snap1 = 'alt'", "  fields = 'tas' 'pr'", "/"])
        if (has_ocean) call w(u, [character(len=400) :: &
            "&term_ocean", "  coef = 'dTo'", "  pattern = 'uniform'", "  fields = 'to'", "/"])
        call field(u, id, "v_t_pd", "t_pd", "K");   call field(u, id, "v_p_pd", "p_pd", "K")
        if (has_terms) then
            call field(u, id, "v_t_alt", "t_alt", "K"); call field(u, id, "v_p_alt", "p_alt", "K")
        end if
        call field(u, id, "v_zs", "zs", "m");       call field(u, id, "v_to", "to", "K")
        call field(u, id, "v_so", "so", "psu")
        ! Yelmo: library defaults, the reference data loads off, masks without named regions; calving "equil"
        ! as the gate configurations (a restart with the default vm-m16 calving is not exact)
        call w(u, [character(len=400) :: &
            "&yelmo", "  nz_aa = "//itoa(nz_aa), "  nml_ytopo = 'ytopo'", "  nml_ycalv = 'ycalv'", &
            "  nml_ydyn = 'ydyn'", "  nml_ytill = 'ytill'", "  nml_ymat = 'ymat'", "  nml_ytrc = 'ytrc'", &
            "  nml_ytherm = 'ytherm'", "  nml_yhyd = 'yhyd'", "  nml_masks = 'yelmo_masks'", &
            "  nml_init_topo = 'yelmo_init_topo'", "  nml_data = 'yelmo_data'", "/", "&ytopo", "/", &
            "&ycalv", "  calv_flt_method = 'equil'", "  calv_grnd_method = 'equil'", "/", "&ydyn", "  solver = 'diva'", "/", &
            "&ytill", "/", "&yhyd", "/", "&ymat", "/", "&ytrc", "/", "&ytherm", "/", &
            "&yelmo_masks", "  basins = 'None'", "  mask_ice_dynamic = 'all'", "  mask_rmse = 'all'", "/", &
            "&yelmo_init_topo", "/", &
            "&yelmo_data", "  pd_topo_load = False", "  pd_tsrf_load = False", "  pd_smb_load = False", &
            "  pd_vel_load = False", "  pd_age_load = False", "/"])
        ! chion and its host (a short spin-up), the orbit of 1950
        call w(u, [character(len=400) :: &
            "&orbit", "  method = 'fixed_orbit'", "/", &
            "&surface_chion", "  sf_a = 0.273", "  sf_b = 273.6", "  sigma_pdd = 5.0", "  time_equil = 2.0", &
            "  dt_days = 1", "  trans_sw = 0.6", "  wind_speed = 5.0", "  rel_hum = 0.8", "/", &
            "&chion", "  model = 'itm'", "  nml_itm = 'itm_chion'", "/", &
            "&itm_chion", "  H_snow_max = 5.0e3", "  Pmaxfrac = 0.6", "  trans_a = 0.46", "  trans_b = 6e-5", &
            "  trans_c = 0.01", "  itm_c = -45.0", "  itm_t = 10.0", "  itm_b = -2.0", "  itm_lat0 = 65.0", &
            "  alb_snow_dry = 0.8", "  alb_snow_wet = 0.65", "  alb_ice = 0.31", "  alb_land = 0.2", &
            "  alb_forest = 0.1", "  alb_ocean = 0.1", "  H_snow_crit_desert = 10.0", "  H_snow_crit_forest = 100.0", &
            "  melt_crit = 0.5", "/"])
        ! marine_shelf (masks without named regions)
        call w(u, [character(len=400) :: &
            "&marine_shelf", "  bmb_method = 'lin'", "  tf_method = 1", "  interp_depth = 'shlf'", &
            "  interp_method = 'interp'", "  find_ocean = True", "  restart = 'none'", "  extrap_shlf = False", &
            "  mask_ocean = 'all'", "  mask_deep_ocean = 'none'", "  mask_pico_deep = 'none'", &
            "  bmb_max = 0.0", "  c_deep = -10.0", &
            "  depth_deep = 500.0", "  depth_const = 0.0", "  depth_min = 1100.0", "  depth_max = 1500.0", &
            "  lambda1 = -0.0575", "  lambda2 = 0.0832", "  lambda3 = 7.59e-4", "  gamma_lin = 630.0", &
            "  gamma_quad = 11100.0", "  gamma_quad_nl = 14500.0", "  gamma_prime = 100.0", "  kappa_grz = 10.0", &
            "  c_grz = 1.0", "  f_grz_shlf = 10.0", "  grz_length = 32.0", "  use_obs = False", &
            "  obs_path = 'none'", "  obs_name = 'bm_eq_reese'", "  obs_scale = 1.0", "  obs_lim = 100.0", "/", &
            "&shelf_corr", "  basins = 'domain'", "  corr_method = 'none'", "  corr_names = ''", "  corr_values = 0.0", &
            "  tf_correction = False", "  tf_path = 'none'", "  tf_name = 'dT_nl'", "/"])
        close(u)
    end subroutine write_domain_par
    subroutine field(u, id, group, name, units)
        integer,          intent(in) :: u, id
        character(len=*), intent(in) :: group, name, units
        write(u,'(a)') "&"//group
        write(u,'(a)') "  filename = 'clim_"//trim(DSFX(id))//".nc'"
        write(u,'(a)') "  name = '"//name//"'"
        write(u,'(a)') "  units = '"//units//"' '"//units//"'"
        write(u,'(a)') "  scaling = 1.0 0.0"
        write(u,'(a)') "  time = 0 0 0 0 1"
        write(u,'(a)') "/"
    end subroutine field
    subroutine w(u, lines)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: lines(:)
        integer :: i
        do i = 1, size(lines)
            if (len_trim(lines(i)) == 0) cycle
            write(u,'(a)') trim(lines(i))
        end do
    end subroutine w

    ! =========================================================================
    ! Helpers of the inputs and processes
    ! =========================================================================

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    logical function child(sub, args, log, quiet)
        ! This executable with `args`, run in the run directory tmp/test_run_domains/<sub> on one thread (bitwise
        ! comparable runs); .true. on exit 0. A failure shows the end of the log unless quiet.
        character(len=*), intent(in) :: sub, args, log
        logical,          intent(in), optional :: quiet
        character(len=512) :: exe
        integer :: es
        logical :: q

        call get_command_argument(0, exe)
        if (exe(1:1) /= "/") exe = trim(repo)//"/"//trim(exe)
        call execute_command_line("cd '"//d//sub//"' && OMP_NUM_THREADS=1 '"//trim(exe)//"' "//args//" > "//log// &
                                  " 2>&1", exitstat=es)
        child = es == 0
        q = .false.
        if (present(quiet)) q = quiet
        if (.not. (child .or. q)) call sh("tail -n 30 '"//d//sub//"/"//log//"'")
    end function child

    logical function in_file(file, text)
        character(len=*), intent(in) :: file, text
        integer :: es
        call execute_command_line("grep -q '"//text//"' '"//file//"'", exitstat=es)
        in_file = es == 0
    end function in_file

    subroutine sh(cmd)
        character(len=*), intent(in) :: cmd
        call execute_command_line(cmd)
    end subroutine sh

    logical function exists(file)
        character(len=*), intent(in) :: file
        inquire(file=file, exist=exists)
    end function exists

    function ftoa(x) result(s)
        real(dp), intent(in) :: x
        character(len=24) :: s
        write(s, '(g0)') x
    end function ftoa

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=12) :: s
        write(s, '(i0)') i
    end function itoa

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        write(*,'(a,a)') merge("  ok    ", "  FAIL  ", ok), what
        if (.not. ok) nfail = nfail + 1
    end subroutine check


end program test_gia_sequence
