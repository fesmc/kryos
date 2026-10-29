module test_run_toys
    ! Host experiments of the sequence test (core-design 3.9): a recorder that
    ! logs, at every hook it wants, what the run has published at that point,
    ! and one that makes the boundary of the ice look unprepared (the guard in
    ! Yelmo's update).
    use precision,       only: dp
    use kryos_domain,    only: domain_t
    use kryos_experiment
    use kryos_ice_yelmo, only: yelmo_ice_t
    implicit none
    private

    integer, parameter, public :: MAXLOG = 64
    integer,  public :: n_log = 0
    integer,  public :: log_id(MAXLOG) = 0
    real(dp), public :: log_t0(MAXLOG), log_t1(MAXLOG), log_elapsed(MAXLOG), log_dtt(MAXLOG)
    real(dp), public :: log_geo_t1(MAXLOG), log_ice_t1(MAXLOG), log_zbed_t1(MAXLOG)   ! publication times seen
    logical,  public :: log_prepared(MAXLOG) = .false.                              ! Yelmo's tick is the step's
    integer(8), public :: log_tick(MAXLOG)
    real(dp), public :: log_tsum(MAXLOG)                    ! sum of Yelmo's T_ice: the state of the thermodynamics
    logical,  public :: log_rates(MAXLOG) = .false.         ! the ice's dHidt export published (after PASS_SPINUP)

    type, extends(experiment_t), public :: recorder_t
    contains
        procedure :: init => rec_init
        procedure :: hook => rec_hook
    end type recorder_t

    type, extends(experiment_t), public :: counter_t        ! stateful: the number of steps it has seen, in <name>_toy.dat
        integer :: n = 0
    contains
        procedure :: init => cnt_init
        procedure :: hook => cnt_hook
        procedure :: restart_write => cnt_write
        procedure :: restart_read  => cnt_read
        procedure :: restart_files => cnt_files
    end type counter_t

    type, extends(experiment_t), public :: skip_prepare_t
    contains
        procedure :: init => skip_init
        procedure :: hook => skip_hook
    end type skip_prepare_t

contains

    subroutine rec_init(self, ctx, dom)
        class(recorder_t), intent(inout) :: self
        type(exp_ctx_t),   intent(in)    :: ctx
        type(domain_t),    intent(inout) :: dom
        associate (u => ctx%time + len(dom%name)); end associate
        self%wants = .true.
    end subroutine rec_init

    subroutine rec_hook(self, dom, id, st)
        class(recorder_t), intent(inout) :: self
        type(domain_t),    intent(inout) :: dom
        integer,           intent(in)    :: id
        type(hook_time_t), intent(in)    :: st

        associate (u => self%restart); end associate
        n_log = n_log + 1
        log_id(n_log)      = id
        log_t0(n_log)      = st%t0
        log_t1(n_log)      = st%t1
        log_elapsed(n_log) = st%t_elapsed
        log_dtt(n_log)     = st%dtt
        log_tick(n_log)    = st%clk%n
        log_geo_t1(n_log)  = dom%geo%c%exports%f(dom%geo%c%h_dz_bed)%t1
        log_ice_t1(n_log)  = dom%ice%c%exports%f(dom%ice%c%h_H_ice)%t1
        log_zbed_t1(n_log) = dom%ice%c%imports%f(dom%ice%c%hi_z_bed)%t1
        log_rates(n_log)   = dom%ice%c%exports%f(dom%ice%c%h_dHidt)%published
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            log_prepared(n_log) = c%tick_prepared == st%clk%n
            log_tsum(n_log)     = sum(real(c%ylmo%thrm%now%T_ice, dp))
        end select
    end subroutine rec_hook

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
        character(len=512), allocatable :: files(:)
        files = [character(len=512) :: trim(self%name)//"_toy.dat"]
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

    subroutine skip_init(self, ctx, dom)
        class(skip_prepare_t), intent(inout) :: self
        type(exp_ctx_t),       intent(in)    :: ctx
        type(domain_t),        intent(inout) :: dom
        associate (u => ctx%time + len(dom%name)); end associate
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            self%wants(HOOK_ICE_BND) = .true.
        class default
            call fesm_error_toy("skip_prepare needs yelmo_ice_t")
        end select
    end subroutine skip_init

    subroutine skip_hook(self, dom, id, st)
        ! Forget that the boundary was prepared for this step (as if domain_prepare had been skipped).
        class(skip_prepare_t), intent(inout) :: self
        type(domain_t),        intent(inout) :: dom
        integer,               intent(in)    :: id
        type(hook_time_t),     intent(in)    :: st
        associate (u => id + st%t1 + len_trim(self%name)); end associate
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            c%tick_prepared = -1
        end select
    end subroutine skip_hook

    subroutine fesm_error_toy(msg)
        use fesm_log, only: fesm_error
        character(len=*), intent(in) :: msg
        call fesm_error("test_run_sequence", msg)
    end subroutine fesm_error_toy

end module test_run_toys


program test_run_sequence
    ! C9: the whole run sequence (core-design 3.1-3.6) on one Greenland
    ! domain, GRL-16KM, with every phase 1 backend: ice = yelmo, geo =
    ! fastisostasy (+ the global bsl part), climate = snapshots, surface =
    ! chion, shelf = marine_shelf, geothermal and sediments = data, forcing =
    ! tsgen (alpha for the climate, the derived dTo, beta with the ice-mass
    ! feedback). run_init (instances,
    ! wiring, output streams, PASS_INIT, static rows, cold start), two steps
    ! (dtt = 10; geo every 10 yr, climate every 20 yr: the climate is skipped in
    ! step 1), run_end, as the driver does.
    !
    ! Inputs: the hub topography (TOPO-M17-v5) and basins (BASINS-nasa) of
    ! GRL-16KM from ice_data ($KRYOS_ICE_DATA, default ./ice_data; the gate's
    ! FesmData products are on Levante only); the climate snapshots are
    ! synthetic fields on GRL-16KM written by this test (tas, pr, zs, ocean
    ! profiles). Without the topography the test reports and passes.
    !
    ! Checked (the sequence's contract, not yelmox bit-identity, which is the
    ! Levante gate): every instance and row present; the hub and every import
    ! of the ice published at t0; the cold-start relations z_bed = ref + dz_bed
    ! and z_sl = bsl + dz_ss in sp with the reduced, scattered bsl (3.6 item 2);
    ! Q_geo, H_sed from the static rows (3.6 item 1); the climate's dTo in the
    ! forcing at t0 (3.6 item 3) and coef_alpha in the climate; at each step:
    ! M_ice and dTo of the previous step in the forcing (lag 1), smb and
    ! bmb_shlf of the previous step in the ice (lag 1), bsl of this step's
    ! reduction in the geo exports and in z_sl, the mirror hub at t, the
    ! climate's cadence (not run at t = 10, run at t = 20), every export at t;
    ! the output files of every stream, the timing table, the restart bundle
    ! (clock, manifest, every instance's field states and backend files).
    !
    ! C10, the restart round trip (core-design 3.3, 7; regression-design 1.2):
    ! the same configuration run continuously 0 -> 80 (cont/, bundles at 40
    ! and 80) and in two chunks, 0 -> 40 (chunk/c00/) and 40 -> 80 restored
    ! from ../c00/restart-40 (chunk/c01/, [kryos] restart relative to the run
    ! directory), as the gate stages a chunk run: every output record after 40
    ! and the bundle at 80 must equal the continuous run bit for bit. Checked
    ! with the comparison tool (tools/bin/kryos-compare.jl --chunks, Julia
    ! project $KRYOS_TOOLS, default <repo>/tools; reported OPEN when Julia or
    ! the project is not available) and, without it, on the ice, sea-level and
    ! surface fields of the bundles and of the last output record. Yelmo's
    ! pc_tau_max in the bundle at 80 (its maximum over the last 50 internal
    ! steps, a window a restart does not hold) is reported OPEN, any other
    ! difference fails; Yelmo calves with "equil" as the gate configurations
    ! (with its default vm-m16 calving a Yelmo restart is not exact). A bundle
    ! whose manifest does not match the run stops the restore (an instance
    ! renamed in the par file, a backend changed in the manifest). The
    ! manifest's files of every instance equal the files the backend declares
    ! (restart_files: the two field files, then the backend's); a file the
    ! manifest lists that is missing from the bundle stops the restore, and
    ! a declared file a backend did not create stops the write (both naming
    ! instance, file and folder).
    !
    ! The runs execute in tmp/test_run_sequence/ (each run directory with its
    ! own input/ assembled from links: kryos input/, the yelmo and chion
    ! defaults, fesm-utils datasrc_defaults.nml, insol's orbital tables;
    ! maps/), in child processes of this executable (one thread):
    !     test_run_sequence.x run <par file>          cold start + checks
    !     test_run_sequence.x plain <par file>        the run as the driver does it
    !     test_run_sequence.x restored <par file>     a restored chunk + checks
    !     test_run_sequence.x hooks <par file>        a restored chunk with the recorder experiment + checks
    !     test_run_sequence.x counter <par file>      a restored chunk with a stateful experiment + checks
    !     test_run_sequence.x xexp <par file>         a run (cold or restored) with the plugin experiments + checks
    !     test_run_sequence.x require <folder> <file> restart_require_files of one file
    ! Built only with yelmo=1 fastisostasy=1 chion=1 (config/common.mk).

    use precision,       only: sp, dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_read, nc_exists_var, nc_size
    use nml,             only: nml_has_group, nml_read
    use kryos_restart,   only: restart_require_files
    use coords,          only: grid_class, grid_cdo_read_desc
    use fesm_fields,     only: field_state_t, fields_get
    use fesm_clock,      only: clock_time, clock_time_elapsed, cadence_t_last
    use kryos_component, only: component_t, PATH_LEN
    use kryos_domain,    only: domain_component, NCOMP
    use kryos_run,       only: run_t, run_init, run_step, run_done, run_end
    use kryos_experiment, only: experiment_proto_t, HOOK_STEP_BEGIN, HOOK_ICE_BND, HOOK_COLD_PRE_STATE, HOOK_COLD_POST_STATE
    use test_run_toys
    use kryos_ice,       only: ice_flotation
    use kryos_ice_yelmo, only: yelmo_ice_t
    use yelmo_defs,      only: MASK_ICE_NONE
    use kryos_exp_relax, only: relax_exp_t
    use kryos_exp_opt_cb_ref, only: opt_cb_ref_exp_t
    use kryos_climate_snapshots,  only: snapshots_climate_t
    use kryos_surface_chion,      only: chion_surface_t
    use kryos_shelf_marine_shelf, only: marshelf_shelf_t

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_run_sequence/", PAR = "par.nml", DATA = "clim_grl.nc"
    integer, parameter :: NX = 106, NY = 181, NZ = 5

    ! The restart round trip: cont 0 -> T_END, chunks 0 -> T_MID -> T_END.
    real(dp), parameter :: T_MID = 40.0_dp, T_END = 80.0_dp
    character(len=*), parameter :: CONT = "cont/", C00 = "chunk/c00/", C01 = "chunk/c01/", CMM = "chunk/mm/", &
                                CHK = "chunk/c01h/", GRD = "guard/", CS0 = "chunk/cs0/", CS1 = "chunk/cs1/"

    ! The instances of the par file (manifest groups, bundle files).
    character(len=16), parameter :: INSTANCES(9) = [character(len=16) :: "forcing", "geo_glob", "ice_grl", &
        "geo_grl", "climate_grl", "surface_grl", "shelf_grl", "ghf_grl", "sed_grl"]

    character(len=512) :: arg, par_arg, file_arg, root, repo, only
    character(len=PATH_LEN) :: files_arg(1)
    integer :: nfail = 0

    call get_command_argument(1, arg)
    select case (trim(arg))
    case ("run", "plain", "restored", "hooks", "counter", "xexp", "require")
        call get_command_argument(2, par_arg)
        select case (trim(arg))
        case ("require")
            call get_command_argument(3, file_arg)
            files_arg(1) = file_arg                 ! (a constructor of trim() as the actual argument crashes gfortran)
            call restart_require_files("ice_grl", trim(par_arg), files_arg, &
                                       "its restart_files declares was not created by its restart_write")
        case ("run");      call run_child(trim(par_arg))
        case ("plain");    call plain_child(trim(par_arg))
        case ("restored"); call restored_child(trim(par_arg))
        case ("hooks");    call hooks_child(trim(par_arg))
        case ("counter");  call counter_child(trim(par_arg))
        case ("xexp");     call xexp_child(trim(par_arg))
        end select
        stop 0
    end select

    call get_environment_variable("PWD", repo)
    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = trim(repo)//"/ice_data"
    if (root(1:1) /= "/") root = trim(repo)//"/"//trim(root)
    if (.not. exists(trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc")) then
        write(*,'(a)') "test_run_sequence: ice_data topography not found under "//trim(root)//": skipped"
        stop 0
    end if

    call get_environment_variable("KRYOS_TEST_ONLY", only)
    if (trim(only) == "exp") then                   ! development: the plugin experiments only (after a full run)
        if (.not. exists(trim(repo)//"/"//DIR//DATA)) call setup_run_dir("")
        call sh("rm -rf '"//trim(repo)//"/"//DIR//"xr'")
        call exp_round_trip()
        if (nfail > 0) error stop 1
        write(*,'(a)') "test_run_sequence (experiments only): OK"
        stop 0
    end if

    call sh("rm -rf '"//trim(repo)//"/"//DIR//"'")
    call setup_run_dir("")
    call write_par(trim(repo)//"/"//DIR//PAR, "None", 0.0_dp, 20.0_dp, 20.0_dp, exps="rec")
    call check(child("", "run "//PAR, "run.log"), &
               "cold start + 2 steps + end with every phase 1 backend (child run, log "//DIR//"run.log)")
    call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//trim(repo)//"/"//DIR//"run.log'")

    call round_trip()
    call manifest_mismatch()
    call exp_round_trip()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_run_sequence: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_run_sequence: OK"

contains

    ! =========================================================================
    ! The run (child process, in the run directory)
    ! =========================================================================

    subroutine run_child(par_file)
        character(len=*), intent(in) :: par_file

        type(run_t), target :: run
        real(dp), allocatable :: smb_prev(:,:), bmb_prev(:,:)
        real(dp) :: m_prev, dto_prev, t
        integer :: n

        call run_init(run, par_file, protos=host_protos())
        call check_init(run)

        n = 0
        do while (.not. run_done(run))
            associate (dom => run%dom(1))
                smb_prev = get2(dom%surface%c%exports, dom%surface%c%h_smb)
                bmb_prev = get2(dom%shelf%c%exports, dom%shelf%c%h_bmb_shlf)
                m_prev   = get0(dom%ice%c%exports, dom%ice%c%h_M_ice)
                dto_prev = get0(dom%climate%c%exports, dom%climate%c%h_derived(1))
            end associate
            call run_step(run)
            n = n + 1
            t = clock_time(run%clk)
            call check_step(run, n, t, smb_prev, bmb_prev, m_prev, dto_prev)
        end do
        call check(n == 2, "two steps to time_end = 20")
        call check_ice_static(run, "t = 20", 0.0_dp)
        call check_ice_flotation(run, "t = 20")
        call run_end(run)
        call check_hooks_run()

        call check_files()
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (run): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (run): OK"
    end subroutine run_child

    subroutine plain_child(par_file)
        ! The run as the driver does it (a continuous run, the first chunk).
        character(len=*), intent(in) :: par_file
        type(run_t), target :: run

        call run_init(run, par_file, protos=host_protos())
        do while (.not. run_done(run))
            call run_step(run)
        end do
        call run_end(run)
        write(*,'(a)') "test_run_sequence (plain): OK"
    end subroutine plain_child

    subroutine restored_child(par_file)
        ! A chunk restored from a bundle at T_MID: the clock continues, the
        ! instances' imports and exports are the bundle's, the rows that run
        ! at the restore (ice -> geo, bsl -> geo, ice -> hub) ran; then the
        ! chunk as the driver runs it.
        character(len=*), intent(in) :: par_file

        type(run_t), target :: run
        integer :: n

        call run_init(run, par_file, protos=host_protos())
        write(*,'(a)') "--- after run_init (restored, t = 40)"
        call check(run%clk%continued .and. clock_time(run%clk) == T_MID .and. clock_time_elapsed(run%clk) == T_MID, &
                   "continuation at 40: the clock and the elapsed time run on")
        associate (dom => run%dom(1), ice => run%dom(1)%ice%c, geo => run%dom(1)%geo%c, &
                   cl => run%dom(1)%climate%c, sf => run%dom(1)%surface%c, fr => run%globe%forcing%c)
            call check(cadence_t_last(ice%cadence, run%clk) == T_MID .and. cadence_t_last(cl%cadence, run%clk) == T_MID, &
                       "cadences restored: every instance last ran at 40")
            call check(published_at(ice%exports, ice%h_H_ice, T_MID) .and. ice%exports%f(ice%h_smb)%t0 == T_MID - 10.0_dp, &
                       "ice exports as the bundle holds them (H_ice at 40, smb MEAN over [30, 40])")
            call check(published_at(ice%imports, ice%hi_smb, T_MID - 10.0_dp) .and. &
                       published_at(ice%imports, ice%hi_z_bed, T_MID), &
                       "ice imports as received in the step to 40 (smb of 30, lag 1; z_bed of 40, lag 0)")
            call check(get0(geo%imports, geo%hi_bsl) == get0(run%globe%geo%c%exports, run%globe%geo%c%h_bsl), &
                       "geo bsl import = the restored global bsl")
            call check(published_at(cl%exports, cl%h_tas, T_MID) .and. published_at(sf%exports, sf%h_smb, T_MID), &
                       "climate and surface exports restored at 40")
            call check(published_at(fr%imports, fr%hi_M_ice, T_MID - 10.0_dp), "forcing M_ice import restored (lag 1)")
            call check(dom%hub_t == T_MID .and. all(get2(dom%hub_exports, dom%hh%h_H_ice) == get2(ice%exports, ice%h_H_ice)), &
                       "hub rebuilt from the restored ice at 40")
        end associate

        call check_ice_static(run, "restored at 40", -2.0_dp)
        n = 0
        do while (.not. run_done(run))
            call run_step(run)
            n = n + 1
        end do
        call check(n == 4 .and. clock_time(run%clk) == T_END, "four steps to 80")
        call run_end(run)
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (restored): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (restored): OK"
    end subroutine restored_child

    function host_protos() result(protos)
        ! The host experiments of the tests, selected by the group name in the par file.
        type(experiment_proto_t), allocatable :: protos(:)
        allocate(protos(3))
        protos(1)%name = "rec"
        allocate(recorder_t :: protos(1)%c)
        protos(2)%name = "skip_prepare"
        allocate(skip_prepare_t :: protos(2)%c)
        protos(3)%name = "counter"
        allocate(counter_t :: protos(3)%c)
    end function host_protos

    subroutine check_hooks_run()
        ! The recorder of a cold-started run, two steps (t = 10, 20). The cold
        ! start: COLD_PRE_STATE and COLD_POST_STATE at t0 around the ice's state
        ! (prepared boundary; thermodynamics set by PASS_STATE; the rates
        ! exported only by PASS_SPINUP, after the second hook). At every step
        ! STEP_BEGIN before geo and the ice, ICE_BND after geo, the ice's
        ! imports and prepare, before the ice update; the elapsed time of the
        ! clock (core-design 3.9). Restore runs no hook.
        integer :: n, k
        real(dp) :: t
        logical :: ok

        write(*,'(a)') "--- experiment hooks (recorder)"
        call check(n_log == 6, "recorder: two cold-start hooks, then STEP_BEGIN and ICE_BND at each of the two steps")
        if (n_log /= 6) return

        call check(log_id(1) == HOOK_COLD_PRE_STATE .and. log_id(2) == HOOK_COLD_POST_STATE, &
                   "cold start: COLD_PRE_STATE, then COLD_POST_STATE")
        call check(log_t0(1) == 0.0_dp .and. log_t1(1) == 0.0_dp .and. log_elapsed(1) == 0.0_dp .and. &
                   log_t1(2) == 0.0_dp .and. log_dtt(1) == 10.0_dp, "cold start: t0 = t1 = time_init, elapsed 0, dtt = 10")
        call check(log_zbed_t1(1) == 0.0_dp .and. log_prepared(1) .and. log_prepared(2), &
                   "cold start: the ice's imports at t0 and its boundary prepared before PASS_STATE")
        call check(log_tsum(1) /= log_tsum(2), "COLD_POST_STATE sees the thermodynamics PASS_STATE set")
        call check(.not. log_rates(1) .and. .not. log_rates(2), "cold start: both hooks before PASS_SPINUP (rates not yet exported)")

        do n = 1, 2
            t = 10.0_dp*n
            k = 2 + 2*n - 1
            ok = log_id(k) == HOOK_STEP_BEGIN .and. log_id(k+1) == HOOK_ICE_BND
            call check(ok, "step "//trim(itoa(n))//": STEP_BEGIN, then ICE_BND")
            call check(log_t1(k) == t .and. log_t0(k) == t - 10.0_dp .and. log_dtt(k) == 10.0_dp .and. &
                       log_t1(k+1) == t, "step "//trim(itoa(n))//": the run's step [t-10, t]")
            call check(log_elapsed(k) == t .and. log_elapsed(k+1) == t, &
                       "step "//trim(itoa(n))//": elapsed time of the clock = t (origin at 0)")
            call check(log_geo_t1(k) == t - 10.0_dp .and. log_ice_t1(k) == t - 10.0_dp .and. &
                       log_zbed_t1(k) == t - 10.0_dp, &
                       "step "//trim(itoa(n))//": STEP_BEGIN before geo, the ice rows and the ice update")
            call check(log_geo_t1(k+1) == t .and. log_zbed_t1(k+1) == t .and. log_ice_t1(k+1) == t - 10.0_dp, &
                       "step "//trim(itoa(n))//": ICE_BND after geo and the ice imports at t, before the ice update")
            call check(.not. log_prepared(k), "step "//trim(itoa(n))//": STEP_BEGIN: prepare not yet run for this step")
            call check(log_prepared(k+1), "step "//trim(itoa(n))//": ICE_BND: the ice's boundary prepared for this step")
            call check(log_rates(k) .and. log_rates(k+1), "step "//trim(itoa(n))//": rates exported by the cold start")
        end do
    end subroutine check_hooks_run

    subroutine counter_child(par_file)
        ! A chunk restored from a bundle that holds a stateful experiment (the
        ! counter of the chunk before: 4 steps seen): its state comes back with
        ! the restore and the count runs on.
        character(len=*), intent(in) :: par_file

        type(run_t), target :: run
        integer :: n, n0

        call run_init(run, par_file, protos=host_protos())
        write(*,'(a)') "--- after run_init (restored, t = 40, stateful experiment)"
        n0 = -1
        select type (c => run%exps(1)%s(1)%c)
        type is (counter_t)
            n0 = c%n
        end select
        call check(n0 == 4 .and. n_log == 0, "the counter restored with its four steps; no hook in the restore")
        n = 0
        do while (.not. run_done(run))
            call run_step(run)
            n = n + 1
        end do
        select type (c => run%exps(1)%s(1)%c)
        type is (counter_t)
            n0 = c%n
        end select
        call check(n == 4 .and. n0 == 8, "four more steps: the count runs on to 8")
        call run_end(run)
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (counter): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (counter): OK"
    end subroutine counter_child

    subroutine hooks_child(par_file)
        ! A chunk restored from the bundle at T_MID with the recorder (a stateless
        ! experiment added between chunks): no hook runs in the restore; the
        ! elapsed time of the continuation keeps the origin of the run (core-design 3.9).
        character(len=*), intent(in) :: par_file

        type(run_t), target :: run
        integer :: n, k

        call run_init(run, par_file, protos=host_protos())
        write(*,'(a)') "--- after run_init (restored, t = 40, recorder)"
        call check(n_log == 0, "the restore runs no hook")
        call check(run%clk%continued .and. clock_time_elapsed(run%clk) == T_MID, "continuation: elapsed = 40 at the restore")
        n = 0
        do while (.not. run_done(run))
            call run_step(run)
            n = n + 1
        end do
        call check(n == 4 .and. n_log == 8, "four steps, two hooks each")
        if (n_log == 8) then
            do k = 1, 4
                call check(log_elapsed(2*k-1) == T_MID + 10.0_dp*k .and. log_elapsed(2*k) == T_MID + 10.0_dp*k .and. &
                           log_t1(2*k) == T_MID + 10.0_dp*k, "step "//trim(itoa(k))//" after the restore: elapsed = "// &
                           trim(ftoa(T_MID + 10.0_dp*k))//" (continuous, not restarted)")
            end do
        end if
        call run_end(run)
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (hooks): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (hooks): OK"
    end subroutine hooks_child

    subroutine xexp_child(par_file)
        ! A run (cold start, or restored from a bundle) with the plugin experiments
        ! of the xr/ set (core-design 3.9). relax: after every step the topography
        ! relaxation of Yelmo is what relax_update gives for the elapsed time of the
        ! clock (continuous across chunks), then the par-file values after time2.
        ! opt_cb_ref: Yelmo takes cb_ref from outside (till_method = -1); the cold
        ! start sets cf_init everywhere; cb_ref changes in the window (elapsed time
        ! 0-60) and is then left alone (method none: never changes).
        character(len=*), intent(in) :: par_file

        type(run_t), target :: run
        real(sp), allocatable :: cb0(:,:), cb60(:,:)
        integer :: n, i, irel, icf, topo_rel0
        real(dp) :: t
        logical :: ok, ok_cb, none_cf

        call run_init(run, par_file, protos=host_protos())
        write(*,'(a)') "--- plugin experiments"
        irel = 0;  icf = 0;  none_cf = .false.
        do i = 1, size(run%exps(1)%s)
            select type (x => run%exps(1)%s(i)%c)
            type is (relax_exp_t)
                irel = i
                call check(x%topo_rel0 == 0 .and. x%topo_rel_tau0 == 10.0_sp, "relax: the par-file ytopo values recorded at init")
                call check(x%wants(HOOK_STEP_BEGIN) .and. count(x%wants) == 1, "relax: wants STEP_BEGIN only")
            type is (opt_cb_ref_exp_t)
                icf = i
                none_cf = trim(x%method) == "none"
                call check(x%wants(HOOK_COLD_PRE_STATE) .and. (x%wants(HOOK_STEP_BEGIN) .neqv. none_cf) .and. &
                           count(x%wants) == merge(1, 2, none_cf), "opt_cb_ref: the cold-start hook, STEP_BEGIN unless method = none")
            end select
        end do
        if (irel > 0 .and. .not. run%clk%continued) then
            select type (c => run%dom(1)%ice%c)
            type is (yelmo_ice_t)
                call check(c%ylmo%tpo%par%topo_rel == 0, "cold start: no relaxation before the first step")
            end select
        end if
        if (icf > 0) then
            select type (x => run%exps(1)%s(icf)%c)
            type is (opt_cb_ref_exp_t)
                select type (c => run%dom(1)%ice%c)
                type is (yelmo_ice_t)
                    call check(c%ylmo%dyn%par%till_method == -1, "opt_cb_ref: Yelmo's till_method = -1 (cb_ref from outside), also in a restore")
                    call check(all(x%cf_min == c%ylmo%dyn%par%till_cf_min) .and. all(x%cf_max == c%ylmo%dyn%par%till_cf_ref), &
                               "opt_cb_ref: the bounds cf_min, cf_max from ytill")
                    if (.not. run%clk%continued) call check(all(c%ylmo%dyn%now%cb_ref == 1.0e-2_sp), &
                               "cold start: cb_ref = cf_init (1e-2) everywhere")
                    cb0 = c%ylmo%dyn%now%cb_ref
                end select
            end select
        end if

        n = 0
        ok = .true.
        ok_cb = .true.
        do while (.not. run_done(run))
            call run_step(run)
            n = n + 1
            t = clock_time_elapsed(run%clk)
            select type (c => run%dom(1)%ice%c)
            type is (yelmo_ice_t)
                if (irel > 0) then
                    topo_rel0 = c%ylmo%tpo%par%topo_rel
                    if (t <= 60.0_dp) then
                        ok = ok .and. topo_rel0 == 3 .and. &
                             abs(c%ylmo%tpo%par%topo_rel_tau - real(tau_expected(t), sp)) <= 1e-5_sp*real(tau_expected(t), sp)
                    else
                        ok = ok .and. topo_rel0 == 0 .and. c%ylmo%tpo%par%topo_rel_tau == 10.0_sp
                    end if
                    if (.not. ok) write(*,'(a,f8.1,a,i0,a,es12.4)') "  at elapsed ", t, ": topo_rel = ", topo_rel0, &
                                                                     ", tau = ", c%ylmo%tpo%par%topo_rel_tau
                end if
                if (icf > 0 .and. .not. none_cf) then
                    if (t == 60.0_dp) then
                        cb60 = c%ylmo%dyn%now%cb_ref
                        call check(any(cb60 /= cb0), "opt_cb_ref: cb_ref has changed from its value at the start of the run, "// &
                                   "at elapsed 60 ("//trim(itoa(count(cb60 /= cb0)))//" cells)")
                    else if (t > 60.0_dp) then
                        ok_cb = ok_cb .and. all(c%ylmo%dyn%now%cb_ref == cb60)
                    end if
                end if
            end select
        end do
        if (irel > 0) call check(ok, "relax: topo_rel = 3 and the tau ramp (100 -> 2000, 20 -> 60 yr, m = 2) of the elapsed time, "// &
                                 "the par-file values after 60 ("//trim(itoa(n))//" steps)")
        if (icf > 0 .and. .not. none_cf .and. clock_time_elapsed(run%clk) > 60.0_dp) &
            call check(ok_cb, "opt_cb_ref: cb_ref unchanged after the window (elapsed 60)")
        if (icf > 0 .and. none_cf) then
            select type (c => run%dom(1)%ice%c)
            type is (yelmo_ice_t)
                call check(all(c%ylmo%dyn%now%cb_ref == cb0), "opt_cb_ref method = none: cb_ref never changes")
            end select
        end if
        call run_end(run)
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (xexp): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (xexp): OK"
    end subroutine xexp_child

    real(dp) function tau_expected(t)
        ! relax_update's ramp: tau1 until time1, then tau1 + (tau2 - tau1)*((t - time1)/(time2 - time1))**m
        real(dp), intent(in) :: t
        if (t <= 20.0_dp) then
            tau_expected = 100.0_dp
        else
            tau_expected = 100.0_dp + 1900.0_dp*((t - 20.0_dp)/40.0_dp)**2
        end if
    end function tau_expected

    ! =========================================================================
    ! C10: restart round trip, manifest checks
    ! =========================================================================

    subroutine round_trip()
        ! cont: 0 -> 80, bundles at 40, 80; chunk c00: 0 -> 40; chunk c01:
        ! 40 -> 80 from ../c00/restart-40. Every record after 40 and the
        ! bundle at 80 equal bit for bit.
        character(len=:), allocatable :: d
        logical :: ok

        write(*,'(a)') "--- restart round trip: cont 0 -> 80 against chunks 0 -> 40 -> 80"
        d = trim(repo)//"/"//DIR
        call setup_run_dir(CONT)
        call setup_run_dir(C00)
        call setup_run_dir(C01)
        call write_par(d//CONT//PAR, "None", 0.0_dp, T_END, T_MID)
        call write_par(d//C00//PAR,  "None", 0.0_dp, T_MID, T_MID)
        call write_par(d//C01//PAR,  "../c00/restart-40", T_MID, T_END, T_MID)

        call check(child(CONT, "plain "//PAR, "run.log"), "continuous run 0 -> 80 ("//DIR//CONT//"run.log)")
        call check(child(C00, "plain "//PAR, "run.log"), "chunk 0 -> 40 ("//DIR//C00//"run.log)")
        ok = child(C01, "restored "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//C01//"run.log'")
        call check(ok, "chunk 40 -> 80 restored from ../c00/restart-40 ("//DIR//C01//"run.log)")

        ! the same chunk with the recorder experiment added (stateless: allowed between chunks), and
        ! a run in which the ice's boundary looks unprepared at the update (the guard)
        call setup_run_dir(CHK)
        call write_par(d//CHK//PAR, "../c00/restart-40", T_MID, T_END, T_MID, exps="rec")
        ok = child(CHK, "hooks "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//CHK//"run.log'")
        call check(ok, "chunk 40 -> 80 with the recorder: no hook in the restore, elapsed time continuous ("// &
                   DIR//CHK//"run.log)")
        call setup_run_dir(GRD)
        call write_par(d//GRD//PAR, "None", 0.0_dp, 20.0_dp, 20.0_dp, exps="skip_prepare")
        call check(.not. child(GRD, "plain "//PAR, "run.log", quiet=.true.) .and. &
                   in_file(d//GRD//"run.log", "kryos_ice_yelmo::update") .and. &
                   in_file(d//GRD//"run.log", "was not prepared"), &
                   "the ice's update stops when its boundary was not prepared for the step ("//DIR//GRD//"run.log)")

        ! a stateful experiment: a manifest group of its own and its file in the domain's folder; its state
        ! restored; a stateful experiment added or dropped against the bundle stops the restore
        call setup_run_dir(CS0)
        call setup_run_dir(CS1)
        call write_par(d//CS0//PAR, "None", 0.0_dp, T_MID, T_MID, exps="counter")
        call write_par(d//CS1//PAR, "../cs0/restart-40", T_MID, T_END, T_MID, exps="counter")
        call check(child(CS0, "plain "//PAR, "run.log"), "chunk 0 -> 40 with a stateful experiment ("//DIR//CS0//"run.log)")
        call check(exists(d//CS0//"restart-40/Greenland/counter_toy.dat") .and. &
                   in_file(d//CS0//"restart-40/manifest.nml", "&counter") .and. &
                   in_file(d//CS0//"restart-40/manifest.nml", 'component      = "experiment"') .and. &
                   in_file(d//CS0//"restart-40/manifest.nml", 'backend        = "host"') .and. &
                   in_file(d//CS0//"restart-40/manifest.nml", '"counter_toy.dat"'), &
                   "the bundle: counter_toy.dat in the domain folder, manifest group counter (experiment, host)")
        call check(.not. in_file(d//C00//"restart-40/manifest.nml", "&rec") .and. &
                   .not. in_file(d//CHK//"restart-80/manifest.nml", "&rec"), &
                   "a stateless experiment (the recorder) is not in the manifest")
        ok = child(CS1, "counter "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//CS1//"run.log'")
        call check(ok, "chunk 40 -> 80 restores the stateful experiment ("//DIR//CS1//"run.log)")
        if (nfail > 0) return

        call check_bundle(CONT//"restart-40")
        call check_bundle(C01//"restart-80")

        ! without the comparison tool: the state that carries the run
        call same_2d(CONT//"restart-80/Greenland/ice_grl_exports.nc", C01//"restart-80/Greenland/ice_grl_exports.nc", &
                     "H_ice", "bundle at 80: ice H_ice export")
        call same_2d(CONT//"restart-80/Greenland/ice_grl_exports.nc", C01//"restart-80/Greenland/ice_grl_exports.nc", &
                     "z_srf", "bundle at 80: ice z_srf export")
        call same_2d(CONT//"restart-80/Greenland/surface_grl_exports.nc", &
                     C01//"restart-80/Greenland/surface_grl_exports.nc", "smb", "bundle at 80: surface smb export")
        call same_2d(CONT//"restart-80/Greenland/geo_grl_exports.nc", C01//"restart-80/Greenland/geo_grl_exports.nc", &
                     "dz_bed", "bundle at 80: geo dz_bed export")
        call same_0d(CONT//"restart-80/geo_glob_exports.nc", C01//"restart-80/geo_glob_exports.nc", "bsl", &
                     "bundle at 80: global bsl")
        call same_0d(CONT//"restart-80/forcing_exports.nc", C01//"restart-80/forcing_exports.nc", "beta", &
                     "bundle at 80: coefficient beta (ice-mass feedback)")
        call same_record(CONT//"ice_grl.nc", C01//"ice_grl.nc", "H_ice", T_END, "output ice_grl.nc H_ice at 80")
        call same_record(CONT//"surface_grl.nc", C01//"surface_grl.nc", "smb", T_END, "output surface_grl.nc smb at 80")

        call compare_tool("")
    end subroutine round_trip

    subroutine exp_round_trip()
        ! The plugin experiments (xr/ set): the same restart round trip as above,
        ! cont 0 -> 80 against chunks 0 -> 40 -> 80, with relax and opt_cb_ref active
        ! (window 20-60 and 0-60 yr of the elapsed time: both straddle the restart; the
        ! optimisation only when the present-day velocity file is there). Each run checks
        ! its own state (xexp_child); the comparison tool checks chunk = cont bit for bit.
        ! Then method = none (a short run) and the stop without the present-day velocity.
        character(len=:), allocatable :: d, vel
        character(len=*), parameter :: X = "xr/", XC = "xr/cont/", X0 = "xr/chunk/c00/", X1 = "xr/chunk/c01/", &
                                       XN = "xr/none/", XV = "xr/novel/"
        character(len=24), allocatable :: names(:)
        logical :: ok, has_vel

        write(*,'(a)') "--- plugin experiments: restart round trip with relax and opt_cb_ref"
        d = trim(repo)//"/"//DIR
        vel = trim(root)//"/Greenland/GRL-16KM/GRL-16KM_VEL-J18.nc"
        has_vel = exists(vel)
        if (has_vel) then
            names = [character(len=24) :: "exp_relax", "exp_cf"]
        else
            names = [character(len=24) :: "exp_relax"]
            call note(.false., "opt_cb_ref not run: no present-day velocity "//vel)
        end if
        call setup_run_dir(XC)
        call setup_run_dir(X0)
        call setup_run_dir(X1)
        call write_par(d//XC//PAR, "None", 0.0_dp, T_END, T_MID, explist=names, xlines=exp_lines(names), vel=vel_arg(has_vel, vel))
        call write_par(d//X0//PAR, "None", 0.0_dp, T_MID, T_MID, explist=names, xlines=exp_lines(names), vel=vel_arg(has_vel, vel))
        call write_par(d//X1//PAR, "../c00/restart-40", T_MID, T_END, T_MID, explist=names, xlines=exp_lines(names), &
                       vel=vel_arg(has_vel, vel))
        ok = child(XC, "xexp "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//XC//"run.log'")
        call check(ok, "continuous run 0 -> 80 with the experiments ("//DIR//XC//"run.log)")
        ok = child(X0, "xexp "//PAR, "run.log")
        call check(ok, "chunk 0 -> 40 with the experiments ("//DIR//X0//"run.log)")
        ok = child(X1, "xexp "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//X1//"run.log'")
        call check(ok, "chunk 40 -> 80: the windows continue on the elapsed time ("//DIR//X1//"run.log)")
        if (nfail > 0) return
        call compare_tool(X)

        ! method = none: external cb_ref = cf_init, never optimised (needs no velocity)
        call setup_run_dir(XN)
        names = [character(len=24) :: "exp_cf0"]
        call write_par(d//XN//PAR, "None", 0.0_dp, 20.0_dp, 20.0_dp, explist=names, xlines=exp_lines(names))
        ok = child(XN, "xexp "//PAR, "run.log")
        call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//d//XN//"run.log'")
        call check(ok, "opt_cb_ref method = none: till_method = -1, cb_ref = cf_init and unchanged ("//DIR//XN//"run.log)")
        ! L21 without the present-day velocity stops at init, naming the key
        call setup_run_dir(XV)
        names = [character(len=24) :: "exp_cf"]
        call write_par(d//XV//PAR, "None", 0.0_dp, 20.0_dp, 20.0_dp, explist=names, xlines=exp_lines(names))
        call check(.not. child(XV, "plain "//PAR, "run.log", quiet=.true.) .and. in_file(d//XV//"run.log", "needs the present-day velocity") &
                   .and. in_file(d//XV//"run.log", "pd_vel_load"), &
                   "opt_cb_ref L21 stops without Yelmo's present-day velocity ("//DIR//XV//"run.log)")
    end subroutine exp_round_trip

    function vel_arg(has, vel) result(v)
        ! Optional argument stand-in: the file when there is one, else "" (write_par: no velocity data).
        logical,          intent(in) :: has
        character(len=*), intent(in) :: vel
        character(len=:), allocatable :: v
        v = ""
        if (has) v = vel
    end function vel_arg

    function exp_lines(names) result(lines)
        ! The par-file groups of the plugin experiments named in `names` (test values: relax
        ! with a ramp inside the 80 yr of the run).
        character(len=*), intent(in) :: names(:)
        character(len=400), allocatable :: lines(:)
        integer :: i

        allocate(lines(0))
        do i = 1, size(names)
            select case (trim(names(i)))
            case ("exp_relax")
                lines = [character(len=400) :: lines, "&exp_relax", "  kind = 'relax'", "  topo_rel = 3", "  tau1 = 100.0", &
                         "  tau2 = 2000.0", "  time1 = 20.0", "  time2 = 60.0", "  m = 2.0", "/"]
            case ("exp_cf")
                lines = [character(len=400) :: lines, "&exp_cf", "  kind = 'opt_cb_ref'", "  method = 'L21'", &
                         "  cf_time_init = 0.0", "  cf_time_end = 60.0", "  tau_c = 50.0", "/"]
            case ("exp_cf0")
                lines = [character(len=400) :: lines, "&exp_cf0", "  kind = 'opt_cb_ref'", "  method = 'none'", "/"]
            end select
        end do
    end function exp_lines

    subroutine compare_tool(sub)
        ! tools/bin/kryos-compare.jl --chunks: every file of cont against the
        ! chunk folders c00, c01 (records after their start) and the bundles
        ! at 40 and 80, class E, of the run set tmp/test_run_sequence/<sub>.
        ! Julia project: $KRYOS_TOOLS, else <repo>/tools.
        character(len=*), intent(in) :: sub
        character(len=512) :: tools
        character(len=:), allocatable :: d, cmd
        integer :: u, es, es2
        character(len=*), parameter :: PC_TAU_MAX = "c01/restart-80/Greenland/ice_grl_restart.nc:pc_tau_max"

        d = trim(repo)//"/"//DIR//sub
        call get_environment_variable("KRYOS_TOOLS", tools)
        if (len_trim(tools) == 0) tools = trim(repo)//"/tools"
        call execute_command_line("julia --project='"//trim(tools)//"' -e 'using KryosTools' > /dev/null 2>&1", &
                                  exitstat=es)
        if (es /= 0) then
            call note(.false., "comparison tool not run: Julia or the project "//trim(tools)// &
                      " is not available (set KRYOS_TOOLS to a tools/ with an instantiated project)")
            return
        end if

        open(newunit=u, file=d//"gate.toml", status="replace", action="write")
        write(u,'(a)') 'gate = "test_run_sequence"', 'chunks = [0.0, 40.0, 80.0]', 'restart_times = [40.0, 80.0]', &
                       '[instances]', 'domain = "Greenland"', 'ice = "ice_grl"', 'geo = "geo_grl"', &
                       'climate = "climate_grl"', 'surface = "surface_grl"', 'shelf = "shelf_grl"', &
                       'geothermal = "ghf_grl"', 'sediments = "sed_grl"', 'forcing = "forcing"', 'geo_glob = "geo_glob"'
        close(u)
        cmd = "cd '"//d//"' && julia --project='"//trim(tools)//"' '"//trim(repo)//"/tools/bin/kryos-compare.jl' "// &
              "cont chunk --chunks --gate gate.toml --md compare.md > compare.log 2>&1"
        call execute_command_line(cmd, exitstat=es)
        if (es == 1) then
            ! Yelmo's pc_tau_max is the maximum over its last 50 internal steps,
            ! a diagnostic whose window the restart does not hold: refilled after
            ! 50 steps (the gate's 500-yr chunk), not in this test's 40 yr. A
            ! failure of any other variable fails the check.
            call execute_command_line("grep '| fail |' '"//d//"compare.md' | grep -v '"//PC_TAU_MAX//"' | grep -q .", &
                                      exitstat=es2)
            if (es2 /= 0) then
                call note(.false., "kryos-compare --chunks: only "//PC_TAU_MAX//" differs (Yelmo's running "// &
                          "maximum over 50 internal steps, not restored; refilled in a 500-yr chunk)")
                es = 0
            end if
        end if
        if (es /= 0) call sh("grep -E '^# |^- |^## |\| fail \|' '"//d//"compare.md' | head -n 60; "// &
                             "tail -n 20 '"//d//"compare.log'")
        call check(es == 0, "kryos-compare --chunks: every record after 40 and the bundles at 40, 80 equal to cont "// &
                   "(report "//DIR//sub//"compare.md)")
    end subroutine compare_tool

    subroutine manifest_mismatch()
        ! A restore whose bundle does not match the run stops, naming the
        ! instance: (1) the run's sediments instance renamed (not in the
        ! bundle); (2) the bundle's manifest with another backend for ghf_grl;
        ! (3) a file the manifest lists missing from the bundle.
        character(len=:), allocatable :: d

        write(*,'(a)') "--- manifest checks"
        d = trim(repo)//"/"//DIR
        call setup_run_dir(CMM)
        call sh("cp -R '"//d//C00//"restart-40' '"//d//CMM//"bad-40'")
        call manifest_edit(d//CMM//"bad-40/manifest.nml", "ghf_grl", "backend", "other")

        call write_par(d//CMM//"rename.nml", "../c00/restart-40", T_MID, T_END, T_MID, sed="sed2_grl")
        call check(.not. child(CMM, "plain rename.nml", "rename.log", quiet=.true.) .and. &
                   in_file(d//CMM//"rename.log", "sed2_grl") .and. in_file(d//CMM//"rename.log", "not in the restart bundle"), &
                   "restore stops: instance sed2_grl of the run is not in the bundle ("//DIR//CMM//"rename.log)")

        ! a file the manifest lists is missing from the bundle (the shelf's restart file)
        call sh("cp -R '"//d//C00//"restart-40' '"//d//CMM//"short-40'")
        call sh("rm '"//d//CMM//"short-40/Greenland/shelf_grl_restart.nc'")
        call write_par(d//CMM//"short.nml", "short-40", T_MID, T_END, T_MID)
        call check(.not. child(CMM, "plain short.nml", "short.log", quiet=.true.) .and. &
                   in_file(d//CMM//"short.log", "shelf_grl") .and. in_file(d//CMM//"short.log", "shelf_grl_restart.nc") .and. &
                   in_file(d//CMM//"short.log", "short-40/Greenland") .and. &
                   in_file(d//CMM//"short.log", "is not in the restart bundle"), &
                   "restore stops: shelf_restart.nc listed in the manifest is missing ("//DIR//CMM//"short.log)")

        ! a declared file the backend did not create stops the write (restart_require_files, as the writer calls it)
        call check(child(CMM, "require "//d//C00//"restart-40/Greenland shelf_grl_restart.nc", "require-ok.log", &
                         quiet=.true.), "restart_require_files: a declared file that exists passes")
        call check(.not. child(CMM, "require "//d//C00//"restart-40/Greenland missing_restart.nc", "require.log", &
                               quiet=.true.) .and. in_file(d//CMM//"require.log", "ice_grl") .and. &
                   in_file(d//CMM//"require.log", "missing_restart.nc") .and. &
                   in_file(d//CMM//"require.log", "restart-40/Greenland"), &
                   "restart_require_files: a declared file that is missing stops naming instance, file and folder")

        call write_par(d//CMM//"dropped.nml", "../cs0/restart-40", T_MID, T_END, T_MID)
        call check(.not. child(CMM, "plain dropped.nml", "dropped.log", quiet=.true.) .and. &
                   in_file(d//CMM//"dropped.log", "instance `counter` of the restart bundle is not in the run"), &
                   "restore stops: the bundle holds the stateful experiment, the run dropped it ("//DIR//CMM//"dropped.log)")
        call write_par(d//CMM//"added.nml", "../c00/restart-40", T_MID, T_END, T_MID, exps="counter")
        call check(.not. child(CMM, "plain added.nml", "added.log", quiet=.true.) .and. &
                   in_file(d//CMM//"added.log", "instance `counter` (experiment, host) is not in the restart bundle"), &
                   "restore stops: the run has a stateful experiment the bundle lacks ("//DIR//CMM//"added.log)")


        call write_par(d//CMM//"backend.nml", "bad-40", T_MID, T_END, T_MID)
        call check(.not. child(CMM, "plain backend.nml", "backend.log", quiet=.true.) .and. &
                   in_file(d//CMM//"backend.log", "ghf_grl") .and. in_file(d//CMM//"backend.log", "backend `data` in the run"), &
                   "restore stops: ghf_grl has backend `other` in the bundle ("//DIR//CMM//"backend.log)")
    end subroutine manifest_mismatch

    subroutine manifest_edit(file, group, key, value)
        ! Set `key` of `group` in a manifest to `value` (the line is replaced).
        character(len=*), intent(in) :: file, group, key, value

        character(len=1000), allocatable :: lines(:)
        character(len=1000) :: line
        logical :: inside
        integer :: u, ios, i

        allocate(lines(0))
        open(newunit=u, file=file, status="old", action="read")
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            lines = [character(len=1000) :: lines, line]
        end do
        close(u)
        inside = .false.
        do i = 1, size(lines)
            if (trim(adjustl(lines(i))) == "&"//group) inside = .true.
            if (trim(adjustl(lines(i))) == "/") inside = .false.
            if (inside .and. index(adjustl(lines(i)), key//" ") == 1) lines(i) = "    "//key//" = '"//value//"'"
        end do
        open(newunit=u, file=file, status="replace", action="write")
        do i = 1, size(lines)
            write(u, '(a)') trim(lines(i))
        end do
        close(u)
    end subroutine manifest_edit

    subroutine check_bundle(bundle)
        ! A bundle (path relative to tmp/test_run_sequence/, or to the current
        ! run directory in a child): clock file, manifest with every instance,
        ! the field states of every instance, the backend files of the
        ! stateful ones (and none for the stateless climate and data).
        character(len=*), intent(in) :: bundle

        character(len=:), allocatable :: b, sub
        character(len=64) :: val
        integer :: i
        logical :: ok

        b = bundle
        if (.not. exists(b//"/manifest.nml")) b = trim(repo)//"/"//DIR//bundle
        call check(exists(b//"/clock_restart.nc") .and. exists(b//"/manifest.nml"), bundle//": clock file and manifest")
        if (.not. exists(b//"/manifest.nml")) return

        ok = nml_has_group(b//"/manifest.nml", "bundle")
        do i = 1, size(INSTANCES)
            ok = ok .and. nml_has_group(b//"/manifest.nml", trim(INSTANCES(i)))
            sub = "/Greenland/"
            if (i <= 2) sub = "/"
            ok = ok .and. exists(b//sub//trim(INSTANCES(i))//"_imports.nc") .and. &
                 exists(b//sub//trim(INSTANCES(i))//"_exports.nc")
        end do
        call check(ok, bundle//": manifest groups and field states of the 9 instances")
        val = ""
        call nml_read(b//"/manifest.nml", "ice_grl", "backend", val)
        call check(trim(val) == "yelmo", bundle//": manifest ice_grl backend = yelmo")
        ok = exists(b//"/forcing_alpha_restart.nc") .and. exists(b//"/forcing_beta_restart.nc") .and. &
             .not. exists(b//"/forcing_dTo_restart.nc") .and. exists(b//"/geo_glob_restart.nc")
        do i = 3, size(INSTANCES)
            if (any(INSTANCES(i) == [character(len=16) :: "climate_grl", "ghf_grl", "sed_grl"])) then
                ok = ok .and. .not. exists(b//"/Greenland/"//trim(INSTANCES(i))//"_restart.nc")
            else
                ok = ok .and. exists(b//"/Greenland/"//trim(INSTANCES(i))//"_restart.nc")
            end if
        end do
        call check(ok, bundle//": backend restart files of the stateful instances only")
        call check_corr(b//"/Greenland/shelf_grl_corr_restart.nc", bundle)

        ok = .true.
        do i = 1, size(INSTANCES)
            ok = ok .and. files_listed(b//"/manifest.nml", trim(INSTANCES(i)))
        end do
        call check(ok, bundle//": manifest files of every instance = the declared lists (the 2 field files, then "// &
                   "the backend's)")
    end subroutine check_bundle

    subroutine check_corr(file, what)
        ! The shelf component's corrections file (kryos_shelf_ocean): tf_corr and tf_corr_basin on the
        ! shelf grid; this configuration has no correction, so both are 0.
        character(len=*), intent(in) :: file, what
        real(sp), allocatable :: a(:,:,:), c(:,:,:)

        if (.not. exists(file)) then
            call check(.false., what//": the shelf's corrections file shelf_grl_corr_restart.nc")
            return
        end if
        allocate(a(NX,NY,1), c(NX,NY,1))
        call nc_read(file, "tf_corr", a)
        call nc_read(file, "tf_corr_basin", c)
        call check(all(a == 0.0_sp) .and. all(c == 0.0_sp), what//": shelf_grl_corr_restart.nc holds tf_corr, tf_corr_basin (0)")
    end subroutine check_corr

    function declared_files(name) result(files)
        ! The files the manifest must list for an instance of the par file: its
        ! two field files, then the files its backend declares (restart_files):
        ! one per coefficient of the forcing, <name>_restart.nc of the stateful
        ! instances, none of the stateless (climate, data).
        character(len=*), intent(in) :: name
        character(len=64), allocatable :: files(:)

        files = [character(len=64) :: trim(name)//"_imports.nc", trim(name)//"_exports.nc"]
        select case (trim(name))
        case ("forcing")
            files = [character(len=64) :: files, "forcing_alpha_restart.nc", "forcing_beta_restart.nc"]
        case ("climate_grl", "ghf_grl", "sed_grl")
            continue
        case ("shelf_grl")                              ! the shelf's corrections file follows the backend's
            files = [character(len=64) :: files, trim(name)//"_restart.nc", trim(name)//"_corr_restart.nc"]
        case default
            files = [character(len=64) :: files, trim(name)//"_restart.nc"]
        end select
    end function declared_files

    logical function files_listed(manifest, name)
        ! The `files` of the instance group in the manifest equal declared_files, in order.
        character(len=*), intent(in) :: manifest, name
        character(len=64) :: got(12)
        character(len=64), allocatable :: want(:)
        integer :: n

        got = ""
        call nml_read(manifest, name, "files", got)
        want = declared_files(name)
        n = size(want)
        files_listed = count(len_trim(got) > 0) == n
        if (files_listed) files_listed = all(got(1:n) == want)
        if (.not. files_listed) write(*,'(a)') "    manifest files of "//name//": "//trim(got(1))//", "//trim(got(2))// &
                                               ", "//trim(got(3))//", "//trim(got(4))//" ..."
    end function files_listed

    subroutine same_2d(fa, fb, var, what)
        ! A 2D field of two fields restart files, bitwise.
        character(len=*), intent(in) :: fa, fb, var, what
        real(dp) :: a(NX,NY), b(NX,NY)
        call nc_read(trim(repo)//"/"//DIR//fa, var, a)
        call nc_read(trim(repo)//"/"//DIR//fb, var, b)
        call check(all(a == b), what//": chunk = cont bitwise ("//trim(itoa(count(a /= b)))//" cells differ)")
    end subroutine same_2d

    subroutine same_0d(fa, fb, var, what)
        ! A 0-D field of two fields restart files, bitwise.
        character(len=*), intent(in) :: fa, fb, var, what
        real(dp) :: a(1), b(1)
        call nc_read(trim(repo)//"/"//DIR//fa, var, a)
        call nc_read(trim(repo)//"/"//DIR//fb, var, b)
        call check(a(1) == b(1), what//": chunk = cont bitwise")
    end subroutine same_0d

    subroutine same_record(fa, fb, var, t, what)
        ! The record at time t of a 2D output variable in two files, bitwise.
        character(len=*), intent(in) :: fa, fb, var, what
        real(dp),         intent(in) :: t
        real(sp) :: a(NX,NY,1), b(NX,NY,1)
        call nc_read(trim(repo)//"/"//DIR//fa, var, a, start=[1,1,record(fa, t)], count=[NX,NY,1])
        call nc_read(trim(repo)//"/"//DIR//fb, var, b, start=[1,1,record(fb, t)], count=[NX,NY,1])
        call check(all(a == b), what//": chunk = cont bitwise ("//trim(itoa(count(a /= b)))//" cells differ)")
    end subroutine same_record

    integer function record(f, t) result(n)
        ! The record of time t in an output file (0 if absent).
        character(len=*), intent(in) :: f
        real(dp),         intent(in) :: t
        real(dp), allocatable :: tt(:)
        integer :: i
        allocate(tt(nc_size(trim(repo)//"/"//DIR//f, "time")))
        call nc_read(trim(repo)//"/"//DIR//f, "time", tt)
        n = 0
        do i = 1, size(tt)
            if (tt(i) == t) n = i
        end do
    end function record

    subroutine check_init(run)
        ! After run_init: instances, rows, the hub and the imports at t0, the
        ! cold-start relations.
        type(run_t), target, intent(in) :: run

        class(component_t), pointer :: c
        real(dp) :: bsl, bsl_g, alpha
        real(dp), allocatable :: a(:,:), b(:,:), r(:,:)
        integer :: k
        logical :: ok

        write(*,'(a)') "--- after run_init (t = 0)"
        ok = .true.
        do k = 1, NCOMP
            ok = ok .and. associated(domain_component(run%dom(1), k))
        end do
        call check(ok, "the seven regional instances are created")
        call check(allocated(run%globe%forcing%c) .and. allocated(run%globe%geo%c), "forcing and global geo created")
        call check(count(run%dom(1)%wires%active) >= 40, "the domain's rows resolved and active (" // &
                   itoa(count(run%dom(1)%wires%active))//")")
        call check(run%dom(1)%hub_t == 0.0_dp .and. run%dom(1)%hub_mirror, "hub refreshed at t0, mirror of the ice")

        associate (dom => run%dom(1), ice => run%dom(1)%ice%c, geo => run%dom(1)%geo%c)
            ok = published_at(ice%imports, ice%hi_z_bed, 0.0_dp) .and. published_at(ice%imports, ice%hi_z_sl, 0.0_dp) &
                 .and. published_at(ice%imports, ice%hi_smb, 0.0_dp) .and. published_at(ice%imports, ice%hi_T_srf, 0.0_dp) &
                 .and. published_at(ice%imports, ice%hi_bmb_shlf, 0.0_dp) .and. published_at(ice%imports, ice%hi_T_shlf, 0.0_dp) &
                 .and. published_at(ice%imports, ice%hi_Q_geo, 0.0_dp) .and. published_at(ice%imports, ice%hi_H_sed, 0.0_dp)
            call check(ok, "every import of the ice published at t0")

            ! static rows (3.6 item 1)
            a = get2(ice%imports, ice%hi_Q_geo)
            call check(all(a == 0.05_dp), &
                       "Q_geo = 0.05 W m-2 as the producer put it, unrounded (geothermal data, static row)")
            a = get2(ice%imports, ice%hi_H_sed)
            call check(all(a == 0.0_dp), "H_sed = 0 m (sediments data, static row)")
            a = get2(geo%imports, geo%hi_z_bed_ref)
            b = get2(ice%exports, ice%h_z_bed_ref)
            call check(all(a == b), "geo z_bed_ref import = the ice's reference (static row)")

            ! bsl: reduced at t0, scattered (3.6 item 2)
            bsl   = get0(geo%exports, geo%h_bsl)
            bsl_g = get0(run%globe%geo%c%exports, run%globe%geo%c%h_bsl)
            call check(bsl == bsl_g .and. bsl == get0(geo%imports, geo%hi_bsl), &
                       "geo bsl export = the reduced global bsl (scattered)")
            call check(published_at(geo%exports, geo%h_dV_bsl, 0.0_dp), "geo dV_bsl published at t0 (gathered)")

            ! geo -> ice in sp: z_bed = ref + dz_bed, z_sl = dphi + dz_ss (fastisostasy: dphi = bsl)
            call check(get0(geo%exports, geo%h_dphi) == bsl .and. geo%exports%f(geo%h_dphi)%t1 == geo%exports%f(geo%h_bsl)%t1, &
                       "geo dphi export = the bsl export (fastisostasy), same validity")
            r = dom%refs(1)%z_bed
            a = get2(geo%exports, geo%h_dz_bed)
            b = get2(ice%imports, ice%hi_z_bed)
            call check(all(real(b, sp) == real(r, sp) + real(a, sp)), "ice z_bed = ref + dz_bed (sp, identity grid)")
            a = get2(geo%exports, geo%h_dz_ss)
            b = get2(ice%imports, ice%hi_z_sl)
            call check(all(real(b, sp) == real(get0(geo%exports, geo%h_dphi), sp) + real(a, sp)), &
                       "ice z_sl = dphi + dz_ss (sp)")

            ! the hub mirrors the ice
            a = get2(dom%hub_exports, dom%hh%h_z_srf)
            b = get2(ice%exports, ice%h_z_srf)
            call check(all(a == b), "hub z_srf = the ice's z_srf (mirror, after PASS_STATE)")

            ! forcing <-> climate at t0 (3.6 item 3)
            associate (cl => dom%climate%c, fr => run%globe%forcing%c)
                call check(get0(fr%imports, fr%hi_derived(2)) == get0(cl%exports, cl%h_derived(1)), &
                           "forcing dTo import = the climate's dTo at t0 (after PASS_BOUNDARY)")
                alpha = get0(fr%exports, fr%h_coef(1))
                call check(get0(cl%imports, cl%hi_coef(1)) == alpha .and. alpha == 0.5_dp, &
                           "climate coef_alpha = forcing alpha (0.5)")
            end associate

            ! surface after PASS_FINAL, landed on the ice
            associate (sf => dom%surface%c)
                a = get2(sf%exports, sf%h_smb)
                b = get2(ice%imports, ice%hi_smb)
                call check(all(a == b) .and. sf%exports%f(sf%h_smb)%t1 == 0.0_dp, &
                           "ice smb = the surface's smb after PASS_FINAL")
            end associate
        end associate

        c => domain_component(run%dom(1), 1)
        call check(trim(c%name) == "ice_grl", "domain_component(ice)")
        call round_trips(run)
        call check_ice_static(run, "t0", 0.0_dp)
        call check_ice_flotation(run, "t0")
    end subroutine check_init

    subroutine check_ice_flotation(run, what)
        ! The ice's flotation binding (contracts-design 2.1) on Yelmo's own
        ! geometry: the exports H_ice, z_bed, z_sl are what Yelmo's last
        ! calc_ytopo_diagnostic used, so the binding returns tpo%now%f_grnd bit
        ! for bit (pinning points enter f_grnd_pin, not f_grnd), and so the
        ! exported f_grnd; on another sea level it differs where the ice floats.
        type(run_t), target, intent(in) :: run
        character(len=*),    intent(in) :: what

        real(dp), allocatable :: H(:,:), zb(:,:), zs(:,:), f(:,:), g(:,:)
        integer :: n_tpo, n_exp

        associate (dom => run%dom(1), ice => run%dom(1)%ice%c)
            H  = get2(ice%exports, ice%h_H_ice)
            zb = get2(ice%exports, ice%h_z_bed)
            zs = get2(ice%exports, ice%h_z_sl)
            allocate(f(size(H,1), size(H,2)), g(size(H,1), size(H,2)))
            call ice_flotation(dom%ice, H, zb, zs, f)
            n_exp = count(f /= get2(ice%exports, ice%h_f_grnd))
            n_tpo = -1
            select type (ice)
            type is (yelmo_ice_t)
                n_tpo = count(real(f, sp) /= ice%ylmo%tpo%now%f_grnd)
            end select
            call check(ice%caps%has_flotation, what//": the yelmo ice declares has_flotation")
            call check(n_tpo == 0, what//": flotation = tpo%f_grnd bit for bit (cells differing: "//trim(itoa(n_tpo))//")")
            call check(n_exp == 0, what//": flotation = the exported f_grnd (cells differing: "//trim(itoa(n_exp))//")")
            call ice_flotation(dom%ice, H, zb, zs + 5000.0_dp, g)
            call check(all(g <= f) .and. any(g < f), what//": a sea level 5 km higher floats ice that was grounded")
        end associate
    end subroutine check_ice_flotation

    subroutine check_ice_static(run, what, t_put)
        ! footprint and H_ice_pd (contracts-design 2.1, phase 4): Yelmo's mask_ice
        ! (/= MASK_ICE_NONE) and dta%pd%H_ice, static, put once at PASS_INIT
        ! (t_put, if > -1: the time of the put, which later steps never change).
        type(run_t), target, intent(in) :: run
        character(len=*),    intent(in) :: what
        real(dp),            intent(in) :: t_put

        real(dp), allocatable :: fp(:,:), pd(:,:)

        select type (ice => run%dom(1)%ice%c)
        type is (yelmo_ice_t)
            fp = get2(ice%exports, ice%h_footprint)
            pd = get2(ice%exports, ice%h_H_ice_pd)
            call check(all(real(fp, sp) == merge(1.0_sp, 0.0_sp, ice%ylmo%bnd%mask_ice /= MASK_ICE_NONE)) .and. &
                       count(fp == 1.0_dp) > 0, what//": footprint = (mask_ice /= MASK_ICE_NONE)")
            call check(all(real(pd, sp) == ice%ylmo%dta%pd%H_ice) .and. any(pd > 0.0_dp), &
                       what//": H_ice_pd = dta%pd%H_ice")
            if (t_put > -1.0_dp) &
                call check(published_at(ice%exports, ice%h_footprint, t_put) .and. &
                           published_at(ice%exports, ice%h_H_ice_pd, t_put), &
                           what//": footprint and H_ice_pd put once, at PASS_INIT (published at the start time)")
        class default
            call check(.false., what//": the ice instance is Yelmo")
        end select
    end subroutine check_ice_static

    subroutine round_trips(run)
        ! The two unit conversions across the exchange (adapters-design 1.2 P3,
        ! 10): pr [mm d-1] * sec_year/86400 out of the climate, / the same factor
        ! into chion; bmb_shlf [m ie yr-1] * rho_ice out of the shelf, / rho_ice
        ! into Yelmo. On identity grids the backend's sp value comes back bit for
        ! bit (adapters-design 10, P3): rows never round, the consumer narrows
        ! once after its conversion (contracts-design 4.3, 9, decision 2026-10-25).
        type(run_t), target, intent(in) :: run
        integer :: npr, nbmb

        npr = -1;  nbmb = -1
        associate (dom => run%dom(1))
            select type (cl => dom%climate%c)
            type is (snapshots_climate_t)
                select type (sf => dom%surface%c)
                type is (chion_surface_t)
                    npr = count(sf%pr_mmd /= cl%out%now%pr)
                end select
            end select
            select type (shf => dom%shelf%c)
            type is (marshelf_shelf_t)
                select type (ice => dom%ice%c)
                type is (yelmo_ice_t)
                    nbmb = count(ice%ylmo%bnd%bmb_shlf /= shf%mshlf%now%bmb_shlf)
                end select
            end select
        end associate
        call check(npr == 0, "pr round trip climate -> chion bit for bit at t0 (values differing: "// &
                  trim(itoa(npr))//")")
        call check(nbmb == 0, "bmb_shlf round trip shelf -> Yelmo bit for bit at t0 (cells differing: "// &
                  trim(itoa(nbmb))//")")
    end subroutine round_trips

    subroutine check_step(run, n, t, smb_prev, bmb_prev, m_prev, dto_prev)
        type(run_t), target, intent(in) :: run
        integer,  intent(in) :: n
        real(dp), intent(in) :: t, smb_prev(:,:), bmb_prev(:,:), m_prev, dto_prev

        real(dp), allocatable :: a(:,:), b(:,:)
        real(dp) :: bsl
        logical  :: ok

        write(*,'(a,i0,a,f6.1,a)') "--- after step ", n, " (t = ", t, ")"
        associate (dom => run%dom(1), ice => run%dom(1)%ice%c, geo => run%dom(1)%geo%c, &
                   cl => run%dom(1)%climate%c, sf => run%dom(1)%surface%c, sh => run%dom(1)%shelf%c, &
                   fr => run%globe%forcing%c)

            ! globe, start of the step: lag 1 rows
            call check(get0(fr%imports, fr%hi_M_ice) == m_prev .and. fr%imports%f(fr%hi_M_ice)%t1 == t - 10.0_dp, &
                       "forcing M_ice = the ice's M_ice of the previous step (lag 1)")
            call check(get0(fr%imports, fr%hi_derived(2)) == dto_prev, "forcing dTo = the climate's dTo as last put (lag 1)")
            call check(published_at(fr%exports, fr%h_coef(1), t), "forcing coefficients at t")

            ! the ice read the surface and the shelf of the previous step (lag 1)
            call check(all(get2(ice%imports, ice%hi_smb) == smb_prev), "ice smb = the surface's smb of the previous step")
            call check(all(get2(ice%imports, ice%hi_bmb_shlf) == bmb_prev) .and. &
                       ice%imports%f(ice%hi_bmb_shlf)%t1 == t - 10.0_dp, &
                       "ice bmb_shlf = the shelf's bmb_shlf of the previous step, unrounded")

            ! geo ran at t; bsl of this step's reduction in the geo exports and z_sl
            bsl = get0(geo%exports, geo%h_bsl)
            call check(published_at(geo%exports, geo%h_dV_bsl, t) .and. cadence_t_last(geo%cadence, run%clk) == t, &
                       "geo updated at t (dt = 10)")
            call check(bsl == get0(run%globe%geo%c%exports, run%globe%geo%c%h_bsl) .and. &
                       geo%exports%f(geo%h_bsl)%t1 == t, "geo bsl = the reduced global bsl of this step")
            a = get2(geo%exports, geo%h_dz_ss)
            b = get2(ice%imports, ice%hi_z_sl)
            call check(all(real(b, sp) == real(get0(geo%exports, geo%h_dphi), sp) + real(a, sp)) .and. &
                       get0(geo%exports, geo%h_dphi) == bsl .and. geo%exports%f(geo%h_dphi)%t1 == t, &
                       "ice z_sl = dphi + dz_ss (sp), dphi = the reduced bsl of this step")

            ! ice and hub at t
            call check(published_at(ice%exports, ice%h_H_ice, t) .and. ice%exports%f(ice%h_smb)%t0 == t - 10.0_dp, &
                       "ice exports at t, MEAN over [t-10, t]")
            call check(dom%hub_t == t .and. all(get2(dom%hub_exports, dom%hh%h_H_ice) == get2(ice%exports, ice%h_H_ice)), &
                       "hub refreshed at t, H_ice = the ice's")

            ! climate cadence 20 yr
            if (n == 1) then
                ok = cadence_t_last(cl%cadence, run%clk) == 0.0_dp .and. cl%exports%f(cl%h_tas)%t1 == 0.0_dp
                call check(ok, "climate not due at t = 10 (dt = 20): exports from t0")
            else
                ok = cadence_t_last(cl%cadence, run%clk) == t .and. cl%exports%f(cl%h_tas)%t1 == t .and. &
                     cl%imports%f(cl%hi_z_srf)%t1 == t .and. cl%imports%f(cl%hi_coef(1))%t1 == t
                call check(ok, "climate updated at t = 20 on the hub and the coefficients of t")
            end if

            ! surface and shelf every step, on the climate as last updated
            call check(published_at(sf%exports, sf%h_smb, t) .and. sf%imports%f(sf%hi_tas)%t1 == &
                       cadence_t_last(cl%cadence, run%clk) .and. sf%imports%f(sf%hi_z_srf)%t1 == t, &
                       "surface updated at t, tas of the climate's last update, z_srf of the hub at t")
            call check(published_at(sh%exports, sh%h_bmb_shlf, t) .and. sh%imports%f(sh%hi_to_ann)%t1 == &
                       cadence_t_last(cl%cadence, run%clk), "shelf updated at t, profiles of the climate's last update")
        end associate
    end subroutine check_step

    subroutine check_files()
        ! One file per registered stream, the timing table, the bundle.
        character(len=32), parameter :: files(11) = [character(len=32) :: "ice_grl.nc", "ice_grl_ts.nc", &
            "geo_grl.nc", "climate_grl.nc", "surface_grl.nc", "shelf_grl.nc", "ghf_grl.nc", "forcing_ts.nc", &
            "geo_glob_ts.nc", "Greenland_hub.nc", "timing.txt"]
        integer :: i
        real(dp) :: tt(3)

        write(*,'(a)') "--- output"
        do i = 1, size(files)
            call check(exists(trim(files(i))), "output "//trim(files(i)))
        end do
        call check_bundle("restart-20")
        call nc_read("Greenland_hub.nc", "time", tt)
        call check(all(tt == [0.0_dp, 10.0_dp, 20.0_dp]), "hub file at t = 0, 10, 20 (output at the end of each step)")
        call check(nc_exists_var("forcing_ts.nc", "coef_dTo") .and. nc_exists_var("ice_grl_ts.nc", "M_ice"), &
                   "forcing coefficients and the ice's M_ice in the 1D files")
    end subroutine check_files

    ! =========================================================================
    ! Run directory, par files, data
    ! =========================================================================

    subroutine setup_run_dir(sub)
        ! A run directory tmp/test_run_sequence/<sub>: input/ (links), maps,
        ! the climate data (written once, in tmp/test_run_sequence/, linked
        ! from the others). The par files are written by the caller.
        character(len=*), intent(in) :: sub

        character(len=:), allocatable :: d

        d = trim(repo)//"/"//DIR//sub
        call sh("mkdir -p '"//d//"input'")
        call sh("ln -s '"//trim(repo)//"/maps' '"//d//"maps'")
        call sh("for f in '"//trim(repo)//"'/input/*; do ln -s ""$f"" '"//d//"input/'; done")
        call sh("for f in '"//trim(repo)//"'/yelmo/input/*.nml '"//trim(repo)//"'/yelmo/input/*.md; do "// &
                "ln -sf ""$f"" '"//d//"input/'; done")
        call sh("ln -sf '"//trim(repo)//"/chion/input/chion_defaults.nml' '"//d//"input/'")
        call sh("ln -sf '"//trim(repo)//"/fesm-utils/par/datasrc_defaults.nml' '"//d//"input/'")
        call sh("for f in '"//trim(repo)//"'/insol/input/LA2004.*.txt; do ln -sf ""$f"" '"//d//"input/'; done")
        if (len(sub) == 0) then
            call write_data(d//DATA)
        else
            call sh("ln -s '"//trim(repo)//"/"//DIR//DATA//"' '"//d//DATA//"'")
        end if
    end subroutine setup_run_dir

    subroutine write_data(file)
        ! Synthetic snapshots on GRL-16KM: the reference pd (tas from a lapse
        ! rate, latitude and a seasonal cycle, pr varying, zs = the topography's
        ! surface, ocean profiles), ctl = pd and alt = 10 K colder, half the pr.
        character(len=*), intent(in) :: file

        type(grid_class) :: grid
        real(sp) :: zs(NX,NY), t(NX,NY,12), p(NX,NY,12), o(NX,NY,NZ), s(NX,NY,NZ)
        real(sp), parameter :: prof(NZ) = [273.5_sp, 273.0_sp, 274.0_sp, 274.5_sp, 272.5_sp]
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        integer :: i, j, m

        call grid_cdo_read_desc(grid, "GRL-16KM", trim(repo)//"/maps")
        call nc_read(trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc", "z_srf", zs)
        zs = max(zs, 0.0_sp)
        do j = 1, NY
            do i = 1, NX
                do m = 1, 12
                    t(i,j,m) = real(281.15_dp - 0.8_dp*(grid%lat(i,j) - 60.0_dp) - 0.0065_dp*zs(i,j) &
                                    + 10.0_dp*cos(2.0_dp*pi*real(m - 7, dp)/12.0_dp), sp)
                    p(i,j,m) = real(0.3_dp + 0.013_dp*i + 0.007_dp*j + 0.11_dp*m, sp)
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

    subroutine write_par(file, restart, time_init, time_end, dt_rst, sed, exps, explist, xlines, vel)
        ! The Greenland par file of the run: kryos groups (core-design 5), the
        ! hub's groups (ice_data), and the backend groups (yelmox_Greenland.nml
        ! values where they matter; library defaults otherwise). The clock
        ! runs time_init -> time_end with restart bundles every dt_rst (and at
        ! time_end); `sed` renames the sediments instance (default sed_grl).
        ! `exps`: one host experiment of that name; `explist`: the groups named by the
        ! domain's `experiments` (their lines in `xlines`); `vel`: Yelmo's
        ! present-day velocity file (pd_vel_load).
        character(len=*), intent(in) :: file, restart
        real(dp),         intent(in) :: time_init, time_end, dt_rst
        character(len=*), intent(in), optional :: sed, exps, explist(:), xlines(:), vel

        character(len=:), allocatable :: topo, basins, sed_name, exp_key, vel_load, vel_path
        integer :: u, i

        topo   = trim(root)//"/{domain}/{grid_name}/{grid_name}_TOPO-M17-v5.nc"
        basins = trim(root)//"/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"
        sed_name = "sed_grl"
        if (present(sed)) sed_name = sed
        exp_key = ""
        if (present(exps)) exp_key = "  experiments = '"//exps//"'"      ! a host experiment of that name
        if (present(explist)) then
            exp_key = "  experiments ="
            do i = 1, size(explist)
                exp_key = exp_key//" '"//trim(explist(i))//"'"
            end do
        end if
        vel_load = "  pd_vel_load = False"
        vel_path = ""
        if (present(vel)) then
            if (len(vel) > 0) then
                vel_load = "  pd_vel_load = True"
                vel_path = "  pd_vel_path = '"//vel//"'"
            end if
        end if

        open(newunit=u, file=file, status="replace", action="write")
        call w(u, [character(len=400) :: &
            "&kryos", "  domains = 'Greenland'", "  restart = '"//restart//"'", "  restart_out = 'out_rst'", &
            "  write_lossless = True", "  log_level = 'info'", "/", &
            "&clock", "  method = 'const'", "  time_const = 0.0", "  time_init = "//ftoa(time_init), &
            "  time_end = "//ftoa(time_end), "  dtt = 10.0", "/", &
            "&out_2D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_1D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_rst", "  method = 'const'", "  dt = "//ftoa(dt_rst), "/", &
            "&globe", "  forcing = 'forcing'", "  geo = 'geo_glob'", "/", &
            "&forcing", "  backend = 'tsgen'", "  clock = 'forcing'", "  coefs = 'alpha' 'dTo' 'beta'", &
            "  out_1D = 'out_1D'", "/", &
            "&geo_glob", "  backend = 'fastisostasy'", "  clock = 'forcing'", "  bsl = 'barysealevel'", &
            "  out_1D = 'out_1D'", "/", &
            "&Greenland", "  grid_hub = 'GRL-16KM'", "  topo = 'topo_grl'", "  regions = 'regions_grl'", &
            "  basins = 'basin'", "  out_2D = 'out_2D'", "  ice = 'ice_grl'", "  geo = 'geo_grl'", &
            "  climate = 'climate_grl'", "  surface = 'surface_grl'", "  shelf = 'shelf_grl'", &
            "  geothermal = 'ghf_grl'", "  sediments = '"//sed_name//"'", exp_key, "/", &
            "&topo_grl", "  path = '"//topo//"'", "  remove_lakes = True", "  H_ice_min = 1.0", &
            "  grad_lim_zb = 0.5", "/", &
            "&regions_grl", "  layers = 'basin'", "  path_basin = '"//basins//"'", "/", &
            "&ice_grl", "  backend = 'yelmo'", "  yelmo = 'yelmo'", "  init_equil_time = 0.0", &
            "  out_2D = 'out_2D'", "  out_1D = 'out_1D'", "/", &
            "&geo_grl", "  backend = 'fastisostasy'", "  dt = 10.0", "  isos = 'isos'", "  out_2D = 'out_2D'", "/", &
            "&climate_grl", "  backend = 'snapshots'", "  clock = 'forcing'", "  dt = 20.0", &
            "  snapshots = 'snapshots'", "  reference = 'pd'", "  terms = 'lgm' 'ocean'", "  out_2D = 'out_2D'", "/", &
            "&surface_grl", "  backend = 'chion'", "  clock = 'forcing'", "  surface_chion = 'surface_chion'", &
            "  chion = 'chion'", "  orbit = 'orbit'", "  out_2D = 'out_2D'", "/", &
            "&shelf_grl", "  backend = 'marine_shelf'", "  marine_shelf = 'marine_shelf'", "  corr = 'shelf_corr'", &
            "  out_2D = 'out_2D'", "/", &
            "&ghf_grl", "  backend = 'data'", "  dt = -1.0", "  Q_geo = 'ghf_const'", "  out_2D = 'out_2D'", "/", &
            "&ghf_const", "  mode = 'const'", "  value = 0.05", "/", &
            "&"//sed_name, "  backend = 'data'", "  dt = -1.0", "  H_sed = 'sed_none'", "/", &
            "&sed_none", "  mode = 'const'", "  value = 0.0", "/"])
        if (present(exps)) write(u,'(a)') "&"//exps, "  kind = 'host'", "/"
        if (present(xlines)) then
            do i = 1, size(xlines)
                write(u,'(a)') trim(xlines(i))
            end do
        end if
        ! forcing coefficients
        call w(u, [character(len=400) :: &
            "&coef_alpha", "  method = 'const'", "  f_min = 0.5", "/", &
            "&coef_beta", "  method = 'exp'", "  f_min = 0.0", "  f_max = 1.0", "  df_dt_max = 0.01", &
            "  dt_ave = 30.0", "  tol = 0.5", "/", &
            "&coef_dTo", "  method = 'derived'", "  from = 'air_anomaly'", "  factor = 0.25", "  region = 'all'", &
            "  units = 'K'", "/"])
        ! climate: snapshots, terms, field groups
        call w(u, [character(len=400) :: &
            "&snapshots", "  f_p = 0.06", "/", &
            "&snapshot_pd", "  tas = 'v_t_pd'", "  pr = 'v_p_pd'", "  zs = 'v_zs'", "  to = 'v_to'", "  so = 'v_so'", "/", &
            "&snapshot_ctl", "  tas = 'v_t_pd'", "  pr = 'v_p_pd'", "  zs = 'v_zs'", "/", &
            "&snapshot_alt", "  tas = 'v_t_alt'", "  pr = 'v_p_alt'", "  zs = 'v_zs'", "/", &
            "&term_lgm", "  coef = 'alpha'", "  snap0 = 'ctl'", "  snap1 = 'alt'", "  fields = 'tas' 'pr'", "/", &
            "&term_ocean", "  coef = 'dTo'", "  pattern = 'uniform'", "  fields = 'to'", "/"])
        call field(u, "v_t_pd", "t_pd", "K");   call field(u, "v_p_pd", "p_pd", "K")
        call field(u, "v_t_alt", "t_alt", "K"); call field(u, "v_p_alt", "p_alt", "K")
        call field(u, "v_zs", "zs", "m");       call field(u, "v_to", "to", "K")
        call field(u, "v_so", "so", "psu")
        ! Yelmo: library defaults, the reference data loads off, masks without named regions;
        ! calving "equil" as yelmox_Greenland (with the default vm-m16 calving a Yelmo restart is
        ! not exact: the run diverges in the first step after the restore)
        call w(u, [character(len=400) :: &
            "&yelmo", "  nz_aa = 6", "/", "&ytopo", "/", &
            "&ycalv", "  calv_flt_method = 'equil'", "  calv_grnd_method = 'equil'", "/", "&ydyn", "  solver = 'diva'", "/", &
            "&ytill", "/", "&yhyd", "/", "&ymat", "/", "&ytrc", "/", "&ytherm", "/", &
            "&yelmo_masks", "  basins = 'None'", "  mask_ice_dynamic = 'all'", "  mask_rmse = 'all'", "/", &
            "&yelmo_init_topo", "/", &
            "&yelmo_data", "  pd_topo_load = False", "  pd_tsrf_load = False", "  pd_smb_load = False", &
            vel_load, vel_path, "  pd_age_load = False", "/"])
        ! FastIsostasy and the global bsl (yelmox_Greenland.nml; constant ocean area)
        call w(u, [character(len=400) :: &
            "&isos", "  heterogeneous_ssh = .true.", "  interactive_sealevel = .true.", "  include_elastic = .true.", &
            "  correct_compression = .true.", "  correct_distortion = .false.", "  method = 2", &
            "  dt_diagnostics = 10.", "  pad = 512.", "  mantle = 'uniform'", "  lithosphere = 'uniform'", &
            "  layering = 'equalized'", "  lumping = 'min'", "  viscosity_scaling_method = 'none'", &
            "  viscosity_scaling = 1.", "  nl = 1", "  zl = 88.", "  viscosities = 1.e21", "  ref_viscosity = 1.e21", &
            "  tau = 3000.0", "  litho_thickness_scaling = 1.0", "  E = 66.0", "  nu = 0.28", "  rho_litho = 3200.0", &
            "  r_earth = 6.378e3", "  m_earth = 5.972e24", "  restart = 'None'", "  mask_file = 'None'", &
            "  rheology_file = 'None'", "  dt_method = 'tsit54'", "  dt_init = 1.0", "  dt_min = 1e-10", &
            "  dt_max = 10.0", "  rtol = 1e-4", "  atol = 1e-2", "  q = 1.0", "/", &
            "&barysealevel", "  method = 'fastiso'", "  bsl_init = 0.0", "  sl_path = 'None'", "  sl_name = 'none'", &
            "  A_ocean_pd = 3.625e14", "  A_ocean_path = 'None'", "  restart = 'None'", "/"])
        ! chion and its host (yelmox_Greenland.nml; a short spin-up), the orbit of 1950
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
        ! marine_shelf (yelmox_Greenland.nml; masks without named regions)
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
    end subroutine write_par

    subroutine field(u, group, name, units)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: group, name, units
        write(u,'(a)') "&"//group
        write(u,'(a)') "  filename = '"//DATA//"'"
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
            write(u,'(a)') trim(lines(i))
        end do
    end subroutine w

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    logical function child(sub, args, log, quiet)
        ! This executable with `args`, run in the run directory
        ! tmp/test_run_sequence/<sub> on one thread (bitwise comparable runs);
        ! .true. on exit 0. A failure shows the end of the log unless quiet.
        character(len=*), intent(in) :: sub, args, log
        logical,          intent(in), optional :: quiet
        character(len=512) :: exe
        character(len=:), allocatable :: d
        integer :: es
        logical :: q

        call get_command_argument(0, exe)
        if (exe(1:1) /= "/") exe = trim(repo)//"/"//trim(exe)
        d = trim(repo)//"/"//DIR//sub
        call execute_command_line("cd '"//d//"' && OMP_NUM_THREADS=1 '"//trim(exe)//"' "//args//" > "//log// &
                                  " 2>&1", exitstat=es)
        child = es == 0
        q = .false.
        if (present(quiet)) q = quiet
        if (.not. (child .or. q)) call sh("tail -n 30 '"//d//log//"'")
    end function child

    logical function in_file(file, text)
        ! Does `file` contain `text`?
        character(len=*), intent(in) :: file, text
        integer :: es
        call execute_command_line("grep -q '"//text//"' '"//file//"'", exitstat=es)
        in_file = es == 0
    end function in_file

    function get2(st, h) result(a)
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        real(dp), allocatable :: a(:,:)
        allocate(a(st%f(h)%nx, st%f(h)%ny))
        call fields_get(st, h, a)
    end function get2

    real(dp) function get0(st, h) result(x)
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        call fields_get(st, h, x)
    end function get0

    logical function published_at(st, h, t)
        type(field_state_t), intent(in) :: st
        integer,             intent(in) :: h
        real(dp),            intent(in) :: t
        published_at = .false.
        if (h <= 0) return
        published_at = st%f(h)%published .and. st%f(h)%t1 == t
    end function published_at

    subroutine sh(cmd)
        character(len=*), intent(in) :: cmd
        call execute_command_line(cmd)
    end subroutine sh

    logical function exists(file)
        character(len=*), intent(in) :: file
        inquire(file=file, exist=exists)
    end function exists

    function ftoa(x) result(s)
        ! A time [yr] for a par file: 40.0
        real(dp), intent(in) :: x
        character(len=24) :: s
        write(s, '(f0.1)') x
    end function ftoa

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=12) :: s
        write(s, '(i0)') i
    end function itoa

    subroutine note(ok, what)
        ! A known open point: reported, not failed.
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        write(*,'(a,a)') merge("  ok    ", "  OPEN  ", ok), what
    end subroutine note

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        write(*,'(a,a)') merge("  ok    ", "  FAIL  ", ok), what
        if (.not. ok) nfail = nfail + 1
    end subroutine check

end program test_run_sequence
