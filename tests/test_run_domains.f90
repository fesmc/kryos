program test_run_domains
    ! Phase 3 WP1: the run with two domains (core-design 3.1-3.4 for N = 2). Greenland (GRL-16KM) and Antarctica
    ! (ANT-32KM) with every phase 1 backend, climate snapshots written here (synthetic fields on each grid), as
    ! the science case NH + ANT in small: the par file of the two-domain run is made from the two single-domain
    ! files by the merge step (tools/bin/kryos-merge.jl, KryosTools.Convert.merge_par); the comparisons are the
    ! kryos-against-kryos checks of the gates (tools/bin/kryos-gate-kk.jl, KryosTools.Gate), phase3-plan 3.2:
    !
    !   const   bsl = const: the domains are independent, so each domain of the two-domain run equals the
    !           single-domain run of its own configuration bit for bit (every instance file, the bundles at 40
    !           and 80, the globe files except the aggregates over the domains); also with an experiment in each
    !           domain (recon in Greenland, relax in Antarctica: the experiment lists are per domain);
    !   sci     bsl = fastiso: the domain order swapped in [kryos] domains gives identical files and bundles
    !           (the reduction is order-independent); chunks 0 -> 40 -> 80 equal the continuous run (restarts);
    !           the bsl budget closes: bsl(t) - bsl(t-dt) = -sum_d dV_bsl_d / A_ocn, recomputed from the exports
    !           (here and by KryosTools.Gate.bsl_budget);
    !   stop    two climates with the derived coefficient dTo stop at init, naming the coefficient (contracts 7).
    !
    ! Inputs: the hub topographies and basins of GRL-16KM and ANT-32KM from ice_data ($KRYOS_ICE_DATA, default
    ! ./ice_data); without them the test reports and passes. Julia with the project of tools/ ($KRYOS_TOOLS) is
    ! needed for the merge and the comparisons; without it the test reports OPEN and stops there.
    ! The runs execute in tmp/test_run_domains/<run>/ in child processes of this executable (one thread):
    !     test_run_domains.x plain <par file>     the run as the driver does it
    ! Built only with yelmo=1 fastisostasy=1 chion=1 (config/common.mk).

    use precision,       only: sp, dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_read, nc_size
    use coords,          only: grid_class, grid_cdo_read_desc
    use kryos_run,       only: run_t, run_init, run_step, run_done, run_end

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_run_domains/"
    real(dp), parameter :: T_MID = 40.0_dp, T_END = 80.0_dp
    integer, parameter :: NZ = 5, NDOM = 2
    ! the two domains: name, instance suffix, hub grid, topography file
    character(len=16), parameter :: DNAME(NDOM) = [character(len=16) :: "Greenland", "Antarctica"]
    character(len=4),  parameter :: DSFX(NDOM)  = [character(len=4)  :: "grl", "ant"]
    character(len=16), parameter :: DGRID(NDOM) = [character(len=16) :: "GRL-16KM", "ANT-32KM"]
    character(len=32), parameter :: DTOPO(NDOM) = [character(len=32) :: "TOPO-M17-v5", "TOPO-BedMachine"]
    ! Variables the checks do not compare. AGG: the aggregates over the domains (the forcing's ice mass, the gathered
    ! dV_bsl; globe files in the bundles). PCT: Yelmo's pc_tau_max, a running maximum over 50 internal steps that a
    ! restart does not hold (refilled in the 500-yr chunk of the gates, not in this test's 40 yr).
    character(len=*), parameter :: AGGV = '"*forcing_imports.nc:M_ice", "*geo_glob_imports.nc:dV_bsl"', PCTV = '"*:pc_tau_max"', &
                                   AGG = 'skip = ['//AGGV//']', PCT = 'skip = ['//PCTV//']'
    integer, parameter :: DNX(NDOM) = [106, 191], DNY(NDOM) = [181, 191]

    character(len=512) :: arg, par_arg, repo, root, tools
    character(len=:), allocatable :: d
    integer :: nfail = 0, es
    logical :: ok

    call get_command_argument(1, arg)
    if (trim(arg) == "plain") then
        call get_command_argument(2, par_arg)
        call plain_child(trim(par_arg))
        stop 0
    end if

    call get_environment_variable("PWD", repo)
    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = trim(repo)//"/ice_data"
    if (root(1:1) /= "/") root = trim(repo)//"/"//trim(root)
    if (.not. (exists(topo_file(1)) .and. exists(topo_file(2)))) then
        write(*,'(a)') "test_run_domains: ice_data topographies of GRL-16KM and ANT-32KM not found under "// &
                       trim(root)//": skipped"
        stop 0
    end if
    call get_environment_variable("KRYOS_TOOLS", tools)
    if (len_trim(tools) == 0) tools = trim(repo)//"/tools"
    call execute_command_line("julia --project='"//trim(tools)//"' -e 'using KryosTools' > /dev/null 2>&1", exitstat=es)
    if (es /= 0) then
        write(*,'(a)') "test_run_domains: OPEN: Julia or the project "//trim(tools)//" is not available (the merge "// &
                       "step and the comparisons need it; set KRYOS_TOOLS): skipped"
        stop 0
    end if

    d = trim(repo)//"/"//DIR
    call sh("rm -rf '"//d//"'")
    call sh("mkdir -p '"//d//"'")
    call write_data(1, d//"clim_grl.nc")
    call write_data(2, d//"clim_ant.nc")

    call const_set()
    call sci_set()
    call experiment_set()
    call producers_stop()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_run_domains: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_run_domains: OK"

contains

    ! =========================================================================
    ! The run sets
    ! =========================================================================

    subroutine const_set()
        ! bsl = const: each domain of the two-domain run equals its single-domain run.
        write(*,'(a)') "--- bsl = const: the domains are independent"
        call run_single("g_const", 1, "const", 0.0_dp, T_END, "None")
        call run_single("a_const", 2, "const", 0.0_dp, T_END, "None")
        call run_merged("s_const", "const", 0.0_dp, T_END, "None", swap=.false.)
        call write_gate("gate_const.toml", checks= &
            [character(len=400) :: &
             '[[check]]', 'type = "domain"', 'domain = "Greenland"', 'ref = "g_const"', 'run = "s_const"', &
             '[[check]]', 'type = "domain"', 'domain = "Antarctica"', 'ref = "a_const"', 'run = "s_const"'], &
            skip=AGG)
        call gate_check("gate_const.toml", "const", "each domain of the two-domain run equals its single-domain run "// &
                        "(every file, the bundles at 40 and 80, the globe files except the aggregates)")
    end subroutine const_set

    subroutine sci_set()
        ! bsl = fastiso: order swap, chunks, budget.
        write(*,'(a)') "--- bsl = fastiso: order swap, chunks, budget"
        call run_merged("cont", "fastiso", 0.0_dp, T_END, "None", swap=.false.)
        call run_merged("swap", "fastiso", 0.0_dp, T_END, "None", swap=.true.)
        call run_merged("chunk/c00", "fastiso", 0.0_dp, T_MID, "None", swap=.false.)
        call run_merged("chunk/c01", "fastiso", T_MID, T_END, "../c00/restart-40", swap=.false.)
        call write_gate("gate_sci.toml", checks= &
            [character(len=400) :: &
             '[[check]]', 'type = "identical"', 'ref = "cont"', 'run = "swap"', &
             '[[check]]', 'type = "chunks"', 'ref = "cont"', 'run = "chunk"', &
             '[[check]]', 'type = "budget"', 'run = "cont"'], skip=PCT)
        call gate_check("gate_sci.toml", "sci", "domain order swap identical, chunks == cont, the bsl budget closes "// &
                        "(KryosTools.Gate)")
        call check_budget("cont")
        call check_budget("swap")
    end subroutine sci_set

    subroutine experiment_set()
        ! One experiment per domain (recon in Greenland, relax in Antarctica, windows across the restart):
        ! the experiment lists are per domain, so the const equality holds with them too.
        write(*,'(a)') "--- an experiment in each domain (bsl = const)"
        call run_single("xg", 1, "const", 0.0_dp, T_END, "None", exps="exp_recon")
        call run_single("xa", 2, "const", 0.0_dp, T_END, "None", exps="exp_relax")
        call run_merged("x", "const", 0.0_dp, T_END, "None", swap=.false., expg="exp_recon", expa="exp_relax")
        call run_merged("xchunk/c00", "const", 0.0_dp, T_MID, "None", swap=.false., expg="exp_recon", expa="exp_relax")
        call run_merged("xchunk/c01", "const", T_MID, T_END, "../c00/restart-40", swap=.false., expg="exp_recon", &
                        expa="exp_relax")
        call write_gate("gate_x.toml", checks= &
            [character(len=400) :: &
             '[[check]]', 'type = "domain"', 'domain = "Greenland"', 'ref = "xg"', 'run = "x"', &
             '[[check]]', 'type = "domain"', 'domain = "Antarctica"', 'ref = "xa"', 'run = "x"', &
             '[[check]]', 'type = "chunks"', 'ref = "x"', 'run = "xchunk"'], skip='skip = ['//AGGV//', '//PCTV//']')
        call gate_check("gate_x.toml", "exp", "recon in Greenland and relax in Antarctica: each domain equals its "// &
                        "single-domain run, and the chunks equal the continuous run")
    end subroutine experiment_set

    subroutine producers_stop()
        ! Two climates with a term on the derived coefficient dTo: the run stops at init, naming the coefficient.
        character(len=:), allocatable :: s
        write(*,'(a)') "--- two producers of one derived coefficient"
        s = "dto"
        call sh("mkdir -p '"//d//s//"'")
        call setup_run_dir(s)
        call write_domain_par(d//s//"/grl.nml", 1, "fastiso", 0.0_dp, 20.0_dp, "None", dto=.true.)
        call write_domain_par(d//s//"/ant.nml", 2, "fastiso", 0.0_dp, 20.0_dp, "None", dto=.true.)
        call merge_pars(s, swap=.false.)
        ok = child(s, "plain par.nml", "run.log", quiet=.true.)
        call check(.not. ok .and. in_file(d//s//"/run.log", "dTo"), "the run stops at init, naming dTo ("//DIR//s//"/run.log)")
        call check(in_file(d//s//"/run.log", "coefficient of its own"), "the message suggests a coefficient of its own per domain")
    end subroutine producers_stop

    ! =========================================================================
    ! Runs and comparisons
    ! =========================================================================

    subroutine run_single(sub, id, bsl, t0, t1, restart, exps)
        character(len=*), intent(in) :: sub, bsl, restart
        integer,          intent(in) :: id
        real(dp),         intent(in) :: t0, t1
        character(len=*), intent(in), optional :: exps
        call setup_run_dir(sub)
        call write_domain_par(d//sub//"/par.nml", id, bsl, t0, t1, restart, exps=exps)
        call check(child(sub, "plain par.nml", "run.log"), "single-domain run "//sub//" ("//DIR//sub//"/run.log)")
    end subroutine run_single

    subroutine run_merged(sub, bsl, t0, t1, restart, swap, expg, expa)
        character(len=*), intent(in) :: sub, bsl, restart
        real(dp),         intent(in) :: t0, t1
        logical,          intent(in) :: swap
        character(len=*), intent(in), optional :: expg, expa
        call setup_run_dir(sub)
        call write_domain_par(d//sub//"/grl.nml", 1, bsl, t0, t1, restart, exps=expg)
        call write_domain_par(d//sub//"/ant.nml", 2, bsl, t0, t1, restart, exps=expa)
        call merge_pars(sub, swap)
        call check(child(sub, "plain par.nml", "run.log"), "two-domain run "//sub// &
                   trim(merge(" (domains swapped)", "                  ", swap))//" ("//DIR//sub//"/run.log)")
    end subroutine run_merged

    subroutine merge_pars(sub, swap)
        ! par.nml of the run directory: the merge of grl.nml and ant.nml (the order of the files is the order of the domains).
        character(len=*), intent(in) :: sub
        logical,          intent(in) :: swap
        integer :: es
        call execute_command_line("cd '"//d//sub//"' && julia --project='"//trim(tools)//"' '"//trim(repo)// &
                                  "/tools/bin/kryos-merge.jl' grl.nml ant.nml --out par.nml --report merge.md"// &
                                  trim(merge(" --swap", "       ", swap))//" > merge.log 2>&1", exitstat=es)
        if (es /= 0) call sh("tail -n 20 '"//d//sub//"/merge.log'")
        call check(es == 0, "merge of the two single-domain par files ("//DIR//sub//"/merge.log)")
    end subroutine merge_pars

    subroutine write_gate(file, checks, skip)
        ! A kryos-against-kryos gate file for the runs of this test. `skip`: the variables not compared.
        character(len=*),  intent(in) :: file, skip
        character(len=*),  intent(in) :: checks(:)
        integer :: u, i
        open(newunit=u, file=d//file, status="replace", action="write")
        write(u,'(a)') 'gate = "test_run_domains"', 'kind = "kryos"', 'chunks = [0.0, 40.0, 80.0]', &
                       'restart_times = [40.0, 80.0]', 'domain_order = ["Greenland", "Antarctica"]', &
                       skip, &
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
        if (es /= 0) call sh("grep -E '^# |^\| ' '"//d//"gate_"//tag//".md' | head -n 30; "// &
                             "grep -E '\| fail \|' '"//d//"gate_"//tag//".md' | head -n 30; tail -n 10 '"//d//"gate_"//tag//".log'")
        call check(es == 0, what//" (report "//DIR//"gate_"//tag//".md)")
    end subroutine gate_check

    subroutine check_budget(sub)
        ! The bsl budget of the run `sub` from its exports: bsl(t) - bsl(t - dt) = -sum_d dV_bsl_d(t) / A_ocn(t) in the
        ! round-off of the sp sums. Also: the domains' dV_bsl are not all zero (a budget of nothing proves nothing),
        ! and the reduced bsl moved.
        character(len=*), intent(in) :: sub
        real(sp), allocatable :: bsl(:), a(:), dv(:,:)
        real(dp) :: tot, resid, scale, worst
        integer :: nt, k, i

        nt = nc_size(d//sub//"/geo_glob_ts.nc", "time")
        allocate(bsl(nt), a(nt), dv(nt, NDOM))
        call nc_read(d//sub//"/geo_glob_ts.nc", "bsl", bsl)
        call nc_read(d//sub//"/geo_glob_ts.nc", "A_ocn", a)
        do i = 1, NDOM
            call nc_read(d//sub//"/geo_"//trim(DSFX(i))//"_ts.nc", "dV_bsl", dv(:,i))
        end do
        worst = 0.0_dp
        do k = 2, nt
            tot = sum(real(dv(k,:), dp))
            resid = real(bsl(k), dp) - real(bsl(k-1), dp) + tot/real(a(k), dp)
            scale = abs(real(bsl(k), dp)) + abs(real(bsl(k-1), dp)) + sum(abs(real(dv(k,:), dp)))/real(a(k), dp)
            if (scale > 0.0_dp) worst = max(worst, abs(resid)/scale)
        end do
        call check(any(dv /= 0.0_sp) .and. maxval(abs(dv(:,1))) > 0.0_sp .and. maxval(abs(dv(:,2))) > 0.0_sp, &
                   sub//": dV_bsl of both domains is not zero")
        call check(maxval(abs(bsl)) > 0.0_sp, sub//": the reduced bsl moved")
        call check(worst <= 4.0_dp*real(epsilon(1.0_sp), dp), sub//": bsl budget closes: max |residual| / scale = "// &
                   trim(adjustl(num(worst)))//" (bound 4 eps of sp)")
    end subroutine check_budget

    function num(x) result(s)
        real(dp), intent(in) :: x
        character(len=24) :: s
        write(s, '(es10.3)') x
    end function num

    ! =========================================================================
    ! The run (child process, in the run directory)
    ! =========================================================================

    subroutine plain_child(par_file)
        character(len=*), intent(in) :: par_file
        type(run_t), target :: run
        call run_init(run, par_file)
        do while (.not. run_done(run))
            call run_step(run)
        end do
        call run_end(run)
        write(*,'(a)') "test_run_domains (plain): OK"
    end subroutine plain_child

    ! =========================================================================
    ! Inputs
    ! =========================================================================

    function topo_file(id) result(f)
        integer, intent(in) :: id
        character(len=:), allocatable :: f
        f = trim(root)//"/"//trim(DNAME(id))//"/"//trim(DGRID(id))//"/"//trim(DGRID(id))//"_"//trim(DTOPO(id))//".nc"
    end function topo_file

    subroutine setup_run_dir(sub)
        ! tmp/test_run_domains/<sub>: input/ (links), maps, the climate data (linked).
        character(len=*), intent(in) :: sub
        call sh("mkdir -p '"//d//sub//"/input'")
        call sh("ln -s '"//trim(repo)//"/maps' '"//d//sub//"/maps'")
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

    subroutine write_domain_par(file, id, bsl, time_init, time_end, restart, exps, dto)
        ! The single-domain par file of domain id as the converter writes one: the domain's instance groups carry
        ! its suffix, the physics and snapshot groups do not (they collide in a file with both domains: the merge
        ! step renames them). Library defaults except where the two domains differ (values of yelmox_Greenland.nml and
        ! yelmox_Antarctica.nml where they matter). `exps`: the experiment of the domain; `dto`: a term on the
        ! derived coefficient dTo (air anomaly) in the climate.
        character(len=*), intent(in) :: file, bsl, restart
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
            "&kryos", "  domains = '"//nm//"'", "  restart = '"//restart//"'", "  restart_out = 'out_rst'", &
            "  write_lossless = True", "  log_level = 'info'", "/", &
            "&clock", "  method = 'const'", "  time_const = 0.0", "  time_init = "//ftoa(time_init), &
            "  time_end = "//ftoa(time_end), "  dtt = 10.0", "/", &
            "&out_2D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_1D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_rst", "  method = 'const'", "  dt = "//ftoa(T_MID), "/", &
            "&globe", "  forcing = 'forcing'", "  geo = 'geo_glob'", "/", &
            "&forcing", "  backend = 'tsgen'", "  clock = 'forcing'", coefs, "  out_1D = 'out_1D'", "/", &
            "&geo_glob", "  backend = 'fastisostasy'", "  clock = 'forcing'", "  bsl = 'barysealevel'", &
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
            "&geo_"//sx, "  backend = 'fastisostasy'", "  dt = 10.0", "  isos = 'isos'", "  out_2D = 'out_2D'", &
            "  out_1D = 'out_1D'", "/", &
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
        ! FastIsostasy and the global bsl (constant ocean area)
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
            "&barysealevel", "  method = '"//bsl//"'", "  bsl_init = 0.0", "  sl_path = 'None'", "  sl_name = 'none'", &
            "  A_ocean_pd = 3.625e14", "  A_ocean_path = 'None'", "  restart = 'None'", "/"])
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

end program test_run_domains
