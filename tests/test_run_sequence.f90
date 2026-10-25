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
    ! (clock only, C10); and that a restore with instances stops (C10).
    !
    ! The run executes in tmp/test_run_sequence/ (its own input/ assembled from
    ! links: kryos input/, the yelmo and chion defaults, fesm-utils
    ! datasrc_defaults.nml, insol's orbital tables; maps/), in a child process
    ! of this executable:
    !     test_run_sequence.x run <par file>          (from the run directory)
    ! Built only with yelmo=1 fastisostasy=1 chion=1 (config/common.mk).

    use precision,       only: sp, dp
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_read, nc_exists_var
    use coords,          only: grid_class, grid_cdo_read_desc
    use fesm_fields,     only: field_state_t, fields_get
    use fesm_clock,      only: clock_time, cadence_t_last
    use kryos_component, only: component_t
    use kryos_domain,    only: domain_component, NCOMP
    use kryos_run,       only: run_t, run_init, run_step, run_done, run_end
    use kryos_ice_yelmo, only: yelmo_ice_t
    use kryos_climate_snapshots,  only: snapshots_climate_t
    use kryos_surface_chion,      only: chion_surface_t
    use kryos_shelf_marine_shelf, only: marshelf_shelf_t

    implicit none

    character(len=*), parameter :: DIR = "tmp/test_run_sequence/", PAR = "par.nml", DATA = "clim_grl.nc"
    integer, parameter :: NX = 106, NY = 181, NZ = 5

    character(len=512) :: arg, root, repo
    integer :: nfail = 0

    call get_command_argument(1, arg)
    if (trim(arg) == "run") then
        call get_command_argument(2, arg)
        call run_child(trim(arg))
        stop 0
    end if

    call get_environment_variable("PWD", repo)
    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = trim(repo)//"/ice_data"
    if (root(1:1) /= "/") root = trim(repo)//"/"//trim(root)
    if (.not. exists(trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc")) then
        write(*,'(a)') "test_run_sequence: ice_data topography not found under "//trim(root)//": skipped"
        stop 0
    end if

    call setup_run_dir()
    call check(child("run "//PAR, "run.log"), &
               "cold start + 2 steps + end with every phase 1 backend (child run, log "//DIR//"run.log)")
    call sh("grep -E '^  (ok|FAIL|OPEN)|^---' '"//trim(repo)//"/"//DIR//"run.log'")
    call check(.not. child("run par_restore.nml", "restore.log", quiet=.true.) .and. &
               in_file(trim(repo)//"/"//DIR//"restore.log", "milestone C10"), &
               "[kryos] restart with instances stops in run_restore (restore of the instances is C10)")

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

        call run_init(run, par_file)
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
        call run_end(run)

        call check_files()
        if (nfail > 0) then
            write(*,'(a,i0,a)') "test_run_sequence (run): FAILED (", nfail, " checks)"
            error stop 1
        end if
        write(*,'(a)') "test_run_sequence (run): OK"
    end subroutine run_child

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
            call check(all(a == real(real(0.05_dp, sp), dp)), &
                       "Q_geo = 0.05 W m-2 in the ice's precision (geothermal data, static row)")
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

            ! geo -> ice in sp: z_bed = ref + dz_bed, z_sl = bsl + dz_ss
            r = dom%refs(1)%z_bed
            a = get2(geo%exports, geo%h_dz_bed)
            b = get2(ice%imports, ice%hi_z_bed)
            call check(all(real(b, sp) == real(r, sp) + real(a, sp)), "ice z_bed = ref + dz_bed (sp, identity grid)")
            a = get2(geo%exports, geo%h_dz_ss)
            b = get2(ice%imports, ice%hi_z_sl)
            call check(all(real(b, sp) == real(bsl, sp) + real(a, sp)), "ice z_sl = bsl + dz_ss (sp)")

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
    end subroutine check_init

    subroutine round_trips(run)
        ! The two unit conversions across the exchange (adapters-design 1.2 P3,
        ! 10): pr [mm d-1] * sec_year/86400 out of the climate, / the same factor
        ! into chion; bmb_shlf [m ie yr-1] * rho_ice out of the shelf, / rho_ice
        ! into Yelmo. adapters-design 10 expects the backend's sp value back bit
        ! for bit on identity grids (the yelmox gate); the wiring narrows the dp
        ! product to the consumer's sp first (contracts-design 4.3, 9: round_to
        ! on every row, the lapse kernel in sp), so P3 does not hold. Open
        ! (C9 report): reported as OPEN with the count of differing cells, not
        ! failed, until the design decides which rule gives way.
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
        call note(npr == 0, "pr round trip climate -> chion bit for bit at t0 (values differing: "// &
                  trim(itoa(npr))//")")
        call note(nbmb == 0, "bmb_shlf round trip shelf -> Yelmo bit for bit at t0 (cells differing: "// &
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
            call check(all(get2(ice%imports, ice%hi_bmb_shlf) == real(real(bmb_prev, sp), dp)) .and. &
                       ice%imports%f(ice%hi_bmb_shlf)%t1 == t - 10.0_dp, &
                       "ice bmb_shlf = the shelf's bmb_shlf of the previous step (in the ice's sp)")

            ! geo ran at t; bsl of this step's reduction in the geo exports and z_sl
            bsl = get0(geo%exports, geo%h_bsl)
            call check(published_at(geo%exports, geo%h_dV_bsl, t) .and. cadence_t_last(geo%cadence, run%clk) == t, &
                       "geo updated at t (dt = 10)")
            call check(bsl == get0(run%globe%geo%c%exports, run%globe%geo%c%h_bsl) .and. &
                       geo%exports%f(geo%h_bsl)%t1 == t, "geo bsl = the reduced global bsl of this step")
            a = get2(geo%exports, geo%h_dz_ss)
            b = get2(ice%imports, ice%hi_z_sl)
            call check(all(real(b, sp) == real(bsl, sp) + real(a, sp)), "ice z_sl = bsl + dz_ss (sp)")

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
        call check(exists("restart-20/clock_restart.nc"), "restart bundle restart-20 (clock only until C10)")
        call nc_read("Greenland_hub.nc", "time", tt)
        call check(all(tt == [0.0_dp, 10.0_dp, 20.0_dp]), "hub file at t = 0, 10, 20 (output at the end of each step)")
        call check(nc_exists_var("forcing_ts.nc", "coef_dTo") .and. nc_exists_var("ice_grl_ts.nc", "M_ice"), &
                   "forcing coefficients and the ice's M_ice in the 1D files")
    end subroutine check_files

    ! =========================================================================
    ! Run directory, par files, data
    ! =========================================================================

    subroutine setup_run_dir()
        ! tmp/test_run_sequence/: input/ (links), maps, the climate data, par files.
        character(len=:), allocatable :: d

        d = trim(repo)//"/"//DIR
        call sh("rm -rf '"//d//"' && mkdir -p '"//d//"input'")
        call sh("ln -s '"//trim(repo)//"/maps' '"//d//"maps'")
        call sh("for f in '"//trim(repo)//"'/input/*; do ln -s ""$f"" '"//d//"input/'; done")
        call sh("for f in '"//trim(repo)//"'/yelmo/input/*.nml '"//trim(repo)//"'/yelmo/input/*.md; do "// &
                "ln -sf ""$f"" '"//d//"input/'; done")
        call sh("ln -sf '"//trim(repo)//"/chion/input/chion_defaults.nml' '"//d//"input/'")
        call sh("ln -sf '"//trim(repo)//"/fesm-utils/par/datasrc_defaults.nml' '"//d//"input/'")
        call sh("for f in '"//trim(repo)//"'/insol/input/LA2004.*.txt; do ln -sf ""$f"" '"//d//"input/'; done")
        call write_data(d//DATA)
        call write_par(d//PAR, "None")
        call write_par(d//"par_restore.nml", "restart-20")
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

    subroutine write_par(file, restart)
        ! The Greenland par file of the run: kryos groups (core-design 5), the
        ! hub's groups (ice_data), and the backend groups (yelmox_Greenland.nml
        ! values where they matter; library defaults otherwise).
        character(len=*), intent(in) :: file, restart

        character(len=:), allocatable :: topo, basins
        integer :: u

        topo   = trim(root)//"/{domain}/{grid_name}/{grid_name}_TOPO-M17-v5.nc"
        basins = trim(root)//"/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc"

        open(newunit=u, file=file, status="replace", action="write")
        call w(u, [character(len=400) :: &
            "&kryos", "  domains = 'Greenland'", "  restart = '"//restart//"'", "  restart_out = 'out_rst'", &
            "  write_lossless = True", "  log_level = 'info'", "/", &
            "&clock", "  method = 'const'", "  time_const = 0.0", "  time_init = 0.0", "  time_end = 20.0", &
            "  dtt = 10.0", "/", &
            "&out_2D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_1D", "  method = 'const'", "  dt = 10.0", "/", &
            "&out_rst", "  method = 'const'", "  dt = 20.0", "/", &
            "&globe", "  forcing = 'forcing'", "  geo = 'geo_glob'", "/", &
            "&forcing", "  backend = 'tsgen'", "  clock = 'forcing'", "  coefs = 'alpha' 'dTo' 'beta'", &
            "  out_1D = 'out_1D'", "/", &
            "&geo_glob", "  backend = 'fastisostasy'", "  clock = 'forcing'", "  bsl = 'barysealevel'", &
            "  out_1D = 'out_1D'", "/", &
            "&Greenland", "  grid_hub = 'GRL-16KM'", "  topo = 'topo_grl'", "  regions = 'regions_grl'", &
            "  basins = 'basin'", "  out_2D = 'out_2D'", "  ice = 'ice_grl'", "  geo = 'geo_grl'", &
            "  climate = 'climate_grl'", "  surface = 'surface_grl'", "  shelf = 'shelf_grl'", &
            "  geothermal = 'ghf_grl'", "  sediments = 'sed_grl'", "/", &
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
            "&shelf_grl", "  backend = 'marine_shelf'", "  marine_shelf = 'marine_shelf'", "  out_2D = 'out_2D'", "/", &
            "&ghf_grl", "  backend = 'data'", "  dt = -1.0", "  Q_geo = 'ghf_const'", "  out_2D = 'out_2D'", "/", &
            "&ghf_const", "  mode = 'const'", "  value = 0.05", "/", &
            "&sed_grl", "  backend = 'data'", "  dt = -1.0", "  H_sed = 'sed_none'", "/", &
            "&sed_none", "  mode = 'const'", "  value = 0.0", "/"])
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
        ! Yelmo: library defaults, the reference data loads off, masks without named regions
        call w(u, [character(len=400) :: &
            "&yelmo", "  nz_aa = 6", "/", "&ytopo", "/", "&ycalv", "/", "&ydyn", "  solver = 'diva'", "/", &
            "&ytill", "/", "&yhyd", "/", "&ymat", "/", "&ytrc", "/", "&ytherm", "/", &
            "&yelmo_masks", "  basins = 'None'", "  mask_ice_dynamic = 'all'", "  mask_rmse = 'all'", "/", &
            "&yelmo_init_topo", "/", &
            "&yelmo_data", "  pd_topo_load = False", "  pd_tsrf_load = False", "  pd_smb_load = False", &
            "  pd_vel_load = False", "  pd_age_load = False", "/"])
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
            "  basins = 'domain'", "  mask_ocean = 'all'", "  mask_deep_ocean = 'none'", "  mask_pico_deep = 'none'", &
            "  corr_method = 'none'", "  corr_names = ''", "  corr_values = 0.0", "  tf_correction = False", &
            "  tf_path = 'none'", "  tf_name = 'dT_nl'", "  bmb_max = 0.0", "  c_deep = -10.0", &
            "  depth_deep = 500.0", "  depth_const = 0.0", "  depth_min = 1100.0", "  depth_max = 1500.0", &
            "  lambda1 = -0.0575", "  lambda2 = 0.0832", "  lambda3 = 7.59e-4", "  gamma_lin = 630.0", &
            "  gamma_quad = 11100.0", "  gamma_quad_nl = 14500.0", "  gamma_prime = 100.0", "  kappa_grz = 10.0", &
            "  c_grz = 1.0", "  f_grz_shlf = 10.0", "  grz_length = 32.0", "  use_obs = False", &
            "  obs_path = 'none'", "  obs_name = 'bm_eq_reese'", "  obs_scale = 1.0", "  obs_lim = 100.0", "/"])
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

    logical function child(args, log, quiet)
        ! This executable with `args`, run in the run directory; .true. on exit 0.
        ! A failure shows the end of the log unless quiet.
        character(len=*), intent(in) :: args, log
        logical,          intent(in), optional :: quiet
        character(len=512) :: exe
        integer :: es
        logical :: q

        call get_command_argument(0, exe)
        if (exe(1:1) /= "/") exe = trim(repo)//"/"//trim(exe)
        call execute_command_line("cd '"//trim(repo)//"/"//DIR//"' && '"//trim(exe)//"' "//args//" > "//log//" 2>&1", &
                                  exitstat=es)
        child = es == 0
        q = .false.
        if (present(quiet)) q = quiet
        if (.not. (child .or. q)) call sh("tail -n 30 '"//trim(repo)//"/"//DIR//log//"'")
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
