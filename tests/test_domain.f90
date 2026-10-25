module test_domain_toys
    ! Minimal instances for the domain tests: they declare the vocabulary of
    ! their component through the real wrappers (ice_init, geo_init, shelf_init)
    ! and do nothing else; the test puts their exports by hand.
    use precision,       only: dp
    use kryos_component, only: init_ctx_t, step_time_t
    use kryos_ice,       only: ice_t
    use kryos_geo,       only: geo_t
    use kryos_surface,   only: surface_t
    implicit none
    private

    type, extends(ice_t), public :: toy_ice_t
    contains
        procedure :: init => toy_ice_init
        procedure :: update => toy_ice_update
        procedure :: restart_write => toy_ice_rw
        procedure :: restart_read => toy_ice_rr
    end type toy_ice_t

    type, extends(geo_t), public :: toy_geo_t
    contains
        procedure :: init => toy_geo_init
        procedure :: update => toy_geo_update
        procedure :: restart_write => toy_geo_rw
        procedure :: restart_read => toy_geo_rr
    end type toy_geo_t

    type, extends(surface_t), public :: toy_surface_t
    contains
        procedure :: init => toy_surface_init
        procedure :: update => toy_surface_update
        procedure :: restart_write => toy_surface_rw
        procedure :: restart_read => toy_surface_rr
    end type toy_surface_t

contains

    subroutine toy_ice_init(self, ctx)
        class(toy_ice_t),  intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
        self%opt%z_bed_ref = .true.
    end subroutine toy_ice_init
    subroutine toy_ice_update(self, st)
        class(toy_ice_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_ice_update
    subroutine toy_ice_rw(self, dir, t)
        class(toy_ice_t),  intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_ice_rw
    subroutine toy_ice_rr(self, dir, t)
        class(toy_ice_t),  intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_ice_rr

    subroutine toy_geo_init(self, ctx)
        class(toy_geo_t),  intent(inout) :: self
        type(init_ctx_t),  intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
    end subroutine toy_geo_init
    subroutine toy_geo_update(self, st)
        class(toy_geo_t),  intent(inout) :: self
        type(step_time_t), intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_geo_update
    subroutine toy_geo_rw(self, dir, t)
        class(toy_geo_t),  intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_geo_rw
    subroutine toy_geo_rr(self, dir, t)
        class(toy_geo_t),  intent(inout) :: self
        character(len=*),  intent(in)    :: dir
        real(dp),          intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_geo_rr

    subroutine toy_surface_init(self, ctx)
        class(toy_surface_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx
        self%nx = ctx%grid%G%nx
        self%ny = ctx%grid%G%ny
    end subroutine toy_surface_init
    subroutine toy_surface_update(self, st)
        class(toy_surface_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        associate (u => st%t1 + self%nx); end associate
    end subroutine toy_surface_update
    subroutine toy_surface_rw(self, dir, t)
        class(toy_surface_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_surface_rw
    subroutine toy_surface_rr(self, dir, t)
        class(toy_surface_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (u => t + self%nx + len(dir)); end associate
    end subroutine toy_surface_rr

end module test_domain_toys


program test_domain
    ! C6: kryos_domain without instances (core-design sections 2, 5, 6;
    ! contracts-design section 5), on real topography from ice_data (the gate
    ! inputs, FesmData, are on Levante only).
    !
    !   * grid chain: the default chain, an explicit grid, absent components;
    !   * the hub: loaded, and prepared by the [topo] keys (z_bed_f_sd,
    !     remove_lakes, H_ice_min, grad_lim_zb) -- unprepared equals the file bit
    !     for bit, each key does what it says;
    !   * refs: refs(1) the prepared hub reference widened; refs on a second grid
    !     equal the coupler's `con` remap of the hub-precision reference, widened,
    !     bit for bit (and the dp remap of the widened field rounded, core R1);
    !     only for grids in use;
    !   * the hub states published at t0; domain_hub_update: the mirror path
    !     (Greenland, hub = ice grid, fractional f_grnd kept), the anomaly path
    !     (Antarctica, ice on ANT-32KM, hub ANT-16KM: bilinear anomalies +
    !     htopo_update, against an independent computation), no ice (static hub);
    !   * exchange: ice -> hub and hub -> surface rows through domain_wire, and the
    !     lag error when the hub was not refreshed;
    !   * domain_check_reference: passes below 1 m, stops above, in stat form;
    !   * the hub output file (<domain>_hub.nc);
    !   * error cases through stat (missing group, unknown key, empty grid_hub,
    !     missing instance group, duplicate instance, missing topo group).
    !
    ! Toy instances (module test_domain_toys) stand for ice, geo and surface: they
    ! declare their vocabulary through ice_init/geo_init/shelf_init and the test
    ! puts the exports. Run from the repo root (as make check does); par files and
    ! outputs go to tmp/test_domain/. ice_data is read from $KRYOS_ICE_DATA,
    ! default ./ice_data; without the files the test reports and passes.
    !
    ! Driver mode for the bit-identity check against yelmox (Levante):
    !     test_domain.x dump <par file> <out dir> [domain]
    ! initialises the first (or the named) domain of the par file and writes
    ! <out dir>/<domain>_hub.nc (the hub at t0, lossless: yelmox htopo.nc layout
    ! plus f_ice, mask) and <out dir>/<domain>_refs.nc (every reference: z_bed,
    ! H_ice, z_srf, z_bed_sd on the hub grid and on each grid in use, as
    ! double, named <var>_<grid> on dims xc_<grid>, yc_<grid>).

    use, intrinsic :: iso_fortran_env, only: int32, int64
    use precision,       only: sp, dp, wp
    use coords,          only: grid_class, grid_cdo_read_desc
    use coupler,         only: coupler_class, coupler_init, coupler_method, remap
    use fesm_clock,      only: clock_t, clock_init, clock_start, clock_advance, clock_time, &
                               cadence_init, cadence_mark
    use fesm_fields,     only: fields_put, fields_has, fields_handle
    use fesm_log,        only: log_set_level, LOG_WARN
    use htopo,           only: htopo_class, htopo_update
    use ncio,            only: nc_read, nc_set_lossless, nc_exists_var, nc_create, nc_write_dim, nc_write
    use phys_constants,  only: phys_const_class, phys_const_set, phys_const_load
    use kryos_component, only: init_ctx_t, census_t, constraints_t, grid_ref_t, kryos_mask
    use kryos_config,    only: config_t, config_init
    use kryos_wiring,    only: wire_t, wire_node_t, wiring_table
    use kryos_ice,       only: ice_t, ice_adopt, ice_init
    use kryos_geo,       only: geo_t, geo_adopt, geo_init
    use kryos_surface,   only: surface_t, surface_adopt, surface_init
    use kryos_domain
    use test_domain_toys

    implicit none

    character(len=*), parameter :: dir = "tmp/test_domain/"
    character(len=*), parameter :: par = dir//"par.nml"
    integer :: nfail = 0
    character(len=256) :: root, arg
    type(phys_const_class) :: cnst
    type(config_t) :: cfg
    type(clock_t)  :: clk

    call get_command_argument(1, arg)
    if (trim(arg) == "dump") then
        call run_dump()
        stop 0
    end if

    call get_environment_variable("KRYOS_ICE_DATA", root)
    if (len_trim(root) == 0) root = "ice_data"
    if (.not. have_data()) then
        write(*,'(a)') "test_domain: ice_data topography not found under "//trim(root)//": skipped"
        stop 0
    end if

    call execute_command_line("mkdir -p "//dir)
    call log_set_level(LOG_WARN)
    call phys_const_set(cnst, label="Earth", rho_ice=910.0_dp, rho_sw=1028.0_dp, rho_w=1000.0_dp, g=9.81_dp)
    call write_par()
    call config_init(cfg, par)

    call test_errors()
    call test_greenland_hub()
    call test_hub_keys()
    call test_antarctica()
    call test_mirror_and_output()
    call test_anomaly_exchange_check()
    call test_no_ice()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_domain: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_domain: OK"

contains

    ! -------------------------------------------------------------------------
    ! Tests
    ! -------------------------------------------------------------------------

    subroutine test_errors()
        ! Par-file and consistency errors, through stat (all before the hub is loaded).
        type(domain_t) :: d
        integer :: stat

        write(*,"(a)") "errors (stat)"
        call domain_init(d, cfg, "NoSuchDomain", cnst, 0.0_dp, stat=stat)
        call check("missing domain group", stat /= 0)
        call domain_init(d, cfg, "BadKey", cnst, 0.0_dp, stat=stat)
        call check("unknown key in the domain group", stat /= 0)
        call domain_init(d, cfg, "NoHub", cnst, 0.0_dp, stat=stat)
        call check("grid_hub not set", stat /= 0)
        call domain_init(d, cfg, "NoInst", cnst, 0.0_dp, stat=stat)
        call check("instance group missing", stat /= 0)
        call domain_init(d, cfg, "DupInst", cnst, 0.0_dp, stat=stat)
        call check("one instance named by two components", stat /= 0)
        call domain_init(d, cfg, "NoTopo", cnst, 0.0_dp, stat=stat)
        call check("topo group missing from the par file", stat /= 0)
        call domain_init(d, cfg, "NoRegions", cnst, 0.0_dp, stat=stat)
        call check("regions group missing from the par file", stat /= 0)
    end subroutine test_errors

    subroutine test_greenland_hub()
        ! All instances on the default chain: one grid, one reference; the hub
        ! loaded, prepared and published at t0.
        type(domain_t) :: d
        real(sp), allocatable :: raw(:,:)
        integer :: k
        logical :: ok

        write(*,"(a)") "Greenland: chain, hub, references"
        call domain_init(d, cfg, "Greenland", cnst, 0.0_dp)

        ok = .true.
        do k = 1, NCOMP
            ok = ok .and. trim(d%grid(k)) == "GRL-16KM" .and. d%grid_default(k)
        end do
        call check("every component on the hub grid by the default chain", ok)
        call check("backends read from the instance groups", trim(d%backend(C_ICE)) == "yelmo" &
                   .and. trim(d%backend(C_SHELF)) == "marine_shelf" .and. trim(d%backend(C_SEDIMENTS)) == "data")
        call check("hub grid, shape", trim(d%grid_hub) == "GRL-16KM" .and. d%hub%nx == 106 .and. d%hub%ny == 181)
        call check("hub mirrors the ice (hub grid = ice grid)", d%hub_mirror)
        call check("one reference, on the hub grid", size(d%refs) == 1 .and. trim(d%refs(1)%grid) == "GRL-16KM")
        call check("refs(1) is the prepared hub reference widened", &
                   same_dp(d%refs(1)%z_bed, real(d%hub%ref%z_bed, dp)) .and. same_dp(d%refs(1)%H_ice, real(d%hub%ref%H_ice, dp)) &
                   .and. same_dp(d%refs(1)%z_srf, real(d%hub%ref%z_srf, dp)) &
                   .and. same_dp(d%refs(1)%z_bed_sd, real(d%hub%ref%z_bed_sd, dp)))
        call check("domain_ref_index", domain_ref_index(d, "GRL-16KM") == 1 .and. domain_ref_index(d, "ANT-16KM") == 0)

        ! the yelmox preparation keys of [topo_grl] (remove_lakes, H_ice_min = 1, grad_lim_zb = 0.5) are
        ! on; on this file at 16 km they change little (see test_hub_keys for stronger settings)
        call check("H_ice_min = 1: no ice thinner than 1 m in the hub reference", &
                   count(d%hub%ref%H_ice > 0.0_wp .and. d%hub%ref%H_ice < 1.0_wp) == 0)
        call read_raw(trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc", "z_srf", raw)
        call check("the surface of the data is kept by the preparation", same_sp(d%hub%ref%z_srf, raw))

        ! the hub at t0: the reference, published
        call check("hub geometry at t0 is the reference", same_sp(d%hub%z_bed, d%hub%ref%z_bed) &
                   .and. same_sp(d%hub%H_ice, d%hub%ref%H_ice) .and. same_sp(d%hub%z_srf, d%hub%ref%z_srf))
        associate (e => d%hub_exports, h => d%hh)
            call check("hub exports published at t0 (widened)", e%f(h%h_z_bed)%published .and. e%f(h%h_z_bed)%t1 == 0.0_dp &
                       .and. same_dp(e%f(h%h_z_bed)%v(:,:,1), real(d%hub%z_bed, dp)) &
                       .and. same_dp(e%f(h%h_H_ice)%v(:,:,1), real(d%hub%H_ice, dp)) &
                       .and. same_dp(e%f(h%h_f_grnd)%v(:,:,1), real(d%hub%f_grnd, dp)) &
                       .and. same_dp(e%f(h%h_z_srf_ref)%v(:,:,1), real(d%hub%ref%z_srf, dp)))
            call check("hub f_ice at t0: ice or no ice", all(e%f(h%h_f_ice)%v(:,:,1) == merge(1.0_dp, 0.0_dp, d%hub%H_ice > 0.0_wp)))
        end associate
        call check("hub_t = t0", d%hub_t == 0.0_dp)
    end subroutine test_greenland_hub

    subroutine test_hub_keys()
        ! The [topo] keys do what they say, on the same file.
        type(domain_t) :: d
        real(sp), allocatable :: zb(:,:), sd(:,:), hi(:,:), zs(:,:)
        character(len=*), parameter :: f = "/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc"
        logical :: okf

        write(*,"(a)") "the [topo] keys of the hub"
        call read_raw(trim(root)//f, "z_bed", zb)
        call read_raw(trim(root)//f, "z_bed_sd", sd)
        call read_raw(trim(root)//f, "H_ice", hi)
        call read_raw(trim(root)//f, "z_srf", zs)

        call domain_init(d, cfg, "GreenlandRaw", cnst, 0.0_dp)
        call check("no keys: the reference is the file (z_bed, H_ice, z_srf, z_bed_sd)", same_sp(d%hub%ref%z_bed, zb) &
                   .and. same_sp(d%hub%ref%H_ice, hi) .and. same_sp(d%hub%ref%z_srf, zs) .and. same_sp(d%hub%ref%z_bed_sd, sd))
        call check("no ice instance: no grids but the hub, no hub imports", size(d%refs) == 1 .and. d%hub_imports%n == 0)

        call domain_init(d, cfg, "GreenlandSD", cnst, 0.0_dp)
        okf = all(d%hub%ref%z_bed == zb + 0.5_sp*sd)
        call check("z_bed_f_sd = 0.5: z_bed + 0.5*z_bed_sd", okf)
        call check("... and nothing else changes (H_ice, z_srf)", same_sp(d%hub%ref%H_ice, hi) .and. same_sp(d%hub%ref%z_srf, zs))

        call domain_init(d, cfg, "GreenlandGrad", cnst, 0.0_dp)
        call check("grad_lim_zb = 0.02: the bed is smoothed (differs from the file), the surface kept", &
                   any(d%hub%ref%z_bed /= zb) .and. same_sp(d%hub%ref%z_srf, zs))

        call domain_init(d, cfg, "GreenlandMin", cnst, 0.0_dp)
        call check("H_ice_min = 5: ice thinner than 5 m removed, the rest kept", &
                   all(d%hub%ref%H_ice == merge(0.0_sp, hi, hi < 5.0_sp)) .and. same_sp(d%hub%ref%z_bed, zb))
        call check("... and the data had some", count(hi > 0.0_sp .and. hi < 5.0_sp) > 0)
    end subroutine test_hub_keys

    subroutine test_antarctica()
        ! Chain with an explicit grid; refs on a second grid.
        type(domain_t) :: d
        type(coupler_class) :: cpl2
        real(sp), allocatable :: s(:,:)
        real(dp), allocatable :: x(:,:), y(:,:)
        logical :: ok

        write(*,"(a)") "Antarctica: chain with an explicit grid, references on two grids"
        call domain_init(d, cfg, "Antarctica", cnst, 0.0_dp)

        call check("ice: its own grid, not default", trim(d%grid(C_ICE)) == "ANT-32KM" .and. .not. d%grid_default(C_ICE))
        call check("geo, climate, surface follow the ice", trim(d%grid(C_GEO)) == "ANT-32KM" .and. &
                   trim(d%grid(C_CLIMATE)) == "ANT-32KM" .and. trim(d%grid(C_SURFACE)) == "ANT-32KM")
        call check("shelf follows the hub", trim(d%grid(C_SHELF)) == "ANT-16KM" .and. d%grid_default(C_SHELF))
        call check("geothermal, sediments follow the ice", trim(d%grid(C_GEOTHERMAL)) == "ANT-32KM" .and. &
                   trim(d%grid(C_SEDIMENTS)) == "ANT-32KM")
        call check("ice grid /= hub grid: anomaly hub", .not. d%hub_mirror)
        call check("anomaly set of hub imports: dz_bed, dH_ice, z_sl", d%hub_imports%n == 3 .and. &
                   fields_has(d%hub_imports, "dz_bed") .and. fields_has(d%hub_imports, "dH_ice") .and. fields_has(d%hub_imports, "z_sl"))
        call check("refs: hub first, then the ice grid", size(d%refs) == 2 .and. trim(d%refs(1)%grid) == "ANT-16KM" &
                   .and. trim(d%refs(2)%grid) == "ANT-32KM")
        call check("domain_ref_index by grid", domain_ref_index(d, "ANT-32KM") == 2 .and. domain_ref_index(d, "GRL-16KM") == 0)

        call check("refs(1) is the hub reference widened", same_dp(d%refs(1)%z_bed, real(d%hub%ref%z_bed, dp)) &
                   .and. same_dp(d%refs(1)%H_ice, real(d%hub%ref%H_ice, dp)))

        ! refs(2) = the coupler's con remap in the hub's precision, widened, bit for bit
        call coupler_init(cpl2, "maps")
        call remap(cpl2, d%hub%ref%z_bed, "ANT-16KM", s, "ANT-32KM", method="con")
        call check("refs(2)%z_bed = con remap of the hub reference (sp path), widened", same_dp(d%refs(2)%z_bed, real(s, dp)))
        call remap(cpl2, d%hub%ref%H_ice, "ANT-16KM", s, "ANT-32KM", method="con")
        call check("refs(2)%H_ice", same_dp(d%refs(2)%H_ice, real(s, dp)))
        call remap(cpl2, d%hub%ref%z_srf, "ANT-16KM", s, "ANT-32KM", method="con")
        call check("refs(2)%z_srf", same_dp(d%refs(2)%z_srf, real(s, dp)))
        call remap(cpl2, d%hub%ref%z_bed_sd, "ANT-16KM", s, "ANT-32KM", method="con")
        call check("refs(2)%z_bed_sd", same_dp(d%refs(2)%z_bed_sd, real(s, dp)))
        ! the same through the dp path on the widened field, rounded (R1)
        call remap(cpl2, real(d%hub%ref%z_bed, dp), "ANT-16KM", x, "ANT-32KM", method="con")
        call remap(cpl2, real(d%hub%ref%H_ice, dp), "ANT-16KM", y, "ANT-32KM", method="con")
        ok = same_dp(d%refs(2)%z_bed, real(real(x, sp), dp)) .and. same_dp(d%refs(2)%H_ice, real(real(y, sp), dp))
        call check("refs(2) = the dp remap of the widened reference, rounded to sp (R1)", ok)
        call check("refs(2) is on the ice grid (shape 191 x 191)", all(shape(d%refs(2)%z_bed) == [191, 191]))
        call check("refs(2)%H_ice >= 0", all(d%refs(2)%H_ice >= 0.0_dp))

        ! the second domain: the chain broken at the climate, an explicit shelf grid, absent components
        call domain_init(d, cfg, "AntExplicit", cnst, 0.0_dp)
        call check("explicit climate grid overrides the chain (geo still follows the ice)", &
                   trim(d%grid(C_CLIMATE)) == "ANT-16KM" .and. .not. d%grid_default(C_CLIMATE) .and. &
                   trim(d%grid(C_GEO)) == "ANT-32KM" .and. d%grid_default(C_GEO))
        call check("the surface follows the climate, wherever it is", trim(d%grid(C_SURFACE)) == "ANT-16KM")
        call check("explicit shelf grid", trim(d%grid(C_SHELF)) == "ANT-32KM" .and. .not. d%grid_default(C_SHELF))
        call check("absent components: no grid, no instance", len_trim(d%grid(C_GEOTHERMAL)) == 0 .and. &
                   len_trim(d%grid(C_SEDIMENTS)) == 0 .and. len_trim(d%inst(C_SEDIMENTS)) == 0)
        call check("grids in use only: hub, then the ice grid (the others repeat them)", size(d%refs) == 2 .and. &
                   trim(d%refs(1)%grid) == "ANT-16KM" .and. trim(d%refs(2)%grid) == "ANT-32KM")
    end subroutine test_antarctica

    subroutine test_mirror_and_output()
        ! Greenland, toy ice on the hub grid: domain_hub_update mirrors the ice
        ! geometry (fractional f_ice and f_grnd kept); the hub output file.
        type(domain_t), target :: d
        type(wire_t), allocatable :: rows(:)
        real(sp), allocatable :: zb(:,:), hh(:,:), zs(:,:), fi(:,:), fg(:,:), zl(:,:)
        real(sp), allocatable :: v3(:,:,:)
        integer, allocatable :: m3(:,:,:)
        integer :: nx, ny, stat, i, j
        character(len=*), parameter :: file = dir//"Greenland_hub.nc"
        logical :: ok

        write(*,"(a)") "hub update, mirror path (hub grid = ice grid); hub output"
        call setup(d, "Greenland", ice=.true.)
        nx = d%hub%nx;  ny = d%hub%ny

        call execute_command_line("rm -f "//file)
        call nc_set_lossless(.true.)
        call domain_hub_write(d, file, 0.0_dp, .true.)

        rows = rows_between("ice", "hub")
        call domain_wire(d, clk, table=rows, check=.false., stat=stat)
        call check("domain_wire: the ice -> hub rows (mirror set + z_sl)", stat == 0 .and. count(d%wires%active) == 6)

        call clock_advance(clk)
        allocate(zb(nx,ny), hh(nx,ny), zs(nx,ny), fi(nx,ny), fg(nx,ny), zl(nx,ny))
        do j = 1, ny
            do i = 1, nx
                zb(i,j) = real(d%refs(1)%z_bed(i,j), sp) + 3.7_sp*sin(0.1_sp*i + 0.2_sp*j)
                hh(i,j) = max(real(d%refs(1)%H_ice(i,j), sp) + 11.3_sp*cos(0.05_sp*i - 0.1_sp*j), 0.0_sp)
                zs(i,j) = real(d%refs(1)%z_srf(i,j), sp) + 1.1_sp*sin(0.3_sp*i)
                fi(i,j) = merge(0.5_sp + 0.25_sp*sin(0.2_sp*i), 0.0_sp, hh(i,j) > 0.0_sp)
                fg(i,j) = 0.25_sp*real(mod(i + j, 5), sp)
                zl(i,j) = -2.5_sp + 0.01_sp*i
            end do
        end do
        call put_ice(d, 10.0_dp, zb, hh, zs, fi, fg, zl, real(d%refs(1)%z_bed, sp))
        call domain_hub_update(d, 10.0_dp, clk, stat=stat)
        call check("domain_hub_update runs (ice -> hub rows resolved and exchanged)", stat == 0)

        call check("mirror: hub z_bed, H_ice, z_srf are the ice's", same_sp(d%hub%z_bed, zb) .and. same_sp(d%hub%H_ice, hh) &
                   .and. same_sp(d%hub%z_srf, zs))
        call check("mirror: fractional f_grnd and z_sl kept", same_sp(d%hub%f_grnd, fg) .and. same_sp(d%hub%z_sl, zl) &
                   .and. any(d%hub%f_grnd > 0.0_sp .and. d%hub%f_grnd < 1.0_sp))
        associate (e => d%hub_exports, h => d%hh)
            ok = same_dp(e%f(h%h_z_bed)%v(:,:,1), real(zb, dp)) .and. same_dp(e%f(h%h_f_ice)%v(:,:,1), real(fi, dp)) &
                 .and. same_dp(e%f(h%h_f_grnd)%v(:,:,1), real(fg, dp)) .and. same_dp(e%f(h%h_z_sl)%v(:,:,1), real(zl, dp))
            call check("mirror: hub exports published at t = 10, incl. the ice's f_ice", ok .and. e%f(h%h_z_bed)%t1 == 10.0_dp)
            call check("z_srf_ref stays the reference (static)", e%f(h%h_z_srf_ref)%t1 == 0.0_dp .and. &
                       same_dp(e%f(h%h_z_srf_ref)%v(:,:,1), real(d%hub%ref%z_srf, dp)))
        end associate
        call check("hub_t = 10", d%hub_t == 10.0_dp)

        ! the ice did not step: the lag check stops the hub update
        call clock_advance(clk)
        call domain_hub_update(d, 20.0_dp, clk, stat=stat)
        call check("ice not updated at t = 20 (due, published at 10): the lag check stops the hub update", stat /= 0)

        ! hub output: two time levels, readable, lossless
        call domain_hub_write(d, file, 10.0_dp, .false.)
        allocate(v3(nx,ny,2), m3(nx,ny,2))
        call nc_read(file, "z_bed", v3)
        call check("hub file: z_bed at t = 0 (reference) and t = 10 (mirror)", same_sp(v3(:,:,1), d%hub%ref%z_bed) .and. same_sp(v3(:,:,2), zb))
        call nc_read(file, "H_ice", v3)
        call check("hub file: H_ice", same_sp(v3(:,:,1), d%hub%ref%H_ice) .and. same_sp(v3(:,:,2), hh))
        call nc_read(file, "f_grnd", v3)
        call check("hub file: f_grnd", same_sp(v3(:,:,2), fg))
        call nc_read(file, "f_ice", v3)
        call check("hub file: f_ice", same_sp(v3(:,:,2), fi))
        call nc_read(file, "mask", m3)
        call check("hub file: four-value mask from f_ice and f_grnd", all(m3(:,:,2) == kryos_mask(real(fi, dp), real(fg, dp))))
        call check("hub file: z_sl, z_srf, basin, time variables", nc_exists_var(file, "z_sl") .and. nc_exists_var(file, "z_srf") &
                   .and. nc_exists_var(file, "basin") .and. nc_exists_var(file, "time"))
        call nc_set_lossless(.false.)
    end subroutine test_mirror_and_output

    subroutine test_anomaly_exchange_check()
        ! AntChain: toy ice and geo on ANT-32KM, toy shelf on ANT-32KM, hub ANT-16KM.
        ! The anomaly path of domain_hub_update against an independent computation,
        ! hub -> shelf rows, and domain_check_reference.
        type(domain_t), target :: d
        type(coupler_class) :: cpl2
        type(htopo_class) :: hub2
        real(sp), allocatable :: zb(:,:), hh(:,:), zl(:,:), z0(:,:), h0(:,:), dz(:,:), dh(:,:)
        real(sp), allocatable :: dzh(:,:), dhh(:,:), zlh(:,:), a(:,:), b(:,:)
        real(sp) :: fi(191,191), fg(191,191), zs(191,191)
        integer :: nx, ny, stat, i, j, k
        logical :: ok
        character(len=8), parameter :: sh(4) = [character(len=8) :: "H_ice", "z_bed", "z_sl", "z_srf"]
        character(len=32) :: mres
        type(wire_t), allocatable :: rows(:)

        write(*,"(a)") "hub update, anomaly path (hub ANT-16KM, ice ANT-32KM); hub -> surface; reference check"
        call setup(d, "AntChain", ice=.true., geo=.true., surface=.true.)
        nx = 191;  ny = 191
        call coupler_init(cpl2, "maps")

        ! the wiring: ice -> hub and hub -> shelf rows (the toys have no other producers)
        rows = [rows_between("ice", "hub"), rows_between("hub", "surface")]
        call domain_wire(d, clk, table=rows, check=.false., stat=stat)
        call check("domain_wire: rows resolved and primed", stat == 0 .and. size(d%wires) == size(rows))
        call check("rows into the hub: anomaly set active, mirror set inactive", count(d%wires%active .and. d%wires%dst_kind == "hub") == 3)
        call domain_wire(d, clk, stat=stat)
        call check("domain_wire with the full table stops on the toy ice's unwired imports", stat /= 0)
        call domain_wire(d, clk, table=rows, check=.false., stat=stat)

        ! step 1 (t = 10): anomalies against refs(2)
        z0 = real(d%refs(2)%z_bed, sp);  h0 = real(d%refs(2)%H_ice, sp)
        allocate(zb(nx,ny), hh(nx,ny), zl(nx,ny))
        do j = 1, ny
            do i = 1, nx
                zb(i,j) = z0(i,j) + 7.3_sp*sin(0.07_sp*i + 0.11_sp*j)
                hh(i,j) = max(h0(i,j) + 13.1_sp*cos(0.05_sp*i - 0.09_sp*j), 0.0_sp)
                zl(i,j) = 1.5_sp + 0.002_sp*j
            end do
        end do
        fi = 1.0_sp;  fg = 1.0_sp;  zs = 0.0_sp
        call clock_advance(clk)
        call put_ice(d, 10.0_dp, zb, hh, zs, fi, fg, zl, z0)
        hub2 = d%hub
        call domain_hub_update(d, 10.0_dp, clk, stat=stat)
        call check("domain_hub_update runs", stat == 0)

        ! expected: anomalies in sp, refined bilinearly (sp path), htopo_update
        dz = zb - z0
        dh = hh - h0
        call remap(cpl2, dz, "ANT-32KM", dzh, "ANT-16KM", method="bilin")
        call remap(cpl2, dh, "ANT-32KM", dhh, "ANT-16KM", method="bilin")
        call remap(cpl2, zl, "ANT-32KM", zlh, "ANT-16KM", method="bilin")
        call htopo_update(hub2, dzh, dhh, zlh)
        call check("anomaly: hub z_bed = reference + bilinear dz_bed (bitwise)", same_sp(d%hub%z_bed, hub2%z_bed))
        call check("anomaly: hub H_ice = max(reference + bilinear dH_ice, 0), z_srf, f_grnd, z_sl", &
                   same_sp(d%hub%H_ice, hub2%H_ice) .and. same_sp(d%hub%z_srf, hub2%z_srf) &
                   .and. same_sp(d%hub%f_grnd, hub2%f_grnd) .and. same_sp(d%hub%z_sl, hub2%z_sl))
        call check("anomaly: the hub is not the ice (finer: reference present)", any(d%hub%z_bed /= d%hub%ref%z_bed))
        call check("anomaly: f_grnd is 0 or 1, H_ice >= 0", all(d%hub%f_grnd == 0.0_sp .or. d%hub%f_grnd == 1.0_sp) &
                   .and. all(d%hub%H_ice >= 0.0_sp))
        associate (e => d%hub_exports, h => d%hh)
            call check("anomaly: f_ice = ice or no ice, exports at t = 10", e%f(h%h_f_ice)%t1 == 10.0_dp .and. &
                       all(e%f(h%h_f_ice)%v(:,:,1) == merge(1.0_dp, 0.0_dp, d%hub%H_ice > 0.0_sp)) .and. &
                       same_dp(e%f(h%h_H_ice)%v(:,:,1), real(d%hub%H_ice, dp)))
        end associate

        ! hub -> surface (auto, ANT-16KM -> ANT-32KM)
        call exchange(d, "hub", "surface", 10.0_dp, clk, stat=stat)
        call check("exchange hub -> surface", stat == 0)
        mres = coupler_method(cpl2, "ANT-16KM", "ANT-32KM", "auto")
        ok = .true.
        do k = 1, size(sh)
            select case (trim(sh(k)))
            case ("H_ice");  a = d%hub%H_ice
            case ("z_bed");  a = d%hub%z_bed
            case ("z_sl");   a = d%hub%z_sl
            case ("z_srf");  a = d%hub%z_srf
            end select
            call remap(cpl2, a, "ANT-16KM", b, "ANT-32KM", method=trim(mres))
            ok = ok .and. same_dp(real(real(d%surface%c%imports%f(fields_handle(d%surface%c%imports, trim(sh(k))))%v(:,:,1), &
                                            sp), dp), real(b, dp))
        end do
        call check("surface imports, narrowed = the sp remap ("//trim(mres)//") of the hub exports (4 rows, bitwise)", ok)

        ! a hub that was not refreshed this step is caught by the lag check
        call clock_advance(clk)
        call put_ice(d, 20.0_dp, zb, hh, zs, fi, fg, zl, z0)
        call exchange(d, "hub", "surface", 20.0_dp, clk, stat=stat)
        call check("hub not refreshed at t = 20: exchange hub -> surface stops (lag check)", stat /= 0)

        ! zero anomalies: the hub is the reference
        call put_ice(d, 20.0_dp, z0, h0, zs, fi, fg, 0.0_sp*zl, z0)
        call domain_hub_update(d, 20.0_dp, clk, stat=stat)
        call check("zero anomalies: hub z_bed and H_ice are the reference", stat == 0 .and. same_sp(d%hub%z_bed, d%hub%ref%z_bed) &
                   .and. same_sp(d%hub%H_ice, d%hub%ref%H_ice))
        call exchange(d, "hub", "surface", 20.0_dp, clk, stat=stat)
        call check("exchange hub -> surface after the refresh", stat == 0)
        call exchange(d, "geo", "ice", 20.0_dp, clk, stat=stat)
        call check("a pair without active rows exchanges nothing", stat == 0)

        ! ---- reference check -------------------------------------------------
        write(*,"(a)") "domain_check_reference"
        call domain_check_reference(d, stat=stat)
        call check("geo z_bed_ref not published: stop", stat /= 0)
        call fields_put(d%geo%c%exports, d%geo%c%h_z_bed_ref, real(d%refs(2)%z_bed, dp), t1=20.0_dp)
        call fields_put(d%ice%c%exports, d%ice%c%h_z_bed_ref, real(d%refs(2)%z_bed, dp), t1=20.0_dp)
        call domain_check_reference(d, stat=stat)
        call check("geo and ice reference = domain reference: ok", stat == 0)
        call fields_put(d%geo%c%exports, d%geo%c%h_z_bed_ref, real(d%refs(2)%z_bed, dp) + 0.9_dp, t1=20.0_dp)
        call domain_check_reference(d, stat=stat)
        call check("geo reference 0.9 m off: within the 1 m tolerance", stat == 0)
        a = real(d%refs(2)%z_bed, sp)
        a(50,60) = a(50,60) + 1.5_sp
        call fields_put(d%geo%c%exports, d%geo%c%h_z_bed_ref, a, t1=20.0_dp)
        call domain_check_reference(d, stat=stat)
        call check("geo reference 1.5 m off in one cell: stop", stat /= 0)
        call fields_put(d%geo%c%exports, d%geo%c%h_z_bed_ref, real(d%refs(2)%z_bed, dp), t1=20.0_dp)
        call fields_put(d%ice%c%exports, d%ice%c%h_z_bed_ref, real(d%refs(2)%z_bed, dp) - 2.0_dp, t1=20.0_dp)
        call domain_check_reference(d, stat=stat)
        call check("ice reference 2 m off: stop", stat /= 0)
    end subroutine test_anomaly_exchange_check

    subroutine test_no_ice()
        ! A domain without an ice instance: the hub stays the reference, republished.
        type(domain_t), target :: d
        integer :: stat

        write(*,"(a)") "hub update without an ice instance"
        call setup(d, "GreenlandRaw", ice=.false.)
        call clock_advance(clk)
        call domain_hub_update(d, 10.0_dp, clk, stat=stat)
        call check("domain_hub_update: static hub", stat == 0 .and. d%hub_t == 10.0_dp .and. &
                   d%hub_exports%f(d%hh%h_z_bed)%t1 == 10.0_dp .and. same_sp(d%hub%z_bed, d%hub%ref%z_bed))
        call exchange(d, "hub", "surface", 10.0_dp, clk, stat=stat)
        call check("exchange before domain_wire stops", stat /= 0)
        call domain_init(d, cfg, "Greenland", cnst, 0.0_dp)
        call domain_hub_update(d, 0.0_dp, clk, stat=stat)
        call check("an ice instance configured but not created is refused", stat /= 0)
    end subroutine test_no_ice

    ! -------------------------------------------------------------------------
    ! Toy domains
    ! -------------------------------------------------------------------------

    subroutine setup(d, name, ice, geo, surface)
        ! domain_init, a new clock (time_init = 0, dtt = 10) and the requested
        ! toy instances created through the adopt/init wrappers, cadences registered.
        type(domain_t),   intent(inout) :: d
        character(len=*), intent(in)    :: name
        logical,          intent(in), optional :: ice, geo, surface

        class(ice_t),     allocatable :: p_ice     ! prototypes (the adopt wrappers take
        class(geo_t),     allocatable :: p_geo     ! polymorphic allocatables)
        class(surface_t), allocatable :: p_surf
        type(census_t) :: census
        type(init_ctx_t) :: ctx

        call domain_init(d, cfg, name, cnst, 0.0_dp)
        call clock_init(clk, 0.0_dp, 100.0_dp, 10.0_dp, "365_day", "tropical")

        if (present(ice)) then
            if (ice) then
                call toy_ctx(ctx, d, C_ICE, "ice")
                allocate(toy_ice_t :: p_ice)
                call ice_adopt(d%ice, p_ice, constraints_t(), ctx, census)
                call ice_init(d%ice, ctx)
                call cadence_init(d%ice%c%cadence, clk, "ice", .true.)
            end if
        end if
        if (present(geo)) then
            if (geo) then
                call toy_ctx(ctx, d, C_GEO, "geo")
                allocate(toy_geo_t :: p_geo)
                call geo_adopt(d%geo, p_geo, constraints_t(), ctx, census)
                call geo_init(d%geo, ctx)
                call cadence_init(d%geo%c%cadence, clk, "geo", 10.0_dp)
            end if
        end if
        if (present(surface)) then
            if (surface) then
                call toy_ctx(ctx, d, C_SURFACE, "surface")
                allocate(toy_surface_t :: p_surf)
                call surface_adopt(d%surface, p_surf, constraints_t(), ctx, census)
                call surface_init(d%surface, ctx)
                call cadence_init(d%surface%c%cadence, clk, "surface", .true.)
            end if
        end if
        call clock_start(clk)
    end subroutine setup

    subroutine toy_ctx(ctx, d, k, comp)
        type(init_ctx_t), intent(out) :: ctx
        type(domain_t),   intent(in)  :: d
        integer,          intent(in)  :: k
        character(len=*), intent(in)  :: comp

        ctx%instance  = d%inst(k)
        ctx%component = comp
        ctx%backend   = "toy"
        ctx%domain    = d%name
        ctx%grid_name = d%grid(k)
        call grid_cdo_read_desc(ctx%grid, trim(ctx%grid_name), "maps")
    end subroutine toy_ctx

    subroutine put_ice(d, t, zb, h, zs, fi, fg, zl, zbref)
        ! The toy ice's exports at t (sp values widened), then the cadence mark
        ! (what ice_update and the sequence do).
        type(domain_t), intent(inout) :: d
        real(dp),       intent(in)    :: t
        real(sp),       intent(in)    :: zb(:,:), h(:,:), zs(:,:), fi(:,:), fg(:,:), zl(:,:), zbref(:,:)

        associate (c => d%ice%c)
            call fields_put(c%exports, c%h_z_bed,  zb,  t1=t)
            call fields_put(c%exports, c%h_H_ice,  h,   t1=t)
            call fields_put(c%exports, c%h_z_srf,  zs,  t1=t)
            call fields_put(c%exports, c%h_f_ice,  fi,  t1=t)
            call fields_put(c%exports, c%h_f_grnd, fg,  t1=t)
            call fields_put(c%exports, c%h_z_sl,   zl,  t1=t)
            call fields_put(c%exports, c%h_z_bed_ref, zbref, t1=t)
            call cadence_mark(c%cadence, clk)
        end associate
    end subroutine put_ice

    function rows_between(src, dst) result(r)
        character(len=*), intent(in) :: src, dst
        type(wire_t), allocatable :: r(:)
        type(wire_t), allocatable :: t(:)
        integer :: i
        t = wiring_table()
        r = pack(t, [(trim(t(i)%src_kind) == trim(src) .and. trim(t(i)%dst_kind) == trim(dst), i = 1, size(t))])
    end function rows_between

    ! -------------------------------------------------------------------------
    ! Par file
    ! -------------------------------------------------------------------------

    subroutine write_par()
        integer :: u
        character(len=:), allocatable :: g, a

        g = trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc"
        a = trim(root)//"/Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine.nc"

        open(newunit=u, file=par, status="replace", action="write")
        write(u,'(a)') "&kryos"
        write(u,'(a)') "  domains = 'Greenland' 'Antarctica'"
        write(u,'(a)') "  log_level = 'warn'"
        write(u,'(a)') "/"
        ! --- Greenland: every component on the default chain; hub prepared as in yelmox
        write(u,'(a)') "&Greenland"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_grl'"
        write(u,'(a)') "  regions = 'regions_grl'"
        write(u,'(a)') "  basins = 'basin'"
        write(u,'(a)') "  out_2D = 'out_2D'"
        call write_instances(u, "ice_grl", "geo_grl", "climate_grl", "surface_grl", "shelf_grl", "ghf_grl", "sed_grl")
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_grl"
        write(u,'(a)') "  path = '"//trim(root)//"/{domain}/{grid_name}/{grid_name}_TOPO-M17-v5.nc'"
        write(u,'(a)') "  z_bed_f_sd = 0.0"
        write(u,'(a)') "  remove_lakes = .true."
        write(u,'(a)') "  H_ice_min = 1.0"
        write(u,'(a)') "  grad_lim_zb = 0.5"
        write(u,'(a)') "/"
        write(u,'(a)') "&regions_grl"
        write(u,'(a)') "  layers = 'basin'"
        write(u,'(a)') "  path_basin = '"//trim(root)//"/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc'"
        write(u,'(a)') "/"
        ! --- Greenland variants for the [topo] keys; no ice instance
        call write_hub_domain(u, "GreenlandRaw", "GRL-16KM", "topo_raw", "regions_lit", "")
        call write_hub_domain(u, "GreenlandSD",  "GRL-16KM", "topo_sd",  "regions_lit", "")
        call write_hub_domain(u, "GreenlandMin", "GRL-16KM", "topo_min", "regions_lit", "")
        call write_hub_domain(u, "GreenlandGrad", "GRL-16KM", "topo_grad", "regions_lit", "")
        write(u,'(a)') "&regions_lit"
        write(u,'(a)') "  layers = 'basin'"
        write(u,'(a)') "  path_basin = '"//trim(root)//"/Greenland/GRL-16KM/GRL-16KM_BASINS-nasa.nc'"
        write(u,'(a)') "/"
        write(u,'(a)') "&regions_antlit"
        write(u,'(a)') "  layers = 'basin'"
        write(u,'(a)') "  path_basin = '"//trim(root)//"/Antarctica/ANT-16KM/ANT-16KM_BASINS-nasa.nc'"
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_raw"
        write(u,'(a)') "  path = '"//g//"'"
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_sd"
        write(u,'(a)') "  path = '"//g//"'"
        write(u,'(a)') "  z_bed_f_sd = 0.5"
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_grad"
        write(u,'(a)') "  path = '"//g//"'"
        write(u,'(a)') "  grad_lim_zb = 0.02"
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_min"
        write(u,'(a)') "  path = '"//g//"'"
        write(u,'(a)') "  H_ice_min = 5.0"
        write(u,'(a)') "/"
        ! --- Antarctica: ice on ANT-32KM, the rest on the chain, shelf on the hub
        write(u,'(a)') "&Antarctica"
        write(u,'(a)') "  grid_hub = 'ANT-16KM'"
        write(u,'(a)') "  topo = 'topo_ant'"
        write(u,'(a)') "  regions = 'regions_ant'"
        write(u,'(a)') "  basins = 'basin'"
        call write_instances(u, "ice_ant", "geo_ant", "climate_ant", "surface_ant", "shelf_ant", "ghf_ant", "sed_ant")
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_ant"
        write(u,'(a)') "  path = '"//trim(root)//"/{domain}/{grid_name}/{grid_name}_TOPO-BedMachine.nc'"
        write(u,'(a)') "  z_bed_f_sd = 0.0"
        write(u,'(a)') "  remove_lakes = .true."
        write(u,'(a)') "  H_ice_min = 1.0"
        write(u,'(a)') "  grad_lim_zb = 0.5"
        write(u,'(a)') "/"
        write(u,'(a)') "&regions_ant"
        write(u,'(a)') "  layers = 'basin'"
        write(u,'(a)') "  path_basin = '"//trim(root)//"/{domain}/{grid_name}/{grid_name}_BASINS-nasa.nc'"
        write(u,'(a)') "/"
        ! --- AntChain: explicit surface and shelf grids, no geothermal and sediments
        write(u,'(a)') "&AntChain"
        write(u,'(a)') "  grid_hub = 'ANT-16KM'"
        write(u,'(a)') "  topo = 'topo_chain'"
        write(u,'(a)') "  regions = 'regions_antlit'"
        write(u,'(a)') "  ice = 'ice_ant'"
        write(u,'(a)') "  geo = 'geo_ant'"
        write(u,'(a)') "  climate = 'climate_ant'"
        write(u,'(a)') "  surface = 'surface_ant'"
        write(u,'(a)') "/"
        ! --- AntExplicit: the chain broken at the climate, explicit shelf grid, no geothermal and sediments
        write(u,'(a)') "&AntExplicit"
        write(u,'(a)') "  grid_hub = 'ANT-16KM'"
        write(u,'(a)') "  topo = 'topo_chain'"
        write(u,'(a)') "  regions = 'regions_antlit'"
        write(u,'(a)') "  ice = 'ice_ant'"
        write(u,'(a)') "  geo = 'geo_ant'"
        write(u,'(a)') "  climate = 'climate_ant16'"
        write(u,'(a)') "  surface = 'surface_ant'"
        write(u,'(a)') "  shelf = 'shelf_ant32'"
        write(u,'(a)') "/"
        write(u,'(a)') "&topo_chain"
        write(u,'(a)') "  path = '"//a//"'"
        write(u,'(a)') "  H_ice_min = 1.0"
        write(u,'(a)') "/"
        ! --- instance groups
        call write_instance(u, "ice_grl", "yelmo")
        call write_instance(u, "geo_grl", "fastisostasy")
        call write_instance(u, "climate_grl", "snapshots")
        call write_instance(u, "surface_grl", "chion")
        call write_instance(u, "shelf_grl", "marine_shelf")
        call write_instance(u, "ghf_grl", "data")
        call write_instance(u, "sed_grl", "data")
        call write_instance(u, "ice_ant", "yelmo", "ANT-32KM")
        call write_instance(u, "geo_ant", "fastisostasy")
        call write_instance(u, "climate_ant", "snapshots")
        call write_instance(u, "surface_ant", "chion")
        call write_instance(u, "shelf_ant", "marine_shelf")
        call write_instance(u, "ghf_ant", "data")
        call write_instance(u, "sed_ant", "data")
        call write_instance(u, "climate_ant16", "snapshots", "ANT-16KM")
        call write_instance(u, "shelf_ant32", "marine_shelf", "ANT-32KM")
        ! --- error cases
        write(u,'(a)') "&BadKey"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_raw'"
        write(u,'(a)') "  bogus = 1"
        write(u,'(a)') "/"
        write(u,'(a)') "&NoHub"
        write(u,'(a)') "  topo = 'topo_raw'"
        write(u,'(a)') "/"
        write(u,'(a)') "&NoInst"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_raw'"
        write(u,'(a)') "  regions = 'regions_grl'"
        write(u,'(a)') "  ice = 'ice_nope'"
        write(u,'(a)') "/"
        write(u,'(a)') "&DupInst"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_raw'"
        write(u,'(a)') "  regions = 'regions_grl'"
        write(u,'(a)') "  ice = 'ice_grl'"
        write(u,'(a)') "  geo = 'ice_grl'"
        write(u,'(a)') "/"
        write(u,'(a)') "&NoTopo"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_missing'"
        write(u,'(a)') "  regions = 'regions_grl'"
        write(u,'(a)') "/"
        write(u,'(a)') "&NoRegions"
        write(u,'(a)') "  grid_hub = 'GRL-16KM'"
        write(u,'(a)') "  topo = 'topo_raw'"
        write(u,'(a)') "  regions = 'regions_missing'"
        write(u,'(a)') "/"
        close(u)
    end subroutine write_par

    subroutine write_hub_domain(u, name, grid, topo, regions, basins)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, grid, topo, regions, basins
        write(u,'(a)') "&"//name
        write(u,'(a)') "  grid_hub = '"//grid//"'"
        write(u,'(a)') "  topo = '"//topo//"'"
        write(u,'(a)') "  regions = '"//regions//"'"
        if (len_trim(basins) > 0) write(u,'(a)') "  basins = '"//basins//"'"
        write(u,'(a)') "/"
    end subroutine write_hub_domain

    subroutine write_instances(u, ice, geo, climate, surface, shelf, geothermal, sediments)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: ice, geo, climate, surface, shelf, geothermal, sediments
        write(u,'(a)') "  ice = '"//ice//"'"
        write(u,'(a)') "  geo = '"//geo//"'"
        write(u,'(a)') "  climate = '"//climate//"'"
        write(u,'(a)') "  surface = '"//surface//"'"
        write(u,'(a)') "  shelf = '"//shelf//"'"
        write(u,'(a)') "  geothermal = '"//geothermal//"'"
        write(u,'(a)') "  sediments = '"//sediments//"'"
    end subroutine write_instances

    subroutine write_instance(u, name, backend, grid)
        integer,          intent(in) :: u
        character(len=*), intent(in) :: name, backend
        character(len=*), intent(in), optional :: grid
        write(u,'(a)') "&"//name
        write(u,'(a)') "  backend = '"//backend//"'"
        if (present(grid)) write(u,'(a)') "  grid = '"//grid//"'"
        write(u,'(a)') "/"
    end subroutine write_instance

    ! -------------------------------------------------------------------------
    ! Driver mode
    ! -------------------------------------------------------------------------

    subroutine run_dump()
        ! test_domain.x dump <par file> <out dir> [domain]
        character(len=512) :: parf, outd, dname
        type(config_t) :: c
        type(phys_const_class) :: cn
        type(domain_t) :: d
        type(grid_class) :: g
        character(len=:), allocatable :: f, gn
        integer :: i

        call get_command_argument(2, parf)
        call get_command_argument(3, outd)
        call get_command_argument(4, dname)
        if (len_trim(parf) == 0 .or. len_trim(outd) == 0) then
            write(*,'(a)') "usage: test_domain.x dump <par file> <out dir> [domain]"
            error stop 2
        end if

        call config_init(c, trim(parf))
        call log_set_level(c%log_level)
        if (len_trim(dname) == 0) then
            if (size(c%domains) == 0) error stop "test_domain dump: no domain in [kryos] domains"
            dname = c%domains(1)
        end if
        call phys_const_load(cn, trim(c%phys_const_file), group=trim(c%phys_const), &
                             defaults_file="fesm-utils/par/phys_const_earth.nml", defaults_group="Earth")
        call domain_init(d, c, trim(dname), cn, 0.0_dp)

        call execute_command_line("mkdir -p "//trim(outd))
        call nc_set_lossless(.true.)
        f = trim(outd)//"/"//trim(dname)//"_hub.nc"
        call execute_command_line("rm -f "//f)
        call domain_hub_write(d, f, 0.0_dp, .true.)

        f = trim(outd)//"/"//trim(dname)//"_refs.nc"
        call execute_command_line("rm -f "//f)
        call nc_create(f)
        do i = 1, size(d%refs)
            gn = sanitize(d%refs(i)%grid)
            call grid_cdo_read_desc(g, trim(d%refs(i)%grid), trim(c%map_fldr))
            call nc_write_dim(f, "xc_"//gn, x=g%G%x, units="km")
            call nc_write_dim(f, "yc_"//gn, x=g%G%y, units="km")
            call nc_write(f, "z_bed_"//gn,    d%refs(i)%z_bed,    dim1="xc_"//gn, dim2="yc_"//gn, units="m")
            call nc_write(f, "H_ice_"//gn,    d%refs(i)%H_ice,    dim1="xc_"//gn, dim2="yc_"//gn, units="m")
            call nc_write(f, "z_srf_"//gn,    d%refs(i)%z_srf,    dim1="xc_"//gn, dim2="yc_"//gn, units="m")
            call nc_write(f, "z_bed_sd_"//gn, d%refs(i)%z_bed_sd, dim1="xc_"//gn, dim2="yc_"//gn, units="m")
        end do
        write(*,'(a)') "test_domain dump: wrote "//trim(outd)//"/"//trim(dname)//"_hub.nc and _refs.nc"
    end subroutine run_dump

    ! -------------------------------------------------------------------------
    ! Helpers
    ! -------------------------------------------------------------------------

    logical function have_data()
        logical :: a, b, c, d, e
        inquire(file=trim(root)//"/Greenland/GRL-16KM/GRL-16KM_TOPO-M17-v5.nc", exist=a)
        inquire(file=trim(root)//"/Greenland/GRL-16KM/GRL-16KM_BASINS-nasa.nc", exist=b)
        inquire(file=trim(root)//"/Antarctica/ANT-16KM/ANT-16KM_TOPO-BedMachine.nc", exist=c)
        inquire(file=trim(root)//"/Antarctica/ANT-16KM/ANT-16KM_BASINS-nasa.nc", exist=d)
        inquire(file="maps/grid_ANT-32KM.txt", exist=e)
        have_data = a .and. b .and. c .and. d .and. e
    end function have_data

    subroutine read_raw(file, var, v)
        character(len=*),      intent(in)  :: file, var
        real(sp), allocatable, intent(out) :: v(:,:)
        allocate(v(106, 181))
        call nc_read(file, var, v)
    end subroutine read_raw

    function sanitize(name) result(s)
        ! A grid name as a NetCDF name suffix: "ANT-32KM" -> "ANT_32KM".
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: s
        integer :: i
        s = trim(name)
        do i = 1, len(s)
            if (s(i:i) == "-") s(i:i) = "_"
        end do
    end function sanitize

    logical function same_sp(a, b)
        ! Equal as 32-bit patterns.
        real(sp), intent(in) :: a(:,:), b(:,:)
        same_sp = .false.
        if (any(shape(a) /= shape(b))) return
        same_sp = all(transfer(a, 1_int32, size(a)) == transfer(b, 1_int32, size(b)))
    end function same_sp

    logical function same_dp(a, b)
        ! Equal as 64-bit patterns.
        real(dp), intent(in) :: a(:,:), b(:,:)
        same_dp = .false.
        if (any(shape(a) /= shape(b))) return
        same_dp = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
    end function same_dp

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

end program test_domain
