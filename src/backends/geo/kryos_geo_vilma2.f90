module kryos_geo_vilma2
    ! geo backend VILMA2 (phase 4 plan section 3), built with vilma=1, else the
    ! stub: vilma_earth_t, the global part (scope global: owns the solid_earth on
    ! its Gauss grid GAUSS-L<lmax>, resolved on the host's loads) and vilma_geo_t,
    ! the regional port (stateless: it holds the anomalies the core scatters to
    ! the domain's geo grid). The factory picks by scope.
    !
    ! VILMA2 computes in real64 (prec = dp): nothing is narrowed. The host terms
    ! (w, I_g, C, F, the external ice, the reference) are the imports the core
    ! fills (gather, globe_ice); the adapter only hands them to solid_earth_* and
    ! puts what it returns. The densities are the run's (phys_const), not keys.

    use precision,       only: dp
    use coordinates,     only: grid_class, grid_init_gaussian
    use coords,          only: grid_cdo_read_desc, grid_cdo_write_desc_short
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_has_param
    use ncio,            only: nc_open, nc_write, nc_close
    use variable_io,     only: var_io_type, var_io_write
    use kryos_outfile,   only: outfile_dims_lonlat, outfile_open, outfile_row, outfile_table_file
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, NAME_LEN, PATH_LEN, &
                               PASS_GEOMETRY, PASS_STATE, component_write_output
    use kryos_geo,       only: geo_t, geo_opt_t
    use vilma_params,    only: vilma_par_load
    use vilma_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_spinup, &
                               solid_earth_finalize, solid_earth_host_check
    use vilma_io,        only: vilma_restart_write, vilma_restart_read, vilma_io_set_table, vilma_write_step
    use vilma_sht,       only: sht_grid_surface_integral

    implicit none

    private
    public :: vilma_earth_t, vilma_geo_t, vilma2_earth_constraints, vilma2_geo_constraints
    public :: vilma2_gauss_dims, VILMA_DEFAULTS_FILE, VILMA_TABLE_FILE

    ! The files VILMA2 reads at run time: copies in kryos' input/ (input/sync.txt).
    character(len=*), parameter :: VILMA_DEFAULTS_FILE = "input/vilma_defaults.nml"
    character(len=*), parameter :: VILMA_TABLE_FILE    = "input/vilma-variables.md"
    character(len=*), parameter :: VILMA_GROUP         = "vilma"       ! fixed in vilma_par_load (defaults group = group)

    ! Largest relative difference of the grid's planet radius and VILMA2's r_earth that is accepted: the
    ! grid is built with radius = r_earth (written without loss in the planet name), so the two agree to
    ! round-off; a description file of another radius (an older run's) must differ by more than this to stop.
    real(dp), parameter :: RADIUS_TOL = 1.0e-12_dp

    type, extends(geo_t) :: vilma_earth_t
        type(solid_earth) :: se                          ! the model; built at PASS_GEOMETRY
        character(len=NAME_LEN) :: group  = VILMA_GROUP  ! the &vilma group of the par file
        character(len=NAME_LEN) :: ice_bg = "globe_ice"  ! the global ice field group (read by globe_ice)
        logical :: built = .false.                       ! solid_earth_init done
        real(dp) :: g_run = 0.0_dp                       ! the run's g (phys_const), for the log
    contains
        procedure :: init          => earth_init_
        procedure :: init_state    => earth_init_state_
        procedure :: update        => earth_update_
        procedure :: restart_write => earth_restart_write_
        procedure :: restart_read  => earth_restart_read_
        procedure :: write_output  => earth_write_output_
        procedure :: end           => earth_end_
    end type vilma_earth_t

    type, extends(geo_t) :: vilma_geo_t
    contains
        procedure :: init          => port_init_
        procedure :: init_state    => port_init_state_
        procedure :: update        => port_update_
        procedure :: restart_write => port_restart_write_
        procedure :: restart_read  => port_restart_read_
    end type vilma_geo_t

contains

    function vilma2_earth_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(build_flag="vilma=1", regional=.false., global=.true., max_instances=1)
    end function vilma2_earth_constraints

    function vilma2_geo_constraints() result(con)
        ! Regional only (the port), any number of domains; no backend restart file.
        type(constraints_t) :: con
        con = constraints_t(build_flag="vilma=1", restart=.false.)
    end function vilma2_geo_constraints

    subroutine vilma2_gauss_dims(lmax, nlat_par, nphi_par, nlat, nphi, name)
        ! The Gauss grid of a VILMA2 with this lmax: nlat = 2 lmax + 2, nphi = 4 lmax
        ! unless [vilma] sets nlat, nphi (>0): what build_sht of VILMA2 builds
        ! (vilma_coupling.f90). The name is GAUSS-L<lmax>, with the size appended
        ! when it is not the default.
        integer,          intent(in)  :: lmax, nlat_par, nphi_par
        integer,          intent(out) :: nlat, nphi
        character(len=*), intent(out) :: name

        nlat = nlat_par
        if (nlat <= 0) nlat = 2*lmax + 2
        nphi = nphi_par
        if (nphi <= 0) nphi = 4*lmax
        name = "GAUSS-L"//to_str(lmax)
        if (nlat /= 2*lmax + 2 .or. nphi /= 4*lmax) name = trim(name)//"-"//to_str(nphi)//"x"//to_str(nlat)
    end subroutine vilma2_gauss_dims

    ! =========================================================================
    ! Global part
    ! =========================================================================

    subroutine earth_init_(self, ctx)
        ! Plan 3.1 init: the &vilma group over input/vilma_defaults.nml, the run's
        ! densities, the Gauss grid. The model itself is built at PASS_GEOMETRY, when
        ! the reference is known.
        class(vilma_earth_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        if (nml_has_param(ctx%par_file, ctx%instance, "vilma"))  &
            call nml_read(ctx%par_file, ctx%instance, "vilma", self%group)
        if (trim(self%group) /= VILMA_GROUP) then       ! vilma_par_load looks the defaults up under the same group name
            call fesm_error("kryos_geo_vilma2::init", "key `vilma` must name group `"//VILMA_GROUP//"`", &
                            "instance "//trim(self%name)//": VILMA2 reads input/vilma_defaults.nml under that group "// &
                            "name only (got `"//trim(self%group)//"`)")
        end if
        if (nml_has_param(ctx%par_file, ctx%instance, "ice_bg")) &
            call nml_read(ctx%par_file, ctx%instance, "ice_bg", self%ice_bg)

        call vilma_io_set_table(VILMA_TABLE_FILE)
        call vilma_par_load(self%se%par, trim(ctx%par_file), defaults_file=VILMA_DEFAULTS_FILE, group=trim(self%group))

        self%prec = dp
        ! the spin-up that holds the t0 load before the seed (plan 7.4): VILMA2's own, pre_spinup_1d or equil_time_max
        self%opt  = geo_opt_t(host=.true., host_spinup=(self%se%par%pre_spinup_1d .or. self%se%par%equil_time_max > 0.0_dp))

        self%g_run = ctx%cnst%g
        call params_check(self, ctx)
        call gauss_build(self, ctx)
    end subroutine earth_init_

    subroutine params_check(self, ctx)
        ! The keys of &vilma the adapter needs, and the densities of the run.
        class(vilma_earth_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        real(dp) :: rho_ice0, rho_w0

        associate (p => self%se%par)
            if (trim(p%solver) /= "v2") then
                call fesm_error("kryos_geo_vilma2::init", "`solver` of group `"//trim(self%group)//"` must be v2", &
                                "instance "//trim(self%name)//": the host terms need the native solver (got `"// &
                                trim(p%solver)//"`)")
            end if
            if (p%sle_fixed_ocean) then
                call fesm_error("kryos_geo_vilma2::init", "`sle_fixed_ocean` of group `"//trim(self%group)// &
                                "` must be false", "instance "//trim(self%name)// &
                                ": the host terms need a migrating coastline")
            end if
            rho_ice0 = p%rho_ice
            rho_w0   = p%rho_water
            p%rho_ice   = ctx%cnst%rho_ice
            p%rho_water = ctx%cnst%rho_sw
            call log_msg(LOG_INFO, "kryos_geo_vilma2", "densities from the run: rho_ice "//to_str(rho_ice0)//" -> "// &
                         to_str(p%rho_ice)//", rho_water "//to_str(rho_w0)//" -> "//to_str(p%rho_water)//" kg m-3", &
                         "instance "//trim(self%name))
        end associate
    end subroutine params_check

    subroutine gauss_build(self, ctx)
        ! plan 3.5: GAUSS-L<lmax> as VILMA2 builds it (grid_is_gauss compares to 1e-6 degree),
        ! longitudes from 0, a sphere of VILMA2's radius r_earth (the cell areas are the quadrature
        ! areas, so integrals over the grid are VILMA2's own).
        ! The description goes to <map_fldr>/grid_<name>.txt for the couplers; an existing file
        ! must agree (a read-only map folder of a gate stack is fine).
        class(vilma_earth_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx

        integer :: nlat, nphi
        character(len=NAME_LEN) :: name
        real(dp) :: dr

        call vilma2_gauss_dims(self%se%par%lmax, self%se%par%nlat, self%se%par%nphi, nlat, nphi, name)
        if (allocated(self%grid_obj)) deallocate(self%grid_obj)
        allocate(self%grid_obj)
        call grid_init_gaussian(self%grid_obj, trim(name), nlon=nphi, nlat=nlat, lon180=.false., x0=0.0_dp, &
                                radius=self%se%par%r_earth)

        dr = abs(self%grid_obj%cs%planet%R - self%se%par%r_earth) / self%se%par%r_earth
        if (dr > RADIUS_TOL) then
            call fesm_error("kryos_geo_vilma2::init", "the grid's planet radius differs from `r_earth` of group `"// &
                            trim(self%group)//"`", "instance "//trim(self%name)//": grid "//trim(name)//" R = "// &
                            to_str(self%grid_obj%cs%planet%R)//" m, r_earth = "//to_str(self%se%par%r_earth)// &
                            " m (relative difference "//to_str(dr)//", tolerance "//to_str(RADIUS_TOL)//")")
        end if
        call grid_file(self, ctx%map_fldr)

        self%grid = trim(name)
        self%nx   = nphi
        self%ny   = nlat
        call log_msg(LOG_INFO, "kryos_geo_vilma2", "Gauss grid `"//trim(name)//"`: "//to_str(nphi)//" x "//to_str(nlat)// &
                     ", lmax "//to_str(self%se%par%lmax)//"; planet radius "//to_str(self%grid_obj%cs%planet%R)// &
                     " m, VILMA2 r_earth "//to_str(self%se%par%r_earth)//" m (relative difference "//to_str(dr)//")", &
                     "instance "//trim(self%name))
    end subroutine gauss_build

    subroutine grid_file(self, fldr)
        ! The grid description for the couplers: written when absent, else compared.
        class(vilma_earth_t), intent(in) :: self
        character(len=*),     intent(in) :: fldr

        character(len=*), parameter :: proc = "kryos_geo_vilma2::init"
        type(grid_class) :: g0
        character(len=:), allocatable :: file
        logical :: found

        file = trim(fldr)//"/grid_"//trim(self%grid_obj%name)//".txt"
        inquire(file=file, exist=found)
        if (.not. found) then
            call execute_command_line("mkdir -p '"//trim(fldr)//"'")
            call grid_cdo_write_desc_short(self%grid_obj, trim(fldr))
            call log_msg(LOG_INFO, "kryos_geo_vilma2", "wrote the grid description "//file, "instance "//trim(self%name))
            return
        end if

        call grid_cdo_read_desc(g0, trim(self%grid_obj%name), trim(fldr))
        if (g0%G%nx /= self%grid_obj%G%nx .or. g0%G%ny /= self%grid_obj%G%ny) then
            call fesm_error(proc, "grid description "//file//" has another size than the Gauss grid of VILMA2", &
                            "instance "//trim(self%name)//": file "//to_str(g0%G%nx)//" x "//to_str(g0%G%ny)// &
                            ", VILMA2 "//to_str(self%grid_obj%G%nx)//" x "//to_str(self%grid_obj%G%ny))
        end if
        if (abs(g0%cs%planet%R - self%grid_obj%cs%planet%R) > RADIUS_TOL*self%grid_obj%cs%planet%R) then
            call fesm_error(proc, "grid description "//file//" has another planet radius than the Gauss grid of VILMA2", &
                            "instance "//trim(self%name)//": file R = "//to_str(g0%cs%planet%R)//" m, r_earth = "// &
                            to_str(self%grid_obj%cs%planet%R)//" m; remove the file to have it rewritten")
        end if
        if (maxval(abs(g0%G%x - self%grid_obj%G%x)) > 1.0e-9_dp .or. maxval(abs(g0%G%y - self%grid_obj%G%y)) > 1.0e-9_dp) then
            call fesm_error(proc, "grid description "//file//" differs from the Gauss grid of VILMA2", &
                            "instance "//trim(self%name)//": its axes differ by more than 1e-9 degree; remove the file "// &
                            "to have it rewritten")
        end if
    end subroutine grid_file

    subroutine earth_init_state_(self, st, pass)
        ! PASS_GEOMETRY: solid_earth_init at the reference, the exports at it. PASS_STATE
        ! (cold start): VILMA2's own spin-up under the t0 load, then the seed (dt = 0).
        ! PASS_INIT and the others: nothing (PASS_RESTART: the bundle's exports stand,
        ! restart_read ended with solid_earth_sync_host).
        class(vilma_earth_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        integer,              intent(in)    :: pass

        select case (pass)
        case (PASS_GEOMETRY)
            call build_model(self, self%g_run)
            call exports_state(self, st%t1)             ! the reference: nothing has moved
        case (PASS_STATE)
            call spin_up_and_seed(self)
            call exports_state(self, st%t1)
        end select
    end subroutine earth_init_state_

    subroutine build_model(self, g_run)
        ! plan 3.1: no h_ice_init; the reference is the background X on the Gauss grid (the 1 - w
        ! part) and the reference host terms. g_run: the run's g, logged beside VILMA2's own (its
        ! gravity comes from the Earth model; there is no key to set).
        class(vilma_earth_t), intent(inout) :: self
        real(dp),             intent(in)    :: g_run

        real(dp), allocatable :: z_bed_ref(:,:), H_ice_ref(:,:), w(:,:), Ig(:,:), C(:,:), F0(:,:)

        call require_published(self, self%hi_z_bed_ref, "z_bed_ref", "init_state")
        call require_published(self, self%hi_H_ice_ref, "H_ice_ref", "init_state")
        call require_published(self, self%hi_host_w,      "host_w",      "init_state")
        call require_published(self, self%hi_host_Ig_ref, "host_Ig_ref", "init_state")
        call require_published(self, self%hi_host_C_ref,  "host_C_ref",  "init_state")
        allocate(z_bed_ref(self%nx, self%ny), H_ice_ref(self%nx, self%ny), w(self%nx, self%ny), &
                 Ig(self%nx, self%ny), C(self%nx, self%ny), F0(self%nx, self%ny))
        call fields_get(self%imports, self%hi_z_bed_ref,   z_bed_ref)
        call fields_get(self%imports, self%hi_H_ice_ref,   H_ice_ref)
        call fields_get(self%imports, self%hi_host_w,      w)
        call fields_get(self%imports, self%hi_host_Ig_ref, Ig)
        call fields_get(self%imports, self%hi_host_C_ref,  C)
        F0 = 0.0_dp
        call host_check(self, w, Ig, C, F0, "the reference")

        call solid_earth_init(self%se, z_bed_ref, H_ice_ref, grid=self%grid_obj, w=w, I_g=Ig, C=C)
        self%built = .true.
        call log_msg(LOG_INFO, "kryos_geo_vilma2", "Earth built at the reference; surface gravity of its Earth model "// &
                     to_str(self%se%resp%g)//" m s-2 (the run's g "//to_str(g_run)//")", "instance "//trim(self%name))
        if (self%se%remap) then
            call fesm_error("kryos_geo_vilma2::init_state", "VILMA2 does not take grid `"//trim(self%grid)// &
                            "` as its own Gauss grid", "instance "//trim(self%name)// &
                            ": it built a remap (grid_is_gauss compares every axis value to 1e-6 degree)")
        end if
    end subroutine build_model

    subroutine spin_up_and_seed(self)
        ! plan 7.4: the spin-up holds the t0 load (a no-op unless pre_spinup_1d or equil_time_max > 0),
        ! then the seed solves under the same load (the dt = 0 branch uses the load of the last call).
        class(vilma_earth_t), intent(inout) :: self

        real(dp), allocatable :: H_ext(:,:), w(:,:), Ig(:,:), C(:,:), F(:,:)

        call require_built(self, "init_state")
        call load_terms(self, H_ext, w, Ig, C, F, "init_state")
        call host_check(self, w, Ig, C, F, "the t0 load")
        call solid_earth_spinup(self%se, H_ext, w, Ig, C, F)
        call solid_earth_update(self%se, H_ext, 0.0_dp, w, Ig, C, F)
        call log_msg(LOG_INFO, "kryos_geo_vilma2", "spin-up and seed done; VILMA2 clock "//to_str(self%se%time)//" yr", &
                     "instance "//trim(self%name))
    end subroutine spin_up_and_seed

    subroutine earth_update_(self, st)
        ! plan 3.1 update: dt = t1 - t0 of the instance's cadence; the load is the one the core gathered
        ! (lag 1). VILMA2 ramps from the previous call's load to this one.
        class(vilma_earth_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st

        real(dp), allocatable :: H_ext(:,:), w(:,:), Ig(:,:), C(:,:), F(:,:)

        call require_built(self, "update")
        call load_terms(self, H_ext, w, Ig, C, F, "update")
        call host_check(self, w, Ig, C, F, "the step load")
        call solid_earth_update(self%se, H_ext, st%t1 - st%t0, w, Ig, C, F)
        call exports_state(self, st%t1)
    end subroutine earth_update_

    subroutine load_terms(self, H_ext, w, Ig, C, F, what)
        ! The imports of one solve: the external ice and the host terms of the gathered load.
        class(vilma_earth_t), intent(inout) :: self
        real(dp), allocatable, intent(out)  :: H_ext(:,:), w(:,:), Ig(:,:), C(:,:), F(:,:)
        character(len=*),     intent(in)    :: what

        call require_published(self, self%hi_H_ice_ext, "H_ice_ext", what)
        call require_published(self, self%hi_host_w,    "host_w",    what)
        call require_published(self, self%hi_host_Ig,   "host_Ig",   what)
        call require_published(self, self%hi_host_C,    "host_C",    what)
        call require_published(self, self%hi_host_F,    "host_F",    what)
        allocate(H_ext(self%nx, self%ny), w(self%nx, self%ny), Ig(self%nx, self%ny), C(self%nx, self%ny), &
                 F(self%nx, self%ny))
        call fields_get(self%imports, self%hi_H_ice_ext, H_ext)
        call fields_get(self%imports, self%hi_host_w,    w)
        call fields_get(self%imports, self%hi_host_Ig,   Ig)
        call fields_get(self%imports, self%hi_host_C,    C)
        call fields_get(self%imports, self%hi_host_F,    F)
    end subroutine load_terms

    subroutine host_check(self, w, Ig, C, F, what)
        ! VILMA2 stops on inadmissible host terms; the same check here names the instance and the load.
        class(vilma_earth_t), intent(in) :: self
        real(dp),             intent(in) :: w(:,:), Ig(:,:), C(:,:), F(:,:)
        character(len=*),     intent(in) :: what

        character(len=160) :: msg

        if (.not. solid_earth_host_check(w, Ig, C, F, msg)) then
            call fesm_error("kryos_geo_vilma2::host_check", "host terms of "//what//" are not admissible", &
                            "instance "//trim(self%name)//": "//trim(msg))
        end if
    end subroutine host_check

    subroutine earth_restart_write_(self, dir, t)
        ! plan 3.1: the restart file of VILMA2, labelled with its own clock (se%time, not the model time),
        ! always lossless; the folder exists (the bundle's).
        class(vilma_earth_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t

        associate (unused => t); end associate
        call require_built(self, "restart_write")
        call vilma_restart_write(self%se, self%se%time, filename=trim(self%name)//"_restart.nc", folder=trim(dir))
    end subroutine earth_restart_write_

    subroutine earth_restart_read_(self, dir, t)
        ! The last slice of the file, into the model PASS_GEOMETRY built at the reference: VILMA2 checks
        ! the reference against it and ends with solid_earth_sync_host.
        class(vilma_earth_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t

        character(len=:), allocatable :: file
        logical :: found

        associate (unused => t); end associate
        call require_built(self, "restart_read")
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
        inquire(file=file, exist=found)
        if (.not. found) then
            call fesm_error("kryos_geo_vilma2::restart_read", "restart file not found: "//file, &
                            "instance "//trim(self%name))
        end if
        call vilma_restart_read(self%se, file)
    end subroutine earth_restart_read_

    subroutine earth_write_output_(self, file, t, stream, init)
        ! Streams 2D and 1D: the default writer (the exports on lon, lat; rounded by the table `geo`). `gia`: the global
        ! ice field VILMA2 got (the host imports, table `geo-vilma2`). `vilma`: VILMA2's own diagnostics file
        ! (vilma_write_step, rounded by its own table; the ice, sea level, bed, ocean function and the
        ! scalar series: not the viscoelastic memory).
        class(vilma_earth_t), intent(inout) :: self
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init

        character(len=*), parameter :: DIMS_LL(3) = [character(len=32) :: "lon", "lat", "time"]
        character(len=*), parameter :: VARS_VILMA(12) = [character(len=12) :: "h_ice", "rsl", "z_bed", "C_ocean", &
            "dz_bed", "dz_ss", "bsl", "dphi", "bsl_ive", "bsl_ive_host", "fill_int", "fill_int_own"]
        type(var_io_type) :: row
        real(dp), allocatable :: a(:,:)
        integer :: hs(5), i, ncid, n
        logical :: found
        character(len=12), parameter :: NAMES(5) = [character(len=12) :: "H_ice_ext", "host_w", "host_Ig", "host_C", "host_F"]

        select case (trim(stream))
        case ("gia")
            call require_built(self, "write_output")
            hs = [self%hi_H_ice_ext, self%hi_host_w, self%hi_host_Ig, self%hi_host_C, self%hi_host_F]
            if (.not. self%imports%f(hs(1))%published) return                ! nothing gathered yet
            if (init) call outfile_dims_lonlat(file, self%grid_obj%G%x, self%grid_obj%G%y, t)
            call outfile_open(file, t, ncid, n)
            call nc_write(file, "t_load", self%imports%f(hs(1))%t1, dim1="time", start=[n], count=[1], ncid=ncid, &
                          units="year", long_name="time of the ice the host terms stand for")
            allocate(a(self%nx, self%ny))
            do i = 1, size(hs)
                call fields_get(self%imports, hs(i), a)
                call outfile_row("geo-vilma2", trim(NAMES(i)), row, found)
                if (.not. found) then
                    call fesm_error("kryos_geo_vilma2::write_output", "variable `"//trim(NAMES(i))//"` has no row in "// &
                                    outfile_table_file("geo-vilma2"), "instance "//trim(self%name))
                end if
                row%units = self%imports%f(hs(i))%units
                call var_io_write(file, row, a, n, ncid=ncid, dims=DIMS_LL)
            end do
            call nc_close(ncid)
        case ("vilma")
            call require_built(self, "write_output")
            call vilma_write_step(self%se, file, t, nms=VARS_VILMA, init=init)
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine earth_write_output_

    subroutine earth_end_(self)
        class(vilma_earth_t), intent(inout) :: self
        call solid_earth_finalize(self%se)
        self%built = .false.
    end subroutine earth_end_

    subroutine exports_state(self, t1)
        ! INSTANT at t1 (plan 3.1): dz_bed = u, dz_ss = N (VILMA2's dz_ss minus the uniform shift dphi),
        ! the diagnostics. A_ocn = integral of C over VILMA2's sphere. At the reference (PASS_GEOMETRY)
        ! this puts zeros and the reference ocean function.
        class(vilma_earth_t), intent(inout) :: self
        real(dp),             intent(in)    :: t1

        real(dp) :: A_ocn

        A_ocn = self%se%par%r_earth**2 * sht_grid_surface_integral(self%se%sht, self%se%gg%C)
        associate (se => self%se)
            call fields_put(self%exports, self%h_dz_bed,     se%gg%dz_bed,               t1=t1)
            call fields_put(self%exports, self%h_dz_ss,      se%gg%dz_ss - se%dphi,      t1=t1)
            call fields_put(self%exports, self%h_C_ocn,      se%gg%C,                    t1=t1)
            call fields_put(self%exports, self%h_dphi,       se%dphi,                    t1=t1)
            call fields_put(self%exports, self%h_bsl,        se%bsl,                     t1=t1)
            call fields_put(self%exports, self%h_bsl_ive,    se%bsl_ive,                 t1=t1)
            call fields_put(self%exports, self%h_bsl_ive_host, se%bsl_ive_host,          t1=t1)
            call fields_put(self%exports, self%h_A_ocn,      A_ocn,                      t1=t1)
            call fields_put(self%exports, self%h_mass_resid, se%worst_mass_resid,        t1=t1)
        end associate
    end subroutine exports_state

    subroutine require_built(self, what)
        class(vilma_earth_t), intent(in) :: self
        character(len=*),     intent(in) :: what
        if (.not. self%built) then
            call fesm_error("kryos_geo_vilma2::"//what, "the Earth has not been built", &
                            "instance "//trim(self%name)//": PASS_GEOMETRY (solid_earth_init at the reference) comes first")
        end if
    end subroutine require_built

    subroutine require_published(self, h, name, what)
        class(vilma_earth_t), intent(in) :: self
        integer,              intent(in) :: h
        character(len=*),     intent(in) :: name, what
        if (h <= 0) then
            call fesm_error("kryos_geo_vilma2::"//what, "import "//name//" is not declared", &
                            "instance "//trim(self%name))
        else if (.not. self%imports%f(h)%published) then
            call fesm_error("kryos_geo_vilma2::"//what, "import "//name//" has not been published", &
                            "instance "//trim(self%name)//": the core fills it (gather, globe_ice) before this call")
        end if
    end subroutine require_published

    ! =========================================================================
    ! Regional port
    ! =========================================================================

    subroutine port_init_(self, ctx)
        ! The port is on the grid the chain gives it (the ice grid by default); all it holds is the
        ! scattered anomalies. Its cadence is the global's (checked by the run, kryos_run).
        class(vilma_geo_t), intent(inout) :: self
        type(init_ctx_t),   intent(in)    :: ctx

        self%prec = dp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = geo_opt_t(ref=.true., port=.true.)
    end subroutine port_init_

    subroutine port_init_state_(self, st, pass)
        ! PASS_GEOMETRY: the reference (zero anomalies, z_bed_ref from the ice). PASS_STATE: the
        ! imports the core scattered. Other passes: nothing (PASS_RESTART: the bundle's exports stand).
        class(vilma_geo_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        integer,            intent(in)    :: pass

        real(dp) :: zero2(self%nx, self%ny)

        select case (pass)
        case (PASS_GEOMETRY)
            zero2 = 0.0_dp
            call fields_put(self%exports, self%h_dz_bed, zero2,    t1=st%t1)
            call fields_put(self%exports, self%h_dz_ss,  zero2,    t1=st%t1)
            call fields_put(self%exports, self%h_dphi,   0.0_dp,   t1=st%t1)
            call fields_put(self%exports, self%h_bsl,    0.0_dp,   t1=st%t1)
            call put_z_bed_ref(self, "init_state")
        case (PASS_STATE)
            call copy_imports(self, "init_state")
        end select
    end subroutine port_init_state_

    subroutine port_update_(self, st)
        class(vilma_geo_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        associate (unused => st%t1); end associate
        call copy_imports(self, "update")
    end subroutine port_update_

    subroutine copy_imports(self, what)
        ! dz_bed, dz_ss, dphi, bsl as scattered, with their validity; z_bed_ref from the ice.
        class(vilma_geo_t), intent(inout) :: self
        character(len=*),   intent(in)    :: what

        real(dp) :: a(self%nx, self%ny), s, t0, t1

        call port_published(self, self%hi_dz_bed, "dz_bed", what)
        call fields_get(self%imports, self%hi_dz_bed, a, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_dz_bed, a, t1=t1, t0=t0)
        call port_published(self, self%hi_dz_ss, "dz_ss", what)
        call fields_get(self%imports, self%hi_dz_ss, a, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_dz_ss, a, t1=t1, t0=t0)
        call port_published(self, self%hi_dphi, "dphi", what)
        call fields_get(self%imports, self%hi_dphi, s, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_dphi, s, t1=t1, t0=t0)
        call port_published(self, self%hi_bsl, "bsl", what)
        call fields_get(self%imports, self%hi_bsl, s, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_bsl, s, t1=t1, t0=t0)
        call put_z_bed_ref(self, what)
    end subroutine copy_imports

    subroutine put_z_bed_ref(self, what)
        ! The ice's static reference bed, passed through (domain_check_reference compares it).
        class(vilma_geo_t), intent(inout) :: self
        character(len=*),   intent(in)    :: what

        real(dp) :: a(self%nx, self%ny), t0, t1

        call port_published(self, self%hi_z_bed_ref, "z_bed_ref", what)
        call fields_get(self%imports, self%hi_z_bed_ref, a, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_z_bed_ref, a, t1=t1, t0=t0)
    end subroutine put_z_bed_ref

    subroutine port_published(self, h, name, what)
        class(vilma_geo_t), intent(in) :: self
        integer,            intent(in) :: h
        character(len=*),   intent(in) :: name, what
        if (h <= 0) then
            call fesm_error("kryos_geo_vilma2::"//what, "import "//name//" is not declared", &
                            "instance "//trim(self%name))
        else if (.not. self%imports%f(h)%published) then
            call fesm_error("kryos_geo_vilma2::"//what, "import "//name//" has not been published", &
                            "instance "//trim(self%name)//": the core exchanges it (scatter, ice) before this call")
        end if
    end subroutine port_published

    subroutine port_restart_write_(self, dir, t)
        ! Stateless: constraints restart = .false., so the core never calls this.
        class(vilma_geo_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine port_restart_write_

    subroutine port_restart_read_(self, dir, t)
        class(vilma_geo_t), intent(inout) :: self
        character(len=*),   intent(in)    :: dir
        real(dp),           intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine port_restart_read_

end module kryos_geo_vilma2
