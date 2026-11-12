module kryos_climate_paleoemu
    ! climate backend paleoemu, built with paleoemu=1, else the stub: the
    ! CESM-3Ma emulator of PaleoEmulator (fesmc/PaleoEmulator, module paleoemu of
    ! its fortran/) as a two-level climate, as geo = vilma2 is a two-level geo:
    !
    !   paleoemu_climate_t  the global part (scope global, one per run) on its own
    !                       lon-lat grid (key `grid`, the background grid of
    !                       [globe_ice]): the ice geometry and its reference come
    !                       from the core (globe_ice: the background with the
    !                       domains merged in), CO2 from the forcing (coef_co2), the
    !                       orbit from the forcing (orbit_ecc, orbit_obl,
    !                       orbit_per: the run's orbit, also the surface's
    !                       insolation; without them the emulator's own table at
    !                       the evaluation time). Each update binds the geometry to the
    !                       emulator's ice grid, steps the emulator and applies its
    !                       anomalies to the pre-industrial reference climate on
    !                       the grid (key `reference`): tas at the reference surface
    !                       z_srf_clim, pr, lapse = Gamma, beta_p. The temperature
    !                       adjustments (paleoemu_downscale): scale_dtas times the
    !                       anomaly (a colder or warmer climate with its pattern),
    !                       the bias correction of the file `bias` (e.g. PMIP4's
    !                       LGM minus the emulator's) times bias_scale and the NH
    !                       ice-volume weight; beta_adj carries them into the
    !                       precipitation, exp(beta_adj dT). beta_p, the lapse
    !                       rows' precipitation per K of elevation cooling, stays 0
    !                       by default: the emulator's dlogpr already holds the
    !                       response to the ice sheets (only dtas is at the
    !                       reference surface). The ocean: the annual ocean
    !                       temperature change dto at the emulator's levels
    !                       (axis `level`, 4-1475 m), scale_dto (dto - dto0) +
    !                       dto_nsd dto_sd, its land and below-floor cells filled
    !                       per level (ocean_fill), interpolated bilinearly onto the
    !                       grid: the export dto_lev [K] (opt%ocn_lev). dto0, the
    !                       present-day anchor (key dto_anchor, default on): dto of
    !                       a copy of the emulator at present-day inputs (0 ka, its
    !                       own orbit table, CO2 dto_anchor_co2, the reference ice,
    !                       the filters in equilibrium), once the reference is set;
    !                       it removes the offset of the emulator's reference
    !                       inputs (CESM-3Ma's 0 ka, whose filters hold the
    !                       deglaciation) from a run at present day.
    !   paleoemu_port_t     the regional port (stateless): the global part's
    !                       atmosphere remapped onto the domain's grid by the rows
    !                       climate_g -> climate, passed on; the lapse rows
    !                       climate -> surface take tas and pr to the ice surface.
    !                       With a present-day ocean on its grid (keys to_ref,
    !                       so_ref: var_defs groups, read as the snapshots read
    !                       theirs, on the file's depth axis) it exports the
    !                       profiles of set A: dto_ann = the imported dto_lev
    !                       interpolated linearly in depth onto the reference's
    !                       depths (constant above the top and below the bottom
    !                       level), to_ann = to_ref + dto_ann, so_ann = so_ref.
    !                       With basins (key dto_basins, a basin set of the
    !                       domain's regions, extended over the whole domain by
    !                       nearest neighbour: the sets end at the continental
    !                       shelf) dto_ann is the mean of each basin's sample (key
    !                       dto_sample, a regions selection: by default the ocean
    !                       cells of the shelf-break buffer, FesmData's zone
    !                       between the continental shelf and the open ocean) in
    !                       all of its cells: a basin-wide change from the ocean
    !                       the shelf exchanges with, not the emulator's filled
    !                       values on the shelf itself. A basin without a sampled
    !                       cell samples its ocean cells instead (the reference
    !                       topography of the init context: bed below sea level,
    !                       not grounded).
    !
    ! The factory picks by scope. The emulator works in dp; nothing is narrowed.
    ! Precipitation: the emulator's mm d-1, exported as kg m-2 yr-1 (days of the
    ! run's year, sec_year / 86400).

    use precision,       only: dp
    use coordinates,     only: grid_class
    use coords,          only: grid_cdo_read_desc
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use interp2D,        only: fill_poisson, fill_nearest
    use nml,             only: nml_read, nml_has_param
    use varslice,        only: varslice_class, varslice_init_nml, varslice_update, varslice_end
    use kryos_snap_nml,  only: group_file
    use regions,         only: regions_basin_ids, regions_select
    use ncio,            only: nc_read, nc_read_attr, nc_size, nc_create, nc_write_dim, nc_write
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, axis_t, axis_find, NAME_LEN, PATH_LEN, &
                               PASS_INIT, PASS_BOUNDARY, PASS_RESTART, component_write_output
    use kryos_climate,   only: climate_t, climate_opt_t
    use paleoemu,        only: paleoemu_class, paleoemu_target_class, paleoemu_init, paleoemu_set_reference, &
                               paleoemu_update, paleoemu_end, paleoemu_target_init, paleoemu_downscale, paleoemu_bias_init, &
                               paleoemu_target_bilinear, paleoemu_write_init, paleoemu_write_step

    implicit none

    private
    public :: paleoemu_climate_t, paleoemu_port_t, paleoemu_global_constraints, paleoemu_port_constraints
    public :: ocean_fill

    integer,  parameter :: NMON = 12
    real(dp), parameter :: MV = -9999.0_dp                 ! missing value of ocean_fill

    type, extends(climate_t) :: paleoemu_climate_t
        character(len=PATH_LEN) :: file      = ""          ! the emulator file
        character(len=PATH_LEN) :: reference = ""          ! the reference climate on the instance grid
        real(dp) :: beta_p   = 0.0_dp                      ! [K-1] lapse rows: precipitation per K of elevation cooling
        real(dp) :: beta_adj = 0.05_dp                     ! [K-1] precipitation per K of the adjustments (scale_dtas, bias)
        real(dp) :: scale_dtas = 1.0_dp                    ! the emulator's temperature anomaly times this
        character(len=PATH_LEN) :: bias = ""               ! the bias correction file ("" none)
        real(dp) :: bias_scale = 1.0_dp                    ! the bias correction times this (and the ice-volume weight)
        real(dp) :: scale_dto = 1.0_dp                     ! the ocean: scale_dto dto + dto_nsd dto_sd
        real(dp) :: dto_nsd   = 0.0_dp
        logical  :: dto_anchor = .true.                    ! dto0: the emulator's dto at present-day inputs, subtracted
        real(dp) :: dto_anchor_co2 = 280.0_dp              ! [ppm] the CO2 of the anchor's inputs
        real(dp), allocatable :: dto0(:,:,:)               ! [lon, lat, level] the anchor on the CESM grid (0 without it)
        real(dp) :: day_year = 365.0_dp                    ! days of the run's year (pr: mm d-1 -> kg m-2 yr-1)
        type(paleoemu_class)        :: emu
        type(paleoemu_target_class) :: tgt                 ! the instance grid with its reference climate
        integer,  allocatable :: bin(:,:)                  ! emulator ice-grid cell of each instance cell (0: outside)
        real(dp), allocatable :: wsum(:,:)                 ! area of the instance cells in each ice-grid cell [m2]
        logical  :: ref_set = .false.
        real(dp), allocatable :: filt_rst(:,:)             ! the filter state of a restart, applied at PASS_RESTART
        real(dp) :: time_rst = 0.0_dp
        logical  :: has_rst = .false.
    contains
        procedure :: init          => global_init_
        procedure :: init_state    => global_init_state_
        procedure :: update        => global_update_
        procedure :: restart_write => global_restart_write_
        procedure :: restart_read  => global_restart_read_
        procedure :: write_output  => global_write_output_
        procedure :: end           => global_end_
    end type paleoemu_climate_t

    type, extends(climate_t) :: paleoemu_port_t
        character(len=NAME_LEN) :: to_ref = "", so_ref = ""     ! var_defs groups of the present-day ocean ("" none)
        real(dp), allocatable :: to0(:,:,:), so0(:,:,:)         ! [x, y, depth] the reference ocean [K], [psu]
        real(dp), allocatable :: W(:,:)                         ! [depth, level] linear-in-depth weights
        character(len=NAME_LEN) :: dto_basins = ""               ! basin set of the basin means ("" none)
        character(len=512)      :: dto_sample = "zone:shelf_break_buffer"   ! the cells each basin's mean samples
        integer, allocatable :: basin(:,:)                      ! [x, y] basin id, extended over the domain
        logical, allocatable :: sample(:,:)                     ! [x, y] sampled cells (with the fallback's)
        integer, allocatable :: nsample(:)                      ! [basin] sampled cells of each basin
    contains
        procedure :: init          => port_init_
        procedure :: init_state    => port_init_state_
        procedure :: update        => port_update_
        procedure :: restart_write => port_restart_write_
        procedure :: restart_read  => port_restart_read_
    end type paleoemu_port_t

contains

    function paleoemu_global_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(build_flag="paleoemu=1", regional=.false., global=.true., max_instances=1)
    end function paleoemu_global_constraints

    function paleoemu_port_constraints() result(con)
        ! Regional only (the port), any number of domains; no backend restart file.
        type(constraints_t) :: con
        con = constraints_t(build_flag="paleoemu=1", restart=.false.)
    end function paleoemu_port_constraints

    ! =========================================================================
    ! Global part
    ! =========================================================================

    subroutine global_init_(self, ctx)
        ! Keys, the emulator file, the instance grid (a grid_<name>.txt of the map folder), the reference
        ! climate on it, the binning onto the emulator's ice grid; the coefficient co2, the orbit (optional).
        class(paleoemu_climate_t), intent(inout) :: self
        type(init_ctx_t),          intent(in)    :: ctx

        character(len=*), parameter :: proc = "kryos_climate_paleoemu::init"
        character(len=NAME_LEN) :: grid_name
        character(len=PATH_LEN+32) :: bias_txt

        call nml_read(ctx%par_file, ctx%instance, "file", self%file)
        call nml_read(ctx%par_file, ctx%instance, "reference", self%reference)
        if (nml_has_param(ctx%par_file, ctx%instance, "beta_p")) &
            call nml_read(ctx%par_file, ctx%instance, "beta_p", self%beta_p)
        if (nml_has_param(ctx%par_file, ctx%instance, "beta_adj")) &
            call nml_read(ctx%par_file, ctx%instance, "beta_adj", self%beta_adj)
        if (nml_has_param(ctx%par_file, ctx%instance, "scale_dtas")) &
            call nml_read(ctx%par_file, ctx%instance, "scale_dtas", self%scale_dtas)
        if (nml_has_param(ctx%par_file, ctx%instance, "bias")) call nml_read(ctx%par_file, ctx%instance, "bias", self%bias)
        if (nml_has_param(ctx%par_file, ctx%instance, "bias_scale")) &
            call nml_read(ctx%par_file, ctx%instance, "bias_scale", self%bias_scale)
        if (nml_has_param(ctx%par_file, ctx%instance, "scale_dto")) &
            call nml_read(ctx%par_file, ctx%instance, "scale_dto", self%scale_dto)
        if (nml_has_param(ctx%par_file, ctx%instance, "dto_nsd")) &
            call nml_read(ctx%par_file, ctx%instance, "dto_nsd", self%dto_nsd)
        if (nml_has_param(ctx%par_file, ctx%instance, "dto_anchor")) &
            call nml_read(ctx%par_file, ctx%instance, "dto_anchor", self%dto_anchor)
        if (nml_has_param(ctx%par_file, ctx%instance, "dto_anchor_co2")) &
            call nml_read(ctx%par_file, ctx%instance, "dto_anchor_co2", self%dto_anchor_co2)
        grid_name = ""
        if (nml_has_param(ctx%par_file, ctx%instance, "grid")) call nml_read(ctx%par_file, ctx%instance, "grid", grid_name)
        if (len_trim(self%file) == 0 .or. len_trim(self%reference) == 0 .or. len_trim(grid_name) == 0) then
            call fesm_error(proc, "keys `file`, `reference` and `grid` are required", &
                            "instance "//trim(self%name)//": the emulator file, the reference climate on the grid, "// &
                            "the grid (the background grid of [globe_ice], a grid_<name>.txt of the map folder)")
        end if

        self%prec = dp
        self%opt  = climate_opt_t(host=.true., orbit=.true., ocn_lev=.true.)
        self%day_year = ctx%sec_year/86400.0_dp

        allocate(self%grid_obj)
        call grid_cdo_read_desc(self%grid_obj, trim(grid_name), trim(ctx%map_fldr))
        self%grid = trim(grid_name)
        self%nx   = self%grid_obj%G%nx
        self%ny   = self%grid_obj%G%ny

        call paleoemu_init(self%emu, trim(self%file))
        if (len_trim(self%bias) > 0) call paleoemu_bias_init(self%emu, trim(self%bias))
        call reference_load(self)
        call bins_init(self)
        call level_axis(self)

        allocate(self%coef_names(1), self%coef_units(1))
        self%coef_names(1) = "co2"
        self%coef_units(1) = "ppm"

        bias_txt = "none"
        if (len_trim(self%bias) > 0) bias_txt = trim(self%bias)//" x "//to_str(self%bias_scale)
        call log_msg(LOG_INFO, "kryos_climate_paleoemu", "emulator "//trim(self%file)//" ("//trim(self%emu%temperature)// &
                     "), grid "//trim(grid_name)//" "//to_str(self%nx)//" x "//to_str(self%ny)//", beta_p "// &
                     to_str(self%beta_p)//", beta_adj "//to_str(self%beta_adj)//", scale_dtas "//to_str(self%scale_dtas)//", bias "//trim(bias_txt)// &
                     ", ocean: "//to_str(self%emu%nlev)//" levels "//to_str(self%emu%ocean_depth(1))//"-"// &
                     to_str(self%emu%ocean_depth(self%emu%nlev))//" m, scale_dto "//to_str(self%scale_dto)//", dto_nsd "// &
                     to_str(self%dto_nsd)//", dto_anchor "//merge("on ", "off", self%dto_anchor), "instance "//trim(self%name))
    end subroutine global_init_

    subroutine level_axis(self)
        ! The ocean axis `level` [m, positive down] of dto_lev: the emulator's levels (its ports declare
        ! their import on it; domain_create passes it to them).
        class(paleoemu_climate_t), intent(inout) :: self
        type(axis_t) :: ax

        ax%name  = "level"
        ax%units = "m"
        ax%v     = self%emu%ocean_depth
        if (.not. allocated(self%axes)) allocate(self%axes(0))
        self%axes = [self%axes, ax]
    end subroutine level_axis

    subroutine reference_load(self)
        ! The reference climate (tas_ref [K], pr_ref [mm d-1] by month, z_srf_ref [m]) on the instance grid,
        ! with its bilinear weights from the CESM grid (paleoemu_target_init).
        class(paleoemu_climate_t), intent(inout) :: self

        character(len=*), parameter :: proc = "kryos_climate_paleoemu::init"
        real(dp), allocatable :: tas(:,:,:), pr(:,:,:), zs(:,:), lon(:), lat(:)
        character(len=NAME_LEN) :: gname
        integer :: nx, ny

        nx = nc_size(trim(self%reference), "lon")
        ny = nc_size(trim(self%reference), "lat")
        if (nx /= self%nx .or. ny /= self%ny) then
            call fesm_error(proc, "the reference climate is not on the instance grid", "instance "//trim(self%name)// &
                            ": "//trim(self%reference)//" is "//to_str(nx)//" x "//to_str(ny)//", grid "//trim(self%grid)// &
                            " "//to_str(self%nx)//" x "//to_str(self%ny))
        end if
        allocate(lon(nx), lat(ny), tas(nx, ny, NMON), pr(nx, ny, NMON), zs(nx, ny))
        call nc_read(trim(self%reference), "lon", lon)
        call nc_read(trim(self%reference), "lat", lat)
        if (maxval(abs(lon - self%grid_obj%lon(:, 1))) > 1.0e-6_dp .or. maxval(abs(lat - self%grid_obj%lat(1, :))) > 1.0e-6_dp) then
            call fesm_error(proc, "the axes of the reference climate differ from the instance grid", &
                            "instance "//trim(self%name)//": "//trim(self%reference)//" against grid "//trim(self%grid))
        end if
        gname = ""
        call nc_read_attr(trim(self%reference), "grid_name", gname)
        if (trim(gname) /= trim(self%grid)) then
            call fesm_error(proc, "the reference climate names another grid", "instance "//trim(self%name)//": "// &
                            trim(self%reference)//" has grid_name `"//trim(gname)//"`, the instance grid is "//trim(self%grid))
        end if
        call nc_read(trim(self%reference), "tas_ref", tas)
        call nc_read(trim(self%reference), "pr_ref", pr)
        call nc_read(trim(self%reference), "z_srf_ref", zs)
        call paleoemu_target_init(self%tgt, self%emu, self%grid_obj%lon, self%grid_obj%lat, tas, pr, zs)
    end subroutine reference_load

    subroutine bins_init(self)
        ! Each instance cell in the emulator's ice-grid cell that holds its centre (0 outside the ice grid);
        ! every ice-grid cell must hold at least one (an instance grid at least as fine as the ice grid).
        class(paleoemu_climate_t), intent(inout) :: self

        integer :: i, j, ii, jj, nix, niy
        real(dp) :: dlon, dlat, lon0, lat0

        nix = size(self%emu%ice_lon); niy = size(self%emu%ice_lat)
        dlon = self%emu%ice_lon(2) - self%emu%ice_lon(1)
        dlat = self%emu%ice_lat(2) - self%emu%ice_lat(1)
        lon0 = self%emu%ice_lon(1) - 0.5_dp*dlon
        lat0 = self%emu%ice_lat(1) - 0.5_dp*dlat
        allocate(self%bin(self%nx, self%ny), self%wsum(nix, niy))
        self%bin = 0
        self%wsum = 0.0_dp
        do j = 1, self%ny
            do i = 1, self%nx
                ii = floor(modulo(self%grid_obj%lon(i, j) - lon0, 360.0_dp)/dlon) + 1
                jj = floor((self%grid_obj%lat(i, j) - lat0)/dlat) + 1
                if (ii < 1 .or. ii > nix .or. jj < 1 .or. jj > niy) cycle
                self%bin(i, j) = ii + (jj - 1)*nix
                self%wsum(ii, jj) = self%wsum(ii, jj) + self%grid_obj%area(i, j)
            end do
        end do
        if (any(self%wsum <= 0.0_dp)) then
            call fesm_error("kryos_climate_paleoemu::init", "the instance grid is coarser than the emulator's ice grid", &
                            "instance "//trim(self%name)//": "//to_str(count(self%wsum <= 0.0_dp))//" ice-grid cells "// &
                            "hold no cell centre of grid "//trim(self%grid))
        end if
    end subroutine bins_init

    subroutine to_ice_grid(self, X, Xi)
        ! The area-weighted mean over the instance cells of each ice-grid cell.
        class(paleoemu_climate_t), intent(in)  :: self
        real(dp),                  intent(in)  :: X(:,:)
        real(dp),                  intent(out) :: Xi(:,:)

        real(dp), allocatable :: acc(:)
        integer  :: i, j

        allocate(acc(size(Xi)))
        acc = 0.0_dp
        do j = 1, self%ny
            do i = 1, self%nx
                if (self%bin(i, j) > 0) acc(self%bin(i, j)) = acc(self%bin(i, j)) + self%grid_obj%area(i, j)*X(i, j)
            end do
        end do
        Xi = reshape(acc, shape(Xi))/self%wsum
    end subroutine to_ice_grid

    subroutine global_init_state_(self, st, pass)
        ! PASS_INIT: the reference climate (no anomaly). PASS_BOUNDARY (cold start): the reference ice from
        ! the core, the first update (its filters in equilibrium) on the geometry the core put (the background
        ! at t0). PASS_RESTART: the reference again, then the filter state of the bundle; the exports stand.
        class(paleoemu_climate_t), intent(inout) :: self
        type(step_time_t),         intent(in)    :: st
        integer,                   intent(in)    :: pass

        select case (pass)
        case (PASS_INIT)
            call exports_reference(self, st%t1)
        case (PASS_BOUNDARY)
            call reference_set(self, "init_state")
            call step(self, st, "init_state")
        case (PASS_RESTART)
            call reference_set(self, "init_state")
            if (.not. self%has_rst) call fesm_error("kryos_climate_paleoemu::init_state", "no restart file was read", &
                                                    "instance "//trim(self%name))
            self%emu%filt = self%filt_rst
            self%emu%now%time = self%time_rst
            self%emu%now%active = .true.
        end select
    end subroutine global_init_state_

    subroutine global_update_(self, st)
        class(paleoemu_climate_t), intent(inout) :: self
        type(step_time_t),         intent(in)    :: st
        call step(self, st, "update")
    end subroutine global_update_

    subroutine reference_set(self, what)
        ! The reference ice (the spliced present day of globe_ice) on the emulator's ice grid; the anchor on it.
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*),          intent(in)    :: what

        real(dp), allocatable :: H(:,:), zb(:,:), zs(:,:)

        call geometry_get(self, self%hi_H_ice_ref, self%hi_z_bed_ref, self%hi_z_srf_ref, H, zb, zs, what)
        call paleoemu_set_reference(self%emu, H, zb, zs)
        self%ref_set = .true.
        call anchor_set(self, H, zb, zs)
    end subroutine reference_set

    subroutine anchor_set(self, H, zb, zs)
        ! dto0: the dto of one update of a copy of the emulator (its filters, just restarted by the reference,
        ! start in equilibrium; the run's emulator is untouched) at present-day inputs: time 0, the emulator's
        ! own orbit table at 0 ka (the forcing's elements are not available at init), CO2 dto_anchor_co2, the
        ! reference ice H, zb, zs. Deterministic: a restart recomputes it bitwise. 0 without the anchor.
        class(paleoemu_climate_t), intent(inout) :: self
        real(dp),                  intent(in)    :: H(:,:), zb(:,:), zs(:,:)

        type(paleoemu_class), allocatable :: e0            ! on the heap (the emulator's matrices)
        integer :: k
        character(len=:), allocatable :: txt
        character(len=32) :: buf

        if (allocated(self%dto0)) deallocate(self%dto0)
        allocate(self%dto0(self%emu%nlon, self%emu%nlat, self%emu%nlev))
        self%dto0 = 0.0_dp
        if (.not. self%dto_anchor) return
        allocate(e0)
        e0 = self%emu
        call paleoemu_update(e0, 0.0_dp, self%dto_anchor_co2, H, zb, zs)
        self%dto0 = e0%now%dto
        deallocate(e0)
        txt = ""
        do k = 1, self%emu%nlev
            write(buf, "(f7.2,a,f6.2)") minval(self%dto0(:, :, k), mask=self%emu%ocean_mask(:, :, k) == 1), " to", &
                                        maxval(self%dto0(:, :, k), mask=self%emu%ocean_mask(:, :, k) == 1)
            txt = txt//" "//to_str(nint(self%emu%ocean_depth(k)))//" m"//trim(buf)//";"
        end do
        call log_msg(LOG_INFO, "kryos_climate_paleoemu", "dto anchor at present-day inputs (0 ka, CO2 "// &
                     to_str(self%dto_anchor_co2)//" ppm, the reference ice, filters in equilibrium), dto0 [K] over the "// &
                     "ocean:"//txt, "instance "//trim(self%name))
    end subroutine anchor_set

    subroutine step(self, st, what)
        ! One emulator update at the evaluation time (the forcing time with clock = "forcing") with the
        ! run's orbit (the forcing's elements; without them the emulator's own table at that time), the
        ! coefficient co2 and the ice geometry the core put; then the
        ! exports at t1.
        class(paleoemu_climate_t), intent(inout) :: self
        type(step_time_t),         intent(in)    :: st
        character(len=*),          intent(in)    :: what

        real(dp), allocatable :: H(:,:), zb(:,:), zs(:,:)
        real(dp) :: co2, orb(3)

        if (.not. self%ref_set) call fesm_error("kryos_climate_paleoemu::"//what, "no reference ice", &
                                                "instance "//trim(self%name)//": PASS_BOUNDARY sets it")
        call require_published(self, self%hi_coef(1), "coef_co2", what)
        call fields_get(self%imports, self%hi_coef(1), co2)
        call geometry_get(self, self%hi_H_ice, self%hi_z_bed, self%hi_z_srf, H, zb, zs, what)
        if (self%imports%f(self%hi_orbit(1))%published) then
            call fields_get(self%imports, self%hi_orbit(1), orb(1))
            call fields_get(self%imports, self%hi_orbit(2), orb(2))
            call fields_get(self%imports, self%hi_orbit(3), orb(3))
            call paleoemu_update(self%emu, st%t_eval, co2, H, zb, zs, orbit=orb)
        else
            call paleoemu_update(self%emu, st%t_eval, co2, H, zb, zs)
        end if
        call paleoemu_downscale(self%tgt, self%emu, scale=[self%scale_dtas, 1.0_dp], bias_scale=self%bias_scale, &
                                beta=self%beta_adj)
        call exports_put(self, st%t1)
        call ocean_put(self, st%t1)
    end subroutine step

    subroutine ocean_put(self, t1)
        ! dto_lev: scale_dto (dto - dto0) + dto_nsd dto_sd on the CESM grid, filled (ocean_fill), each level
        ! interpolated bilinearly onto the instance grid.
        class(paleoemu_climate_t), intent(inout) :: self
        real(dp),                  intent(in)    :: t1

        real(dp), allocatable :: d(:,:,:), x(:,:,:)
        integer :: k

        allocate(d(self%emu%nlon, self%emu%nlat, self%emu%nlev), x(self%nx, self%ny, self%emu%nlev))
        d = self%scale_dto*(self%emu%now%dto - self%dto0) + self%dto_nsd*self%emu%now%dto_sd
        call ocean_fill(self%emu%ocean_mask, d)
        do k = 1, self%emu%nlev
            x(:, :, k) = paleoemu_target_bilinear(self%tgt, d(:, :, k))
        end do
        call fields_put(self%exports, self%h_dto_lev, x, t1=t1)
    end subroutine ocean_put

    subroutine ocean_fill(mask, d)
        ! The cells of d [lon, lat, level] outside the ocean (mask /= 1: land, below the floor; the
        ! emulator's values there are 0) filled level by level from the ocean: a Poisson fill, periodic
        ! in longitude, from the zonal mean of the ocean as first guess (fill_poisson, method 1). A level
        ! without ocean takes the level above. The ocean cells keep their values.
        integer,  intent(in)    :: mask(:,:,:)
        real(dp), intent(inout) :: d(:,:,:)
        integer :: k

        do k = 1, size(d, 3)
            if (count(mask(:, :, k) == 1) == 0) then
                if (k == 1) call fesm_error("kryos_climate_paleoemu::ocean_fill", "the top level holds no ocean")
                d(:, :, k) = d(:, :, k - 1)
                cycle
            end if
            where (mask(:, :, k) /= 1) d(:, :, k) = MV
            call fill_poisson(d(:, :, k), MV, 1, .true.)
        end do
    end subroutine ocean_fill

    subroutine geometry_get(self, hH, hb, hs, H, zb, zs, what)
        ! Three imports (ice thickness, bed, surface on the instance grid) on the emulator's ice grid.
        class(paleoemu_climate_t), intent(in)  :: self
        integer,                   intent(in)  :: hH, hb, hs
        real(dp), allocatable,     intent(out) :: H(:,:), zb(:,:), zs(:,:)
        character(len=*),          intent(in)  :: what

        real(dp), allocatable :: a(:,:)
        integer  :: nix, niy

        nix = size(self%emu%ice_lon); niy = size(self%emu%ice_lat)
        allocate(H(nix, niy), zb(nix, niy), zs(nix, niy), a(self%nx, self%ny))
        call require_published(self, hH, self%imports%f(hH)%name, what)
        call require_published(self, hb, self%imports%f(hb)%name, what)
        call require_published(self, hs, self%imports%f(hs)%name, what)
        call fields_get(self%imports, hH, a);  call to_ice_grid(self, a, H)
        call fields_get(self%imports, hb, a);  call to_ice_grid(self, a, zb)
        call fields_get(self%imports, hs, a);  call to_ice_grid(self, a, zs)
    end subroutine geometry_get

    subroutine exports_put(self, t1)
        ! tas at the reference surface z_srf_clim, pr [kg m-2 yr-1], lapse = Gamma, beta_p.
        class(paleoemu_climate_t), intent(inout) :: self
        real(dp),                  intent(in)    :: t1

        real(dp), allocatable :: pr(:,:,:), lapse(:,:,:), bp(:,:)   ! on the heap, not as expression temporaries
                                                                    ! (25 MB on GLOBAL-0.5DEG: ifx puts those on the stack)

        allocate(pr(self%nx, self%ny, NMON), lapse(self%nx, self%ny, NMON), bp(self%nx, self%ny))
        pr = self%tgt%pr*self%day_year
        lapse = self%tgt%lapse
        bp = self%beta_p
        call fields_put(self%exports, self%h_tas,        self%tgt%tas,                t1=t1)
        call fields_put(self%exports, self%h_pr,         pr,                          t1=t1)
        call fields_put(self%exports, self%h_z_srf_clim, self%tgt%z_srf_ref,          t1=t1)
        call fields_put(self%exports, self%h_lapse,      lapse,                       t1=t1)
        call fields_put(self%exports, self%h_beta_p,     bp,                          t1=t1)
    end subroutine exports_put

    subroutine exports_reference(self, t1)
        ! The reference climate itself (PASS_INIT: no update yet), no ocean temperature change.
        class(paleoemu_climate_t), intent(inout) :: self
        real(dp),                  intent(in)    :: t1
        real(dp), allocatable :: x(:,:,:)
        self%tgt%tas = self%tgt%tas_ref
        self%tgt%pr  = self%tgt%pr_ref
        call exports_put(self, t1)
        allocate(x(self%nx, self%ny, self%emu%nlev))
        x = 0.0_dp
        call fields_put(self%exports, self%h_dto_lev, x, t1=t1)
    end subroutine exports_reference

    subroutine global_restart_write_(self, dir, t)
        ! The emulator's memory: its filter states and the time of its last update.
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*),          intent(in)    :: dir
        real(dp),                  intent(in)    :: t

        character(len=:), allocatable :: file

        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
        call nc_create(file)
        call nc_write_dim(file, "slow", x=1, nx=size(self%emu%filt, 1))
        call nc_write_dim(file, "tau",  x=1, nx=size(self%emu%filt, 2))
        call nc_write_dim(file, "one",  x=1, nx=1)
        call nc_write(file, "filt", self%emu%filt, dim1="slow", dim2="tau", &
                      long_name="filter states of log(CO2/CO2_ref) and the regional ice volume changes")
        call nc_write(file, "time_last", [self%emu%now%time], dim1="one", units="yr", &
                      long_name="time of the last update")
        call nc_write(file, "time_bundle", [t], dim1="one", units="yr")
    end subroutine global_restart_write_

    subroutine global_restart_read_(self, dir, t)
        ! The filter state; applied at PASS_RESTART, after the reference (paleoemu_set_reference restarts the filters).
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*),          intent(in)    :: dir
        real(dp),                  intent(in)    :: t

        character(len=:), allocatable :: file
        real(dp) :: tl(1)
        logical  :: found

        associate (unused => t); end associate
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
        inquire(file=file, exist=found)
        if (.not. found) call fesm_error("kryos_climate_paleoemu::restart_read", "restart file not found: "//file, &
                                         "instance "//trim(self%name))
        allocate(self%filt_rst(size(self%emu%filt, 1), size(self%emu%filt, 2)))
        call nc_read(file, "filt", self%filt_rst)
        call nc_read(file, "time_last", tl)
        self%time_rst = tl(1)
        self%has_rst = .true.
    end subroutine global_restart_read_

    subroutine global_write_output_(self, file, t, stream, init)
        ! Streams 2D and 1D: the default writer (the exports on lon, lat). `emu`: the emulator's state on the
        ! CESM grid (paleoemu_write_*: inputs, dV, dz, dtas, dlogpr, dto and their standard deviations).
        class(paleoemu_climate_t), intent(inout) :: self
        character(len=*),          intent(in)    :: file, stream
        real(dp),                  intent(in)    :: t
        logical,                   intent(in)    :: init

        select case (trim(stream))
        case ("emu")
            if (.not. self%emu%now%active) return
            if (init) call paleoemu_write_init(file, self%emu, t)
            call paleoemu_write_step(file, self%emu)
        case default
            call component_write_output(self, file, t, stream, init)
        end select
    end subroutine global_write_output_

    subroutine global_end_(self)
        class(paleoemu_climate_t), intent(inout) :: self
        call paleoemu_end(self%emu)
    end subroutine global_end_

    subroutine require_published(self, h, name, what)
        class(climate_t), intent(in) :: self
        integer,          intent(in) :: h
        character(len=*), intent(in) :: name, what
        if (h <= 0) then
            call fesm_error("kryos_climate_paleoemu::"//what, "import "//trim(name)//" is not declared", &
                            "instance "//trim(self%name))
        else if (.not. self%imports%f(h)%published) then
            call fesm_error("kryos_climate_paleoemu::"//what, "import "//trim(name)//" has not been published", &
                            "instance "//trim(self%name)//": the core (globe_ice), the forcing or the global part "// &
                            "puts it before this call")
        end if
    end subroutine require_published

    ! =========================================================================
    ! Regional port
    ! =========================================================================

    subroutine port_init_(self, ctx)
        ! On the grid the chain gives it; it holds the global part's atmosphere, remapped. Its cadence is the
        ! global part's (kryos_run). With keys to_ref and so_ref: the present-day ocean on its grid, its
        ! axis `depth`, the weights from the global part's axis `level` (in ctx%axes: domain_create).
        class(paleoemu_port_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        self%prec = dp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = climate_opt_t(port=.true.)
        if (nml_has_param(ctx%par_file, ctx%instance, "to_ref")) call nml_read(ctx%par_file, ctx%instance, "to_ref", self%to_ref)
        if (nml_has_param(ctx%par_file, ctx%instance, "so_ref")) call nml_read(ctx%par_file, ctx%instance, "so_ref", self%so_ref)
        if ((len_trim(self%to_ref) == 0) .neqv. (len_trim(self%so_ref) == 0)) then
            call fesm_error("kryos_climate_paleoemu::init", "keys `to_ref` and `so_ref` go together", &
                            "instance "//trim(self%name)//": both name the present-day ocean (var_defs groups), "// &
                            "or both are empty (no ocean)")
        end if
        if (len_trim(self%to_ref) == 0) then
            if (nml_has_param(ctx%par_file, ctx%instance, "dto_basins")) then
                call nml_read(ctx%par_file, ctx%instance, "dto_basins", self%dto_basins)
                if (len_trim(self%dto_basins) > 0) &
                    call fesm_error("kryos_climate_paleoemu::init", "key `dto_basins` needs the port's ocean", &
                                    "instance "//trim(self%name)//": set to_ref and so_ref")
            end if
            return
        end if

        call ocean_ref_load(self, ctx)
        self%caps%has_ocn = .true.
        self%opt%ocn_lev  = .true.
        if (nml_has_param(ctx%par_file, ctx%instance, "dto_basins")) &
            call nml_read(ctx%par_file, ctx%instance, "dto_basins", self%dto_basins)
        if (nml_has_param(ctx%par_file, ctx%instance, "dto_sample")) &
            call nml_read(ctx%par_file, ctx%instance, "dto_sample", self%dto_sample)
        if (len_trim(self%dto_basins) > 0) call basins_init(self, ctx)
        call log_msg(LOG_INFO, "kryos_climate_paleoemu", "port ocean: present day "//trim(self%to_ref)//", "// &
                     trim(self%so_ref)//" on "//to_str(size(self%to0, 3))//" depths "//to_str(self%axes(1)%v(1))//"-"// &
                     to_str(self%axes(1)%v(size(self%axes(1)%v)))//" m; + dto_lev of "//to_str(size(self%W, 2))// &
                     " levels", "instance "//trim(self%name))
    end subroutine port_init_

    subroutine ocean_ref_load(self, ctx)
        ! The present-day ocean of the groups to_ref, so_ref (the par file, else the first library file
        ! holding them; the keys of a snapshot's to and so: units, scaling, fill), its annual mean, on the
        ! file's depths (positive down: the axis `depth` of set A); the weights W [depth, level] from the
        ! global part's levels, linear in depth, constant above the top and below the bottom level.
        class(paleoemu_port_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        character(len=*), parameter :: proc = "kryos_climate_paleoemu::init"
        real(dp), allocatable :: zt(:), zs(:), lev(:)
        type(axis_t) :: ax
        integer :: k, l, nz, nl

        call ocean_read(self, ctx, self%to_ref, self%to0, zt)
        call ocean_read(self, ctx, self%so_ref, self%so0, zs)
        if (size(zt) /= size(zs) .or. any(zt /= zs)) then
            call fesm_error(proc, "to_ref and so_ref are on different depths", "instance "//trim(self%name)// &
                            ": groups "//trim(self%to_ref)//", "//trim(self%so_ref))
        end if

        ax%name  = "depth"
        ax%units = "m"
        ax%v     = zt
        if (.not. allocated(self%axes)) allocate(self%axes(0))
        self%axes = [self%axes, ax]

        k = 0
        if (allocated(ctx%axes)) k = axis_find(ctx%axes, "level")
        if (k == 0) then
            call fesm_error(proc, "the port's ocean needs the global part's axis `level`", "instance "// &
                            trim(self%name)//": a global paleoemu climate ([globe] climate) passes it to the domains")
        end if
        lev = ctx%axes(k)%v
        nz = size(zt);  nl = size(lev)
        allocate(self%W(nz, nl))
        self%W = 0.0_dp
        do k = 1, nz
            if (zt(k) <= lev(1)) then
                self%W(k, 1) = 1.0_dp
            else if (zt(k) >= lev(nl)) then
                self%W(k, nl) = 1.0_dp
            else
                l = count(lev <= zt(k))                        ! lev(l) <= zt(k) < lev(l+1)
                self%W(k, l + 1) = (zt(k) - lev(l))/(lev(l + 1) - lev(l))
                self%W(k, l)     = 1.0_dp - self%W(k, l + 1)
            end if
        end do
    end subroutine ocean_ref_load

    subroutine basins_init(self, ctx)
        ! The basins of dto_basins on the port's grid (regions_basin_ids), extended over the whole domain:
        ! a cell outside every basin (the ocean beyond the continental shelf, where the basin sets end)
        ! takes the id of its nearest basin cell (fill_nearest), so that the sample at the shelf break and
        ! the shelves the ice may advance onto belong to a basin. The sampled cells (regions_select of
        ! dto_sample); a basin with none samples its ocean cells instead (the fallback: the port imports no
        ! geometry, so the reference topography of the init context, ctx%ref: bed below sea level and not
        ! grounded, rho_ice H < -rho_sw z_bed); each basin's count of them. A basin without ocean either
        ! keeps its cells' values.
        class(paleoemu_port_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        character(len=:), allocatable :: txt, txt_fb
        logical, allocatable :: ocean(:,:)
        integer :: b

        self%basin  = regions_basin_ids(ctx%reg, trim(self%dto_basins))
        if (all(self%basin <= 0)) then
            call fesm_error("kryos_climate_paleoemu::init", "dto_basins `"//trim(self%dto_basins)//"` holds no basin", &
                            "instance "//trim(self%name)//" on grid "//trim(ctx%grid_name))
        end if
        call fill_nearest(self%basin, 0)
        self%sample = regions_select(ctx%reg, trim(self%dto_sample))
        allocate(self%nsample(max(maxval(self%basin), 0)))
        txt = "";  txt_fb = ""
        do b = 1, size(self%nsample)
            self%nsample(b) = count(self%basin == b .and. self%sample)
            if (count(self%basin == b) == 0 .or. self%nsample(b) > 0) cycle
            if (.not. allocated(ocean)) ocean = ref_ocean(self, ctx)
            self%nsample(b) = count(self%basin == b .and. ocean)
            if (self%nsample(b) > 0) then
                self%sample = self%sample .or. (self%basin == b .and. ocean)
                txt_fb = txt_fb//" "//to_str(b)//" ("//to_str(self%nsample(b))//")"
            else
                txt = txt//" "//to_str(b)
            end if
        end do
        call log_msg(LOG_INFO, "kryos_climate_paleoemu", "port ocean: dto_ann the basin means of "// &
                     trim(self%dto_basins)//" ("//to_str(size(self%nsample))//" basins) over `"//trim(self%dto_sample)// &
                     "` ("//to_str(sum(self%nsample))//" cells); basins without a sample, over their ocean cells of the "// &
                     "reference topography (cells):"//merge(txt_fb, " none", len(txt_fb) > 0)// &
                     "; basins without ocean either keep their cells' values:"//merge(txt, " none", len(txt) > 0), &
                     "instance "//trim(self%name))
    end subroutine basins_init

    function ref_ocean(self, ctx) result(ocean)
        ! The ocean of the reference topography on the port's grid (ctx%ref, sea level 0): bed below sea level
        ! and not grounded, rho_ice H_ice < -rho_sw z_bed (H_ice = 0 included).
        class(paleoemu_port_t), intent(in) :: self
        type(init_ctx_t),       intent(in) :: ctx
        logical, allocatable :: ocean(:,:)

        if (.not. (allocated(ctx%ref%z_bed) .and. allocated(ctx%ref%H_ice))) then
            call fesm_error("kryos_climate_paleoemu::init", "the fallback of an unsampled basin needs the reference "// &
                            "topography", "instance "//trim(self%name)//" on grid "//trim(ctx%grid_name)// &
                            ": ctx%ref (z_bed, H_ice) is not set")
        end if
        ocean = ctx%ref%z_bed < 0.0_dp .and. ctx%cnst%rho_ice*ctx%ref%H_ice < -ctx%cnst%rho_sw*ctx%ref%z_bed
    end function ref_ocean

    subroutine ocean_read(self, ctx, group, x, z)
        ! One ocean field of the var_defs group `group` on the port's grid: its annual mean (a monthly
        ! climatology is averaged), widened, and its depths |z| [m].
        class(paleoemu_port_t), intent(in)  :: self
        type(init_ctx_t),       intent(in)  :: ctx
        character(len=*),       intent(in)  :: group
        real(dp), allocatable,  intent(out) :: x(:,:,:), z(:)

        type(varslice_class) :: vs
        character(len=512)   :: vdb
        integer :: nz, nt

        vdb = group_file(trim(ctx%par_file), ctx%library, trim(group), "&"//trim(self%name)//" "//trim(group))
        call varslice_init_nml(vs, trim(vdb), trim(group), domain=trim(ctx%domain), grid_name=trim(ctx%grid_name), &
                               grid=ctx%grid)
        if (vs%par%with_time) then
            call varslice_update(vs, [ctx%time_rel], method="extrap")
        else
            call varslice_update(vs)
        end if
        nz = size(vs%var, 3);  nt = size(vs%var, 4)
        if (size(vs%var, 1) /= self%nx .or. size(vs%var, 2) /= self%ny .or. nz < 2 .or. (nt /= 1 .and. nt /= NMON)) then
            call fesm_error("kryos_climate_paleoemu::init", "group "//trim(group)//" is not an ocean field on the grid", &
                            "instance "//trim(self%name)//": [x, y, depth] (and 1 or 12 months) on "//trim(ctx%grid_name)// &
                            " expected")
        end if
        x = real(sum(vs%var, dim=4)/real(nt, kind(vs%var)), dp)
        z = abs(real(vs%z(1:nz), dp))
        call varslice_end(vs)
    end subroutine ocean_read

    subroutine port_init_state_(self, st, pass)
        ! PASS_BOUNDARY: the imports the core passed (cold start). Other passes: nothing (PASS_RESTART: the
        ! bundle's exports stand).
        class(paleoemu_port_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st
        integer,                intent(in)    :: pass
        associate (unused => st%t1); end associate
        if (pass == PASS_BOUNDARY) call copy_imports(self, "init_state")
    end subroutine port_init_state_

    subroutine port_update_(self, st)
        class(paleoemu_port_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st
        associate (unused => st%t1); end associate
        call copy_imports(self, "update")
    end subroutine port_update_

    subroutine copy_imports(self, what)
        ! tas, pr, z_srf_clim, lapse, beta_p as remapped, with their validity; with an ocean, the profiles.
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*),       intent(in)    :: what

        real(dp), allocatable :: a3(:,:,:), a2(:,:)
        real(dp) :: t0, t1

        allocate(a3(self%nx, self%ny, NMON), a2(self%nx, self%ny))

        call require_published(self, self%hi_tas, "tas", what)
        call fields_get(self%imports, self%hi_tas, a3, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_tas, a3, t1=t1, t0=t0)
        call require_published(self, self%hi_pr, "pr", what)
        call fields_get(self%imports, self%hi_pr, a3, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_pr, a3, t1=t1, t0=t0)
        call require_published(self, self%hi_lapse, "lapse", what)
        call fields_get(self%imports, self%hi_lapse, a3, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_lapse, a3, t1=t1, t0=t0)
        call require_published(self, self%hi_z_srf_clim, "z_srf_clim", what)
        call fields_get(self%imports, self%hi_z_srf_clim, a2, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_z_srf_clim, a2, t1=t1, t0=t0)
        call require_published(self, self%hi_beta_p, "beta_p", what)
        call fields_get(self%imports, self%hi_beta_p, a2, t0=t0, t1=t1)
        call fields_put(self%exports, self%h_beta_p, a2, t1=t1, t0=t0)
        if (self%opt%ocn_lev) call ocean_profiles(self, what)
    end subroutine copy_imports

    subroutine ocean_profiles(self, what)
        ! dto_ann = W dto_lev in each cell, to_ann = to_ref + dto_ann, so_ann = so_ref, with the validity
        ! of the import.
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*),       intent(in)    :: what

        real(dp), allocatable :: dl(:,:,:), dz(:,:,:)
        real(dp) :: t0, t1
        integer  :: i, j

        allocate(dl(self%nx, self%ny, size(self%W, 2)), dz(self%nx, self%ny, size(self%W, 1)))
        call require_published(self, self%hi_dto_lev, "dto_lev", what)
        call fields_get(self%imports, self%hi_dto_lev, dl, t0=t0, t1=t1)
        if (allocated(self%basin)) call basin_means(self, dl)
        do j = 1, self%ny
            do i = 1, self%nx
                dz(i, j, :) = matmul(self%W, dl(i, j, :))
            end do
        end do
        call fields_put(self%exports, self%h_dto_ann, dz, t1=t1, t0=t0)
        call fields_put(self%exports, self%h_to_ann, self%to0 + dz, t1=t1, t0=t0)
        call fields_put(self%exports, self%h_so_ann, self%so0, t1=t1, t0=t0)
    end subroutine ocean_profiles

    subroutine basin_means(self, dl)
        ! dl [x, y, level] replaced in each basin by the mean over its sampled cells (the fallback's included),
        ! level by level (W is linear: the same as the mean of dto_ann); basins without any keep their values.
        class(paleoemu_port_t), intent(in)    :: self
        real(dp),               intent(inout) :: dl(:,:,:)

        real(dp) :: m
        integer  :: b, k

        do b = 1, size(self%nsample)
            if (self%nsample(b) == 0) cycle
            do k = 1, size(dl, 3)
                m = sum(dl(:, :, k), mask=self%basin == b .and. self%sample)/real(self%nsample(b), dp)
                where (self%basin == b) dl(:, :, k) = m
            end do
        end do
    end subroutine basin_means

    subroutine port_restart_write_(self, dir, t)
        ! Stateless: constraints restart = .false., so the core never calls this.
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine port_restart_write_

    subroutine port_restart_read_(self, dir, t)
        class(paleoemu_port_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t
        associate (unused => t + len(dir) + self%nx); end associate
    end subroutine port_restart_read_

end module kryos_climate_paleoemu
