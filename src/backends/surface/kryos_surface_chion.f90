module kryos_chion_host
    ! chion as the yelmox surface model (surface_method = "chion").
    !
    ! chion (fesmc/chion) steps independent snowpack columns one step at a time
    ! and leaves the time loop to the host. This wrapper is that host: it packs
    ! the (nx,ny) fields of grid_surface into chion's column list, drives one
    ! annual cycle of daily steps from monthly climate, and aggregates the result
    ! to the annual fields yelmox lands on Yelmo (dom%schn%ann), in smbpal's
    ! units so couple_surface_to_yelmo treats both models alike.
    !
    ! Inputs chion does not own, supplied here (as smbpal does internally):
    !   - daily top-of-atmosphere insolation (ITM, BESSI; fesmc/insol, on the orbit the
    !     host passes: the run's, from the forcing) from a per-day latitude table (insol_dlat)
    !     interpolated to the columns, as shortwave_down: as is
    !     for ITM (which applies its own transmissivity); for BESSI (surface
    !     shortwave) times the clear-sky transmissivity of chion's cloud proxy
    !     at the surface elevation and one cloud fraction,
    !     tau_clear(z)*(1 - cloud); and as the TOA itself (toa_shortwave) for
    !     BESSI's cloud-proxy longwave, so chion does not divide by its own
    !     fixed-orbit 365-day TOA;
    !   - the annual positive degree days (ITM only: its vegetation proxy for the
    !     critical snow depth and the snow-free land albedo), with chion's
    !     Calov-Greve integral on every pdd_dday-th day of the daily temperature,
    !     as smbpal does, fixed for the year;
    !   - the snowfall/rainfall split (smbpal's calc_snowfrac);
    !   - the wind speed from the surface elevation and slope (below), a constant
    !     relative humidity, the air pressure from the
    !     surface elevation (barometric, with the annual-mean air temperature, as
    !     Chion.jl's air_pressure_from_surface_height) and the solar longitude of
    !     the day (BESSI's energy balance; the climate supplies none of them).
    ! Monthly -> daily is chion_forcing_monthly (mean-preserving).
    !
    ! Daily temperature variability (optional, tvar_sigma_max > 0): the smooth daily
    ! air temperature plus sigma*xi. xi is a red-noise sequence (AR(1), e-folding
    ! tvar_tau days) of unit variance and zero mean in every month, the same for all
    ! columns and every annual cycle (seed tvar_seed; a portable Park-Miller generator),
    ! so the monthly means are unchanged and the spin-up equilibrates. sigma follows
    ! the column's monthly mean temperature, lower over melting surfaces (Seguinot and
    ! Rogozhina, 2014, J. Glaciol. 60, 603):
    !     sigma = min(max(tvar_sigma0 + tvar_sigma_slope*(T_mon - T0), tvar_sigma_min), tvar_sigma_max).
    ! ERA5 1981-2010 over the Greenland and Antarctic ice: sigma = 2.0 - 0.14 (T - T0) in [2, 6] K.
    !
    ! Wind (BESSI's turbulent fluxes; the climate supplies none), constant through
    ! the year and set each annual cycle from the geometry:
    !     U = wind_speed + wind_lapse*max(z_srf,0)[km] + wind_slope*|grad z_srf|[m/km],
    ! the slope a centred difference over a fixed baseline of about
    ! wind_slope_length (2n cells, n = max(1, nint(L/(2 dx)))), so that it does not
    ! depend on the grid's resolution. A first-order stand-in for the katabatic
    ! wind; the coefficients are per domain (ERA5 1981-2010 summer 10-m wind over
    ! the ice: Greenland 2.5 + 1.27 z[km], Antarctica 5.5 + 0.2 z[km]).
    !
    ! Active columns: land or ice, not open ocean (H_ice > 0 or z_bed > z_sl),
    ! refreshed every annual cycle. A column switched on starts from chion's
    ! cold state. Inactive columns: smb = 0, tsrf = annual-mean air temperature.
    !
    ! smb is the surface mass balance, sf + rf - runoff - subl (smbpal's smb),
    ! not chion_get_smb's ice-facing flux. tsrf: for ITM and PDD chion's
    ! chion_tsrf_annual on the annual means: over ice the annual-mean air
    ! temperature warmed by the gross annual refreezing, capped at freezing,
    ! min(T0, t_air + firn_fac*max(refrz,0)), elsewhere the air temperature.
    ! This is intended: smbpal PDD did the same, while smbpal ITM warmed by
    ! the net refreezing max(refrz - melt, 0), i.e. less where it melts.
    ! BESSI: the annual mean of its per-step surface temperature
    ! (chion_get_surface). alb_s is the annual mean of
    ! chion's per-step albedo (missing for PDD, which has none).
    !
    ! The per-column host work (forcing, accumulation, PDDs, insolation lookup)
    ! is OpenMP-parallel over the active columns, like chion's own step.
    !
    ! Supported chion models: "itm", "bessi" and "pdd". BESSI's downwelling
    ! longwave is chion's own (longwave_scheme; default "cloud_proxy", whose
    ! cloudiness 1 - SWdn/(TOA*tau_clear(z)) is then the key cloud everywhere).
    ! The cloud fraction is per domain: MAR v3.14 1981-2010 over Greenland's
    ! ice gives SWdn/TOA = 0.61 + 0.05 z[km] in summer, cloud = 0.28.
    ! PDD takes the air temperature and the snowfall/rainfall only.
    !
    ! Optionally (BESSI, surface_chion_seb_enable) the monthly means of the
    ! last annual cycle of BESSI's surface energy balance, albedo, surface
    ! temperature and melt (sc%seb): chion_get_surface_fluxes after every
    ! step, by the month of the step's first day.

    use nml,        only : nml_read, nml_has_param
    use ncio,       only : nc_open, nc_close, nc_write, nc_read, nc_exists_var
    use phys_constants, only : phys_const_class, sec_day
    use insolation, only : insol_class, insol_init, calc_insol_days, insol_wp => wp, &
                           insol_stencil_class, insol_stencil_init, insol_stencil_point
    use insolation, only : orbit_pars_class
    use chion,      only : wp, dp, wp_acc, MV, chion_class, chion_init, chion_init_state, &
                           chion_update, chion_set_active_mask, &
                           chion_get_surface_flux_totals, chion_get_surface, &
                           chion_get_surface_fluxes, &
                           chion_restart_write, chion_restart_read, &
                           monthly_to_daily_class, monthly_to_daily_init, &
                           monthly_to_daily_controls, interp_monthly_to_day, &
                           pdd_expected_positive_temperature, chion_tsrf_annual, &
                           clear_sky_transmissivity

    implicit none

    private

    integer,  parameter :: nmon       = 12
    integer,  parameter :: nday_mon   = 30       ! 360-day year, as smbpal and the climate
    real(dp), parameter :: insol_dlat = 0.1_dp   ! [deg] latitude spacing of the insolation table
    integer,  parameter :: pdd_dday   = 10       ! [d] sampling of the annual PDD sum (as smbpal)
    real(wp), parameter :: p_sl       = 101325.0_wp  ! [Pa] sea-level air pressure (Chion.jl)
    real(wp), parameter :: R_dry      = 287.058_wp   ! [J kg-1 K-1] dry-air gas constant (Chion.jl R/M_air)
    real(wp), parameter :: days_trop  = 365.2422_wp  ! [d] tropical year, for the solar longitude

    type surface_chion_param_class
        real(wp)           :: sf_a, sf_b    ! snowfall fraction, -0.5*tanh(sf_a*(T-sf_b))+0.5
        real(wp)           :: sigma_pdd     ! [K] temperature standard deviation for the PDDs
        real(wp)           :: time_equil    ! [yr] cold-start snowpack spin-up
        integer            :: dt_days       ! [d] chion step; divides the 360-day year
        real(wp)           :: cloud         ! [1] cloud fraction (BESSI): SWdn = TOA*tau_clear(z)*(1 - cloud)
        real(wp)           :: wind_speed    ! [m s-1] wind speed at sea level on flat ground
        real(wp)           :: wind_lapse        = 0.0_wp   ! [m s-1 km-1] its change with the surface elevation
        real(wp)           :: wind_slope        = 0.0_wp   ! [m s-1 per m/km] its change with the surface slope
        real(wp)           :: wind_slope_length = 50.0_wp  ! [km] the baseline of the slope
        real(wp)           :: rel_hum       ! [1] constant relative humidity
        real(wp)           :: swe_init_ice = 0.0_wp   ! [kg m-2] cold start: the snowpack of ice columns (0: chion's own)
        real(wp)           :: tvar_sigma0      = 0.0_wp   ! [K] daily temperature variability at T_mon = T0
        real(wp)           :: tvar_sigma_slope = 0.0_wp   ! [K K-1] its change with the monthly mean temperature
        real(wp)           :: tvar_sigma_min   = 0.0_wp   ! [K] its lower bound
        real(wp)           :: tvar_sigma_max   = 0.0_wp   ! [K] its upper bound (0: no variability)
        real(wp)           :: tvar_tau         = 3.0_wp   ! [d] e-folding time of the red noise
        integer            :: tvar_seed        = 1        ! seed of the noise sequence
    end type surface_chion_param_class

    type surface_chion_ann_class
        ! Annual means on grid_surface. Mass fluxes in [mm w.e./yr].
        real(wp), allocatable :: mask(:,:)      ! [1] chion active (land or ice)
        real(wp), allocatable :: t2m(:,:)       ! [K]
        real(wp), allocatable :: pr(:,:)        ! [mm w.e./yr]
        real(wp), allocatable :: sf(:,:)        ! [mm w.e./yr]
        real(wp), allocatable :: S(:,:)         ! [W m-2] TOA insolation
        real(wp), allocatable :: PDDs(:,:)      ! [K d]
        real(wp), allocatable :: melt(:,:)      ! [mm w.e./yr]
        real(wp), allocatable :: runoff(:,:)    ! [mm w.e./yr]
        real(wp), allocatable :: refrz(:,:)     ! [mm w.e./yr]
        real(wp), allocatable :: smb(:,:)       ! [mm w.e./yr] sf + rf - runoff - subl
        real(wp), allocatable :: tsrf(:,:)      ! [K] surface temperature
        real(wp), allocatable :: alb_s(:,:)     ! [1] surface albedo
        real(wp), allocatable :: H_snow(:,:)    ! [mm w.e.] end-of-year snowpack (itm, pdd)
    end type surface_chion_ann_class

    ! The monthly surface energy balance (sc%seb%f(:,:,month,k)): BESSI's terms
    ! [W m-2], positive into the surface, the albedo, the surface temperature
    ! and the melt rate.
    integer, parameter, public :: NSEB = 10
    character(len=4), parameter, public :: SEB_NAME(NSEB) = &
        [character(len=4) :: "swd", "swn", "lwd", "lwn", "shf", "lhf", "prh", "alb", "tsrf", "melt"]
    character(len=10), parameter, public :: SEB_UNITS(NSEB) = &
        [character(len=10) :: "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", "W m-2", "1", "K", "mm w.e./d"]
    character(len=40), parameter, public :: SEB_LONG(NSEB) = [character(len=40) :: &
        "Shortwave down", "Net shortwave", "Longwave down", "Net longwave", "Sensible heat flux", &
        "Turbulent latent heat flux", "Heat of the precipitation", "Surface albedo", "Surface temperature", "Melt"]

    type surface_chion_seb_class
        logical               :: on = .false.
        real(wp), allocatable :: f(:,:,:,:)            ! (nx,ny,nmon,NSEB), MV off the active columns
    end type surface_chion_seb_class

    type surface_chion_insol_class
        ! Daily TOA insolation on a regular latitude table, rebuilt only when
        ! the orbital elements change, and the fixed linear-interpolation
        ! stencil of every column into it (insol's latitude stencil). The state
        ! has no orbital tables: the host passes the elements.
        type(insol_class)      :: is                   ! insol state on the table latitudes (no tables)
        type(orbit_pars_class) :: orbit                ! the elements of the table
        logical               :: set = .false.         ! the table holds the insolation of orbit
        real(dp), allocatable :: lat(:)                ! (nlat) [deg N]
        real(wp), allocatable :: S(:,:)                ! (nlat,nday) [W m-2]
        type(insol_stencil_class) :: st                ! the table latitudes onto the columns
    end type surface_chion_insol_class

    type surface_chion_class
        type(surface_chion_param_class) :: par
        type(chion_class)               :: chn
        type(monthly_to_daily_class)    :: md
        type(surface_chion_insol_class) :: ins
        type(surface_chion_ann_class)   :: ann
        type(surface_chion_seb_class)   :: seb
        integer                         :: nx, ny, ncol
        real(wp)                        :: dx       ! [m] grid spacing (0: not a metric grid, no slope)
        real(wp), allocatable           :: xi(:)    ! [1] the daily temperature noise, one value per step of the year
        logical                         :: itm      ! model = "itm" (ITM-only forcing)
        logical                         :: insol    ! the model takes insolation (itm, bessi)
        logical                         :: bessi    ! model = "bessi" (surface shortwave from the cloud fraction)
    end type surface_chion_class

    public :: surface_chion_class
    public :: surface_chion_init
    public :: surface_chion_spinup
    public :: surface_chion_seb_enable
    public :: surface_chion_update
    public :: surface_chion_restart_write
    public :: surface_chion_restart_read

contains

    subroutine surface_chion_init(sc, filename, lats, dx, group, chion_group, cnst)
        ! Parameters, chion (cold state) and the annual fields on an (nx,ny) grid.
        ! chion takes its shared physical constants from the domain's record.
        type(surface_chion_class), intent(inout) :: sc
        character(len=*),          intent(in)    :: filename
        real(wp),                  intent(in)    :: lats(:,:)      ! [deg N]
        real(wp),                  intent(in)    :: dx             ! [m] grid spacing; 0 if not a metric grid
        character(len=*),          intent(in)    :: group          ! &surface_chion
        character(len=*),          intent(in)    :: chion_group    ! &chion
        type(phys_const_class),    intent(in)    :: cnst

        call surface_chion_par_load(sc%par, filename, group)

        sc%nx   = size(lats,1)
        sc%ny   = size(lats,2)
        sc%ncol = sc%nx*sc%ny
        sc%dx   = dx
        if (sc%par%wind_slope /= 0.0_wp .and. .not. sc%dx > 0.0_wp) then
            write(*,*) "surface_chion_init:: error: wind_slope needs a grid in metres or km; dx = ", dx
            stop 1
        end if

        call chion_init(sc%chn, filename, sc%ncol, group=chion_group, cnst=cnst)

        select case(trim(sc%chn%par%model))
            case("itm")
                sc%itm   = .true.
                sc%bessi = .false.
                sc%insol = .true.
            case("bessi")
                sc%itm   = .false.
                sc%bessi = .true.
                sc%insol = .true.
            case("pdd")
                sc%itm   = .false.
                sc%bessi = .false.
                sc%insol = .false.
            case default
                write(*,*) "surface_chion_init:: error: chion model must be itm, bessi or pdd &
                           &(the host supplies no forcing for others); got "//trim(sc%chn%par%model)
                stop 1
        end select

        call monthly_to_daily_init(sc%md, nmon, nday_mon)

        if (sc%par%dt_days < 1 .or. mod(sc%md%nday_year, sc%par%dt_days) /= 0) then
            write(*,*) "surface_chion_init:: error: dt_days must divide the year; dt_days, nday_year = ", &
                       sc%par%dt_days, sc%md%nday_year
            stop 1
        end if

        allocate(sc%xi(sc%md%nday_year/sc%par%dt_days))
        sc%xi = 0.0_wp
        if (sc%par%tvar_sigma_max > 0.0_wp) &
            sc%xi = red_noise_monthly(sc%md%nday_year, sc%par%dt_days, sc%par%tvar_tau, sc%par%tvar_seed)

        ! Column icol = i + (j-1)*nx: reshape order.
        sc%chn%forc%latitude_deg = reshape(lats, [sc%ncol])
        if (sc%insol) call insol_table_init(sc%ins, sc%chn%forc%latitude_deg, sc%md%nday_year)

        ! Constant for the run; the wind until the first annual cycle sets it from the geometry.
        sc%chn%forc%wind_speed            = sc%par%wind_speed
        sc%chn%forc%relative_humidity     = sc%par%rel_hum
        sc%chn%forc%has_relative_humidity = .true.
        sc%chn%forc%has_toa_shortwave     = .true.

        call chion_init_state(sc%chn)

        call ann_alloc(sc%ann, sc%nx, sc%ny)
    end subroutine surface_chion_init

    subroutine surface_chion_seb_enable(sc)
        ! Keep the monthly surface energy balance of every annual cycle (BESSI only: the
        ! other models have none).
        type(surface_chion_class), intent(inout) :: sc

        if (trim(sc%chn%par%model) /= "bessi") then
            write(*,*) "surface_chion_seb_enable:: error: the surface energy balance is BESSI's; model = "// &
                       trim(sc%chn%par%model)
            stop 1
        end if
        sc%seb%on = .true.
        allocate(sc%seb%f(sc%nx, sc%ny, nmon, NSEB))
        sc%seb%f = MV
    end subroutine surface_chion_seb_enable

    subroutine surface_chion_spinup(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, orbit)
        ! Cold-start snowpack: time_equil annual cycles under the initial climate
        ! and geometry (smbpal's 100-yr ITM equilibration). With swe_init_ice > 0
        ! the cold start is rebuilt first, now that the climate and the ice are
        ! known: ice columns (H_ice > 0) start with that much snow and firn, land
        ! columns bare, everything (snow, surface, BESSI's ice substrate) at the
        ! column's annual-mean air temperature (chion_init_state swe, t_init).
        ! An ice sheet is never snow-free; started bare, a very dry column can
        ! settle into a thin, sublimating state (chion D46).
        type(surface_chion_class), intent(inout) :: sc
        real(wp), intent(in) :: t2m(:,:,:), pr(:,:,:)       ! [K], [mm w.e./d] monthly
        real(wp), intent(in) :: z_srf(:,:), z_bed(:,:), H_ice(:,:), z_sl(:,:)
        type(orbit_pars_class), intent(in), optional :: orbit   ! the orbital elements (required with insolation)

        integer :: n

        if (sc%par%swe_init_ice > 0.0_wp) then
            call chion_init_state(sc%chn, swe=reshape(merge(sc%par%swe_init_ice, 0.0_wp, H_ice > 0.0_wp), [sc%ncol]), &
                                  t_init=reshape(sum(t2m, dim=3)/real(size(t2m, 3), wp), [sc%ncol]))
        end if

        do n = 1, nint(sc%par%time_equil)
            call surface_chion_update(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, orbit)
        end do
    end subroutine surface_chion_spinup

    subroutine surface_chion_update(sc, t2m, pr, z_srf, z_bed, H_ice, z_sl, orbit)
        ! One annual cycle of daily chion steps from monthly climate, aggregated
        ! to sc%ann.
        type(surface_chion_class), intent(inout) :: sc
        real(wp), intent(in) :: t2m(:,:,:)          ! [K] monthly near-surface air temperature
        real(wp), intent(in) :: pr(:,:,:)           ! [mm w.e./d] monthly precipitation
        real(wp), intent(in) :: z_srf(:,:)          ! [m]
        real(wp), intent(in) :: z_bed(:,:)          ! [m]
        real(wp), intent(in) :: H_ice(:,:)          ! [m]
        real(wp), intent(in) :: z_sl(:,:)           ! [m]
        type(orbit_pars_class), intent(in), optional :: orbit   ! the orbital elements (required with insolation)

        real(wp), allocatable :: t_ctl(:,:), p_ctl(:,:), mon(:,:), ctl(:,:), col(:), t_ann(:)
        real(wp), allocatable :: t_sum(:), pr_sum(:), sf_sum(:), S_sum(:), PDDs(:)
        real(wp), allocatable :: tsrf_c(:), alb_c(:), ts_sum(:), alb_sum(:)
        real(wp), allocatable :: t_air_c(:), refrz_c(:), tsrf_ann_c(:)
        real(wp_acc), allocatable :: melt0(:), runoff0(:), refrz0(:), subl0(:)
        real(wp_acc), allocatable :: melt1(:), runoff1(:), refrz1(:), subl1(:)
        real(wp), allocatable :: melt_a(:), runoff_a(:), refrz_a(:), subl_a(:)
        real(wp), allocatable :: fx(:,:), seb_sum(:,:,:), sw_fac(:), sig(:,:)
        integer,  allocatable :: idx(:)
        integer  :: na, i, icol, m, k, day, nday
        real(wp) :: dt, spd, T0, t_d, p_d, sf_d, S_d

        nday = sc%md%nday_year
        dt   = real(sc%par%dt_days, wp)
        spd  = real(sec_day, wp)
        T0   = sc%chn%c%T0

        ! Active columns, from the current geometry: land or ice.
        call chion_set_active_mask(sc%chn, reshape(H_ice > 0.0_wp .or. z_bed > z_sl, [sc%ncol]))
        na  = sc%chn%grd%n_active
        idx = sc%chn%grd%active_idx(1:na)

        if (sc%insol) then
            if (.not. present(orbit)) then
                write(*,*) "surface_chion_update:: Error: the insolation (itm, bessi) needs the orbital elements"
                error stop 1
            end if
            call insol_table_update(sc%ins, orbit, sc%par%dt_days)
        end if

        ! Mean-preserving daily control values, (nmon,na) so a column is contiguous.
        allocate(mon(na,nmon), ctl(na,nmon))
        do m = 1, nmon
            col = reshape(t2m(:,:,m), [sc%ncol])
            mon(:,m) = col(idx)
        end do
        ! The daily temperature variability of each column and month, from its monthly mean.
        sig = min(max(sc%par%tvar_sigma0 + sc%par%tvar_sigma_slope*(mon - T0), sc%par%tvar_sigma_min), &
                  sc%par%tvar_sigma_max)
        call monthly_to_daily_controls(sc%md, mon, ctl)
        t_ctl = transpose(ctl)
        do m = 1, nmon
            col = reshape(pr(:,:,m), [sc%ncol])
            mon(:,m) = col(idx)
        end do
        call monthly_to_daily_controls(sc%md, mon, ctl)
        p_ctl = transpose(ctl)
        deallocate(mon, ctl)

        ! Fixed for the year: geometry, the shortwave transmissivity and the air
        ! pressure from the surface elevation (the latter with the annual-mean air
        ! temperature), the ice thickness
        ! (H_ice = 0 is a land column: no ice albedo, substrate or ablation)
        ! and (ITM) the annual positive degree days.
        col   = reshape(sum(t2m, dim=3)/real(nmon, wp), [sc%ncol])
        t_ann = col(idx)
        col = reshape(z_srf, [sc%ncol])
        sc%chn%forc%surface_height(idx) = col(idx)
        ! shortwave_down = sw_fac*TOA: ITM applies its own transmissivity; BESSI the
        ! cloud proxy's clear sky at the surface and the cloud fraction.
        allocate(sw_fac(na))
        if (sc%bessi) then
            sw_fac = clear_sky_transmissivity(sc%chn%c, col(idx))*(1.0_wp - sc%par%cloud)
        else
            sw_fac = 1.0_wp
        end if
        sc%chn%forc%air_pressure(idx)   = p_sl*exp(-sc%chn%c%grav*col(idx)/(R_dry*t_ann))
        col = reshape(wind_speed_field(sc%par, z_srf, sc%dx), [sc%ncol])
        sc%chn%forc%wind_speed(idx) = col(idx)
        col = reshape(H_ice, [sc%ncol])
        sc%chn%forc%H_ice(idx) = col(idx)

        allocate(PDDs(na))
        PDDs = 0.0_wp
        if (sc%itm) then
            !$omp parallel do default(shared) private(i,day,t_d)
            do i = 1, na
                do day = 1, nday, pdd_dday
                    call interp_monthly_to_day(sc%md, t_ctl(:,i), day, t_d)
                    PDDs(i) = PDDs(i) + real(pdd_dday, wp) &
                            * real(pdd_expected_positive_temperature(real(t_d - T0, dp), &
                                                                     real(sc%par%sigma_pdd, dp)), wp)
                end do
                sc%chn%forc%PDDs(idx(i)) = PDDs(i)
            end do
            !$omp end parallel do
        end if

        ! Annual flux totals are the difference of chion's running totals.
        allocate(melt0(sc%ncol), runoff0(sc%ncol), refrz0(sc%ncol), subl0(sc%ncol))
        allocate(melt1(sc%ncol), runoff1(sc%ncol), refrz1(sc%ncol), subl1(sc%ncol))
        call chion_get_surface_flux_totals(sc%chn, melt=melt0, runoff=runoff0, refrz=refrz0, subl=subl0)

        allocate(t_sum(na), pr_sum(na), sf_sum(na), S_sum(na), ts_sum(na), alb_sum(na))
        t_sum = 0.0_wp;  pr_sum = 0.0_wp;  sf_sum = 0.0_wp;  S_sum = 0.0_wp
        ts_sum = 0.0_wp; alb_sum = 0.0_wp
        allocate(tsrf_c(sc%ncol), alb_c(sc%ncol))
        if (sc%seb%on) then
            allocate(fx(sc%ncol, NSEB), seb_sum(na, nmon, NSEB))
            seb_sum = 0.0_wp
        end if

        do day = 1, nday, sc%par%dt_days
            sc%chn%forc%day_of_year         = real(day, wp)
            sc%chn%forc%solar_longitude_deg = solar_longitude_deg(day, nday)

            m = (day - 1)/nday_mon + 1
            k = (day - 1)/sc%par%dt_days + 1
            !$omp parallel do default(shared) private(i,icol,t_d,p_d,sf_d,S_d)
            do i = 1, na
                icol = idx(i)
                call interp_monthly_to_day(sc%md, t_ctl(:,i), day, t_d)
                t_d = t_d + sig(i,m)*sc%xi(k)
                call interp_monthly_to_day(sc%md, p_ctl(:,i), day, p_d)
                ! The mean-preserving controls can undershoot zero for a peaked field.
                p_d  = max(p_d, 0.0_wp)
                sf_d = p_d*calc_snowfrac(t_d, sc%par%sf_a, sc%par%sf_b)
                S_d  = 0.0_wp
                if (sc%insol) then
                    S_d = insol_stencil_point(sc%ins%st, icol, sc%ins%S(:,day))
                end if

                sc%chn%forc%air_temperature(icol) = t_d
                sc%chn%forc%snowfall_rate(icol)   = sf_d/spd             ! [mm/d] -> [kg m-2 s-1]
                sc%chn%forc%rainfall_rate(icol)   = (p_d - sf_d)/spd
                sc%chn%forc%shortwave_down(icol)  = sw_fac(i)*S_d
                sc%chn%forc%toa_shortwave(icol)   = S_d

                t_sum(i)  = t_sum(i)  + t_d*dt
                pr_sum(i) = pr_sum(i) + p_d*dt
                sf_sum(i) = sf_sum(i) + sf_d*dt
                S_sum(i)  = S_sum(i)  + S_d*dt
            end do
            !$omp end parallel do

            call chion_update(sc%chn, dt)

            ! Surface state after the step, into the annual means.
            call chion_get_surface(sc%chn, t_srf=tsrf_c, albedo=alb_c)
            !$omp parallel do default(shared) private(i,icol)
            do i = 1, na
                icol = idx(i)
                ts_sum(i)  = ts_sum(i)  + tsrf_c(icol)*dt
                alb_sum(i) = alb_sum(i) + alb_c(icol)*dt
            end do
            !$omp end parallel do

            if (sc%seb%on) then
                call chion_get_surface_fluxes(sc%chn, melt=fx(:,10), sw_down=fx(:,1), sw_net=fx(:,2), &
                                              lw_down=fx(:,3), lw_net=fx(:,4), sh=fx(:,5), lh=fx(:,6), prh=fx(:,7))
                fx(:,8)  = alb_c
                fx(:,9)  = tsrf_c
                fx(:,10) = fx(:,10)*spd                                ! [kg m-2 s-1] -> [mm w.e./d]
                !$omp parallel do default(shared) private(i)
                do i = 1, na
                    seb_sum(i,m,:) = seb_sum(i,m,:) + fx(idx(i),:)*dt
                end do
                !$omp end parallel do
            end if
        end do

        if (sc%seb%on) then
            do k = 1, NSEB
                do m = 1, nmon
                    sc%seb%f(:,:,m,k) = unpack_col(seb_sum(:,m,k)/real(nday_mon, wp), idx, sc, MV)
                end do
            end do
        end if

        call chion_get_surface_flux_totals(sc%chn, melt=melt1, runoff=runoff1, refrz=refrz1, subl=subl1)

        ! Annual totals [mm w.e./yr] == [kg m-2/yr]; MV (unresolved) -> 0.
        melt_a   = annual_total(melt1(idx),   melt0(idx))
        runoff_a = annual_total(runoff1(idx), runoff0(idx))
        refrz_a  = annual_total(refrz1(idx),  refrz0(idx))
        subl_a   = annual_total(subl1(idx),   subl0(idx))

        sc%ann%mask = unpack_col(spread(1.0_wp, 1, na), idx, sc, 0.0_wp)
        sc%ann%t2m  = sum(t2m, dim=3)/real(nmon, wp)
        call scatter(sc%ann%t2m,    t_sum/real(nday, wp), idx, sc)
        sc%ann%pr   = sum(pr, dim=3)/real(nmon, wp)*real(nday, wp)
        call scatter(sc%ann%pr,     pr_sum, idx, sc)
        sc%ann%sf     = unpack_col(sf_sum,                       idx, sc, 0.0_wp)
        sc%ann%S      = unpack_col(S_sum/real(nday, wp),         idx, sc, 0.0_wp)
        sc%ann%PDDs   = unpack_col(PDDs,                         idx, sc, 0.0_wp)
        sc%ann%melt   = unpack_col(melt_a,                       idx, sc, 0.0_wp)
        sc%ann%runoff = unpack_col(runoff_a,                     idx, sc, 0.0_wp)
        sc%ann%refrz  = unpack_col(refrz_a,                      idx, sc, 0.0_wp)
        sc%ann%smb    = unpack_col(pr_sum - runoff_a - subl_a,   idx, sc, 0.0_wp)
        sc%ann%alb_s  = unpack_col(alb_sum/real(nday, wp),       idx, sc, MV)

        ! Surface temperature (inactive columns keep t2m): ITM and PDD, the
        ! annual-mean air temperature warmed by the gross annual refreezing
        ! (chion_tsrf_annual; H_ice is chion's forcing of the year; smbpal ITM
        ! used the net refreezing refrz - melt instead); BESSI, the annual mean
        ! of its per-step surface temperature.
        sc%ann%tsrf = sc%ann%t2m
        if (trim(sc%chn%par%model) == "bessi") then
            call scatter(sc%ann%tsrf, ts_sum/real(nday, wp), idx, sc)
        else
            allocate(t_air_c(sc%ncol), refrz_c(sc%ncol), tsrf_ann_c(sc%ncol))
            t_air_c = MV
            refrz_c = MV
            t_air_c(idx) = t_sum/real(nday, wp)
            refrz_c(idx) = refrz_a
            call chion_tsrf_annual(sc%chn, t_air_c, refrz_c, tsrf_ann_c)
            call scatter(sc%ann%tsrf, tsrf_ann_c(idx), idx, sc)
        end if

        sc%ann%H_snow = MV
        select case(trim(sc%chn%par%model))
            case("itm")
                call scatter(sc%ann%H_snow, sc%chn%itm%now%H_snow(idx), idx, sc)
            case("pdd")
                call scatter(sc%ann%H_snow, sc%chn%pdd%now%snowpack_swe(idx), idx, sc)
        end select
    end subroutine surface_chion_update

    subroutine surface_chion_restart_write(sc, filename, time)
        ! chion's prognostic state (chion_restart_write), plus the annual fields
        ! of the last annual cycle (sc%ann, on chion's `column` dimension as
        ! ann_<name>): they are the surface output Yelmo receives at the next
        ! step (the one-step coupling lag) and cannot be rebuilt from the state
        ! without stepping it another year.
        type(surface_chion_class), intent(in) :: sc
        character(len=*),          intent(in) :: filename
        real(dp),                  intent(in) :: time

        integer :: ncid

        call chion_restart_write(sc%chn, filename, time)

        call nc_open(filename, ncid, writable=.true.)
        call write_ann("mask",   sc%ann%mask,   "1",          "chion active (land or ice)")
        call write_ann("t2m",    sc%ann%t2m,    "K",          "Near-surface air temperature")
        call write_ann("pr",     sc%ann%pr,     "mm w.e./yr", "Precipitation")
        call write_ann("sf",     sc%ann%sf,     "mm w.e./yr", "Snowfall")
        call write_ann("S",      sc%ann%S,      "W m-2",      "TOA insolation")
        call write_ann("PDDs",   sc%ann%PDDs,   "K d",        "Positive degree days")
        call write_ann("melt",   sc%ann%melt,   "mm w.e./yr", "Melt")
        call write_ann("runoff", sc%ann%runoff, "mm w.e./yr", "Runoff")
        call write_ann("refrz",  sc%ann%refrz,  "mm w.e./yr", "Refreezing")
        call write_ann("smb",    sc%ann%smb,    "mm w.e./yr", "Surface mass balance (sf+rf-runoff-subl)")
        call write_ann("tsrf",   sc%ann%tsrf,   "K",          "Surface temperature")
        call write_ann("alb_s",  sc%ann%alb_s,  "1",          "Surface albedo")
        call write_ann("H_snow", sc%ann%H_snow, "mm w.e.",    "Snowpack thickness (end of year)")
        call nc_close(ncid)

    contains

        subroutine write_ann(name, var, units, long_name)
            character(len=*), intent(in) :: name, units, long_name
            real(wp),         intent(in) :: var(:,:)

            call nc_write(filename, "ann_"//name, reshape(var, [sc%ncol]), dim1="column", dim2="time", &
                          start=[1,1], units=units, long_name="Annual "//long_name, &
                          grid_mapping="", ncid=ncid)
        end subroutine write_ann

    end subroutine surface_chion_restart_write

    subroutine surface_chion_restart_read(sc, filename, ann_restored)
        ! chion's prognostic state (and active mask) and the annual fields of
        ! the last annual cycle from a restart file written by
        ! surface_chion_restart_write. ann_restored is false for a file without
        ! the annual fields (written before they were added): sc%ann is then
        ! unset and the host must rebuild it with an annual update, which
        ! steps the snowpack one extra year (the restart is not exact).
        type(surface_chion_class), intent(inout) :: sc
        character(len=*),          intent(in)    :: filename
        logical,                   intent(out)   :: ann_restored

        real(dp) :: time_file
        integer  :: ncid

        call chion_restart_read(sc%chn, filename, time_file)

        ann_restored = nc_exists_var(filename, "ann_smb")
        if (.not. ann_restored) then
            write(*,*) "surface_chion_restart_read:: warning: "//trim(filename)//" has no annual &
                       &fields (ann_*); they are rebuilt by an annual update, so the restart is &
                       &not exact."
            return
        end if

        call nc_open(filename, ncid, writable=.false.)
        call read_ann("mask",   sc%ann%mask)
        call read_ann("t2m",    sc%ann%t2m)
        call read_ann("pr",     sc%ann%pr)
        call read_ann("sf",     sc%ann%sf)
        call read_ann("S",      sc%ann%S)
        call read_ann("PDDs",   sc%ann%PDDs)
        call read_ann("melt",   sc%ann%melt)
        call read_ann("runoff", sc%ann%runoff)
        call read_ann("refrz",  sc%ann%refrz)
        call read_ann("smb",    sc%ann%smb)
        call read_ann("tsrf",   sc%ann%tsrf)
        call read_ann("alb_s",  sc%ann%alb_s)
        call read_ann("H_snow", sc%ann%H_snow)
        call nc_close(ncid)

    contains

        subroutine read_ann(name, var)
            character(len=*), intent(in)    :: name
            real(wp),         intent(inout) :: var(:,:)

            real(wp) :: col(sc%ncol)

            call nc_read(filename, "ann_"//name, col, start=[1,1], count=[sc%ncol,1], ncid=ncid)
            var = reshape(col, [sc%nx, sc%ny])
        end subroutine read_ann

    end subroutine surface_chion_restart_read

    subroutine surface_chion_par_load(par, filename, group)
        type(surface_chion_param_class), intent(out) :: par
        character(len=*),                intent(in)  :: filename
        character(len=*),                intent(in)  :: group

        call nml_read(filename, group, "sf_a",        par%sf_a)
        call nml_read(filename, group, "sf_b",        par%sf_b)
        call nml_read(filename, group, "sigma_pdd",   par%sigma_pdd)
        call nml_read(filename, group, "time_equil",  par%time_equil)
        call nml_read(filename, group, "dt_days",     par%dt_days)
        call nml_read(filename, group, "cloud",       par%cloud)
        call nml_read(filename, group, "wind_speed",  par%wind_speed)
        call nml_read(filename, group, "rel_hum",     par%rel_hum)
        if (nml_has_param(filename, group, "wind_lapse"))        call nml_read(filename, group, "wind_lapse", par%wind_lapse)
        if (nml_has_param(filename, group, "wind_slope"))        call nml_read(filename, group, "wind_slope", par%wind_slope)
        if (nml_has_param(filename, group, "wind_slope_length")) &
            call nml_read(filename, group, "wind_slope_length", par%wind_slope_length)
        if (nml_has_param(filename, group, "tvar_sigma0"))      call nml_read(filename, group, "tvar_sigma0", par%tvar_sigma0)
        if (nml_has_param(filename, group, "tvar_sigma_slope")) &
            call nml_read(filename, group, "tvar_sigma_slope", par%tvar_sigma_slope)
        if (nml_has_param(filename, group, "tvar_sigma_min"))   call nml_read(filename, group, "tvar_sigma_min", par%tvar_sigma_min)
        if (nml_has_param(filename, group, "tvar_sigma_max"))   call nml_read(filename, group, "tvar_sigma_max", par%tvar_sigma_max)
        if (nml_has_param(filename, group, "tvar_tau"))         call nml_read(filename, group, "tvar_tau", par%tvar_tau)
        if (nml_has_param(filename, group, "tvar_seed"))        call nml_read(filename, group, "tvar_seed", par%tvar_seed)
        if (nml_has_param(filename, group, "swe_init_ice")) call nml_read(filename, group, "swe_init_ice", par%swe_init_ice)
    end subroutine surface_chion_par_load

    subroutine ann_alloc(ann, nx, ny)
        type(surface_chion_ann_class), intent(inout) :: ann
        integer,                       intent(in)    :: nx, ny

        allocate(ann%mask(nx,ny), ann%t2m(nx,ny), ann%pr(nx,ny), ann%sf(nx,ny))
        allocate(ann%S(nx,ny), ann%PDDs(nx,ny), ann%melt(nx,ny), ann%runoff(nx,ny))
        allocate(ann%refrz(nx,ny), ann%smb(nx,ny), ann%tsrf(nx,ny), ann%alb_s(nx,ny))
        allocate(ann%H_snow(nx,ny))

        ann%mask   = 0.0_wp
        ann%t2m    = 0.0_wp
        ann%pr     = 0.0_wp
        ann%sf     = 0.0_wp
        ann%S      = 0.0_wp
        ann%PDDs   = 0.0_wp
        ann%melt   = 0.0_wp
        ann%runoff = 0.0_wp
        ann%refrz  = 0.0_wp
        ann%smb    = 0.0_wp
        ann%tsrf   = 0.0_wp
        ann%alb_s  = MV
        ann%H_snow = MV
    end subroutine ann_alloc

    function unpack_col(val, idx, sc, fill) result(var)
        ! Active-column values onto the (nx,ny) grid, fill elsewhere.
        real(wp),                  intent(in) :: val(:)
        integer,                   intent(in) :: idx(:)
        type(surface_chion_class), intent(in) :: sc
        real(wp),                  intent(in) :: fill
        real(wp) :: var(sc%nx,sc%ny)

        var = fill
        call scatter(var, val, idx, sc)
    end function unpack_col

    function red_noise_monthly(nday, dt_days, tau, seed) result(xi)
        ! One year of red noise, a value per step of dt_days: AR(1) with the lag-one correlation
        ! exp(-dt_days/tau) driven by standard normals (Park-Miller uniforms, Box-Muller), then
        ! its mean removed in every 30-day month and the whole scaled to unit variance.
        integer,  intent(in) :: nday, dt_days, seed
        real(wp), intent(in) :: tau        ! [d]
        real(wp) :: xi(nday/dt_days)

        integer, parameter :: i8 = selected_int_kind(18)
        integer(i8) :: x
        integer  :: n, k, m, nm
        real(wp) :: r, u1, u2, g

        n = nday/dt_days
        x = int(mod(abs(seed), 2147483646), i8) + 1_i8
        r = exp(-real(dt_days, wp)/max(tau, 1.0e-6_wp))
        do k = 1, n
            x  = mod(16807_i8*x, 2147483647_i8); u1 = real(x, wp)/2147483647.0_wp
            x  = mod(16807_i8*x, 2147483647_i8); u2 = real(x, wp)/2147483647.0_wp
            g  = sqrt(-2.0_wp*log(u1))*cos(8.0_wp*atan(1.0_wp)*u2)
            if (k == 1) then
                xi(k) = g
            else
                xi(k) = r*xi(k-1) + sqrt(1.0_wp - r*r)*g
            end if
        end do
        nm = max(1, nday_mon/dt_days)
        do m = 1, n, nm
            xi(m:min(m+nm-1, n)) = xi(m:min(m+nm-1, n)) - sum(xi(m:min(m+nm-1, n)))/real(min(nm, n-m+1), wp)
        end do
        xi = xi/sqrt(max(sum(xi*xi)/real(n, wp), tiny(1.0_wp)))
    end function red_noise_monthly

    function wind_speed_field(par, z_srf, dx) result(u)
        ! U = wind_speed + wind_lapse*max(z,0)[km] + wind_slope*|grad z|[m/km], at least 0; the
        ! gradient a centred difference over 2n cells, n = max(1, nint(L/(2 dx))), one-sided
        ! over the cells available at the edges.
        type(surface_chion_param_class), intent(in) :: par
        real(wp), intent(in) :: z_srf(:,:)          ! [m]
        real(wp), intent(in) :: dx                  ! [m]
        real(wp) :: u(size(z_srf,1), size(z_srf,2))

        integer  :: i, j, n, nx, ny, i0, i1, j0, j1
        real(wp) :: gx, gy

        u = par%wind_speed + par%wind_lapse*max(z_srf, 0.0_wp)*1.0e-3_wp
        if (par%wind_slope /= 0.0_wp) then
            nx = size(z_srf,1); ny = size(z_srf,2)
            n  = max(1, nint(par%wind_slope_length*1.0e3_wp/(2.0_wp*dx)))
            do j = 1, ny
                j0 = max(j-n, 1); j1 = min(j+n, ny)
                do i = 1, nx
                    i0 = max(i-n, 1); i1 = min(i+n, nx)
                    gx = 0.0_wp; gy = 0.0_wp
                    if (i1 > i0) gx = (z_srf(i1,j) - z_srf(i0,j))/(real(i1-i0, wp)*dx)
                    if (j1 > j0) gy = (z_srf(i,j1) - z_srf(i,j0))/(real(j1-j0, wp)*dx)
                    u(i,j) = u(i,j) + par%wind_slope*1.0e3_wp*sqrt(gx*gx + gy*gy)
                end do
            end do
        end if
        u = max(u, 0.0_wp)
    end function wind_speed_field

    subroutine scatter(var, val, idx, sc)
        ! Overwrite the active columns of an (nx,ny) field.
        real(wp),                  intent(inout) :: var(:,:)
        real(wp),                  intent(in)    :: val(:)
        integer,                   intent(in)    :: idx(:)
        type(surface_chion_class), intent(in)    :: sc

        real(wp), allocatable :: col(:)

        col = reshape(var, [sc%ncol])
        col(idx) = val
        var = reshape(col, [sc%nx, sc%ny])
    end subroutine scatter

    elemental function annual_total(tot1, tot0) result(y)
        ! Difference of two running totals; chion reports MV for a flux the
        ! model does not resolve, which counts as zero here.
        real(wp_acc), intent(in) :: tot1, tot0
        real(wp) :: y

        y = 0.0_wp
        if (tot1 /= real(MV, wp_acc) .and. tot0 /= real(MV, wp_acc)) y = real(tot1 - tot0, wp)
    end function annual_total

    subroutine insol_table_init(ins, lats, nday)
        ! The regular latitude table (insol_dlat)
        ! over [-90,90] with its insol state, and each column's
        ! linear-interpolation stencil into it.
        type(surface_chion_insol_class), intent(inout) :: ins
        real(wp),                        intent(in)    :: lats(:)     ! (ncol) [deg N]
        integer,                         intent(in)    :: nday

        integer :: nlat, j

        nlat = nint(180.0_dp/insol_dlat) + 1
        allocate(ins%lat(nlat), ins%S(nlat,nday))
        do j = 1, nlat
            ins%lat(j) = -90.0_dp + real(j-1, dp)*insol_dlat
        end do
        ins%S   = 0.0_wp
        ins%set = .false.

        call insol_init(ins%is, ins%lat, day_year=nday, tables=.false.)

        call insol_stencil_init(ins%st, ins%lat, real(lats, dp))
    end subroutine insol_table_init

    subroutine insol_table_update(ins, orbit, dt_days)
        ! Daily TOA insolation on the latitude table for the stepped days, only
        ! when the orbital elements have changed since the last call.
        type(surface_chion_insol_class), intent(inout) :: ins
        type(orbit_pars_class),          intent(in)    :: orbit
        integer,                         intent(in)    :: dt_days

        integer, allocatable :: days(:)
        real(insol_wp), allocatable :: S(:,:)
        integer :: q

        if (ins%set) then
            if (orbit%ecc == ins%orbit%ecc .and. orbit%obl == ins%orbit%obl .and. orbit%per == ins%orbit%per) return
        end if

        days = [(q, q = 1, size(ins%S,2), dt_days)]
        allocate(S(size(ins%lat), size(days)))
        call calc_insol_days(ins%is, S, days, orbit)
        do q = 1, size(days)
            ins%S(:,days(q)) = real(S(:,q), wp)
        end do
        ins%orbit = orbit
        ins%set   = .true.
    end subroutine insol_table_update


    pure function solar_longitude_deg(day, nday) result(lon)
        ! [deg] Solar longitude of a day of the model year (nday days), with
        ! Chion.jl's calendar-day formula (_solar_longitude_deg_from_calendar_day)
        ! on the day stretched to the tropical year.
        integer, intent(in) :: day, nday
        real(wp) :: lon

        real(dp) :: d, mean_lon, mean_anom
        real(dp), parameter :: deg = acos(-1.0_dp)/180.0_dp

        d         = real(day - 1, dp)*real(days_trop, dp)/real(nday, dp)
        mean_lon  = 280.46646_dp + 0.98564736_dp*d
        mean_anom = 357.52911_dp + 0.98560028_dp*d
        lon = real(modulo(mean_lon + 1.914602_dp*sin(mean_anom*deg) &
                          + 0.019993_dp*sin(2.0_dp*mean_anom*deg), 360.0_dp), wp)
    end function solar_longitude_deg

    elemental function calc_snowfrac(t2m, a, b) result(f)
        ! smbpal calc_snowfrac: snow fraction of total precipitation.
        real(wp), intent(in) :: t2m, a, b
        real(wp) :: f

        f = -0.5_wp*tanh(a*(t2m - b)) + 0.5_wp
    end function calc_snowfrac

end module kryos_chion_host

module kryos_surface_chion
    ! surface backend chion (adapters-design section 5), built with chion=1
    ! (and insol=1: the host's insolation table); otherwise
    ! kryos_surface_chion_stub.f90 provides this module.
    !
    ! The adapter absorbs the yelmox chion host: the module above is
    ! libs/surface_chion.f90 of yelmox 74f3a2b (branch topo-prep), unchanged
    ! except its module name (kryos_chion_host); compiling it here, in the file
    ! of the backend, takes it out of the build when chion=0. The orbit of its
    ! insolation is the run's: the forcing's elements orbit_ecc, orbit_obl,
    ! orbit_per (imports, with itm or bessi), resolved there by kryos_orbit.
    !
    ! chion computes in wp (sp), accumulators in dp; the exchange is dp: imports
    ! are narrowed once, exports widened from chion's values (P1, P5). One
    ! annual cycle per call, whatever t1 - t0 (adapters-design 5.2).

    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use fesm_units,      only: units_conv
    use ncio,            only: nc_close
    use nml,             only: nml_read, nml_has_param
    use kryos_outfile,   only: outfile_dims_2d, outfile_open, outfile_var2d, outfile_axis, outfile_var3d
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, component_write_output, &
                               sec_year_check, NAME_LEN, PASS_BOUNDARY, PASS_FINAL
    use kryos_surface,   only: surface_t, surface_opt_t
    use chion,           only: cwp => wp, chion_end
    use insolation,      only: orbit_pars_class, orbit_pars_fixed
    use kryos_chion_host, only: surface_chion_class, surface_chion_init, surface_chion_spinup, &
                                surface_chion_update, surface_chion_restart_write, &
                                surface_chion_restart_read, surface_chion_seb_enable, &
                                NSEB, SEB_NAME, SEB_UNITS, SEB_LONG

    implicit none

    private
    public :: chion_surface_t, chion_surface_constraints

    integer, parameter :: NMON = 12

    ! Precisions of the 2D stream (input/kryos-variables-<TABLE>.md)
    character(len=*), parameter :: TABLE = "surface-chion"

    type, extends(surface_t) :: chion_surface_t
        type(surface_chion_class) :: sc               ! chion + its yelmox host
        character(len=NAME_LEN) :: host_group  = "surface_chion"   ! host parameters ([surface_chion])
        character(len=NAME_LEN) :: chion_group = "chion"           ! chion's group
        logical :: seb = .false.                      ! stream seb on (out_seb): the monthly energy balance
        real(dp) :: pr_scale = 1.0_dp                 ! pr [kg m-2 yr-1] = pr [mm d-1] * pr_scale
        real(cwp), allocatable :: t2m(:,:,:), pr_mmd(:,:,:)        ! narrowed imports (nx,ny,12)
        real(dp),  allocatable :: pr_in(:,:,:)
        real(cwp), allocatable :: z_srf_s(:,:), z_bed_s(:,:), H_ice_s(:,:), z_sl_s(:,:)
    contains
        procedure :: init          => chion_init_
        procedure :: init_state    => chion_init_state_
        procedure :: update        => chion_update_
        procedure :: restart_write => chion_restart_write_
        procedure :: restart_read  => chion_restart_read_
        procedure :: write_output  => chion_write_output_
        procedure :: end           => chion_end_
    end type chion_surface_t

contains

    function chion_surface_constraints() result(con)
        ! Regional, any grid, several instances.
        type(constraints_t) :: con
        con = constraints_t(build_flag="chion=1")
    end function chion_surface_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine chion_init_(self, ctx)
        ! surface_chion_init (kryos.f90:410-412) on the grid's latitudes.
        class(chion_surface_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        real(dp) :: offset, dx

        call group_read(self, ctx)
        self%prec = cwp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny
        self%opt  = surface_opt_t(melt=.true., runoff=.true., refrz=.true., albedo=.false.)

        ! The grid spacing [m] for the slope of the wind; none on a lon-lat grid.
        dx = 0.0_dp
        if (trim(ctx%grid%cs%units) /= "degrees") dx = ctx%grid%G%dx*ctx%grid%cs%xy_conv
        call surface_chion_init(self%sc, trim(ctx%par_file), real(ctx%grid%lat, cwp), real(dx, cwp), &
                                group=trim(self%host_group), chion_group=trim(self%chion_group), &
                                cnst=ctx%cnst)
        self%opt%orbit = self%sc%insol                       ! the forcing's orbit for the insolation
        if (self%seb) call surface_chion_seb_enable(self%sc)

        select case (trim(self%sc%chn%par%restart))
        case ("none", "None", "")
        case default
            call fesm_error("kryos_surface_chion::init", "["//trim(self%chion_group)//"] restart must be 'None'", &
                            "instance "//trim(self%name)//": restarts come from the bundle, got '"// &
                            trim(self%sc%chn%par%restart)//"'")
        end select

        ! pr: the climate exports kg m-2 yr-1 = mm d-1 * scale (scale = sec_year/86400, the
        ! year of the run calendar); chion takes mm d-1 (P3).
        call sec_year_check(ctx, self%name)
        call units_conv("mm d-1", "kg m-2 yr-1", self%pr_scale, offset, sec_year=ctx%sec_year)

        allocate(self%t2m(self%nx, self%ny, NMON), self%pr_mmd(self%nx, self%ny, NMON), &
                 self%pr_in(self%nx, self%ny, NMON))
        allocate(self%z_srf_s(self%nx, self%ny), self%z_bed_s(self%nx, self%ny), &
                 self%H_ice_s(self%nx, self%ny), self%z_sl_s(self%nx, self%ny))

        call log_msg(LOG_INFO, "kryos_surface_chion", "chion model "//trim(self%sc%chn%par%model)// &
                     ": one 360-day cycle per call; its annual totals are used as per-year rates", &
                     "instance "//trim(self%name)//"; pr: kg m-2 yr-1 / "//to_str(self%pr_scale)//" = mm d-1")
    end subroutine chion_init_

    subroutine chion_init_state_(self, st, pass)
        ! PASS_BOUNDARY: spin-up and one cycle (kryos_coupling.f90:492-495,
        ! kryos_startup.f90:162-169); PASS_FINAL: one cycle on the final initial
        ! geometry (:192). Both step the snowpack at t_init on purpose. Other
        ! passes: nothing (the annual fields come with the restart file).
        class(chion_surface_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st
        integer,                intent(in)    :: pass

        select case (pass)
        case (PASS_BOUNDARY)
            call imports_narrow(self)
            call surface_chion_spinup(self%sc, self%t2m, self%pr_mmd, self%z_srf_s, self%z_bed_s, &
                                      self%H_ice_s, self%z_sl_s, orbit_now(self))
            call surface_chion_update(self%sc, self%t2m, self%pr_mmd, self%z_srf_s, self%z_bed_s, &
                                      self%H_ice_s, self%z_sl_s, orbit_now(self))
            call exports_put(self, st)
        case (PASS_FINAL)
            call imports_narrow(self)
            call surface_chion_update(self%sc, self%t2m, self%pr_mmd, self%z_srf_s, self%z_bed_s, &
                                      self%H_ice_s, self%z_sl_s, orbit_now(self))
            call exports_put(self, st)
        end select
    end subroutine chion_init_state_

    subroutine chion_update_(self, st)
        ! step_surface (kryos_coupling.f90:492-495): one annual cycle at the
        ! forcing time of t1, exported as the rate of [t0, t1].
        class(chion_surface_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st

        call imports_narrow(self)
        call surface_chion_update(self%sc, self%t2m, self%pr_mmd, self%z_srf_s, self%z_bed_s, &
                                  self%H_ice_s, self%z_sl_s, orbit_now(self))
        call exports_put(self, st)
    end subroutine chion_update_

    function orbit_now(self) result(orbit)
        ! The orbital elements the forcing exported (rows forcing -> surface, lag 0); the elements of 1950
        ! without insolation (then unused).
        class(chion_surface_t), intent(in) :: self
        type(orbit_pars_class) :: orbit
        real(dp) :: e(3)
        integer  :: k
        orbit = orbit_pars_fixed(0.016725_dp, 23.4458_dp, 102.063_dp)
        if (.not. self%opt%orbit) return
        do k = 1, 3
            if (.not. self%imports%f(self%hi_orbit(k))%published) then
                call fesm_error("kryos_surface_chion::update", "import "//trim(self%imports%f(self%hi_orbit(k))%name)// &
                                " has not been published", "instance "//trim(self%name)//": the insolation (itm, bessi) "// &
                                "takes the run's orbit from the forcing ([forcing] orbit = 'orbit')")
            end if
            call fields_get(self%imports, self%hi_orbit(k), e(k))
        end do
        orbit = orbit_pars_fixed(e(1), e(2), e(3))
    end function orbit_now

    subroutine chion_restart_write_(self, dir, t)
        ! kryos_startup.f90:431-432: chion state plus the annual fields.
        class(chion_surface_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t
        call surface_chion_restart_write(self%sc, restart_file(self, dir), t)
    end subroutine chion_restart_write_

    subroutine chion_restart_read_(self, dir, t)
        ! kryos_startup.f90:518-519; a file without the annual fields stops
        ! (adapters-design 1.4) instead of re-deriving them by an extra year.
        class(chion_surface_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t

        logical :: ann_restored
        associate (unused => t); end associate

        call surface_chion_restart_read(self%sc, restart_file(self, dir), ann_restored)
        if (.not. ann_restored) then
            call fesm_error("kryos_surface_chion::restart_read", "restart file has no annual fields (ann_*)", &
                            "instance "//trim(self%name)//": file "//restart_file(self, dir))
        end if
    end subroutine chion_restart_read_

    subroutine chion_write_output_(self, file, t, stream, init)
        ! Stream 2D: the 13 annual fields of surface_chion_write_step
        ! (kryos_output.f90:326-361), same units and precision. Stream seb:
        ! BESSI's monthly surface energy balance, albedo, surface temperature
        ! and melt of the last annual cycle (xc, yc, month). Others: default.
        class(chion_surface_t), intent(inout) :: self
        character(len=*),       intent(in)    :: file, stream
        real(dp),               intent(in)    :: t
        logical,                intent(in)    :: init

        integer :: ncid, n, k, m

        select case (trim(stream))
        case ("2D")
        case ("seb")
            if (init) then
                call outfile_dims_2d(file, self%xkm, self%ykm, t)
                call outfile_axis(file, "month", [(real(m, dp), m = 1, NMON)], "month")
            end if
            call outfile_open(file, t, ncid, n)
            do k = 1, NSEB
                call outfile_var3d(file, trim(SEB_NAME(k)), self%sc%seb%f(:,:,:,k), "month", n, ncid, &
                                   trim(SEB_UNITS(k)), "Monthly mean "//trim(SEB_LONG(k)))
            end do
            call nc_close(ncid)
            return
        case default
            call component_write_output(self, file, t, stream, init)
            return
        end select

        if (init) call outfile_dims_2d(file, self%xkm, self%ykm, t)
        call outfile_open(file, t, ncid, n)
        associate (ann => self%sc%ann)
            call outfile_var2d(file, "smb",    ann%smb,    n, ncid, "mm w.e./yr", "Surface mass balance (sf+rf-runoff-subl)", table=TABLE)
            call outfile_var2d(file, "tsrf",   ann%tsrf,   n, ncid, "K",          "Surface temperature", table=TABLE)
            call outfile_var2d(file, "melt",   ann%melt,   n, ncid, "mm w.e./yr", "Melt", table=TABLE)
            call outfile_var2d(file, "runoff", ann%runoff, n, ncid, "mm w.e./yr", "Runoff", table=TABLE)
            call outfile_var2d(file, "refrz",  ann%refrz,  n, ncid, "mm w.e./yr", "Refreezing", table=TABLE)
            call outfile_var2d(file, "pr",     ann%pr,     n, ncid, "mm w.e./yr", "Precipitation", table=TABLE)
            call outfile_var2d(file, "sf",     ann%sf,     n, ncid, "mm w.e./yr", "Snowfall", table=TABLE)
            call outfile_var2d(file, "t2m",    ann%t2m,    n, ncid, "K",          "Near-surface air temperature", table=TABLE)
            call outfile_var2d(file, "S",      ann%S,      n, ncid, "W m-2",      "TOA insolation", table=TABLE)
            call outfile_var2d(file, "PDDs",   ann%PDDs,   n, ncid, "K d",        "Positive degree days", table=TABLE)
            call outfile_var2d(file, "alb_s",  ann%alb_s,  n, ncid, "1",          "Surface albedo", table=TABLE)
            call outfile_var2d(file, "H_snow", ann%H_snow, n, ncid, "mm w.e.",    "Snowpack thickness (end of year)", table=TABLE)
            call outfile_var2d(file, "mask",   ann%mask,   n, ncid, "1",          "chion active (land or ice)", table=TABLE)
        end associate
        call nc_close(ncid)
    end subroutine chion_write_output_

    subroutine chion_end_(self)
        class(chion_surface_t), intent(inout) :: self
        call chion_end(self%sc%chn)
    end subroutine chion_end_

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    subroutine group_read(self, ctx)
        ! Adapter keys of the instance group; the core has validated them.
        class(chion_surface_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        character(len=NAME_LEN) :: sched

        if (nml_has_param(ctx%par_file, ctx%instance, "surface_chion")) &
            call nml_read(ctx%par_file, ctx%instance, "surface_chion", self%host_group)
        if (nml_has_param(ctx%par_file, ctx%instance, "chion")) &
            call nml_read(ctx%par_file, ctx%instance, "chion", self%chion_group)
        ! The stream seb (its schedule out_seb) needs the monthly accumulation, set up at init.
        if (nml_has_param(ctx%par_file, ctx%instance, "out_seb")) then
            call nml_read(ctx%par_file, ctx%instance, "out_seb", sched)
            self%seb = .not. any(trim(sched) == ["none", "None", "    "])
        end if
    end subroutine group_read

    function restart_file(self, dir) result(file)
        class(chion_surface_t), intent(in) :: self
        character(len=*),       intent(in) :: dir
        character(len=:), allocatable      :: file
        file = trim(dir)//"/"//trim(self%name)//"_restart.nc"
    end function restart_file

    subroutine imports_narrow(self)
        ! Every import once: tas, geometry as widened sp values narrowed back
        ! (P1); pr divided by the factor the climate multiplied with, then
        ! narrowed (P3). The re-elevation (A6) is the core's wiring transform.
        class(chion_surface_t), intent(inout) :: self

        call fields_get(self%imports, self%hi_tas, self%t2m)
        call fields_get(self%imports, self%hi_pr, self%pr_in)
        self%pr_mmd = real(self%pr_in / self%pr_scale, cwp)
        call fields_get(self%imports, self%hi_z_srf, self%z_srf_s)
        call fields_get(self%imports, self%hi_z_bed, self%z_bed_s)
        call fields_get(self%imports, self%hi_H_ice, self%H_ice_s)
        call fields_get(self%imports, self%hi_z_sl,  self%z_sl_s)
    end subroutine imports_narrow

    subroutine exports_put(self, st)
        ! MEAN over [t0, t1], one put each: the annual fields of the last cycle
        ! in kg m-2 yr-1 (= mm w.e./yr numerically, adapters-design 5.3).
        class(chion_surface_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st

        associate (ann => self%sc%ann)
            call fields_put(self%exports, self%h_smb,    ann%smb,    t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_T_srf,  ann%tsrf,   t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_melt,   ann%melt,   t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_runoff, ann%runoff, t1=st%t1, t0=st%t0)
            call fields_put(self%exports, self%h_refrz,  ann%refrz,  t1=st%t1, t0=st%t0)
        end associate
    end subroutine exports_put

end module kryos_surface_chion
