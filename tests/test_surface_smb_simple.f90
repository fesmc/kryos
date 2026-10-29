! ---------------------------------------------------------------------------
! The yelmox module below is a verbatim copy of libs/smb_simple.f90 of yelmox 5ac89c8
! (branch topo-prep): the reference of the test. Do not edit it.
! ---------------------------------------------------------------------------
!=====================================================================
! Module: smb_simple_m
!
! Fortran sibling of `smb_simple.jl`. Two surface mass balance schemes
! sharing the same climate snow line:
!
!   * calc_smb_simple_syn
!         Piecewise-linear SMB(z) profile evaluated on a synthetic
!         surface elevation derived from a target ice mask. Production
!         entry point; one call returns SMB (mm w.e./yr) and a
!         surface-temperature field t_srf (K, sea-level temperature +
!         elevation lapse rate).
!
!         Profile per cell (accumulation minus ablation):
!             SMB = acc - abl
!             acc = c0 * f_lat(phi) * f_dist(d_in)    [maritime; decays inland]
!             abl = beta_abl(phi) * max(0, z_SL - z)  [ablation below the snow line]
!         evaluated on z = z_syn(target mask, real elevation). acc applies
!         everywhere and is continuous at the snow line (abl -> 0 there).
!
!   * calc_smb_simple_tg24
!         Original Talento-Ganopolski (2023) "simple" scheme:
!             dz > 0:  smb_cm =  (b1 + b4*|grad z|)*exp(-dz/b5)
!             dz <= 0: smb_cm = -b2*dz**2
!         optionally + b3*(mask - 1) when use_mask = .true.
!         smb_cm is in cm/yr; the kernel converts to mm w.e./yr (× 10).
!
! Shared snow line (phi = |latitude|, so the latitude terms are hemisphere-agnostic):
!     z_SL = a1 + a2*phi + a3*CO2 + a4*(f - fmean) [ + dz_SL ]
!
! Public interface:
!   - smb_params_syn      : parameter container for the synthetic-
!                           elevation scheme (snow line + ablation +
!                           accumulation + profile + clip).
!   - smb_params_tg24     : parameter container for the TG24 scheme.
!   - calc_smb_simple_syn : production entry point for the synthetic
!                           scheme. Computes signed distance, z_syn,
!                           SMB, and the outside-mask floor in one call.
!   - calc_smb_simple_pwl : low-level piecewise-linear SMB(z) kernel,
!                           exposed for diagnostics / testing.
!   - calc_smb_simple_tg24: 2-D field entry point for the TG24 scheme.
!   - compute_signed_distance : signed distance (m) to target-mask
!                               boundary, from 2-D x/y coordinate pairs
!                               (Cartesian or lon/lat per `units`).
!   - compute_z_syn_linear : linear wedge profile inside the mask.
!   - compute_z_syn_plastic: Nye/Vialov plastic profile inside the mask.
!   - apply_smb_min_outside: post-process SMB floor outside the mask.
!   - compute_t_srf_lapse  : surface temperature from t_sl and an
!                            elevation lapse-rate correction.
!=====================================================================
module smb_simple_m

    use iso_fortran_env, only: sp => real32
    use nml
    use ncio
    use phys_constants, only : phys_const_class, phys_const_require, phys_const_get
    implicit none
    private

    public :: smb_params_syn
    public :: smb_params_tg24
    public :: smb_simple_class
    public :: smb_simple_init
    public :: smb_simple_par_load
    public :: smb_simple_set_mask
    public :: smb_simple_update
    public :: calc_smb_simple_syn
    public :: calc_smb_simple_pwl
    public :: calc_smb_simple_tg24
    public :: compute_signed_distance
    public :: compute_z_syn_linear
    public :: compute_z_syn_plastic
    public :: apply_smb_min_outside
    public :: compute_t_srf_lapse

    !-----------------------------------------------------------------
    ! Parameter container for the synthetic-elevation piecewise-linear
    ! scheme. Defaults match the tuned set in smb_simple.jl: LGM
    ! Tarasov-fit τ₀ = 36 kPa, snow line anchored at
    !   z_SL(45 deg N) ~ 1500 m
    !   z_SL(75 deg N) ~    0 m
    !-----------------------------------------------------------------
    type :: smb_params_syn

        ! Climate snow line z_SL = a1 + a2*phi + a3*CO2 + a4*(f - fmean)
        real(sp) :: a1     = 3750.0_sp      ! m (constant)
        real(sp) :: a2     =  -50.0_sp      ! m / deg N (negative: drops poleward)
        real(sp) :: a3     =    0.0_sp      ! m / ppm           (0 in lat-only form)
        real(sp) :: a4     =    0.0_sp      ! m / (W/m^2)       (0 in lat-only form)
        real(sp) :: fmean  =  494.5_sp      ! W/m^2

        ! Ablation lapse rate beta_abl(phi)
        real(sp) :: beta0      = 4.0_sp     ! (mm/yr)/m at phiref
        real(sp) :: beta1      = 0.05_sp    ! (mm/yr)/m per deg
        real(sp) :: phiref     = 60.0_sp    ! deg N
        real(sp) :: beta_floor = 0.0_sp     ! (mm/yr)/m

        ! Accumulation: acc = c0 * f_lat(phi) * f_dist(d_in)
        !   f_dist = exp(-d_in_km / L_acc)   maritime decay inland (1 at margin -> 0 deep inland)
        !   f_lat  = clamp(1 - k_acc_lat*max(0, phi - phi_acc_ref), facc_lat_min, 1)
        ! d_in_km is distance inland from the target-mask margin (km).
        real(sp) :: c0    = 150.0_sp   ! mm w.e./yr, maximum (coastal) accumulation at the margin
        real(sp) :: L_acc = 1000.0_sp  ! km, e-folding length of inland accumulation decay
        real(sp) :: k_acc_lat    = 0.0_sp   ! frac. accum decrease per deg N beyond phi_acc_ref (0 => off)
        real(sp) :: phi_acc_ref  = 0.0_sp   ! deg N, latitude above which accum starts dropping
        real(sp) :: facc_lat_min = 0.3_sp   ! floor on latitude modifier

        ! Synthetic-elevation profile
        logical  :: use_plastic = .true.      ! .true. → plastic Nye/Vialov; .false. → linear
        real(sp) :: slope       = 6.09e-3_sp  ! m/m, linear inside (when use_plastic=.false.)
        real(sp) :: tau0        = 36.0e3_sp   ! Pa, plastic basal yield stress
        real(sp) :: slope_out   = 6.09e-3_sp  ! m/m, linear outside ramp (both kinds)
        real(sp) :: z_max_in    = 2500.0_sp   ! m, inside cap
        real(sp) :: z_max_out   =  200.0_sp   ! m, outside floor

        ! Plastic profile physical constants. Taken from the phys_const_class
        ! passed to smb_simple_init when there is one; these Earth values are
        ! the standalone fallback, overridable from the namelist.
        real(sp) :: rho_ice = 910.0_sp        ! kg/m^3
        real(sp) :: g       =   9.81_sp       ! m/s^2

        ! Surface temperature  t_srf = t_sl - gamma_t*max(z_syn, 0)
        real(sp) :: gamma_t   = 6.5e-3_sp     ! K/m (lapse rate)
        real(sp) :: t_ice_max = 273.15_sp     ! K, clamp over ice-covered cells

        ! Outside-mask SMB floor (must be <= 0)
        real(sp) :: smb_min = -2000.0_sp      ! mm w.e./yr

    end type smb_params_syn

    !-----------------------------------------------------------------
    ! Parameter container for the TG24 scheme.
    !-----------------------------------------------------------------
    type :: smb_params_tg24

        real(sp) :: a1     = 6500.0_sp
        real(sp) :: a2     = -200.0_sp
        real(sp) :: a3     =   19.0_sp
        real(sp) :: a4     =   20.0_sp
        real(sp) :: fmean  =  494.5_sp

        real(sp) :: b1 = 10.08_sp    ! cm/yr
        real(sp) :: b2 = 1.2e-5_sp   ! cm/yr/m^2
        real(sp) :: b3 = 227.6_sp    ! cm/yr (use_mask only)
        real(sp) :: b4 = 3.3e-7_sp   ! cm/yr
        real(sp) :: b5 = 8.0e5_sp    ! m

    end type smb_params_tg24

    !-----------------------------------------------------------------
    ! Container holding the smb_simple settings, the (static) grid, the
    ! target ice mask, and the output fields. Lets a driver program keep
    ! a single object instead of many loose variables.
    !-----------------------------------------------------------------
    type :: smb_simple_class

        ! Settings / parameters
        type(smb_params_syn) :: par         ! synthetic-scheme parameters
        character(len=56)    :: scheme      ! "syn" (tg24 not yet supported)
        real(sp)             :: co2         ! [ppm]  constant CO2
        real(sp)             :: f           ! [W/m2] constant insolation
        character(len=512)   :: mask_file = "" ! target-mask file ("" => H_ice_ref)
        character(len=56)    :: mask_var  = "" ! target-mask variable name
        integer              :: mask_idx = 1 ! index along the 3rd (e.g. time) dim of mask_var
        character(len=16)    :: units       ! coordinate units for distance ("m")

        ! Static grid (stored once so update needs no grid arguments)
        real(sp), allocatable :: x(:,:)      ! x coordinate (per units)
        real(sp), allocatable :: y(:,:)      ! y coordinate (per units)
        real(sp), allocatable :: lat(:,:)    ! latitude [deg N]

        ! State
        logical,  allocatable :: mask(:,:)   ! target ice mask
        real(sp), allocatable :: smb(:,:)    ! [mm w.e./yr]
        real(sp), allocatable :: t_srf(:,:)  ! [K]

    end type smb_simple_class

contains

    !=================================================================
    ! Class-based driver interface (init / set_mask / update)
    !=================================================================

    !-----------------------------------------------------------------
    !> Initialize an smb_simple_class: load the namelist settings, store
    !> the (static) coordinate grid, and allocate the state fields.
    !> Build the target mask separately with smb_simple_set_mask.
    !-----------------------------------------------------------------
    subroutine smb_simple_init(smbs, filename, x, y, lat, group, units, init, cnst)

        type(smb_simple_class), intent(inout) :: smbs
        character(len=*),       intent(in)    :: filename
        real(sp),               intent(in)    :: x(:,:)
        real(sp),               intent(in)    :: y(:,:)
        real(sp),               intent(in)    :: lat(:,:)
        character(len=*),       intent(in), optional :: group
        character(len=*),       intent(in), optional :: units
        logical,                intent(in), optional :: init
        type(phys_const_class), intent(in), optional :: cnst

        integer :: nx, ny

        nx = size(x, 1)
        ny = size(x, 2)

        if (size(y, 1)   /= nx .or. size(y, 2)   /= ny .or. &
            size(lat, 1) /= nx .or. size(lat, 2) /= ny) then
            error stop "smb_simple_init: x, y, lat shape mismatch"
        end if

        ! Load settings + synthetic-scheme parameters
        call smb_simple_par_load(smbs, filename, group=group, init=init, cnst=cnst)

        ! Coordinate units for the signed-distance calculation
        smbs%units = "m"
        if (present(units)) smbs%units = units

        ! Store the static grid
        if (allocated(smbs%x))   deallocate(smbs%x)
        if (allocated(smbs%y))   deallocate(smbs%y)
        if (allocated(smbs%lat)) deallocate(smbs%lat)
        allocate(smbs%x(nx, ny), smbs%y(nx, ny), smbs%lat(nx, ny))
        smbs%x   = x
        smbs%y   = y
        smbs%lat = lat

        ! Allocate state fields
        if (allocated(smbs%mask))  deallocate(smbs%mask)
        if (allocated(smbs%smb))   deallocate(smbs%smb)
        if (allocated(smbs%t_srf)) deallocate(smbs%t_srf)
        allocate(smbs%mask(nx, ny), smbs%smb(nx, ny), smbs%t_srf(nx, ny))
        smbs%mask  = .false.
        smbs%smb   = 0.0_sp
        smbs%t_srf = 0.0_sp

    end subroutine smb_simple_init

    !-----------------------------------------------------------------
    !> Define the target ice mask. If mask_file is set, read slice
    !> mask_idx along the 3rd dimension of mask_var (cells > 0 are ice);
    !> otherwise derive the mask from the reference ice thickness
    !> H_ice_ref (cells > 0 are ice). Call again after H_ice_ref is
    !> updated to keep the mask consistent.
    !-----------------------------------------------------------------
    subroutine smb_simple_set_mask(smbs, H_ice_ref)

        type(smb_simple_class), intent(inout) :: smbs
        real(sp),               intent(in)    :: H_ice_ref(:,:)

        integer :: nx, ny
        real(sp), allocatable :: tmp(:,:)

        if (.not. allocated(smbs%mask)) then
            error stop "smb_simple_set_mask: smbs not initialized"
        end if
        if (size(H_ice_ref, 1) /= size(smbs%mask, 1) .or. &
            size(H_ice_ref, 2) /= size(smbs%mask, 2)) then
            error stop "smb_simple_set_mask: H_ice_ref shape mismatch"
        end if

        nx = size(smbs%mask, 1)
        ny = size(smbs%mask, 2)

        if (len_trim(smbs%mask_file) > 0) then
            ! Read the mask_idx slice along the 3rd (e.g. time/index) dimension
            allocate(tmp(nx, ny))
            call nc_read(smbs%mask_file, smbs%mask_var, tmp, &
                         start=[1, 1, smbs%mask_idx], count=[nx, ny, 1])
            smbs%mask = (tmp > 0.0_sp)
            deallocate(tmp)
        else
            smbs%mask = (H_ice_ref > 0.0_sp)
        end if

    end subroutine smb_simple_set_mask

    !-----------------------------------------------------------------
    !> Compute the SMB (mm w.e./yr) and surface temperature (K) fields
    !> for the active scheme, storing them in smbs%smb and smbs%t_srf.
    !> Only the synthetic-elevation scheme ("syn") is currently wired;
    !> "tg24" is not yet supported.
    !-----------------------------------------------------------------
    subroutine smb_simple_update(smbs, z_srf, t_sl)

        type(smb_simple_class), intent(inout) :: smbs
        real(sp),               intent(in)    :: z_srf(:,:)
        real(sp),               intent(in)    :: t_sl(:,:)

        select case (trim(smbs%scheme))

            case ("syn")
                call calc_smb_simple_syn(smbs%smb, smbs%t_srf, z_srf, smbs%mask, &
                                         smbs%x, smbs%y, smbs%lat, t_sl, &
                                         smbs%co2, smbs%f, smbs%par, units=smbs%units)

            case default
                write(*,*) "smb_simple_update: scheme not supported: ", trim(smbs%scheme)
                write(*,*) "Only scheme='syn' is currently implemented."
                error stop "smb_simple_update: unsupported scheme"

        end select

    end subroutine smb_simple_update

    !=================================================================
    ! Production entry point — synthetic-elevation scheme
    !=================================================================

    !-----------------------------------------------------------------
    !> Compute the SMB field (mm w.e./yr) and surface-temperature
    !> field t_srf (K) of the synthetic-elevation scheme. One call, no
    !> persistent state.
    !>
    !> Pipeline:
    !>   1. signed distance to target-mask boundary (m), from the 2-D
    !>      coordinate pairs x/y interpreted per `units` ("m"/"km"
    !>      Cartesian, "degrees" lon/lat)
    !>   2. z_syn from (d, z_sur, mask_target) per p%use_plastic
    !>      (.true. → plastic Nye/Vialov; .false. → linear wedge)
    !>   3. SMB = acc - abl kernel on z_syn (dz_SL = 0), with acc the maritime
    !>      accumulation field from compute_acc; all tuned params are in
    !>      mm w.e./yr, so no unit conversion is needed.
    !>   4. outside-mask SMB floor at p%smb_min
    !>   5. t_srf = t_sl - p%gamma_t*max(z_syn, 0), capped at p%t_ice_max
    !>      over ice-covered (target-mask) cells; ocean / below-sea-level
    !>      cells (z_syn <= 0) keep the imposed t_sl.
    !>
    !> x(:,:), y(:,:) drive the distance (Cartesian or lon/lat per
    !> `units`); lat(:,:) is the latitude in degrees used by the SMB
    !> climatology and is required independently of `units`. `units` is
    !> optional and defaults to "m".
    !-----------------------------------------------------------------
    subroutine calc_smb_simple_syn(smb, t_srf, z_sur, mask_target, x, y, lat, &
                                   t_sl, CO2, f, p, units)
        real(sp), intent(out) :: smb(:,:)
        real(sp), intent(out) :: t_srf(:,:)
        real(sp), intent(in)  :: z_sur(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(sp), intent(in)  :: x(:,:)
        real(sp), intent(in)  :: y(:,:)
        real(sp), intent(in)  :: lat(:,:)
        real(sp), intent(in)  :: t_sl(:,:)
        real(sp), intent(in)  :: CO2
        real(sp), intent(in)  :: f
        type(smb_params_syn), intent(in) :: p
        character(len=*), intent(in), optional :: units

        integer :: nx, ny
        real(sp), allocatable :: d_m(:,:), z_syn(:,:)
        real(sp), allocatable :: dz_SL(:,:), acc(:,:)
        character(len=16) :: units_use

        nx = size(z_sur, 1)
        ny = size(z_sur, 2)

        units_use = "m"
        if (present(units)) units_use = units

        if (size(smb, 1)         /= nx .or. size(smb, 2)         /= ny .or. &
            size(t_srf, 1)       /= nx .or. size(t_srf, 2)       /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny .or. &
            size(t_sl, 1)        /= nx .or. size(t_sl, 2)        /= ny .or. &
            size(x, 1)           /= nx .or. size(x, 2)           /= ny .or. &
            size(y, 1)           /= nx .or. size(y, 2)           /= ny .or. &
            size(lat, 1)         /= nx .or. size(lat, 2)         /= ny) then
            error stop "calc_smb_simple_syn: shape mismatch among inputs"
        end if
        if (p%smb_min > 0.0_sp) then
            error stop "calc_smb_simple_syn: p%smb_min must be <= 0"
        end if
        if (p%L_acc <= 0.0_sp) then
            error stop "calc_smb_simple_syn: p%L_acc must be > 0"
        end if

        allocate(d_m  (nx, ny), z_syn(nx, ny), dz_SL(nx, ny), &
                 acc  (nx, ny))

        call compute_signed_distance(d_m, mask_target, x, y, units_use)

        if (p%use_plastic) then
            call compute_z_syn_plastic(z_syn, d_m, z_sur, mask_target,     &
                                       p%tau0, p%slope_out,                &
                                       p%z_max_in, p%z_max_out,            &
                                       p%rho_ice, p%g)
        else
            call compute_z_syn_linear(z_syn, d_m, z_sur, mask_target,      &
                                      p%slope, p%z_max_in, p%z_max_out)
        end if

        dz_SL = 0.0_sp
        call compute_acc(acc, d_m, lat, p)

        call calc_smb_simple_pwl(smb, z_syn, dz_SL, acc, lat, CO2, f, p)
        call apply_smb_min_outside(smb, mask_target, p%smb_min)
        call compute_t_srf_lapse(t_srf, z_syn, t_sl, mask_target, &
                                 p%gamma_t, p%t_ice_max)

        deallocate(d_m, z_syn, dz_SL, acc)
    end subroutine calc_smb_simple_syn

    !=================================================================
    ! Piecewise-linear SMB(z) kernel — `_pwl`
    !=================================================================

    !-----------------------------------------------------------------
    !> SMB = accumulation - ablation on a fully-specified elevation field z,
    !> with snow-line offset dz_SL(:,:) and accumulation field acc(:,:).
    !> Accumulation applies everywhere; ablation grows linearly with depth
    !> below the snow line and is zero above it, so SMB is continuous across
    !> the snow line. The kernel is agnostic to whether z is the real surface
    !> or a synthetic field.
    !-----------------------------------------------------------------
    subroutine calc_smb_simple_pwl(smb, z, dz_SL, acc, lat, CO2, f, p)
        real(sp),             intent(out) :: smb(:,:)
        real(sp),             intent(in)  :: z(:,:)
        real(sp),             intent(in)  :: dz_SL(:,:)
        real(sp),             intent(in)  :: acc(:,:)
        real(sp),             intent(in)  :: lat(:,:)
        real(sp),             intent(in)  :: CO2
        real(sp),             intent(in)  :: f
        type(smb_params_syn), intent(in)  :: p

        integer  :: i, j, nx, ny
        real(sp) :: df, phi, zSL, beta, dz_abl

        nx = size(z, 1)
        ny = size(z, 2)

        if (.not. (size(smb,1)   == nx .and. size(smb,2)   == ny .and. &
                   size(dz_SL,1) == nx .and. size(dz_SL,2) == ny .and. &
                   size(acc,1)   == nx .and. size(acc,2)   == ny .and. &
                   size(lat,1)   == nx .and. size(lat,2)   == ny)) then
            error stop "calc_smb_simple_pwl: shape mismatch among inputs"
        end if

        df = f - p%fmean

        do j = 1, ny
            do i = 1, nx
                phi    = abs(lat(i, j))    ! poleward latitude: hemisphere-agnostic snow line / ablation
                zSL    = p%a1 + p%a2 * phi + p%a3 * CO2 + p%a4 * df + dz_SL(i, j)
                beta   = max(p%beta_floor, p%beta0 + p%beta1 * (p%phiref - phi))
                dz_abl = max(0.0_sp, zSL - z(i, j))
                smb(i, j) = acc(i, j) - beta * dz_abl
            end do
        end do
    end subroutine calc_smb_simple_pwl

    !-----------------------------------------------------------------
    !> Maritime accumulation field: maximal at the ice margin and decaying
    !> exponentially inland (continentality), optionally modulated by latitude.
    !>
    !>   acc    = c0 * f_lat(phi) * f_dist(d_in)
    !>   f_lat  = clamp(1 - k_acc_lat*max(0, phi - phi_acc_ref), facc_lat_min, 1)
    !>   f_dist = exp(-d_in_km / L_acc)
    !>
    !> d_m is the signed distance to the mask boundary in metres (positive
    !> inside); only the inland part (d_m > 0) reduces accumulation. Outside the
    !> mask d_in = 0 (acc = c0*f_lat), but those cells are overwritten by the
    !> smb_min floor afterwards. With k_acc_lat = 0 the latitude factor is 1.
    !-----------------------------------------------------------------
    subroutine compute_acc(acc, d_m, lat, p)
        real(sp),             intent(out) :: acc(:,:)
        real(sp),             intent(in)  :: d_m(:,:)
        real(sp),             intent(in)  :: lat(:,:)
        type(smb_params_syn), intent(in)  :: p

        integer  :: i, j, nx, ny
        real(sp) :: phi, d_in_km, f_lat, f_dist

        nx = size(acc, 1)
        ny = size(acc, 2)

        if (size(d_m,1) /= nx .or. size(d_m,2) /= ny .or. &
            size(lat,1) /= nx .or. size(lat,2) /= ny) then
            error stop "compute_acc: shape mismatch among inputs"
        end if

        do j = 1, ny
            do i = 1, nx
                phi   = abs(lat(i, j))    ! poleward latitude (hemisphere-agnostic)
                f_lat = 1.0_sp - p%k_acc_lat * max(0.0_sp, phi - p%phi_acc_ref)
                f_lat = max(p%facc_lat_min, min(1.0_sp, f_lat))

                d_in_km = max(0.0_sp, d_m(i, j)) / 1000.0_sp
                f_dist  = exp(-d_in_km / p%L_acc)

                acc(i, j) = p%c0 * f_lat * f_dist
            end do
        end do
    end subroutine compute_acc

    !=================================================================
    ! TG24 2-D kernel
    !=================================================================

    subroutine calc_smb_simple_tg24(smb, z_sur, dz_SL, dzdx, dzdy, mask, &
                                    lat, CO2, f, p, use_mask)
        real(sp),              intent(out) :: smb(:,:)
        real(sp),              intent(in)  :: z_sur(:,:)
        real(sp),              intent(in)  :: dz_SL(:,:)
        real(sp),              intent(in)  :: dzdx(:,:)
        real(sp),              intent(in)  :: dzdy(:,:)
        real(sp),              intent(in)  :: mask(:,:)
        real(sp),              intent(in)  :: lat(:,:)
        real(sp),              intent(in)  :: CO2
        real(sp),              intent(in)  :: f
        type(smb_params_tg24), intent(in)  :: p
        logical, optional,     intent(in)  :: use_mask

        integer  :: i, j, nx, ny
        real(sp) :: df, phi, z_SL, dz, smb_cm, grad_z_mag
        logical  :: apply_mask

        nx = size(z_sur, 1)
        ny = size(z_sur, 2)

        if (.not. (size(smb,1)   == nx .and. size(smb,2)   == ny .and. &
                   size(dz_SL,1) == nx .and. size(dz_SL,2) == ny .and. &
                   size(dzdx,1)  == nx .and. size(dzdx,2)  == ny .and. &
                   size(dzdy,1)  == nx .and. size(dzdy,2)  == ny .and. &
                   size(mask,1)  == nx .and. size(mask,2)  == ny .and. &
                   size(lat,1)   == nx .and. size(lat,2)   == ny)) then
            error stop "calc_smb_simple_tg24: shape mismatch among inputs"
        end if

        apply_mask = .false.
        if (present(use_mask)) apply_mask = use_mask

        df = f - p%fmean

        do j = 1, ny
            do i = 1, nx
                phi        = abs(lat(i, j))    ! poleward latitude (hemisphere-agnostic)
                z_SL       = p%a1 + p%a2 * phi + p%a3 * CO2 + p%a4 * df + dz_SL(i, j)
                dz         = z_sur(i, j) - z_SL
                grad_z_mag = sqrt(dzdx(i, j)**2 + dzdy(i, j)**2)
                if (dz > 0.0_sp) then
                    smb_cm = (p%b1 + p%b4 * grad_z_mag) * exp(-dz / p%b5)
                else
                    smb_cm = -p%b2 * dz**2
                end if
                if (apply_mask) then
                    smb_cm = smb_cm + p%b3 * (mask(i, j) - 1.0_sp)
                end if
                smb(i, j) = smb_cm * 10.0_sp     ! cm/yr -> mm w.e./yr (density ratio in tuned params)
            end do
        end do
    end subroutine calc_smb_simple_tg24

    !=================================================================
    ! Signed distance to mask boundary (m)
    !=================================================================

    !-----------------------------------------------------------------
    !> Per-cell signed distance to the nearest target-mask boundary, in
    !> metres. Positive inside the mask, negative outside.
    !>
    !> x(:,:), y(:,:) are the (unique) 2-D coordinates of each cell.
    !> `units` selects how they are interpreted (optional, default "m"):
    !>   "m"       Cartesian, already metres: sqrt(dx^2 + dy^2)
    !>   "km"      Cartesian in km; result converted to metres (x1000)
    !>   "degrees" x = lon, y = lat (deg); local flat-Earth metric
    !>             (R = 6.371e6 m, mid-latitude cos(phi) factor for dlon,
    !>             longitude wrap to [-180, 180)).
    !> The result is always in metres, regardless of `units`.
    !-----------------------------------------------------------------
    subroutine compute_signed_distance(d, mask_target, x, y, units)
        real(sp), intent(out) :: d(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(sp), intent(in)  :: x(:,:)
        real(sp), intent(in)  :: y(:,:)
        character(len=*), intent(in), optional :: units

        integer  :: nx, ny, i, j, ib, jb, k, nb, kk
        logical  :: m, is_bnd, latlon
        real(sp) :: xq, yq, xb, yb, dx, dy, dphi, dlam
        real(sp) :: phi_mid, dy_m, dx_m, dist, dmin, scale
        character(len=16) :: units_use
        integer, allocatable :: bi(:), bj(:)
        logical, allocatable :: bm(:)
        real(sp), parameter :: R_earth = 6.371e6_sp   ! m
        real(sp), parameter :: deg2rad = 3.141592653589793_sp / 180.0_sp

        nx = size(mask_target, 1)
        ny = size(mask_target, 2)

        units_use = "m"
        if (present(units)) units_use = units

        if (size(d, 1) /= nx .or. size(d, 2) /= ny) then
            error stop "compute_signed_distance: d shape /= mask_target"
        end if
        if (size(x, 1) /= nx .or. size(x, 2) /= ny .or. &
            size(y, 1) /= nx .or. size(y, 2) /= ny) then
            error stop "compute_signed_distance: x, y shape /= mask_target"
        end if
        if (nx < 2 .or. ny < 2) then
            error stop "compute_signed_distance: nx, ny must be >= 2"
        end if

        select case (trim(units_use))
        case ("degrees")
            latlon = .true.;  scale = 1.0_sp
        case ("km")
            latlon = .false.; scale = 1000.0_sp
        case ("m")
            latlon = .false.; scale = 1.0_sp
        case default
            error stop "compute_signed_distance: units must be 'm', 'km', or 'degrees'"
        end select

        ! Pass 1: count boundary cells (4-connectivity).
        nb = 0
        do j = 1, ny
            do i = 1, nx
                m      = mask_target(i, j)
                is_bnd = .false.
                if (i > 1) then
                    if (mask_target(i - 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. i < nx) then
                    if (mask_target(i + 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j > 1) then
                    if (mask_target(i, j - 1) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j < ny) then
                    if (mask_target(i, j + 1) .neqv. m) is_bnd = .true.
                end if
                if (is_bnd) nb = nb + 1
            end do
        end do

        if (nb == 0) then
            if (mask_target(1, 1)) then
                d = huge(1.0_sp)
            else
                d = -huge(1.0_sp)
            end if
            return
        end if

        allocate(bi(nb), bj(nb), bm(nb))

        ! Pass 2: fill boundary lists.
        k = 0
        do j = 1, ny
            do i = 1, nx
                m      = mask_target(i, j)
                is_bnd = .false.
                if (i > 1) then
                    if (mask_target(i - 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. i < nx) then
                    if (mask_target(i + 1, j) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j > 1) then
                    if (mask_target(i, j - 1) .neqv. m) is_bnd = .true.
                end if
                if (.not. is_bnd .and. j < ny) then
                    if (mask_target(i, j + 1) .neqv. m) is_bnd = .true.
                end if
                if (is_bnd) then
                    k = k + 1
                    bi(k) = i
                    bj(k) = j
                    bm(k) = m
                end if
            end do
        end do

        ! Pass 3: min distance to opposite-mask boundary, per query cell.
        do j = 1, ny
            do i = 1, nx
                m     = mask_target(i, j)
                xq    = x(i, j)
                yq    = y(i, j)
                dmin  = huge(1.0_sp)
                do kk = 1, nb
                    if (bm(kk) .eqv. m) cycle
                    ib    = bi(kk)
                    jb    = bj(kk)
                    xb    = x(ib, jb)
                    yb    = y(ib, jb)
                    if (latlon) then
                        dphi    = yq - yb
                        dlam    = xq - xb
                        dlam    = modulo(dlam + 180.0_sp, 360.0_sp) - 180.0_sp
                        phi_mid = 0.5_sp * (yq + yb)
                        dy_m    = R_earth * dphi * deg2rad
                        dx_m    = R_earth * cos(phi_mid * deg2rad) * dlam * deg2rad
                        dist    = sqrt(dy_m * dy_m + dx_m * dx_m)
                    else
                        dx   = xq - xb
                        dy   = yq - yb
                        dist = scale * sqrt(dx * dx + dy * dy)
                    end if
                    if (dist < dmin) dmin = dist
                end do
                if (m) then
                    d(i, j) =  dmin
                else
                    d(i, j) = -dmin
                end if
            end do
        end do

        deallocate(bi, bj, bm)
    end subroutine compute_signed_distance

    !=================================================================
    ! Synthetic-elevation profiles (linear and plastic)
    !=================================================================

    !-----------------------------------------------------------------
    !> Linear wedge profile both inside and outside the mask. d_m is the
    !> signed distance in metres and slope is in m/m.
    !>     z_syn_raw  = clamp(slope * d_m, -z_max_out, +z_max_in)
    !>     z_syn(i,j) = max(z_syn_raw, z_sur(i,j))  if  mask_target(i,j)
    !>                = min(z_syn_raw, z_sur(i,j))  otherwise
    !-----------------------------------------------------------------
    subroutine compute_z_syn_linear(z_syn, d_m, z_sur, mask_target, &
                                    slope, z_max_in, z_max_out)
        real(sp), intent(out) :: z_syn(:,:)
        real(sp), intent(in)  :: d_m(:,:)
        real(sp), intent(in)  :: z_sur(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(sp), intent(in)  :: slope
        real(sp), intent(in)  :: z_max_in
        real(sp), intent(in)  :: z_max_out

        integer  :: nx, ny, i, j
        real(sp) :: zr

        nx = size(z_syn, 1)
        ny = size(z_syn, 2)

        if (size(d_m, 1)         /= nx .or. size(d_m, 2)         /= ny .or. &
            size(z_sur, 1)       /= nx .or. size(z_sur, 2)       /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "compute_z_syn_linear: shape mismatch among inputs"
        end if
        if (z_max_in < 0.0_sp .or. z_max_out < 0.0_sp) then
            error stop "compute_z_syn_linear: z_max_in and z_max_out must be >= 0"
        end if

        do j = 1, ny
            do i = 1, nx
                zr = slope * d_m(i, j)
                if (zr >  z_max_in)  zr =  z_max_in
                if (zr < -z_max_out) zr = -z_max_out
                if (mask_target(i, j)) then
                    z_syn(i, j) = max(zr, z_sur(i, j))
                else
                    z_syn(i, j) = min(zr, z_sur(i, j))
                end if
            end do
        end do
    end subroutine compute_z_syn_linear

    !-----------------------------------------------------------------
    !> Perfect-plasticity (Nye/Vialov) profile inside; linear outside.
    !> d_m is the signed distance in metres and slope_out is in m/m.
    !>   inside  (d >= 0):  z_syn_raw = min(z_max_in,  C * sqrt(d_m))
    !>                       C = sqrt(2 * tau0 / (rho_ice*g))
    !>   outside (d <  0):  z_syn_raw = max(-z_max_out, slope_out*d_m)
    !-----------------------------------------------------------------
    subroutine compute_z_syn_plastic(z_syn, d_m, z_sur, mask_target,         &
                                     tau0, slope_out, z_max_in, z_max_out,   &
                                     rho_ice, g)
        real(sp), intent(out) :: z_syn(:,:)
        real(sp), intent(in)  :: d_m(:,:)
        real(sp), intent(in)  :: z_sur(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(sp), intent(in)  :: tau0
        real(sp), intent(in)  :: slope_out
        real(sp), intent(in)  :: z_max_in
        real(sp), intent(in)  :: z_max_out
        real(sp), intent(in)  :: rho_ice
        real(sp), intent(in)  :: g

        integer  :: nx, ny, i, j
        real(sp) :: zr, d, C

        nx = size(z_syn, 1)
        ny = size(z_syn, 2)

        if (size(d_m, 1)         /= nx .or. size(d_m, 2)         /= ny .or. &
            size(z_sur, 1)       /= nx .or. size(z_sur, 2)       /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "compute_z_syn_plastic: shape mismatch among inputs"
        end if
        if (z_max_in < 0.0_sp .or. z_max_out < 0.0_sp) then
            error stop "compute_z_syn_plastic: z_max_in and z_max_out must be >= 0"
        end if
        if (tau0 <= 0.0_sp) then
            error stop "compute_z_syn_plastic: tau0 must be > 0"
        end if

        if (rho_ice <= 0.0_sp .or. g <= 0.0_sp) then
            error stop "compute_z_syn_plastic: rho_ice and g must be > 0"
        end if

        C = sqrt(2.0_sp * tau0 / (rho_ice * g))

        do j = 1, ny
            do i = 1, nx
                d = d_m(i, j)
                if (d >= 0.0_sp) then
                    zr = C * sqrt(d)
                    if (zr > z_max_in) zr = z_max_in
                else
                    zr = slope_out * d
                    if (zr < -z_max_out) zr = -z_max_out
                end if
                if (mask_target(i, j)) then
                    z_syn(i, j) = max(zr, z_sur(i, j))
                else
                    z_syn(i, j) = min(zr, z_sur(i, j))
                end if
            end do
        end do
    end subroutine compute_z_syn_plastic

    !=================================================================
    ! Outside-mask SMB floor
    !=================================================================

    !-----------------------------------------------------------------
    !> In-place SMB floor outside the target mask. For each cell with
    !> .not. mask_target(i,j) and smb(i,j) > smb_min, set
    !> smb(i,j) = smb_min. Inside-mask cells are untouched.
    !-----------------------------------------------------------------
    subroutine apply_smb_min_outside(smb, mask_target, smb_min)
        real(sp), intent(inout) :: smb(:,:)
        logical,  intent(in)    :: mask_target(:,:)
        real(sp), intent(in)    :: smb_min

        integer :: i, j, nx, ny

        nx = size(smb, 1)
        ny = size(smb, 2)

        if (size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "apply_smb_min_outside: mask_target shape /= smb"
        end if

        do j = 1, ny
            do i = 1, nx
                if (.not. mask_target(i, j) .and. smb(i, j) > smb_min) then
                    smb(i, j) = smb_min
                end if
            end do
        end do
    end subroutine apply_smb_min_outside

    !=================================================================
    ! Surface temperature from sea-level temperature + lapse rate
    !=================================================================

    !-----------------------------------------------------------------
    !> Surface temperature (K) from a sea-level temperature field and a
    !> lapse-rate correction applied to elevation:
    !>
    !>     t_srf(i,j) = t_sl(i,j) - gamma_t * max(z(i,j), 0)
    !>
    !> Below-sea-level / ocean cells (z <= 0) keep the imposed t_sl.
    !> Over ice-covered cells (mask_target(i,j) .true.) the result is
    !> capped at t_ice_max. z is typically the synthetic elevation z_syn.
    !-----------------------------------------------------------------
    subroutine compute_t_srf_lapse(t_srf, z, t_sl, mask_target, &
                                   gamma_t, t_ice_max)
        real(sp), intent(out) :: t_srf(:,:)
        real(sp), intent(in)  :: z(:,:)
        real(sp), intent(in)  :: t_sl(:,:)
        logical,  intent(in)  :: mask_target(:,:)
        real(sp), intent(in)  :: gamma_t
        real(sp), intent(in)  :: t_ice_max

        integer  :: i, j, nx, ny
        real(sp) :: t

        nx = size(t_srf, 1)
        ny = size(t_srf, 2)

        if (size(z, 1)           /= nx .or. size(z, 2)           /= ny .or. &
            size(t_sl, 1)        /= nx .or. size(t_sl, 2)        /= ny .or. &
            size(mask_target, 1) /= nx .or. size(mask_target, 2) /= ny) then
            error stop "compute_t_srf_lapse: shape mismatch among inputs"
        end if

        do j = 1, ny
            do i = 1, nx
                t = t_sl(i, j) - gamma_t * max(z(i, j), 0.0_sp)
                if (mask_target(i, j) .and. t > t_ice_max) then
                    t = t_ice_max
                end if
                t_srf(i, j) = t
            end do
        end do
    end subroutine compute_t_srf_lapse

    !=================================================================
    ! Namelist loading (yelmox integration)
    !=================================================================

    !-----------------------------------------------------------------
    !> Load smb_simple control settings and the synthetic-scheme (`syn`)
    !> parameters from a yelmox-style namelist file.
    !>
    !>   group : control fields + the `syn` parameters
    !>           (default "smb_simple")
    !>
    !> The TG24 scheme is not yet wired into yelmox, so its parameter
    !> type is not loaded here.
    !-----------------------------------------------------------------
    subroutine smb_simple_par_load(smbs, filename, group, init, cnst)

        type(smb_simple_class), intent(inout) :: smbs
        character(len=*),       intent(in)    :: filename
        character(len=*),       intent(in), optional :: group
        logical,                intent(in), optional :: init
        ! Physical constants come from here instead of the namelist when
        ! present, so the &smb_simple group need not declare rho_ice and g.
        type(phys_const_class), intent(in), optional :: cnst

        ! Local variables
        logical           :: init_pars
        character(len=32) :: nml_group

        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "smb_simple"
        end if

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        ! Control fields
        call nml_read(filename,nml_group,"scheme",    smbs%scheme,    init=init_pars)
        call nml_read(filename,nml_group,"co2_const", smbs%co2,       init=init_pars)
        call nml_read(filename,nml_group,"f_const",   smbs%f,         init=init_pars)
        call nml_read(filename,nml_group,"mask_file", smbs%mask_file, init=init_pars)
        call nml_read(filename,nml_group,"mask_var",  smbs%mask_var,  init=init_pars)
        call nml_read(filename,nml_group,"mask_idx",  smbs%mask_idx,  init=init_pars)

        ! Synthetic-elevation (syn) scheme parameters
        call nml_read(filename,nml_group,"a1",         smbs%par%a1,         init=init_pars)
        call nml_read(filename,nml_group,"a2",         smbs%par%a2,         init=init_pars)
        call nml_read(filename,nml_group,"a3",         smbs%par%a3,         init=init_pars)
        call nml_read(filename,nml_group,"a4",         smbs%par%a4,         init=init_pars)
        call nml_read(filename,nml_group,"fmean",      smbs%par%fmean,      init=init_pars)
        call nml_read(filename,nml_group,"beta0",      smbs%par%beta0,      init=init_pars)
        call nml_read(filename,nml_group,"beta1",      smbs%par%beta1,      init=init_pars)
        call nml_read(filename,nml_group,"phiref",     smbs%par%phiref,     init=init_pars)
        call nml_read(filename,nml_group,"beta_floor", smbs%par%beta_floor, init=init_pars)
        call nml_read(filename,nml_group,"c0",         smbs%par%c0,         init=init_pars)
        call nml_read(filename,nml_group,"L_acc",      smbs%par%L_acc,      init=init_pars)
        call nml_read(filename,nml_group,"k_acc_lat",    smbs%par%k_acc_lat,    init=init_pars)
        call nml_read(filename,nml_group,"phi_acc_ref",  smbs%par%phi_acc_ref,  init=init_pars)
        call nml_read(filename,nml_group,"facc_lat_min", smbs%par%facc_lat_min, init=init_pars)
        call nml_read(filename,nml_group,"use_plastic",smbs%par%use_plastic,init=init_pars)
        call nml_read(filename,nml_group,"slope",      smbs%par%slope,      init=init_pars)
        call nml_read(filename,nml_group,"tau0",       smbs%par%tau0,       init=init_pars)
        call nml_read(filename,nml_group,"slope_out",  smbs%par%slope_out,  init=init_pars)
        call nml_read(filename,nml_group,"z_max_in",   smbs%par%z_max_in,   init=init_pars)
        call nml_read(filename,nml_group,"z_max_out",  smbs%par%z_max_out,  init=init_pars)
        if (present(cnst)) then
            call phys_const_require(cnst, "smb_simple_init")
            call phys_const_get(cnst, "rho_ice", smbs%par%rho_ice)
            call phys_const_get(cnst, "g",       smbs%par%g)
        else
            call nml_read(filename,nml_group,"rho_ice",smbs%par%rho_ice,init=init_pars)
            call nml_read(filename,nml_group,"g",      smbs%par%g,      init=init_pars)
        end if
        call nml_read(filename,nml_group,"gamma_t",    smbs%par%gamma_t,    init=init_pars)
        call nml_read(filename,nml_group,"t_ice_max",  smbs%par%t_ice_max,  init=init_pars)
        call nml_read(filename,nml_group,"smb_min",    smbs%par%smb_min,    init=init_pars)

    end subroutine smb_simple_par_load

end module smb_simple_m

program test_surface_smb_simple
    ! surface = smb_simple (make check): the backend and its kernels against the literal
    ! yelmox module (above), the same inputs in, bit-identical smb and T_srf out.
    !
    !   kernels   calc_smb_simple_syn of yelmox and of kryos_smb_simple on the same fields
    !             (plastic and linear profile, units m / km / degrees, with and without
    !             accumulation latitude terms, a mask without boundary, isolated cells)
    !   backend   the component (surface_create, surface_init, PASS_BOUNDARY, update) against the
    !             yelmox call chain of step_surface (kryos_coupling.f90:477-484): smb_simple_init,
    !             smb_simple_set_mask(H_ice_ref), smb_simple_update(z_srf, tsl_ann); the mask
    !             follows H_ice_ref when it changes (recon); the scheme with a mask_file (read once)
    !             imports no H_ice_ref; a plastic and a linear parameter set
    !   contract  the declared imports (z_srf, tsl_ann, H_ice_ref; no tas, pr, H_ice, z_bed, z_sl),
    !             the exports (smb, T_srf), no restart file, MEAN publication
    ! On a polar stereographic grid of 32 km cells (NH-32KM like) with synthetic fields.

    use precision,        only: sp, dp
    use ncio,             only: nc_create, nc_write_dim, nc_write
    use coordinates,      only: grid_class, grid_init
    use phys_constants,   only: phys_const_class, phys_const_set, sec_year_tropical
    use fesm_fields,      only: fields_get, fields_put
    use kryos_component,  only: init_ctx_t, step_time_t, census_t, PASS_BOUNDARY
    use kryos_surface,    only: surface_slot_t, surface_create, surface_init, surface_init_state, surface_update, surface_end
    use smb_simple_m,     only: y_class => smb_simple_class, y_init => smb_simple_init, y_mask => smb_simple_set_mask, &
                                y_update => smb_simple_update, y_syn => calc_smb_simple_syn, y_par => smb_params_syn
    use kryos_smb_simple, only: k_syn => calc_smb_simple_syn, k_par => smb_params_syn

    implicit none

    integer, parameter :: nx = 48, ny = 56
    character(len=*), parameter :: fpar = "test_smb_simple.nml", fmask = "test_smb_mask_North_TEST-SMB.nc"

    type(grid_class)       :: grid
    type(phys_const_class) :: cnst
    real(sp), allocatable  :: x(:,:), y(:,:), lat(:,:), z_srf(:,:), tsl(:,:), href(:,:), href2(:,:), tsl2(:,:), z2(:,:)
    integer :: nfail, seed

    nfail = 0
    seed  = 20261029
    call phys_const_set(cnst, label="Earth", rho_ice=910.0_dp, rho_sw=1028.0_dp, rho_w=1000.0_dp, g=9.81_dp)
    call grid_init(grid, name="TEST-SMB", mtype="polar_stereographic", units="kilometers", planet="WGS84", &
                   lambda=-45.0_dp, phi=70.0_dp, x0=-800.0_dp, dx=32.0_dp, nx=nx, y0=-3300.0_dp, dy=32.0_dp, ny=ny)
    allocate(x(nx,ny), y(nx,ny), lat(nx,ny))
    x   = real(grid%x * grid%cs%xy_conv, sp)
    y   = real(grid%y * grid%cs%xy_conv, sp)
    lat = real(grid%lat, sp)
    call check("grid: latitudes of a northern domain", minval(lat) > 40.0_sp .and. maxval(lat) < 90.0_sp)

    call fields(z_srf, tsl, href, 1)
    call fields(z2, tsl2, href2, 2)
    call write_par()
    call write_mask(href2)

    call test_kernels()
    call test_backend("surface_nh", "smb_simple",     .false.)
    call test_backend("surface_lin", "smb_simple_lin", .false.)
    call test_backend("surface_file", "smb_simple_file", .true.)

    call execute_command_line("rm -f "//fpar//" "//fmask)
    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_surface_smb_simple: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_surface_smb_simple: OK"

contains

    ! =========================================================================
    ! The kernels
    ! =========================================================================

    subroutine test_kernels()
        type(y_par) :: yp
        type(k_par) :: kp
        real(sp), allocatable :: ys(:,:), yt(:,:), ks(:,:), kt(:,:)
        logical,  allocatable :: mask(:,:)
        character(len=16) :: units(3) = [character(len=16) :: "m", "km", "degrees"]
        real(sp), allocatable :: xx(:,:), yy(:,:)
        integer :: iu, ip

        allocate(ys(nx,ny), yt(nx,ny), ks(nx,ny), kt(nx,ny), mask(nx,ny), xx(nx,ny), yy(nx,ny))
        mask = href > 0.0_sp
        write(*,'(a)') "--- kernels: calc_smb_simple_syn, yelmox against kryos"
        do ip = 1, 4
            call params(ip, yp, kp)
            do iu = 1, 3
                select case (iu)
                case (1);  xx = x;            yy = y
                case (2);  xx = x*1e-3_sp;    yy = y*1e-3_sp
                case (3);  xx = real(grid%lon, sp); yy = lat
                end select
                call y_syn(ys, yt, z_srf, mask, xx, yy, lat, tsl, 190.0_sp, 494.5_sp, yp, units=trim(units(iu)))
                call k_syn(ks, kt, z_srf, mask, xx, yy, lat, tsl, 190.0_sp, 494.5_sp, kp, units=trim(units(iu)))
                call check("parameter set "//itoa(ip)//", units "//trim(units(iu))//": smb and t_srf bit-identical", &
                           all(ys == ks) .and. all(yt == kt))
            end do
        end do

        ! a mask without a boundary (all land, all ocean) and isolated cells
        call params(1, yp, kp)
        mask = .true.
        call y_syn(ys, yt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, yp)
        call k_syn(ks, kt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, kp)
        call check("mask all ice: bit-identical", all(ys == ks) .and. all(yt == kt))
        mask = .false.
        call y_syn(ys, yt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, yp)
        call k_syn(ks, kt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, kp)
        call check("mask no ice: bit-identical (smb at or below the floor)", all(ys == ks) .and. all(yt == kt) .and. &
               all(ks <= kp%smb_min))
        mask(5,7) = .true.;  mask(20,30) = .true.;  mask(21,30) = .true.
        call y_syn(ys, yt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, yp)
        call k_syn(ks, kt, z_srf, mask, x, y, lat, tsl, 190.0_sp, 494.5_sp, kp)
        call check("isolated ice cells: bit-identical", all(ys == ks) .and. all(yt == kt))
    end subroutine test_kernels

    subroutine params(ip, yp, kp)
        ! Four parameter sets: the module defaults, the NH set of yelmox_North.nml, linear profile
        ! with the accumulation latitude term, a snow line with CO2 and insolation terms.
        integer,     intent(in)  :: ip
        type(y_par), intent(out) :: yp
        type(k_par), intent(out) :: kp

        select case (ip)
        case (2)
            yp%a1 = 3900.0_sp;  yp%beta0 = 5.0_sp;  yp%c0 = 300.0_sp;  yp%L_acc = 300.0_sp
            kp%a1 = 3900.0_sp;  kp%beta0 = 5.0_sp;  kp%c0 = 300.0_sp;  kp%L_acc = 300.0_sp
        case (3)
            yp%use_plastic = .false.;  yp%slope = 5.0e-3_sp;  yp%k_acc_lat = 0.02_sp;  yp%phi_acc_ref = 55.0_sp
            kp%use_plastic = .false.;  kp%slope = 5.0e-3_sp;  kp%k_acc_lat = 0.02_sp;  kp%phi_acc_ref = 55.0_sp
        case (4)
            yp%a3 = 2.0_sp;  yp%a4 = 10.0_sp;  yp%beta_floor = 1.0_sp;  yp%tau0 = 80.0e3_sp;  yp%t_ice_max = 268.0_sp
            kp%a3 = 2.0_sp;  kp%a4 = 10.0_sp;  kp%beta_floor = 1.0_sp;  kp%tau0 = 80.0e3_sp;  kp%t_ice_max = 268.0_sp
        end select
    end subroutine params

    ! =========================================================================
    ! The backend
    ! =========================================================================

    subroutine test_backend(inst, group, with_file)
        character(len=*), intent(in) :: inst, group
        logical,          intent(in) :: with_file

        type(surface_slot_t) :: sf
        type(init_ctx_t)     :: ctx
        type(census_t)       :: census
        type(step_time_t)    :: st
        type(y_class)        :: ys
        real(dp), allocatable :: e(:,:)
        real(sp), allocatable :: hmask(:,:)
        character(len=:), allocatable :: tag

        tag = inst//" ["//group//"]: "
        write(*,'(a)') "--- backend "//inst//" ("//group//")"
        ctx%par_file  = fpar
        ctx%instance  = inst
        ctx%component = "surface"
        ctx%backend   = "smb_simple"
        ctx%scope     = "regional"
        ctx%domain    = "North"
        ctx%grid_name = "TEST-SMB"
        ctx%grid      = grid
        ctx%cnst      = cnst
        ctx%sec_year  = sec_year_tropical
        call surface_create(sf, ctx, census)
        call surface_init(sf, ctx)

        associate (c => sf%c)
            call check(tag//"imports z_srf and tsl_ann, none of tas, pr, H_ice, z_bed, z_sl", &
                       c%hi_z_srf > 0 .and. c%hi_tsl_ann > 0 .and. c%hi_tas == 0 .and. c%hi_pr == 0 .and. &
                       c%hi_H_ice == 0 .and. c%hi_z_bed == 0 .and. c%hi_z_sl == 0)
            call check(tag//"H_ice_ref imported iff there is no mask file", (c%hi_H_ice_ref > 0) .eqv. (.not. with_file))
            call check(tag//"exports smb and T_srf only", c%exports%n == 2 .and. c%h_smb > 0 .and. c%h_T_srf > 0 .and. &
                       c%h_melt == 0 .and. c%h_runoff == 0 .and. c%h_refrz == 0 .and. c%h_albedo == 0)
            call check(tag//"stateless: no restart file", .not. c%con%restart)
        end associate

        ! yelmox: smb_simple_init (kryos.f90:395-406), then step_surface (kryos_coupling.f90:477-484)
        call y_init(ys, fpar, x=x, y=y, lat=lat, group=group, units="m", cnst=cnst)
        if (with_file) then
            call nml_replace_(ys%mask_file, "{domain}", "North")
            call nml_replace_(ys%mask_file, "{grid_name}", "TEST-SMB")
        end if

        ! step 1: PASS_BOUNDARY on (z_srf, tsl, H_ice_ref)
        st = step_time_t(t0=-10.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=10.0_dp, dtt=10.0_dp)
        call put_imports(sf, z_srf, tsl, href, 0.0_dp, with_file)
        call surface_init_state(sf, st, PASS_BOUNDARY)
        call y_mask(ys, href)
        call y_update(ys, z_srf, tsl)
        call compare(sf, ys, tag//"PASS_BOUNDARY")

        ! step 2: a new surface, temperature and reference ice (recon), update
        st = step_time_t(t0=0.0_dp, t1=10.0_dp, t_eval=10.0_dp, dt=10.0_dp, dtt=10.0_dp)
        call put_imports(sf, z2, tsl2, href2, 10.0_dp, with_file)
        call surface_update(sf, st)
        call y_mask(ys, href2)
        call y_update(ys, z2, tsl2)
        call compare(sf, ys, tag//"update")
        call check(tag//"exports MEAN over [0, 10]", sf%c%exports%f(sf%c%h_smb)%t1 == 10.0_dp .and. &
                   sf%c%exports%f(sf%c%h_smb)%t0 == 0.0_dp .and. sf%c%exports%f(sf%c%h_T_srf)%published)

        ! the two masks differ (H_ice_ref changed), so the follow-up above is not trivial
        if (.not. with_file) then
            hmask = href;  call y_mask(ys, href);  call y_update(ys, z2, tsl2)
            allocate(e(nx,ny));  call fields_get(sf%c%exports, sf%c%h_smb, e)
            call check(tag//"smb depends on the mask (the reference ice of step 2 is not that of step 1)", any(e /= real(ys%smb, dp)))
        end if
        call surface_end(sf)
    end subroutine test_backend

    subroutine put_imports(sf, z, t, h, time, with_file)
        type(surface_slot_t), intent(inout) :: sf
        real(sp),             intent(in)    :: z(:,:), t(:,:), h(:,:)
        real(dp),             intent(in)    :: time
        logical,              intent(in)    :: with_file
        call fields_put(sf%c%imports, sf%c%hi_z_srf,   z, t1=time)
        call fields_put(sf%c%imports, sf%c%hi_tsl_ann, t, t1=time)
        if (.not. with_file) call fields_put(sf%c%imports, sf%c%hi_H_ice_ref, h, t1=time)
    end subroutine put_imports

    subroutine compare(sf, ys, what)
        type(surface_slot_t), intent(inout) :: sf
        type(y_class),        intent(in)    :: ys
        character(len=*),     intent(in)    :: what
        real(dp), allocatable :: e(:,:), f(:,:)
        allocate(e(nx,ny), f(nx,ny))
        call fields_get(sf%c%exports, sf%c%h_smb,   e)
        call fields_get(sf%c%exports, sf%c%h_T_srf, f)
        call check(what//": smb bit-identical to yelmox smb_simple_update", all(e == real(ys%smb, dp)))
        call check(what//": T_srf bit-identical to yelmox smb_simple_update", all(f == real(ys%t_srf, dp)))
        call check(what//": the fields are not trivial (smb spans the snow line, T_srf below freezing)", &
                   minval(ys%smb) < -100.0_sp .and. maxval(ys%smb) > 0.0_sp .and. minval(ys%t_srf) < 270.0_sp)
    end subroutine compare

    ! =========================================================================
    ! Inputs
    ! =========================================================================

    subroutine fields(z, t, h, variant)
        ! Synthetic surface elevation, sea-level temperature and reference ice (an ice cap with a
        ! hole and a satellite; variant 2: another cap, as after a reconstruction).
        real(sp), allocatable, intent(out) :: z(:,:), t(:,:), h(:,:)
        integer,               intent(in)  :: variant
        integer  :: i, j
        real(sp) :: r, cx, cy, noise

        allocate(z(nx,ny), t(nx,ny), h(nx,ny))
        cx = merge(22.0_sp, 28.0_sp, variant == 1)
        cy = merge(28.0_sp, 24.0_sp, variant == 1)
        do j = 1, ny
            do i = 1, nx
                r = sqrt((real(i,sp) - cx)**2 + (real(j,sp) - cy)**2)
                h(i,j) = max(0.0_sp, 2800.0_sp*(1.0_sp - (r/16.0_sp)**2))
                if (abs(real(i,sp) - cx - 4.0_sp) < 2.0_sp .and. abs(real(j,sp) - cy) < 2.0_sp) h(i,j) = 0.0_sp
                if (i > 38 .and. j > 44 .and. i < 42 .and. j < 48) h(i,j) = 300.0_sp
                noise = lcg() - 0.5_sp
                z(i,j) = h(i,j) + 150.0_sp*noise + merge(0.0_sp, 400.0_sp*(1.0_sp - r/40.0_sp), r > 40.0_sp)
                t(i,j) = 262.0_sp - 0.45_sp*(real(lat(i,j), sp) - 60.0_sp) + 2.0_sp*(lcg() - 0.5_sp)*real(variant, sp)
            end do
        end do
    end subroutine fields

    real(sp) function lcg()
        ! A deterministic uniform number in [0, 1).
        seed = int(mod(1103515245_8*int(seed, 8) + 12345_8, 2147483648_8))
        lcg = real(seed, sp)/2147483648.0_sp
    end function lcg

    subroutine write_mask(h)
        ! The mask file: three slices on the surface grid, slice 2 the mask of the run (mask_idx = 2).
        real(sp), intent(in) :: h(:,:)
        real(sp) :: m(nx,ny,3)
        m(:,:,1) = 0.0_sp
        m(:,:,2) = merge(1.0_sp, 0.0_sp, h > 0.0_sp)
        m(:,:,3) = 1.0_sp
        call nc_create(fmask)
        call nc_write_dim(fmask, "xc", x=real(grid%G%x, dp), units="km")
        call nc_write_dim(fmask, "yc", x=real(grid%G%y, dp), units="km")
        call nc_write_dim(fmask, "idx", x=[1.0_dp, 2.0_dp, 3.0_dp])
        call nc_write(fmask, "mask", m, dim1="xc", dim2="yc", dim3="idx")
    end subroutine write_mask

    subroutine write_par()
        ! Three instances: the NH set of yelmox_North.nml (mask from H_ice_ref), a linear profile
        ! with accumulation latitude terms, and the NH set with a mask file.
        integer :: u
        open(newunit=u, file=fpar, status="replace", action="write")
        write(u,'(a)') "&surface_nh", "  backend = 'smb_simple'", "/", &
                       "&surface_lin", "  backend = 'smb_simple'", "  smb_simple = 'smb_simple_lin'", "/", &
                       "&surface_file", "  backend = 'smb_simple'", "  smb_simple = 'smb_simple_file'", "/"
        call group(u, "smb_simple", 3900.0, 5.0, 300.0, 300.0, "True", "''", 1)
        call group(u, "smb_simple_lin", 3750.0, 4.0, 150.0, 1000.0, "False", "''", 1)
        call group(u, "smb_simple_file", 3900.0, 5.0, 300.0, 300.0, "True", "'"//fmask//"'", 2)
        close(u)
    end subroutine write_par

    subroutine group(u, name, a1, beta0, c0, L_acc, plastic, mask_file, idx)
        ! The full group of a yelmox par file (yelmox reads every key).
        integer,          intent(in) :: u, idx
        character(len=*), intent(in) :: name, plastic, mask_file
        real,             intent(in) :: a1, beta0, c0, L_acc
        character(len=:), allocatable :: mf
        mf = mask_file
        if (len(mf) > 2) mf = "'test_smb_mask_{domain}_{grid_name}.nc'"
        write(u,'(a)') "&"//name, "  scheme = 'syn'", "  co2_const = 190.0", "  f_const = 494.5", &
                       "  mask_file = "//mf, "  mask_var = 'mask'", "  mask_idx = "//itoa(idx), &
                       "  a1 = "//ftoa(a1), "  a2 = -50.0", "  a3 = 0.0", "  a4 = 0.0", "  fmean = 494.5", &
                       "  beta0 = "//ftoa(beta0), "  beta1 = 0.05", "  phiref = 60.0", "  beta_floor = 0.0", &
                       "  c0 = "//ftoa(c0), "  L_acc = "//ftoa(L_acc), &
                       "  k_acc_lat = "//merge("0.02", "0.0 ", name == "smb_simple_lin"), &
                       "  phi_acc_ref = "//merge("55.0", "0.0 ", name == "smb_simple_lin"), "  facc_lat_min = 0.3", &
                       "  use_plastic = "//plastic, "  slope = 5.0e-3", "  tau0 = 36.0e3", "  slope_out = 6.09e-3", &
                       "  z_max_in = 2500.0", "  z_max_out = 200.0", "  gamma_t = 6.5e-3", "  t_ice_max = 273.15", &
                       "  smb_min = -2000.0", "/"
    end subroutine group

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    subroutine nml_replace_(s, text, rep)
        use nml, only: nml_replace
        character(len=*), intent(inout) :: s
        character(len=*), intent(in)    :: text, rep
        call nml_replace(s, text, rep)
    end subroutine nml_replace_

    function itoa(i) result(s)
        integer, intent(in) :: i
        character(len=:), allocatable :: s
        character(len=16) :: b
        write(b,'(i0)') i
        s = trim(b)
    end function itoa

    function ftoa(v) result(s)
        real, intent(in) :: v
        character(len=:), allocatable :: s
        character(len=32) :: b
        write(b,'(f0.4)') v
        s = trim(b)
    end function ftoa

    subroutine check(what, ok)
        character(len=*), intent(in) :: what
        logical,          intent(in) :: ok
        if (ok) then
            write(*,'(a)') "  ok    "//what
        else
            write(*,'(a)') "  FAIL  "//what
            nfail = nfail + 1
        end if
    end subroutine check

end program test_surface_smb_simple
