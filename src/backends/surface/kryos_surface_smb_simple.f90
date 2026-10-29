module kryos_smb_simple
    ! The kernels of the synthetic-elevation SMB scheme `syn` of yelmox
    ! libs/smb_simple.f90 (commit 5ac89c8, branch topo-prep), the expressions
    ! unchanged: the parameter type smb_params_syn, calc_smb_simple_syn (signed
    ! distance to the mask boundary, synthetic elevation, piecewise-linear SMB,
    ! outside-mask floor, surface temperature), calc_smb_simple_pwl, compute_acc,
    ! compute_signed_distance, compute_z_syn_linear/plastic, apply_smb_min_outside
    ! and compute_t_srf_lapse. Left out: the class (smb_simple_class and its
    ! init/set_mask/update; the backend below owns the grid, the mask and the
    ! fields), the namelist loader, and the TG24 scheme (never wired into yelmox).
    ! tests/test_surface_smb_simple.f90 compares them with the literal yelmox
    ! module.
    !
    ! Snow line, with phi = |latitude| (hemisphere-agnostic):
    !     z_SL = a1 + a2*phi + a3*CO2 + a4*(f - fmean)
    ! SMB = acc - abl, acc = c0 * f_lat(phi) * f_dist(d_in) (maritime, decays
    ! inland), abl = beta_abl(phi) * max(0, z_SL - z), evaluated on the synthetic
    ! elevation z_syn of the target mask. SMB in mm w.e./yr, t_srf in K.

    use precision, only: sp

    implicit none
    private

    public :: smb_params_syn
    public :: calc_smb_simple_syn
    public :: calc_smb_simple_pwl
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

        ! Plastic profile physical constants. Taken from the run's phys_const_class
        ! by the backend; these Earth values are the standalone fallback.
        real(sp) :: rho_ice = 910.0_sp        ! kg/m^3
        real(sp) :: g       =   9.81_sp       ! m/s^2

        ! Surface temperature  t_srf = t_sl - gamma_t*max(z_syn, 0)
        real(sp) :: gamma_t   = 6.5e-3_sp     ! K/m (lapse rate)
        real(sp) :: t_ice_max = 273.15_sp     ! K, clamp over ice-covered cells

        ! Outside-mask SMB floor (must be <= 0)
        real(sp) :: smb_min = -2000.0_sp      ! mm w.e./yr

    end type smb_params_syn

contains


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

end module kryos_smb_simple

module kryos_surface_smb_simple
    ! surface backend smb_simple (adapters-design section 5): the synthetic-elevation
    ! scheme `syn` of yelmox, SMB and surface temperature from the surface elevation, the
    ! sea-level air temperature and a target ice mask alone, no library. Stateless:
    ! every call computes the fields from the imports, nothing is kept and no restart
    ! file is written (constraints restart = .false.). The kernels are kryos_smb_simple
    ! (above), computed in sp: imports are narrowed once, exports widened.
    !
    ! Imports: z_srf (hub), tsl_ann (climate; yelmox clim%now%tsl_ann) and H_ice_ref (ice;
    ! yelmox Yelmo bnd%H_ice_ref, every call: it follows recon). The mask is
    ! H_ice_ref > 0, or, with mask_file set, a slice of a file on the surface grid read
    ! once at init (then H_ice_ref is not imported). The par group named by the instance
    ! key `smb_simple` ([smb_simple]) has the yelmox keys, the schema is the template
    ! &smb_simple of input/kryos_defaults.nml. smb is exported as kg m-2 yr-1 = mm w.e./yr
    ! (the number yelmox lands on Yelmo before its unit factor).

    use precision,       only: sp, dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_validate, nml_has_param, nml_replace
    use ncio,            only: nc_read
    use phys_constants,  only: phys_const_require, phys_const_get
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, NAME_LEN, PATH_LEN, PASS_BOUNDARY, PASS_FINAL
    use kryos_surface,   only: surface_t, surface_opt_t
    use kryos_smb_simple, only: smb_params_syn, calc_smb_simple_syn

    implicit none

    private
    public :: smb_simple_surface_t, smb_simple_surface_constraints

    character(len=*), parameter :: DEF_FILE = "input/kryos_defaults.nml"   ! template &smb_simple

    type, extends(surface_t) :: smb_simple_surface_t
        character(len=NAME_LEN) :: group = "smb_simple"      ! the parameter group ([smb_simple])
        type(smb_params_syn)    :: par
        real(sp) :: co2 = 280.0_sp, f = 494.5_sp             ! [ppm], [W m-2] constants of the snow line
        character(len=PATH_LEN) :: mask_file = ""            ! "" => the mask is H_ice_ref > 0
        character(len=NAME_LEN) :: mask_var  = ""            ! its variable and
        integer                 :: mask_idx  = 1             ! the index along its third dimension
        real(sp), allocatable :: x(:,:), y(:,:), lat(:,:)    ! [m], [m], [deg N] of the surface grid
        logical,  allocatable :: mask(:,:)                   ! the target ice mask
        real(sp), allocatable :: z_srf_s(:,:), tsl_s(:,:), Href_s(:,:), smb(:,:), t_srf(:,:)
    contains
        procedure :: init          => smb_simple_init_
        procedure :: init_state    => smb_simple_init_state_
        procedure :: update        => smb_simple_update_
        procedure :: restart_write => smb_simple_restart_write_
        procedure :: restart_read  => smb_simple_restart_read_
    end type smb_simple_surface_t

contains

    function smb_simple_surface_constraints() result(con)
        ! Regional, any grid, several instances; stateless, no restart file.
        type(constraints_t) :: con
        con = constraints_t(restart=.false.)
    end function smb_simple_surface_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine smb_simple_init_(self, ctx)
        ! smb_simple_init (kryos.f90:395-399): the parameters, the grid in metres
        ! and the latitudes; the mask file, if any, once.
        class(smb_simple_surface_t), intent(inout) :: self
        type(init_ctx_t),            intent(in)    :: ctx

        character(len=NAME_LEN) :: scheme
        logical :: from_file

        if (nml_has_param(ctx%par_file, ctx%instance, "smb_simple")) &
            call nml_read(ctx%par_file, ctx%instance, "smb_simple", self%group)
        self%prec = sp
        self%nx   = ctx%grid%G%nx
        self%ny   = ctx%grid%G%ny

        call par_load(self, ctx, scheme)
        if (trim(scheme) /= "syn") then
            call fesm_error("kryos_surface_smb_simple::init", "["//trim(self%group)//"] scheme must be 'syn'", &
                            "instance "//trim(self%name)//": got '"//trim(scheme)//"' (tg24 is not implemented)")
        end if

        from_file = len_trim(self%mask_file) > 0
        self%opt  = surface_opt_t(atm=.false., geom=.false., tsl_ann=.true., H_ice_ref=.not. from_file)

        allocate(self%x(self%nx, self%ny), self%y(self%nx, self%ny), self%lat(self%nx, self%ny))
        self%x   = real(ctx%grid%x * ctx%grid%cs%xy_conv, sp)
        self%y   = real(ctx%grid%y * ctx%grid%cs%xy_conv, sp)
        self%lat = real(ctx%grid%lat, sp)

        allocate(self%mask(self%nx, self%ny), self%smb(self%nx, self%ny), self%t_srf(self%nx, self%ny))
        allocate(self%z_srf_s(self%nx, self%ny), self%tsl_s(self%nx, self%ny), self%Href_s(self%nx, self%ny))
        self%mask = .false.
        if (from_file) call mask_read(self, ctx)

        if (from_file) then
            call log_msg(LOG_INFO, "kryos_surface_smb_simple", "scheme syn: smb [mm w.e./yr] and T_srf from z_srf, tsl_ann "// &
                         "and the mask of "//trim(self%mask_file), "instance "//trim(self%name)//"; stateless")
        else
            call log_msg(LOG_INFO, "kryos_surface_smb_simple", "scheme syn: smb [mm w.e./yr] and T_srf from z_srf, tsl_ann "// &
                         "and the mask H_ice_ref > 0 of the ice", "instance "//trim(self%name)//"; stateless")
        end if
    end subroutine smb_simple_init_

    subroutine smb_simple_init_state_(self, st, pass)
        ! PASS_BOUNDARY and PASS_FINAL: step_surface (kryos_coupling.f90:477-484) on the
        ! boundary and on the final initial geometry (recon changes H_ice_ref). Stateless:
        ! nothing else, nothing on PASS_RESTART (the exports come with the bundle).
        class(smb_simple_surface_t), intent(inout) :: self
        type(step_time_t),           intent(in)    :: st
        integer,                     intent(in)    :: pass

        select case (pass)
        case (PASS_BOUNDARY, PASS_FINAL)
            call compute(self, st)
        end select
    end subroutine smb_simple_init_state_

    subroutine smb_simple_update_(self, st)
        ! step_surface: smb_simple_set_mask, smb_simple_update on the imports of the step.
        class(smb_simple_surface_t), intent(inout) :: self
        type(step_time_t),           intent(in)    :: st
        call compute(self, st)
    end subroutine smb_simple_update_

    subroutine smb_simple_restart_write_(self, dir, t)
        ! Stateless: no file (constraints restart = .false., the wrapper does not call this).
        class(smb_simple_surface_t), intent(inout) :: self
        character(len=*),            intent(in)    :: dir
        real(dp),                    intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine smb_simple_restart_write_

    subroutine smb_simple_restart_read_(self, dir, t)
        class(smb_simple_surface_t), intent(inout) :: self
        character(len=*),            intent(in)    :: dir
        real(dp),                    intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine smb_simple_restart_read_

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    subroutine compute(self, st)
        ! The imports narrowed once; the mask; calc_smb_simple_syn; MEAN over [t0, t1].
        class(smb_simple_surface_t), intent(inout) :: self
        type(step_time_t),           intent(in)    :: st

        call fields_get(self%imports, self%hi_z_srf,   self%z_srf_s)
        call fields_get(self%imports, self%hi_tsl_ann, self%tsl_s)
        if (self%hi_H_ice_ref > 0) then
            call fields_get(self%imports, self%hi_H_ice_ref, self%Href_s)
            self%mask = (self%Href_s > 0.0_sp)
        end if
        call calc_smb_simple_syn(self%smb, self%t_srf, self%z_srf_s, self%mask, self%x, self%y, self%lat, &
                                 self%tsl_s, self%co2, self%f, self%par, units="m")
        call fields_put(self%exports, self%h_smb,   self%smb,   t1=st%t1, t0=st%t0)
        call fields_put(self%exports, self%h_T_srf, self%t_srf, t1=st%t1, t0=st%t0)
    end subroutine compute

    subroutine par_load(self, ctx, scheme)
        ! smb_simple_par_load: the keys of [smb_simple] (template defaults for the ones a
        ! par file leaves out); rho_ice and g from the run's constants, as yelmox passes cnst.
        class(smb_simple_surface_t), intent(inout) :: self
        type(init_ctx_t),            intent(in)    :: ctx
        character(len=*),            intent(out)   :: scheme

        character(len=:), allocatable :: grp, file

        grp  = trim(self%group)
        file = trim(ctx%par_file)
        call nml_validate(file, DEF_FILE, grp, defaults_group="smb_simple")

        call nml_read(file, grp, "scheme",    scheme,       defaults_file=DEF_FILE, defaults_group="smb_simple")
        call nml_read(file, grp, "co2_const", self%co2,     defaults_file=DEF_FILE, defaults_group="smb_simple")
        call nml_read(file, grp, "f_const",   self%f,       defaults_file=DEF_FILE, defaults_group="smb_simple")
        call nml_read(file, grp, "mask_file", self%mask_file, defaults_file=DEF_FILE, defaults_group="smb_simple")
        call nml_read(file, grp, "mask_var",  self%mask_var,  defaults_file=DEF_FILE, defaults_group="smb_simple")
        call nml_read(file, grp, "mask_idx",  self%mask_idx,  defaults_file=DEF_FILE, defaults_group="smb_simple")

        associate (p => self%par)
            call nml_read(file, grp, "a1",           p%a1,           defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "a2",           p%a2,           defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "a3",           p%a3,           defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "a4",           p%a4,           defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "fmean",        p%fmean,        defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "beta0",        p%beta0,        defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "beta1",        p%beta1,        defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "phiref",       p%phiref,       defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "beta_floor",   p%beta_floor,   defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "c0",           p%c0,           defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "L_acc",        p%L_acc,        defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "k_acc_lat",    p%k_acc_lat,    defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "phi_acc_ref",  p%phi_acc_ref,  defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "facc_lat_min", p%facc_lat_min, defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "use_plastic",  p%use_plastic,  defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "slope",        p%slope,        defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "tau0",         p%tau0,         defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "slope_out",    p%slope_out,    defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "z_max_in",     p%z_max_in,     defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "z_max_out",    p%z_max_out,    defaults_file=DEF_FILE, defaults_group="smb_simple")
            call phys_const_require(ctx%cnst, "kryos_surface_smb_simple")
            call phys_const_get(ctx%cnst, "rho_ice", p%rho_ice)
            call phys_const_get(ctx%cnst, "g",       p%g)
            call nml_read(file, grp, "gamma_t",      p%gamma_t,      defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "t_ice_max",    p%t_ice_max,    defaults_file=DEF_FILE, defaults_group="smb_simple")
            call nml_read(file, grp, "smb_min",      p%smb_min,      defaults_file=DEF_FILE, defaults_group="smb_simple")
        end associate
    end subroutine par_load

    subroutine mask_read(self, ctx)
        ! smb_simple_set_mask with mask_file: the slice mask_idx of mask_var on the surface
        ! grid (cells > 0 are ice). {domain} and {grid_name} (the surface grid) are replaced.
        class(smb_simple_surface_t), intent(inout) :: self
        type(init_ctx_t),            intent(in)    :: ctx

        real(sp), allocatable :: tmp(:,:)

        call nml_replace(self%mask_file, "{domain}",    trim(ctx%domain))
        call nml_replace(self%mask_file, "{grid_name}", trim(ctx%grid_name))
        allocate(tmp(self%nx, self%ny))
        call nc_read(trim(self%mask_file), trim(self%mask_var), tmp, &
                     start=[1, 1, self%mask_idx], count=[self%nx, self%ny, 1])
        self%mask = (tmp > 0.0_sp)
    end subroutine mask_read

end module kryos_surface_smb_simple
