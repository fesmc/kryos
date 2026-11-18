module kryos_plastic_taub
    ! The basal shear stress tau_b [Pa] of the ice backend `plastic` (docs/dev/plastic-backend.md), on
    ! the instance grid at a model time. No library: compiled whatever `plastic=`.
    !
    !   "const"    one value everywhere;
    !   "file"     a static gridded field (x, y);
    !   "domains"  Gowan et al. (2021): a map of domain ids (0 = none) and one value per id, static or a
    !              series per id (linear in time), as the product {grid}_IceHistory-PaleoMIST_taub.nc; a par
    !              table (codes, values) overrides the value of some ids, constant in time, or gives ids of
    !              the map the file has no value for; every id of the map must then have one. Cells without
    !              a domain take `default` (0: such a cell in the ice mask stops the solve, kryos_ice_plastic);
    !   "coulomb"  the Coulomb yield law tau_b = c_f lambda(z_bed - z_sl) N_eff (Robinson et al. 2020), solved
    !              with the geometry by the library (plastic_solve_coulomb: N_eff depends on H): here only its
    !              coefficient mu = c_f lambda (taub_mu), with c_f one value per region of a region map (0 =
    !              none: such a cell in the ice mask stops the solve), constant in time, and lambda the bed
    !              factor of kryos_plastic_bed, either a fixed field (taub_set_lambda_fixed) of a bed relative to
    !              present sea level: the field z_lambda of a file (taub_read_lambda_bed; the backend's
    !              taub_lambda_bed = "file") or the run's relaxed reference bed ("relaxed"), or on the bed and sea
    !              level of the solve ("current"); and the law's p, H_N and floor delta (N >= delta rho_i g H).
    !
    ! "domains" may be scaled. A factor per domain (taub_domain_factor_file, the present-day fit of
    ! tests/plastic_taub_fit.f90): the value of an id, after the override table, times exp(log_factor(id)),
    ! constant in time. Then the scale groups (taub_scale_file, the calibration of
    ! ice-earth-recon/tools/plastic_esmda.jl): each id belongs to a scale group (0 = unscaled) and its value is
    ! multiplied by exp(log_factor(group, t)), linear in t between the nodes of the groups and held outside
    ! them. The default of the cells without a domain is scaled by neither.
    !
    ! Any method may be smoothed (taub_set_smooth; the backend's taub_smooth_km): the field at t, after the
    ! factors, is replaced on the cells with tau_b > 0 by its Gaussian mean (sd, within 3 sd) over those
    ! cells (a normalised convolution: steps between domains become ramps over a few sd, values are kept
    ! where the field is uniform). Cells with tau_b <= 0 stay as they are. "coulomb" smooths mu the same way.

    use precision, only: dp
    use fesm_log,  only: fesm_error, to_str
    use ncio,      only: nc_read, nc_read_attr, nc_dims, nc_exists_var
    use kryos_plastic_margin, only: model_time_axis, time_order
    use kryos_plastic_bed,    only: bed_lambda

    implicit none

    private

    integer, parameter, public :: TAUB_M_CONST = 1, TAUB_M_FILE = 2, TAUB_M_DOMAINS = 3, TAUB_M_COULOMB = 4

    type, public :: taub_t
        integer  :: method = 0
        real(dp) :: const = 0.0_dp                   ! [Pa] "const"
        real(dp) :: default = 0.0_dp                 ! [Pa] "domains": cells without a domain
        logical  :: hold = .false.                   ! "domains" series: hold the end values outside their range
        real(dp), allocatable :: field(:,:)          ! [Pa] "file"
        integer,  allocatable :: dom(:,:)            ! "domains": the domain id of each cell (<= 0: none)
        integer,  allocatable :: idx(:,:)            ! its position in ids (0: none)
        integer,  allocatable :: ids(:)              ! the domain ids with a value
        real(dp), allocatable :: val(:)              ! [Pa] their static value
        real(dp), allocatable :: val_t(:,:)          ! [Pa] (id, time): their series; unallocated: static
        real(dp), allocatable :: time(:)             ! [yr] model time of the series, increasing
        integer,  allocatable :: sgrp(:)             ! the scale group of each id (0: unscaled); unallocated: no scaling
        real(dp), allocatable :: stime(:)            ! [yr] model time of the scale nodes, increasing
        real(dp), allocatable :: slf(:,:)            ! (group, node): natural log of the multiplier
        real(dp), allocatable :: dfac(:)             ! the factor of each id (exp of its log factor); unallocated: none
        integer,  allocatable :: reg(:,:)            ! "coulomb": the region id of each cell (<= 0: none)
        integer,  allocatable :: ridx(:,:)           ! its position in rids (0: none)
        integer,  allocatable :: rids(:)             ! the region ids with a c_f
        real(dp), allocatable :: cf(:)               ! their c_f (exp of the log)
        real(dp), allocatable :: lam(:,:)            ! "coulomb": a fixed lambda (unallocated: on the current bed)
        real(dp) :: p = 1.0_dp                       ! "coulomb": N = rho_i g H max((1 - H_f/H)^p, delta)
        real(dp) :: H_N = 0.0_dp                     ! [m] "coulomb": cap N_max = rho_i g H_N
        real(dp) :: delta = 0.05_dp                  ! "coulomb": floor of N as a fraction of the overburden
        real(dp) :: lam_min = 1.0_dp                 ! "coulomb": lambda at a bed zb0 or lower
        real(dp) :: zb0 = -300.0_dp, zb1 = 200.0_dp  ! [m] "coulomb": bed above sea level of lam_min, of lambda = 1
        character(len=512) :: csrc = ""              ! messages: the c_f file
        real(dp) :: smooth = 0.0_dp                  ! [m] Gaussian sd of the smoothing (0: none)
        real(dp) :: sdx = 0.0_dp, sdy = 0.0_dp       ! [m] the grid spacing of the smoothing
        character(len=512) :: dsrc = ""              ! messages: the domain factor file
        character(len=512) :: ssrc = ""              ! messages: the scale file
        character(len=512) :: src = ""               ! messages
    end type taub_t

    public :: taub_method_code, taub_set_const, taub_set_field, taub_read_field
    public :: taub_set_domains, taub_read_domains, taub_at
    public :: taub_set_scale, taub_read_scale, taub_scale_at
    public :: taub_set_domain_factor, taub_read_domain_factor
    public :: taub_set_smooth
    public :: taub_set_coulomb, taub_read_coulomb, taub_mu, taub_set_lambda_fixed, taub_read_lambda_bed

contains

    integer function taub_method_code(name) result(code)
        ! "const" | "file" | "domains" | "coulomb" -> TAUB_*; 0 for anything else.
        character(len=*), intent(in) :: name
        select case (trim(name))
        case ("const");   code = TAUB_M_CONST
        case ("file");    code = TAUB_M_FILE
        case ("domains"); code = TAUB_M_DOMAINS
        case ("coulomb"); code = TAUB_M_COULOMB
        case default;     code = 0
        end select
    end function taub_method_code

    subroutine taub_set_const(tb, value)
        type(taub_t), intent(out) :: tb
        real(dp),     intent(in)  :: value
        if (.not. value > 0.0_dp) then
            call fesm_error("kryos_plastic_taub::taub_set_const", "taub_const must be positive", &
                            "taub_const = "//to_str(value)//" Pa")
        end if
        tb%method = TAUB_M_CONST
        tb%const  = value
        tb%src    = "taub_const"
    end subroutine taub_set_const

    subroutine taub_set_field(tb, field, src)
        type(taub_t),     intent(out) :: tb
        real(dp),         intent(in)  :: field(:,:)
        character(len=*), intent(in)  :: src
        tb%method = TAUB_M_FILE
        tb%field  = field
        tb%src    = src
    end subroutine taub_set_field

    subroutine taub_read_field(tb, file, var, nx, ny)
        ! The static field `var` (x, y) [Pa] of `file` on the instance grid.
        type(taub_t),     intent(out) :: tb
        character(len=*), intent(in)  :: file, var
        integer,          intent(in)  :: nx, ny

        real(dp), allocatable :: f(:,:)

        call check_xy(file, var, nx, ny)
        allocate(f(nx, ny))
        call nc_read(file, var, f)
        call taub_set_field(tb, f, trim(file)//" ("//trim(var)//")")
    end subroutine taub_read_field

    subroutine taub_set_domains(tb, dom, ids, val, default, src, val_t, time, hold, codes, values, table)
        ! The domain map `dom` (<= 0: none), the values `val` [Pa] of the ids `ids` and, optionally, their
        ! series val_t(id, time) at model times `time` [yr] (any order, distinct); then the override
        ! table `codes`, `values` [Pa] (named `table` in messages).
        type(taub_t),     intent(out) :: tb
        integer,          intent(in)  :: dom(:,:), ids(:)
        real(dp),         intent(in)  :: val(:), default
        character(len=*), intent(in)  :: src
        real(dp),         intent(in), optional :: val_t(:,:), time(:)
        logical,          intent(in), optional :: hold
        integer,          intent(in), optional :: codes(:)
        real(dp),         intent(in), optional :: values(:)
        character(len=*), intent(in), optional :: table

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_domains"
        integer, allocatable :: p(:)
        integer :: k

        if (size(ids) /= size(val)) call fesm_error(proc, "one value per domain id is needed", trim(src))
        do k = 1, size(ids)
            if (ids(k) <= 0) call fesm_error(proc, "domain ids are positive (0 = no domain)", trim(src)// &
                                             ": id "//to_str(ids(k)))
            if (count(ids == ids(k)) > 1) call fesm_error(proc, "domain id "//to_str(ids(k))//" is listed twice", &
                                                         trim(src))
        end do
        if (default < 0.0_dp) call fesm_error(proc, "taub_default must be >= 0", "taub_default = "//to_str(default))

        tb%method  = TAUB_M_DOMAINS
        tb%dom     = dom
        tb%ids     = ids
        tb%val     = val
        tb%default = default
        tb%src     = src
        if (present(hold)) tb%hold = hold
        if (present(val_t)) then
            if (.not. present(time)) call fesm_error(proc, "a series needs its times", trim(src))
            if (size(val_t, 1) /= size(ids) .or. size(val_t, 2) /= size(time) .or. size(time) == 0) then
                call fesm_error(proc, "the series is not (id, time)", trim(src))
            end if
            p = time_order(time)
            do k = 2, size(p)
                if (.not. time(p(k)) > time(p(k-1))) call fesm_error(proc, "two records of the series have the "// &
                                                                     "same time", trim(src))
            end do
            tb%time  = time(p)
            tb%val_t = val_t(:, p)
        end if
        if (present(codes)) then
            if (.not. (present(values) .and. present(table))) call fesm_error(proc, "an override needs its values", src)
            call override(tb, codes, values, table)
        end if
        call index_cells(tb)
    end subroutine taub_set_domains

    subroutine taub_read_domains(tb, file, dom_var, id_var, val_var, ser_var, default, hold, nx, ny, codes, values, table)
        ! The product {grid}_IceHistory-PaleoMIST_taub.nc: dom_var (x, y) the domain id of each cell (<= 0: none),
        ! id_var (n) the ids, val_var (n) their tau_b [Pa]; ser_var (n, time) their series when the
        ! file has it (either dimension order; the time coordinate is the variable named like its
        ! dimension, converted from its units). The override table as taub_set_domains.
        type(taub_t),     intent(out) :: tb
        character(len=*), intent(in)  :: file, dom_var, id_var, val_var, ser_var
        real(dp),         intent(in)  :: default
        logical,          intent(in)  :: hold
        integer,          intent(in)  :: nx, ny
        integer,          intent(in), optional :: codes(:)
        real(dp),         intent(in), optional :: values(:)
        character(len=*), intent(in), optional :: table

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_domains"
        character(len=64), allocatable :: names(:)
        integer,  allocatable :: dims(:), dom(:,:), ids(:)
        real(dp), allocatable :: val(:), ser(:,:), t(:)
        character(len=256) :: units
        character(len=:), allocatable :: src
        integer :: n, it
        logical :: exists

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the tau_b domain file does not exist", "taub_domain_file = "//trim(file))
        src = trim(file)

        call check_xy(file, dom_var, nx, ny)
        allocate(dom(nx, ny))
        call nc_read(file, dom_var, dom)

        call nc_dims(file, id_var, names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`"//trim(id_var)//"` is not a list of ids (1-D)", src)
        n = dims(1)
        call nc_dims(file, val_var, names, dims)
        if (size(dims) /= 1 .or. dims(1) /= n) then
            call fesm_error(proc, "`"//trim(val_var)//"` is not one value per id ("//to_str(n)//", 1-D)", &
                            src//"; taub_value_var names the per-domain values")
        end if
        allocate(ids(n), val(n))
        call nc_read(file, id_var, ids)
        call nc_read(file, val_var, val)

        if (.not. nc_exists_var(file, ser_var)) then
            call taub_set_domains(tb, dom, ids, val, default, src, hold=hold, codes=codes, values=values, table=table)
            return
        end if

        call nc_dims(file, ser_var, names, dims)
        if (size(dims) /= 2) call fesm_error(proc, "`"//trim(ser_var)//"` is not a series (id, time)", src)
        if (dims(1) == n) then
            it = 2
        else if (dims(2) == n) then
            it = 1
        else
            call fesm_error(proc, "`"//trim(ser_var)//"` has no dimension of the "//to_str(n)//" ids", src)
        end if
        allocate(t(dims(it)), ser(dims(1), dims(2)))
        call nc_read(file, trim(names(it)), t)
        units = ""
        call nc_read_attr(file, trim(names(it)), "units", units)
        t = model_time_axis(t, units)
        call nc_read(file, ser_var, ser)
        if (it == 1) ser = transpose(ser)
        call taub_set_domains(tb, dom, ids, val, default, src, val_t=ser, time=t, hold=hold, codes=codes, values=values, &
                              table=table)
    end subroutine taub_read_domains

    subroutine override(tb, codes, values, src)
        ! The par table: the value [Pa] of the domain ids `codes`, constant in time; an id the file has
        ! no value for is added (it must occur in the map).
        type(taub_t),     intent(inout) :: tb
        integer,          intent(in)    :: codes(:)
        real(dp),         intent(in)    :: values(:)
        character(len=*), intent(in)    :: src

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_domains"
        real(dp), allocatable :: w(:,:)
        integer :: k, p

        if (size(codes) /= size(values)) call fesm_error(proc, "one value per code is needed", src)
        do k = 1, size(codes)
            if (.not. values(k) > 0.0_dp) call fesm_error(proc, "tau_b must be positive", src//": code "// &
                                                          to_str(codes(k))//", "//to_str(values(k))//" Pa")
            p = findloc(tb%ids, codes(k), dim=1)
            if (p == 0) then
                if (.not. any(tb%dom == codes(k))) then
                    call fesm_error(proc, "code "//to_str(codes(k))//" is not a domain of the map", src// &
                                    "; map "//trim(tb%src))
                end if
                if (allocated(tb%val_t)) then
                    allocate(w(size(tb%ids) + 1, size(tb%time)))
                    w(1:size(tb%ids), :) = tb%val_t
                    w(size(tb%ids) + 1, :) = values(k)
                    call move_alloc(w, tb%val_t)
                end if
                tb%ids = [tb%ids, codes(k)]
                tb%val = [tb%val, values(k)]
            else
                tb%val(p) = values(k)
                if (allocated(tb%val_t)) tb%val_t(p, :) = values(k)
            end if
        end do
    end subroutine override

    subroutine taub_at(tb, t, tau)
        ! tau_b [Pa] at model time t [yr] on every cell (the solver reads it on the ice mask).
        type(taub_t), intent(in)  :: tb
        real(dp),     intent(in)  :: t
        real(dp),     intent(out) :: tau(:,:)

        real(dp), allocatable :: v(:)
        integer :: i, j

        select case (tb%method)
        case (TAUB_M_CONST)
            tau = tb%const
        case (TAUB_M_FILE)
            tau = tb%field
        case (TAUB_M_DOMAINS)
            v = values_at(tb, t)
            if (allocated(tb%dfac)) v = v*tb%dfac
            if (allocated(tb%sgrp)) call scale_values(tb, t, v)
            do j = 1, size(tau, 2)
                do i = 1, size(tau, 1)
                    if (tb%idx(i,j) > 0) then
                        tau(i,j) = v(tb%idx(i,j))
                    else
                        tau(i,j) = tb%default
                    end if
                end do
            end do
        case (TAUB_M_COULOMB)
            call fesm_error("kryos_plastic_taub::taub_at", "taub_method = 'coulomb' has no tau_b before the solve", &
                            trim(tb%src)//": the solve takes mu (taub_mu) and returns tau_b")
        case default
            call fesm_error("kryos_plastic_taub::taub_at", "tau_b is not set", trim(tb%src))
        end select
        if (tb%smooth > 0.0_dp) call smooth_positive(tau, tb%sdx, tb%sdy, tb%smooth)
    end subroutine taub_at

    subroutine taub_set_coulomb(tb, reg, ids, log_cf, p, H_N, delta, lam_min, zb0, zb1, src, csrc)
        ! "coulomb": the region map `reg` (<= 0: none), c_f = exp(log_cf) of the region ids `ids` (an id not in
        ! the map is allowed; a region of the map without c_f stops), the law's p in [0, 1], H_N > 0 [m], floor
        ! delta in (0, 1] and the bed factor's lam_min in (0, 1], zb0 < zb1 [m]. src, csrc name the map and the
        ! c_f table in messages.
        type(taub_t),     intent(out) :: tb
        integer,          intent(in)  :: reg(:,:), ids(:)
        real(dp),         intent(in)  :: log_cf(:), p, H_N, delta, lam_min, zb0, zb1
        character(len=*), intent(in)  :: src, csrc

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_coulomb"
        integer :: k, i, j, q

        if (size(ids) /= size(log_cf)) call fesm_error(proc, "one log_cf per region id is needed", trim(csrc))
        if (.not. all(abs(log_cf) <= huge(1.0_dp))) call fesm_error(proc, "log_cf is not finite", trim(csrc))
        do k = 1, size(ids)
            if (count(ids == ids(k)) > 1) call fesm_error(proc, "region id "//to_str(ids(k))//" is listed twice", &
                                                         trim(csrc))
        end do
        if (.not. (p >= 0.0_dp .and. p <= 1.0_dp)) call fesm_error(proc, "taub_p must be in [0, 1]", &
                                                                  "taub_p = "//to_str(p))
        if (.not. (H_N > 0.0_dp)) call fesm_error(proc, "taub_HN must be positive", "taub_HN = "//to_str(H_N)//" m")
        if (.not. (delta > 0.0_dp .and. delta <= 1.0_dp)) call fesm_error(proc, "taub_delta must be in (0, 1]", &
                                                                          "taub_delta = "//to_str(delta))
        if (.not. (lam_min > 0.0_dp .and. lam_min <= 1.0_dp)) &
            call fesm_error(proc, "taub_lam_min must be in (0, 1]", "taub_lam_min = "//to_str(lam_min))
        if (lam_min < 1.0_dp .and. .not. zb1 > zb0) &
            call fesm_error(proc, "taub_zb1 must be above taub_zb0", "taub_zb0 = "//to_str(zb0)//", taub_zb1 = "// &
                            to_str(zb1))

        tb%method  = TAUB_M_COULOMB
        tb%reg     = reg
        tb%rids    = ids
        tb%cf      = exp(log_cf)
        tb%p       = p
        tb%H_N     = H_N
        tb%delta   = delta
        tb%lam_min = lam_min
        tb%zb0     = zb0
        tb%zb1     = zb1
        tb%src     = src
        tb%csrc    = csrc
        allocate(tb%ridx(size(reg, 1), size(reg, 2)))
        tb%ridx = 0
        do j = 1, size(reg, 2)
            do i = 1, size(reg, 1)
                if (reg(i,j) <= 0) cycle
                q = findloc(ids, reg(i,j), dim=1)
                if (q == 0) call fesm_error(proc, "region "//to_str(reg(i,j))//" of the map has no c_f", trim(csrc)// &
                                            "; map "//trim(src)//" (cell "//to_str(i)//", "//to_str(j)//")")
                tb%ridx(i,j) = q
            end do
        end do
    end subroutine taub_set_coulomb

    subroutine taub_read_coulomb(tb, region_file, cf_file, nx, ny, p, H_N, delta, lam_min, zb0, zb1)
        ! "coulomb" from the region map file (region (x, y), the ids, 0 = none; ice-earth-recon/tools/plastic_cf_regions.jl)
        ! and the c_f file (dimension region: region_id(region), log_cf(region), natural log), then taub_set_coulomb.
        type(taub_t),     intent(out) :: tb
        character(len=*), intent(in)  :: region_file, cf_file
        integer,          intent(in)  :: nx, ny
        real(dp),         intent(in)  :: p, H_N, delta, lam_min, zb0, zb1

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_coulomb"
        character(len=64), allocatable :: names(:)
        integer,  allocatable :: dims(:), reg(:,:), ids(:)
        real(dp), allocatable :: lcf(:)
        integer :: n
        logical :: exists

        inquire(file=trim(region_file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the region file does not exist", "taub_region_file = "// &
                                          trim(region_file))
        inquire(file=trim(cf_file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the c_f file does not exist", "taub_cf_file = "//trim(cf_file))
        call check_xy(region_file, "region", nx, ny)
        allocate(reg(nx, ny))
        call nc_read(region_file, "region", reg)
        call nc_dims(cf_file, "region_id", names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`region_id` is not a list of ids (1-D)", trim(cf_file))
        n = dims(1)
        call nc_dims(cf_file, "log_cf", names, dims)
        if (size(dims) /= 1 .or. dims(1) /= n) call fesm_error(proc, "`log_cf` is not one value per region_id", &
                                                              trim(cf_file))
        allocate(ids(n), lcf(n))
        call nc_read(cf_file, "region_id", ids)
        call nc_read(cf_file, "log_cf", lcf)
        call taub_set_coulomb(tb, reg, ids, lcf, p, H_N, delta, lam_min, zb0, zb1, trim(region_file), trim(cf_file))
    end subroutine taub_read_coulomb

    function taub_read_lambda_bed(file, nx, ny) result(z)
        ! [m] The bed of a fixed lambda (taub_set_lambda_fixed), relative to present sea level: the field z_lambda
        ! (x, y) of `file` on the nx x ny grid (ice-earth-recon/tools/plastic_lambda_bed.jl), every value finite. The backend
        ! (taub_lambda_bed = "file") and the calibration programs (tests/plastic_coulomb_*.f90) read it with this.
        character(len=*), intent(in) :: file
        integer,          intent(in) :: nx, ny
        real(dp), allocatable :: z(:,:)
        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_lambda_bed"
        logical :: exists

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the lambda bed file does not exist", "taub_lambda_file = "//trim(file))
        if (.not. nc_exists_var(file, "z_lambda")) &
            call fesm_error(proc, "the lambda bed file has no `z_lambda`", trim(file))
        call check_xy(file, "z_lambda", nx, ny)
        allocate(z(nx, ny))
        call nc_read(file, "z_lambda", z)
        if (.not. all(abs(z) <= 1.0e5_dp)) call fesm_error(proc, "`z_lambda` is not finite and within 100 km everywhere", &
                                                           trim(file))
    end function taub_read_lambda_bed

    subroutine taub_set_lambda_fixed(tb, z_rel)
        ! "coulomb" (set before): lambda fixed in time, of the bed z_rel [m] relative to present sea level (the field
        ! of taub_read_lambda_bed, or the relaxed present-day bed of kryos_plastic_bed relaxed_bed); taub_mu then
        ! ignores its z_bed, z_sl.
        type(taub_t), intent(inout) :: tb
        real(dp),     intent(in)    :: z_rel(:,:)
        if (tb%method /= TAUB_M_COULOMB) call fesm_error("kryos_plastic_taub::taub_set_lambda_fixed", &
                                                         "a fixed lambda needs taub_method = 'coulomb'", trim(tb%src))
        tb%lam = bed_lambda(z_rel, 0.0_dp, tb%lam_min, tb%zb0, tb%zb1)
    end subroutine taub_set_lambda_fixed

    subroutine taub_mu(tb, z_bed, z_sl, mu)
        ! "coulomb": mu = c_f(region) lambda on the cells with a region, 0 elsewhere, then the smoothing of
        ! taub_set_smooth over the cells with mu > 0; lambda the fixed field (taub_set_lambda_fixed) or of
        ! z_bed - z_sl.
        type(taub_t), intent(in)  :: tb
        real(dp),     intent(in)  :: z_bed(:,:), z_sl(:,:)
        real(dp),     intent(out) :: mu(:,:)
        integer :: i, j

        if (tb%method /= TAUB_M_COULOMB) call fesm_error("kryos_plastic_taub::taub_mu", "mu needs taub_method = "// &
                                                         "'coulomb'", trim(tb%src))
        do j = 1, size(mu, 2)
            do i = 1, size(mu, 1)
                if (tb%ridx(i,j) > 0 .and. allocated(tb%lam)) then
                    mu(i,j) = tb%cf(tb%ridx(i,j))*tb%lam(i,j)
                else if (tb%ridx(i,j) > 0) then
                    mu(i,j) = tb%cf(tb%ridx(i,j))*bed_lambda(z_bed(i,j), z_sl(i,j), tb%lam_min, tb%zb0, tb%zb1)
                else
                    mu(i,j) = 0.0_dp
                end if
            end do
        end do
        if (tb%smooth > 0.0_dp) call smooth_positive(mu, tb%sdx, tb%sdy, tb%smooth)
    end subroutine taub_mu

    subroutine taub_set_smooth(tb, sd, dx, dy)
        ! The Gaussian smoothing of taub_at: sd [m] (0: none) on a grid of spacing dx, dy [m].
        type(taub_t), intent(inout) :: tb
        real(dp),     intent(in)    :: sd, dx, dy
        if (sd < 0.0_dp) call fesm_error("kryos_plastic_taub::taub_set_smooth", "the smoothing sd must be >= 0", &
                                         "taub_smooth_km = "//to_str(sd/1000.0_dp))
        if (sd > 0.0_dp .and. .not. (dx > 0.0_dp .and. dy > 0.0_dp)) &
            call fesm_error("kryos_plastic_taub::taub_set_smooth", "the grid spacing must be positive", trim(tb%src))
        tb%smooth = sd
        tb%sdx = dx
        tb%sdy = dy
    end subroutine taub_set_smooth

    subroutine smooth_positive(f, dx, dy, sd)
        ! f on its cells with f > 0 replaced by sum w f / sum w over those cells, w = exp(-r^2 / (2 sd^2)) within
        ! 3 sd: separable (rows, then columns) for the numerator and the weight sum alike.
        real(dp), intent(inout) :: f(:,:)
        real(dp), intent(in)    :: dx, dy, sd
        real(dp), allocatable :: num(:,:), den(:,:)
        logical,  allocatable :: v(:,:)

        v = f > 0.0_dp
        num = merge(f, 0.0_dp, v)
        den = merge(1.0_dp, 0.0_dp, v)
        call gauss_1d(num, den, 1, dx)
        call gauss_1d(num, den, 2, dy)
        where (v .and. den > 0.0_dp) f = num/den
    contains
        subroutine gauss_1d(a, b, dim, h)
            real(dp), intent(inout) :: a(:,:), b(:,:)
            integer,  intent(in)    :: dim
            real(dp), intent(in)    :: h
            real(dp), allocatable :: w(:), ta(:,:), tb2(:,:)
            integer :: r, k, i, j, n1, n2
            r = ceiling(3.0_dp*sd/h)
            w = [(exp(-(k*h)**2/(2.0_dp*sd*sd)), k = -r, r)]
            n1 = size(a, 1);  n2 = size(a, 2)
            allocate(ta(n1, n2), tb2(n1, n2))
            ta = 0.0_dp;  tb2 = 0.0_dp
            !$omp parallel do private(i, j, k) schedule(static)
            do j = 1, n2
                do i = 1, n1
                    do k = -r, r
                        if (dim == 1) then
                            if (i + k < 1 .or. i + k > n1) cycle
                            ta(i, j) = ta(i, j) + w(k + r + 1)*a(i + k, j)
                            tb2(i, j) = tb2(i, j) + w(k + r + 1)*b(i + k, j)
                        else
                            if (j + k < 1 .or. j + k > n2) cycle
                            ta(i, j) = ta(i, j) + w(k + r + 1)*a(i, j + k)
                            tb2(i, j) = tb2(i, j) + w(k + r + 1)*b(i, j + k)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
            a = ta;  b = tb2
        end subroutine gauss_1d
    end subroutine smooth_positive

    subroutine taub_set_scale(tb, ids, grp, time, log_factor, src)
        ! The multipliers of "domains" (set before): the scale group `grp` of the domain ids `ids` (0 =
        ! unscaled) and log_factor(group, node) at model times `time` [yr] (any order, distinct). An id of
        ! tb without an entry is unscaled; an id with a group > 0 must be a domain of tb.
        type(taub_t),     intent(inout) :: tb
        integer,          intent(in)    :: ids(:), grp(:)
        real(dp),         intent(in)    :: time(:), log_factor(:,:)
        character(len=*), intent(in)    :: src

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_scale"
        integer, allocatable :: p(:)
        integer :: k, q, ng

        if (tb%method /= TAUB_M_DOMAINS) call fesm_error(proc, "a tau_b scale needs taub_method = 'domains'", trim(src))
        ng = size(log_factor, 1)
        if (size(ids) /= size(grp)) call fesm_error(proc, "one group per domain id is needed", trim(src))
        if (ng == 0 .or. size(time) == 0 .or. size(log_factor, 2) /= size(time)) then
            call fesm_error(proc, "log_factor is not (group, node) with at least one of each", trim(src)// &
                            ": "//to_str(size(log_factor, 1))//" x "//to_str(size(log_factor, 2))//", "// &
                            to_str(size(time))//" nodes")
        end if
        if (any(grp < 0 .or. grp > ng)) call fesm_error(proc, "a group is outside 0 (unscaled) to "//to_str(ng), &
                                                       trim(src))
        if (.not. all(abs(log_factor) <= huge(1.0_dp))) call fesm_error(proc, "log_factor is not finite", trim(src))
        p = time_order(time)
        do k = 2, size(p)
            if (.not. time(p(k)) > time(p(k-1))) call fesm_error(proc, "two nodes have the same time", trim(src))
        end do

        allocate(tb%sgrp(size(tb%ids)))
        tb%sgrp = 0
        do k = 1, size(ids)
            if (count(ids == ids(k)) > 1) call fesm_error(proc, "domain id "//to_str(ids(k))//" is listed twice", &
                                                         trim(src))
            q = findloc(tb%ids, ids(k), dim=1)
            if (q > 0) then
                tb%sgrp(q) = grp(k)
            else if (grp(k) > 0) then
                call fesm_error(proc, "domain id "//to_str(ids(k))//" (group "//to_str(grp(k))//") has no tau_b", &
                                trim(src)//"; tau_b "//trim(tb%src))
            end if
        end do
        tb%stime = time(p)
        tb%slf   = log_factor(:, p)
        tb%ssrc  = src
    end subroutine taub_set_scale

    subroutine taub_read_scale(tb, file)
        ! The multiplier file: group(domain) the scale group of each id (0 = unscaled), domain_id(domain),
        ! time(node) [its units: "years ago", "ka BP", ...], log_factor(scale_group, node).
        type(taub_t),     intent(inout) :: tb
        character(len=*), intent(in)    :: file

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_scale"
        character(len=64), allocatable :: names(:)
        integer,  allocatable :: dims(:), ids(:), grp(:)
        real(dp), allocatable :: t(:), lf(:,:)
        character(len=256) :: units
        integer :: nd
        logical :: exists

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the tau_b scale file does not exist", "taub_scale_file = "//trim(file))
        call nc_dims(file, "domain_id", names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`domain_id` is not a list of ids (1-D)", trim(file))
        nd = dims(1)
        call nc_dims(file, "group", names, dims)
        if (size(dims) /= 1 .or. dims(1) /= nd) call fesm_error(proc, "`group` is not one group per domain_id", trim(file))
        call nc_dims(file, "time", names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`time` is not a list of nodes (1-D)", trim(file))
        allocate(t(dims(1)))
        call nc_dims(file, "log_factor", names, dims)
        if (size(dims) /= 2) call fesm_error(proc, "`log_factor` is not (scale_group, node)", trim(file))
        if (dims(2) /= size(t)) call fesm_error(proc, "`log_factor` is not (scale_group, node): "//to_str(dims(2))// &
                                               " nodes, `time` "//to_str(size(t)), trim(file))
        allocate(ids(nd), grp(nd), lf(dims(1), dims(2)))
        call nc_read(file, "domain_id", ids)
        call nc_read(file, "group", grp)
        call nc_read(file, "time", t)
        units = ""
        call nc_read_attr(file, "time", "units", units)
        call nc_read(file, "log_factor", lf)
        call taub_set_scale(tb, ids, grp, model_time_axis(t, units), lf, trim(file))
    end subroutine taub_read_scale

    subroutine taub_set_domain_factor(tb, ids, log_factor, src)
        ! The factor per domain of "domains" (set before, with its override table): exp(log_factor) of the
        ! domain ids `ids`, constant in time. An id of tb without an entry keeps factor 1; an id of the list
        ! that tb has no tau_b for must have log_factor 0.
        type(taub_t),     intent(inout) :: tb
        integer,          intent(in)    :: ids(:)
        real(dp),         intent(in)    :: log_factor(:)
        character(len=*), intent(in)    :: src

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_set_domain_factor"
        integer :: k, q

        if (tb%method /= TAUB_M_DOMAINS) call fesm_error(proc, "a domain factor needs taub_method = 'domains'", trim(src))
        if (size(ids) /= size(log_factor)) call fesm_error(proc, "one log_factor per domain id is needed", trim(src))
        if (.not. all(abs(log_factor) <= huge(1.0_dp))) call fesm_error(proc, "log_factor is not finite", trim(src))
        if (allocated(tb%dfac)) deallocate(tb%dfac)
        allocate(tb%dfac(size(tb%ids)))
        tb%dfac = 1.0_dp
        do k = 1, size(ids)
            if (count(ids == ids(k)) > 1) call fesm_error(proc, "domain id "//to_str(ids(k))//" is listed twice", &
                                                         trim(src))
            q = findloc(tb%ids, ids(k), dim=1)
            if (q > 0) then
                tb%dfac(q) = exp(log_factor(k))
            else if (log_factor(k) /= 0.0_dp) then
                call fesm_error(proc, "domain id "//to_str(ids(k))//" (log_factor "//to_str(log_factor(k))// &
                                ") has no tau_b", trim(src)//"; tau_b "//trim(tb%src))
            end if
        end do
        tb%dsrc = src
    end subroutine taub_set_domain_factor

    subroutine taub_read_domain_factor(tb, file)
        ! The domain factor file: domain_id(domain), log_factor(domain) (natural log of the factor).
        type(taub_t),     intent(inout) :: tb
        character(len=*), intent(in)    :: file

        character(len=*), parameter :: proc = "kryos_plastic_taub::taub_read_domain_factor"
        character(len=64), allocatable :: names(:)
        integer,  allocatable :: dims(:), ids(:)
        real(dp), allocatable :: lf(:)
        integer :: nd
        logical :: exists

        inquire(file=trim(file), exist=exists)
        if (.not. exists) call fesm_error(proc, "the tau_b domain factor file does not exist", &
                                          "taub_domain_factor_file = "//trim(file))
        call nc_dims(file, "domain_id", names, dims)
        if (size(dims) /= 1) call fesm_error(proc, "`domain_id` is not a list of ids (1-D)", trim(file))
        nd = dims(1)
        call nc_dims(file, "log_factor", names, dims)
        if (size(dims) /= 1 .or. dims(1) /= nd) call fesm_error(proc, "`log_factor` is not one value per domain_id", &
                                                               trim(file))
        allocate(ids(nd), lf(nd))
        call nc_read(file, "domain_id", ids)
        call nc_read(file, "log_factor", lf)
        call taub_set_domain_factor(tb, ids, lf, trim(file))
    end subroutine taub_read_domain_factor

    function taub_scale_at(tb, t) result(f)
        ! log_factor of each group at model time t [yr]: linear between the nodes, held outside.
        type(taub_t), intent(in) :: tb
        real(dp),     intent(in) :: t
        real(dp) :: f(size(tb%slf, 1))

        integer  :: k, n
        real(dp) :: a

        n = size(tb%stime)
        if (t <= tb%stime(1)) then
            f = tb%slf(:, 1)
        else if (t >= tb%stime(n)) then
            f = tb%slf(:, n)
        else
            k = 1
            do while (tb%stime(k+1) <= t)
                k = k + 1
            end do
            a = (t - tb%stime(k)) / (tb%stime(k+1) - tb%stime(k))
            f = (1.0_dp - a)*tb%slf(:, k) + a*tb%slf(:, k+1)
        end if
    end function taub_scale_at

    subroutine scale_values(tb, t, v)
        ! v(id) times exp(log_factor) of its group at t (group 0: unchanged).
        type(taub_t), intent(in)    :: tb
        real(dp),     intent(in)    :: t
        real(dp),     intent(inout) :: v(:)

        real(dp), allocatable :: m(:)
        integer :: k

        m = exp(taub_scale_at(tb, t))
        do k = 1, size(v)
            if (tb%sgrp(k) > 0) v(k) = v(k)*m(tb%sgrp(k))
        end do
    end subroutine scale_values

    function values_at(tb, t) result(v)
        ! The value of each id at t: static, or the series linear in time.
        type(taub_t), intent(in) :: tb
        real(dp),     intent(in) :: t
        real(dp), allocatable :: v(:)

        integer  :: k, n
        real(dp) :: a

        if (.not. allocated(tb%val_t)) then
            v = tb%val
            return
        end if
        n = size(tb%time)
        if (t < tb%time(1) .or. t > tb%time(n)) then
            if (.not. tb%hold) then
                call fesm_error("kryos_plastic_taub::taub_at", "t = "//to_str(t)//" yr is outside the time range "// &
                                "of the tau_b series", trim(tb%src)//": ["//to_str(tb%time(1))//", "// &
                                to_str(tb%time(n))//"] yr; margin_extrapolate = 'hold' holds the end values")
            end if
            v = tb%val_t(:, merge(1, n, t < tb%time(1)))
            return
        end if
        k = 1
        do while (k < n)
            if (tb%time(k+1) > t) exit
            k = k + 1
        end do
        if (k == n) then
            v = tb%val_t(:, n)
        else
            a = (t - tb%time(k)) / (tb%time(k+1) - tb%time(k))
            v = (1.0_dp - a)*tb%val_t(:, k) + a*tb%val_t(:, k+1)
        end if
    end function values_at

    subroutine index_cells(tb)
        ! idx: the position in ids of each cell's domain; a domain of the map without a value stops.
        type(taub_t), intent(inout) :: tb
        integer :: i, j, p

        if (allocated(tb%idx)) deallocate(tb%idx)
        allocate(tb%idx(size(tb%dom, 1), size(tb%dom, 2)))
        tb%idx = 0
        do j = 1, size(tb%dom, 2)
            do i = 1, size(tb%dom, 1)
                if (tb%dom(i,j) <= 0) cycle
                p = findloc(tb%ids, tb%dom(i,j), dim=1)
                if (p == 0) then
                    call fesm_error("kryos_plastic_taub::taub_set_domains", "domain "//to_str(tb%dom(i,j))// &
                                    " of the map has no tau_b", trim(tb%src)//" (cell "//to_str(i)//", "//to_str(j)// &
                                    "); give it one in the override table (taub_override)")
                end if
                tb%idx(i,j) = p
            end do
        end do
    end subroutine index_cells

    subroutine check_xy(file, var, nx, ny)
        ! `var` of `file` is a field (x, y) on the instance grid.
        character(len=*), intent(in) :: file, var
        integer,          intent(in) :: nx, ny
        character(len=64), allocatable :: names(:)
        integer, allocatable :: dims(:)

        call nc_dims(file, var, names, dims)
        if (size(dims) /= 2) then
            call fesm_error("kryos_plastic_taub::read", "`"//trim(var)//"` is not a static field (x, y)", &
                            "file "//trim(file)//", "//to_str(size(dims))//" dimensions")
        end if
        if (dims(1) /= nx .or. dims(2) /= ny) then
            call fesm_error("kryos_plastic_taub::read", "`"//trim(var)//"` is not on the instance grid", "file "// &
                            trim(file)//": "//to_str(dims(1))//" x "//to_str(dims(2))//", the grid "//to_str(nx)// &
                            " x "//to_str(ny))
        end if
    end subroutine check_xy

end module kryos_plastic_taub
