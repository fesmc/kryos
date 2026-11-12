module kryos_forcing_tsgen
    ! forcing backend tsgen (adapters-design section 8, contracts-design 2.7
    ! and 7): the coefficient registry of yelmox kryos_forcing.f90 (coefs_init,
    ! coefs_update, coefs_kill, coefs_restart_write/read, coefs_write_step)
    ! on fesm-utils tsgen, minus its domain coupling; no stub (fesm-utils is
    ! always linked). Global, one instance per run.
    !
    ! A coefficient is a scalar of the forcing time, named in the instance key
    ! `coefs`, defined by a group [coef_<name>] (the par file, else the first
    ! library file that has it; unset keys from the template `coef` of
    ! input/forcing_defaults.nml). tsgen coefficients are advanced once per
    ! step on the forcing time (st%t_eval = time_rel) with the response
    ! variable M_ice [Gt]; derived coefficients (method = "derived") are set
    ! by a producer (the climate) through forcing_set, copied here from their
    ! import at each update. Every coefficient is exported INSTANT at t1, in
    ! sp widened (yelmox coef_class%f is sp), so the exchange is exact (P1).
    !
    ! The orbit (instance key `orbit`, a group [orbit] read by kryos_orbit): the
    ! orbital elements at the forcing time (or fixed), exported as orbit_ecc,
    ! orbit_obl, orbit_per for every component that takes the orbit (the
    ! climate, the surface's insolation), so the run has one orbit. The stream
    ! `insol` (out_insol, <instance>_insol.nc) holds them with the monthly mean
    ! top-of-atmosphere insolation on 1° latitudes, a diagnostic.

    use precision,       only: sp, dp, wp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_validate, nml_has_param
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_read, nc_open, nc_close, nc_time_index
    use tsgen,           only: tsgen_class, tsgen_init, tsgen_update, tsgen_restart_write, tsgen_restart_read
    use insolation,      only: insol_class, insol_init, calc_insol_ave, insol_wp => wp
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, NAME_LEN, PATH_LEN, PASS_INIT
    use kryos_forcing,   only: forcing_t
    use kryos_orbit,     only: orbit_ctl_class, orbit_pars_class, orbit_ctl_load, orbit_pars_now

    implicit none

    private
    public :: tsgen_forcing_t, tsgen_forcing_constraints

    integer, parameter :: COEF_MAX = 32            ! max entries of `coefs`

    ! Schema of the [coef_<name>] groups (the forcing component's, shared with
    ! the climate backends that name coefficients).
    character(len=*), parameter :: COEF_DEFAULTS = "input/forcing_defaults.nml"

    type, extends(forcing_t) :: tsgen_forcing_t
        type(tsgen_class), allocatable :: tsg(:)        ! one generator per coefficient (unused for derived)
        real(wp),          allocatable :: f(:)          ! the current value of each coefficient [sp]
        real(wp),          allocatable :: every_df(:)   ! restart_every_df of each (0: off)
        real(wp),          allocatable :: f_mark(:)     ! its value at the last forcing-increment mark
        character(len=32), allocatable :: from(:)       ! derived: its producer kind (air_anomaly, obm_*)
        character(len=NAME_LEN) :: orbit_group = ""      ! the orbit group ("" none)
        type(orbit_ctl_class)   :: orb_ctl
        type(insol_class)       :: orb_is                ! the orbital tables (orbit_pars_now)
        type(orbit_pars_class)  :: orb                   ! the elements of the last update
        type(insol_class)       :: diag_is               ! 1° latitudes without tables (the stream insol)
        integer                 :: day_year = 365        ! the run's calendar (the stream insol)
    contains
        procedure :: init          => tsgen_init_
        procedure :: init_state    => tsgen_init_state_
        procedure :: update        => tsgen_update_
        procedure :: restart_write => tsgen_restart_write_
        procedure :: restart_read  => tsgen_restart_read_
        procedure :: restart_files => tsgen_restart_files_
        procedure :: write_output  => tsgen_write_output_
        procedure :: kill          => tsgen_kill_
        procedure :: restart_due   => tsgen_restart_due_
    end type tsgen_forcing_t

contains

    function tsgen_forcing_constraints() result(con)
        ! Global only, one instance per run.
        type(constraints_t) :: con
        con = constraints_t(regional=.false., global=.true., max_instances=1)
    end function tsgen_forcing_constraints

    ! =========================================================================
    ! Bindings
    ! =========================================================================

    subroutine tsgen_init_(self, ctx)
        ! coefs_init (kryos_forcing.f90:76-160) at the forcing time ctx%time_rel.
        class(tsgen_forcing_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx

        character(len=NAME_LEN) :: names(COEF_MAX)
        character(len=512) :: gfile
        character(len=128) :: group
        character(len=32)  :: method
        integer  :: k, n

        names = ""
        if (nml_has_param(ctx%par_file, ctx%instance, "coefs")) &
            call nml_read(ctx%par_file, ctx%instance, "coefs", names, init=.true.)
        n = count(len_trim(names) > 0)

        self%prec = sp
        allocate(self%coef_names(n), self%coef_units(n), self%coef_derived(n))
        allocate(self%tsg(n), self%f(n), self%from(n), self%every_df(n), self%f_mark(n))
        self%coef_derived = .false.
        self%every_df = 0.0_wp
        self%f    = 0.0_wp
        self%from = ""

        do k = 1, n
            self%coef_names(k) = trim(names(k))
            if (any(names(1:k-1) == names(k))) then
                call fesm_error("kryos_forcing_tsgen::init", "coefficient '"//trim(names(k))//"' is listed twice", &
                                "instance "//trim(self%name)//", key coefs")
            end if
            group = "coef_"//trim(names(k))
            gfile = coef_file(ctx, trim(group))
            call nml_validate(trim(gfile), COEF_DEFAULTS, trim(group), defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "method", method, init=.true., &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "units", self%coef_units(k), init=.true., &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "restart_every_df", self%every_df(k), &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")

            if (trim(method) == "derived") then
                self%coef_derived(k) = .true.
                if (self%every_df(k) > 0.0_wp) then
                    call fesm_error("kryos_forcing_tsgen::init", "["//trim(group)//"] restart_every_df > 0 on a derived "// &
                                    "coefficient", "instance "//trim(self%name)//": the forcing-increment restarts "// &
                                    "follow tsgen coefficients only")
                end if
                call nml_read(trim(gfile), trim(group), "from", self%from(k), init=.true., &
                              defaults_file=COEF_DEFAULTS, defaults_group="coef")
                select case (trim(self%from(k)))
                case ("air_anomaly", "obm_north", "obm_south")
                case default
                    call fesm_error("kryos_forcing_tsgen::init", "["//trim(group)//"] from = '"//trim(self%from(k))// &
                                    "': derived coefficients come from air_anomaly, obm_north or obm_south", &
                                    "instance "//trim(self%name))
                end select
                self%f(k) = 0.0_wp
                call log_msg(LOG_INFO, "kryos_forcing_tsgen", trim(self%coef_names(k))//" (derived from "// &
                             trim(self%from(k))//"), set each step by the model", "instance "//trim(self%name))
            else
                call tsgen_init(self%tsg(k), trim(gfile), ctx%time_rel, group=trim(group), &
                                defaults_file=COEF_DEFAULTS, defaults_group="coef")
                self%f(k) = self%tsg(k)%f_now
                if (self%tsg(k)%par%feedback .or. self%tsg(k)%par%with_kill) self%opt%M_ice = .true.
                call log_msg(LOG_INFO, "kryos_forcing_tsgen", trim(self%coef_names(k))//" ("// &
                             trim(self%tsg(k)%par%method)//"), f = "//to_str(real(self%f(k), dp)), &
                             "instance "//trim(self%name))
            end if
        end do
        self%f_mark = self%f                           ! coefs_init: the mark of a coefficient is its initial value

        if (nml_has_param(ctx%par_file, ctx%instance, "orbit")) &
            call nml_read(ctx%par_file, ctx%instance, "orbit", self%orbit_group)
        if (len_trim(self%orbit_group) > 0) call orbit_init_(self, ctx)
    end subroutine tsgen_init_

    subroutine orbit_init_(self, ctx)
        ! The orbit group, the tables, the elements at the initial forcing time; the 1° latitude state of the
        ! diagnostic stream on the run's calendar (360, 365 or 366 days).
        class(tsgen_forcing_t), intent(inout) :: self
        type(init_ctx_t),       intent(in)    :: ctx
        real(dp) :: lats(180)
        integer  :: j

        self%opt%orbit = .true.
        call orbit_ctl_load(self%orb_ctl, trim(ctx%par_file), trim(self%orbit_group))
        call insol_init(self%orb_is, fldr=trim(self%orb_ctl%fldr))
        self%orb = orbit_pars_now(self%orb_ctl, self%orb_is, ctx%time_rel)
        self%day_year = nint(ctx%sec_year/86400.0_dp)
        if (self%day_year == 365 .and. ctx%sec_year > 365.5_dp*86400.0_dp) self%day_year = 366
        lats = [(-89.5_dp + real(j - 1, dp), j = 1, 180)]
        call insol_init(self%diag_is, lats, day_year=self%day_year, tables=.false.)
        call log_msg(LOG_INFO, "kryos_forcing_tsgen", "orbit ["//trim(self%orbit_group)//"] "//trim(self%orb_ctl%method)// &
                     ": ecc "//to_str(self%orb%ecc)//", obl "//to_str(self%orb%obl)//" deg, per "//to_str(self%orb%per)// &
                     " deg at the forcing time "//to_str(ctx%time_rel)//" yr", "instance "//trim(self%name))
    end subroutine orbit_init_

    subroutine orbit_put(self, t1)
        class(tsgen_forcing_t), intent(inout) :: self
        real(dp),               intent(in)    :: t1
        call fields_put(self%exports, self%h_orbit(1), self%orb%ecc, t1=t1)
        call fields_put(self%exports, self%h_orbit(2), self%orb%obl, t1=t1)
        call fields_put(self%exports, self%h_orbit(3), self%orb%per, t1=t1)
    end subroutine orbit_put

    subroutine tsgen_init_state_(self, st, pass)
        ! PASS_INIT: every coefficient at the initial forcing time (derived: 0,
        ! kryos_forcing.f90:147); the derived imports start at 0 too, so the
        ! first copy in update finds a value before any producer has set one.
        class(tsgen_forcing_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st
        integer,                intent(in)    :: pass

        integer :: k

        if (pass /= PASS_INIT) return
        do k = 1, size(self%f)
            call fields_put(self%exports, self%h_coef(k), self%f(k), t1=st%t1)
            if (self%coef_derived(k)) call fields_put(self%imports, self%hi_derived(k), self%f(k), t1=st%t1)
        end do
        if (self%opt%orbit) then
            self%orb = orbit_pars_now(self%orb_ctl, self%orb_is, st%t_eval)
            call orbit_put(self, st%t1)
        end if
    end subroutine tsgen_init_state_

    subroutine tsgen_update_(self, st)
        ! coefs_update (kryos_forcing.f90:162-176) on the forcing time
        ! st%t_eval with var = M_ice (yelmox.f90:185): tsgen coefficients
        ! advance, derived ones copy their import (the value last set by the
        ! producer, as yelmox keeps it); all are exported at t1.
        class(tsgen_forcing_t), intent(inout) :: self
        type(step_time_t),      intent(in)    :: st

        real(wp) :: m_ice
        integer  :: k

        m_ice = 0.0_wp
        if (self%opt%M_ice) call fields_get(self%imports, self%hi_M_ice, m_ice)

        do k = 1, size(self%f)
            if (self%coef_derived(k)) then
                call fields_get(self%imports, self%hi_derived(k), self%f(k))
            else
                call tsgen_update(self%tsg(k), st%t_eval, var=m_ice)
                self%f(k) = self%tsg(k)%f_now
            end if
            call fields_put(self%exports, self%h_coef(k), self%f(k), t1=st%t1)
        end do
        if (self%opt%orbit) then
            self%orb = orbit_pars_now(self%orb_ctl, self%orb_is, st%t_eval)
            call orbit_put(self, st%t1)
        end if
    end subroutine tsgen_update_

    subroutine tsgen_restart_write_(self, dir, t)
        ! coefs_restart_write: the state of every tsgen coefficient, one file
        ! each. The derived values are in the export state (the core's fields file).
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t

        integer :: k

        do k = 1, size(self%f)
            if (self%coef_derived(k)) cycle
            call tsgen_restart_write(coef_restart_file(self, dir, k), self%tsg(k), t)
        end do
        if (increments_on(self)) then                ! the bundle counter: a restored run numbers on from it
            call nc_create(counter_file(self, dir))
            call nc_write(counter_file(self, dir), "n_restart", self%n_restart)
        end if
    end subroutine tsgen_restart_write_

    function tsgen_restart_files_(self) result(files)
        ! The files tsgen_restart_write_ creates: <instance>_<coef>_restart.nc
        ! of every coefficient that is not derived, in coefficient order.
        class(tsgen_forcing_t), intent(in) :: self
        character(len=PATH_LEN), allocatable :: files(:)

        integer :: k

        allocate(files(0))
        do k = 1, size(self%f)
            if (self%coef_derived(k)) cycle
            files = [character(len=PATH_LEN) :: files, coef_restart_name(self, k)]
        end do
        if (increments_on(self)) files = [character(len=PATH_LEN) :: files, counter_name(self)]
    end function tsgen_restart_files_

    subroutine tsgen_restart_read_(self, dir, t)
        ! coefs_restart_read: every tsgen coefficient from its file (after init).
        ! A missing file stops (no cold-start fallback, adapters-design 1.4).
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t

        character(len=:), allocatable :: file
        logical :: found
        integer :: k

        associate (unused => t); end associate
        do k = 1, size(self%f)
            if (self%coef_derived(k)) cycle
            file = coef_restart_file(self, dir, k)
            inquire(file=file, exist=found)
            if (.not. found) then
                call fesm_error("kryos_forcing_tsgen::restart_read", "restart file of a coefficient not found", &
                                "instance "//trim(self%name)//", coefficient "//trim(self%coef_names(k))// &
                                ": "//file)
            end if
            call tsgen_restart_read(self%tsg(k), file)
            self%f(k) = self%tsg(k)%f_now
            self%f_mark(k) = self%f(k)                 ! coefs_restart_read: the mark restarts at the restored value
        end do
        ! the bundle counter; a bundle of a run without restart_every_df has none: 0
        self%n_restart = 0
        file = counter_file(self, dir)
        inquire(file=file, exist=found)
        if (found) call nc_read(file, "n_restart", self%n_restart)
    end subroutine tsgen_restart_read_

    subroutine tsgen_write_output_(self, file, t, stream, init)
        ! Stream 1D: <instance>_ts.nc, coefs_write_step (kryos_forcing.f90:390-427):
        ! coef_<name> in its units and, for tsgen coefficients, the rate
        ! coef_<name>_df_dt and, with a response history (feedback methods,
        ! kill switch), the windowed response rate coef_<name>_dv_dt.
        ! A derived coefficient is written as last set (its import), as the
        ! registry value of yelmox.
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*),       intent(in)    :: file, stream
        real(dp),               intent(in)    :: t
        logical,                intent(in)    :: init

        integer :: k, ncid, n
        real(wp) :: f
        character(len=72) :: vnm

        if (trim(stream) == "insol") then
            call insol_write(self, file, t, init)
            return
        end if
        if (trim(stream) /= "1D") return
        if (init) then
            call nc_create(file)
            call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
        end if
        if (size(self%f) == 0) return

        call nc_open(file, ncid, writable=.true.)
        n = nc_time_index(file, "time", t, ncid)
        call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)

        do k = 1, size(self%f)
            vnm = "coef_"//trim(self%coef_names(k))
            f = self%f(k)
            if (self%coef_derived(k)) then
                if (self%imports%f(self%hi_derived(k))%published) &
                    call fields_get(self%imports, self%hi_derived(k), f)
            end if
            call nc_write(file, trim(vnm), f, dim1="time", start=[n], count=[1], ncid=ncid, &
                          units=trim(self%coef_units(k)), long_name="Coefficient "//trim(self%coef_names(k)))
            if (self%coef_derived(k)) cycle
            associate (tsg => self%tsg(k))
            call nc_write(file, trim(vnm)//"_df_dt", tsg%df_dt, dim1="time", start=[n], count=[1], ncid=ncid, &
                          units=trim(self%coef_units(k))//"/yr", &
                          long_name="Rate of change of coefficient "//trim(self%coef_names(k)))
            if (allocated(tsg%dv_dt)) &
                call nc_write(file, trim(vnm)//"_dv_dt", tsg%dv_dt_ave, dim1="time", start=[n], count=[1], ncid=ncid, &
                              units="Gt/yr", long_name="Windowed rate of the response (ice volume) of coefficient "// &
                                                       trim(self%coef_names(k)))
            end associate
        end do

        call nc_close(ncid)
    end subroutine tsgen_write_output_

    subroutine insol_write(self, file, t, init)
        ! Stream insol: the orbital elements and the monthly mean top-of-atmosphere insolation [W m-2] on 1°
        ! latitudes (a diagnostic; it depends on latitude only), months of the run's calendar.
        class(tsgen_forcing_t), intent(inout) :: self
        character(len=*),       intent(in)    :: file
        real(dp),               intent(in)    :: t
        logical,                intent(in)    :: init

        integer, parameter :: ML365(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        integer  :: ml(12), d0, m, ncid, n
        real(insol_wp) :: S(180), Sm(180, 12)
        real(dp) :: lats(180)

        if (.not. self%opt%orbit) return
        select case (self%day_year)
        case (360);  ml = 30
        case (366);  ml = ML365;  ml(2) = 29
        case default; ml = ML365
        end select
        d0 = 0
        do m = 1, 12
            call calc_insol_ave(self%diag_is, S, [d0 + 1, d0 + ml(m)], self%orb)
            Sm(:, m) = S
            d0 = d0 + ml(m)
        end do
        if (init) then
            lats = [(-89.5_dp + real(m - 1, dp), m = 1, 180)]
            call nc_create(file)
            call nc_write_dim(file, "lat", x=lats, units="degrees_north")
            call nc_write_dim(file, "month", x=1.0_dp, dx=1.0_dp, nx=12, units="month")
            call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="year", unlimited=.true.)
        end if
        call nc_open(file, ncid, writable=.true.)
        n = nc_time_index(file, "time", t, ncid)
        call nc_write(file, "time", t, dim1="time", start=[n], count=[1], ncid=ncid)
        call nc_write(file, "orbit_ecc", self%orb%ecc, dim1="time", start=[n], count=[1], ncid=ncid, units="1", &
                      long_name="eccentricity")
        call nc_write(file, "orbit_obl", self%orb%obl, dim1="time", start=[n], count=[1], ncid=ncid, units="degrees", &
                      long_name="obliquity")
        call nc_write(file, "orbit_per", self%orb%per, dim1="time", start=[n], count=[1], ncid=ncid, units="degrees", &
                      long_name="longitude of perihelion (Berger 1978)")
        call nc_write(file, "S", Sm, dim1="lat", dim2="month", dim3="time", start=[1, 1, n], count=[180, 12, 1], &
                      ncid=ncid, units="W m-2", long_name="monthly mean top-of-atmosphere insolation")
        call nc_close(ncid)
    end subroutine insol_write

    logical function tsgen_kill_(self) result(kill)
        ! coefs_kill: the tsgen kill switch of any coefficient (the response
        ! equilibrated at a forcing bound). The run stops at the end of the step.
        class(tsgen_forcing_t), intent(in) :: self
        integer :: k
        kill = .false.
        do k = 1, size(self%f)
            if (self%coef_derived(k)) cycle
            kill = kill .or. self%tsg(k)%kill
        end do
    end function tsgen_kill_

    logical function tsgen_restart_due_(self) result(due)
        ! coefs_restart_due: a coefficient has moved by more than its restart_every_df
        ! (either direction) since its own mark; its mark moves to its value. Called
        ! after the step, so f is the value at t1.
        class(tsgen_forcing_t), intent(inout) :: self
        integer :: k
        due = .false.
        do k = 1, size(self%f)
            if (self%every_df(k) <= 0.0_wp) cycle
            if (abs(self%f(k) - self%f_mark(k)) > self%every_df(k)) then
                self%f_mark(k) = self%f(k)
                due = .true.
            end if
        end do
    end function tsgen_restart_due_

    ! =========================================================================
    ! helpers
    ! =========================================================================

    function coef_file(ctx, group) result(gfile)
        ! The file holding [coef_<name>]: the par file, else the first library file.
        use nml, only: nml_find_group
        type(init_ctx_t), intent(in) :: ctx
        character(len=*), intent(in) :: group
        character(len=len(ctx%par_file)) :: gfile
        character(len=len(ctx%par_file)) :: files(1 + size(ctx%library))
        files(1)  = ctx%par_file
        files(2:) = ctx%library
        gfile = nml_find_group(files, group, "&"//trim(ctx%instance)//" coefs")
    end function coef_file

    logical function increments_on(self)
        ! Some coefficient carries restart_every_df: the bundles hold the counter.
        class(tsgen_forcing_t), intent(in) :: self
        increments_on = any(self%every_df > 0.0_wp)
    end function increments_on

    function counter_name(self) result(name)
        class(tsgen_forcing_t), intent(in) :: self
        character(len=:), allocatable :: name
        name = trim(self%name)//"_increment_restart.nc"
    end function counter_name

    function counter_file(self, dir) result(file)
        class(tsgen_forcing_t), intent(in) :: self
        character(len=*),       intent(in) :: dir
        character(len=:), allocatable :: file
        file = trim(dir)//"/"//counter_name(self)
    end function counter_file

    function coef_restart_name(self, k) result(name)
        class(tsgen_forcing_t), intent(in) :: self
        integer,                intent(in) :: k
        character(len=:), allocatable :: name
        name = trim(self%name)//"_"//trim(self%coef_names(k))//"_restart.nc"
    end function coef_restart_name

    function coef_restart_file(self, dir, k) result(file)
        class(tsgen_forcing_t), intent(in) :: self
        character(len=*),       intent(in) :: dir
        integer,                intent(in) :: k
        character(len=:), allocatable :: file
        file = trim(dir)//"/"//coef_restart_name(self, k)
    end function coef_restart_file

end module kryos_forcing_tsgen
