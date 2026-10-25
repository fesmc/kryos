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

    use precision,       only: sp, dp, wp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: fields_get, fields_put
    use nml,             only: nml_read, nml_validate, nml_has_param
    use ncio,            only: nc_create, nc_write_dim, nc_write, nc_open, nc_close, nc_time_index
    use tsgen,           only: tsgen_class, tsgen_init, tsgen_update, tsgen_restart_write, tsgen_restart_read
    use kryos_component, only: init_ctx_t, step_time_t, constraints_t, NAME_LEN, PATH_LEN, PASS_INIT
    use kryos_forcing,   only: forcing_t

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
        character(len=32), allocatable :: from(:)       ! derived: its producer kind (air_anomaly, obm_*)
    contains
        procedure :: init          => tsgen_init_
        procedure :: init_state    => tsgen_init_state_
        procedure :: update        => tsgen_update_
        procedure :: restart_write => tsgen_restart_write_
        procedure :: restart_read  => tsgen_restart_read_
        procedure :: restart_files => tsgen_restart_files_
        procedure :: write_output  => tsgen_write_output_
        procedure :: kill          => tsgen_kill_
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
        real(wp) :: every_df
        integer  :: k, n

        names = ""
        if (nml_has_param(ctx%par_file, ctx%instance, "coefs")) &
            call nml_read(ctx%par_file, ctx%instance, "coefs", names, init=.true.)
        n = count(len_trim(names) > 0)

        self%prec = sp
        allocate(self%coef_names(n), self%coef_units(n), self%coef_derived(n))
        allocate(self%tsg(n), self%f(n), self%from(n))
        self%coef_derived = .false.
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
            call nml_read(trim(gfile), trim(group), "restart_every_df", every_df, &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            if (every_df > 0.0_wp) then
                call fesm_error("kryos_forcing_tsgen::init", "["//trim(group)//"] restart_every_df > 0 is not supported", &
                                "instance "//trim(self%name)//": the forcing-increment restart bundles come with phase 3 "// &
                                "(core-design 7); only the kill switch is kept")
            end if

            if (trim(method) == "derived") then
                self%coef_derived(k) = .true.
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
    end subroutine tsgen_init_

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
        end do
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
