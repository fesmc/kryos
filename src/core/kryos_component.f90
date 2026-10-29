module kryos_component
    ! The component contract (design section 4, contracts-design section 3.1):
    ! the abstract component_t every instance extends, the plain types the core
    ! hands to the wrappers (init_ctx_t, step_time_t), and the free procedures
    ! the component modules share (backend_admit, component_close, kryos_mask).
    !
    ! Style rule (design decision 1): the deferred bindings exist only here and
    ! on the per-component abstract types; they are called from one place each,
    ! the free-procedure wrappers of the component modules (ice_update(ice, st)
    ! wraps ice%c%update(st)). Driver, core and experiments never write %update.

    use, intrinsic :: iso_fortran_env, only: int64
    use precision,      only: sp, dp
    use coordinates,    only: grid_class
    use regions,        only: regions_class
    use phys_constants, only: phys_const_class
    use fesm_log,       only: LOG_WARN, log_msg, log_set_prefix, fesm_error, to_str
    use fesm_fields,    only: field_state_t, fields_finalize_interval, fields_get, units_equivalent, &
                              FIELD_INSTANT, FIELD_NAME_LEN, FIELD_UNITS_LEN
    use fesm_clock,     only: clock_t, cadence_t, cadence_interval, cadence_dt_years, clock_time, &
                              clock_time_rel, clock_dtt
    use ncio,           only: nc_close, nc_write
    use variable_io,    only: var_io_type, var_io_write
    use kryos_outfile,  only: outfile_row, outfile_table_file, outfile_dims_2d, outfile_dims_ts, &
                              outfile_axis, outfile_open

    implicit none

    private

    ! Initialisation passes of init_state (core-design section 3.2; PASS_RESTART
    ! after a restore, review outcome). PASS_INIT: right after every instance's
    ! init, cold start and restore alike, put the exports of the state as
    ! initialised at t0 (what other instances read before the cold start).
    ! PASS_SPINUP: after PASS_STATE, with the experiment hook between them: the
    ! cold-start equilibration of the instance and the exports of its initial
    ! state (yelmo: init_equil_time, init_time_thrm).
    ! PASS_RESTART: after restart_read and the core's restore of the imports
    ! and exports (core-design 3.3): re-derive the backend's internal state
    ! from them; put no export (the bundle's are the run's, bit for bit).
    ! The values are not stored anywhere.
    integer, parameter, public :: PASS_INIT = 0, PASS_GEOMETRY = 1, PASS_BOUNDARY = 2, PASS_STATE = 3, &
                                  PASS_SPINUP = 4, PASS_FINAL = 5, PASS_RESTART = 6

    ! The canonical four-value mask (contracts-design section 8).
    integer, parameter, public :: MASK_OCEAN = 0, MASK_LAND = 1, MASK_GRND = 2, MASK_FLOAT = 3

    integer, parameter, public :: NAME_LEN = 64, PATH_LEN = 512, NLIB = 8

    type, public :: caps_t                      ! capability flags, set by the backend in init
        logical :: has_ocn = .false., has_ocn_shelf = .false., has_smb = .false., has_Qd = .false.
        logical :: has_tsl_ann = .false.        ! climate: the annual sea-level air temperature tsl_ann (surface smb_simple)
        logical :: has_flotation = .false.      ! ice: the flotation binding (contracts-design 2.1)
    end type caps_t

    type, public :: constraints_t               ! decision 12; one function per backend module
        logical           :: linked = .true.    ! .false. in the stub
        character(len=16) :: build_flag = ""    ! "yelmo=1"
        logical           :: regional = .true., global = .false.   ! allowed scopes
        integer           :: max_instances = 0  ! 0: unlimited
        character(len=NAME_LEN) :: fixed_grid = "", domain = ""    ! "": any
        logical           :: restart = .true.   ! .false.: no backend restart file
    end type constraints_t

    type, public :: axis_t                      ! the third axis of rank-3 fields
        character(len=FIELD_NAME_LEN)  :: name  = ""   ! "month" | "depth"
        character(len=FIELD_UNITS_LEN) :: units = ""
        real(dp), allocatable :: v(:)
    end type axis_t

    type, public :: census_t                    ! admitted instances per "kind:backend:scope" (run)
        character(len=3*NAME_LEN), allocatable :: key(:)
        integer,                   allocatable :: n(:)
        character(len=NAME_LEN),   allocatable :: first(:)   ! instance admitted first
    end type census_t

    type, public :: grid_ref_t                  ! decision 4: reference topography on one grid
        character(len=NAME_LEN) :: grid = ""
        real(dp), allocatable   :: z_bed(:,:), H_ice(:,:), z_srf(:,:), z_bed_sd(:,:)
    end type grid_ref_t

    type, public :: init_ctx_t                  ! one per instance, built by the core, freed after init
        character(len=PATH_LEN) :: par_file = "", library(NLIB) = ""
        character(len=NAME_LEN) :: instance = "", component = "", backend = "", domain = ""
        character(len=NAME_LEN) :: grid_name = "", scope = "regional"
        character(len=PATH_LEN) :: map_fldr = "maps"   ! [kryos] map_fldr: grid descriptions and maps (a global instance writes its own grid there)
        type(grid_class)        :: grid         ! from maps/grid_<name>.txt
        type(regions_class)     :: reg          ! domain regions on this grid
        type(grid_ref_t)        :: ref          ! reference topography on this grid
        type(phys_const_class)  :: cnst
        real(dp)                :: time = 0.0_dp, time_rel = 0.0_dp, dtt = 0.0_dp   ! t0
        real(dp)                :: sec_year = 0.0_dp   ! [s] year of the run calendar (fesm_clock year)
        type(axis_t), allocatable :: axes(:)    ! axes declared by the producers initialised so far
    end type init_ctx_t

    type, public :: step_time_t
        real(dp) :: t0 = 0.0_dp, t1 = 0.0_dp    ! model-time interval; t0 = t1 = time_init during init
        real(dp) :: t_eval = 0.0_dp             ! t1, or the forcing time at t1 (clock = "forcing")
        real(dp) :: dt = 0.0_dp                 ! the instance's cadence step [yr]
        real(dp) :: dtt = 0.0_dp                ! the run's step [yr]
        integer(int64) :: tick = -1             ! the clock's tick at t1 (step_time_of); a backend compares ticks, not years
    end type step_time_t

    type, abstract, public :: component_t
        character(len=NAME_LEN) :: name = "", kind = "", backend = "", grid = "", domain = ""
        character(len=16)       :: scope = "regional"   ! "regional" | "global"
        integer                 :: prec = sp            ! kind of the backend's arrays (section 9)
        integer                 :: nx = 0, ny = 0
        type(cadence_t)         :: cadence
        type(caps_t)            :: caps
        type(constraints_t)     :: con
        type(field_state_t)     :: imports, exports
        type(axis_t), allocatable :: axes(:)            ! axes of the instance's rank-3 fields
        character(len=FIELD_NAME_LEN), allocatable :: opt_imports(:)   ! imports that may stay unwired
        character(len=FIELD_NAME_LEN), allocatable :: core_imports(:)  ! imports the core fills outside the wiring rows
        real(dp), allocatable :: xkm(:), ykm(:)         ! grid axes [km], from the init context (component_identity)
        type(grid_class), allocatable :: grid_obj       ! a global instance's own grid (its init sets it with grid, nx, ny);
                                                        ! the core registers it in the couplers (kryos_globe, kryos_run)
        character(len=FIELD_NAME_LEN), allocatable :: out_vars(:)   ! exports the default writer writes; unset: all
    contains
        procedure(c_init),   deferred :: init           ! (ctx): read group; set prec, nx, ny, caps, opt
        procedure                     :: init_state => component_init_state      ! (st, pass)
        procedure                     :: prepare => component_prepare            ! (st): imports -> backend boundary
        procedure(c_update), deferred :: update         ! (st): integrate t0 -> t1; put every export
        procedure(c_rst),    deferred :: restart_write, restart_read             ! (dir, t)
        procedure                     :: restart_files => component_restart_files ! () -> base names restart_write creates
        procedure                     :: write_output => component_write_output  ! (file, t, stream, init)
        procedure                     :: end => component_end
    end type component_t

    abstract interface
        subroutine c_init(self, ctx)
            import :: component_t, init_ctx_t
            class(component_t), intent(inout) :: self
            type(init_ctx_t),   intent(in)    :: ctx
        end subroutine c_init

        subroutine c_update(self, st)
            import :: component_t, step_time_t
            class(component_t), intent(inout) :: self
            type(step_time_t),  intent(in)    :: st
        end subroutine c_update

        subroutine c_rst(self, dir, t)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: dir
            real(dp),           intent(in)    :: t
        end subroutine c_rst
    end interface

    public :: step_time_of, sec_year_check
    public :: backend_admit, component_identity, component_close, component_enter, component_leave
    public :: component_axis, component_optional, component_core_filled, month_axis, axis_find, axes_equal
    public :: census_count, kryos_mask
    public :: component_init_state, component_restart_files, component_write_output, component_end, component_output
    public :: component_prepare, component_prepare_call

contains

    ! =========================================================================
    ! Default bindings
    ! =========================================================================

    subroutine component_init_state(self, st, pass)
        ! Default: nothing to do in any pass.
        class(component_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        integer,            intent(in)    :: pass
        associate (unused => st%t1 + pass); end associate
        associate (unused => self%nx); end associate
    end subroutine component_init_state

    subroutine component_prepare(self, st)
        ! Default: nothing to prepare (the backend reads its imports in update).
        ! A backend that converts its imports into a boundary of its own (yelmo:
        ! bnd) does it here, so that an experiment hook can change that
        ! boundary before update.
        class(component_t), intent(inout) :: self
        type(step_time_t),  intent(in)    :: st
        associate (unused => st%t1); end associate
        associate (unused => self%nx); end associate
    end subroutine component_prepare

    function component_restart_files(self) result(files)
        ! The base names (no folder) of the files restart_write creates in the
        ! bundle folder it is given; kryos_restart lists them in the manifest
        ! and checks they exist. Default: the one file <name>_restart.nc of a
        ! backend with a restart file (constraints restart), none otherwise.
        ! A backend whose files differ overrides it, in step with its
        ! restart_write.
        class(component_t), intent(in) :: self
        character(len=PATH_LEN), allocatable :: files(:)

        if (self%con%restart) then
            files = [character(len=PATH_LEN) :: trim(self%name)//"_restart.nc"]
        else
            allocate(files(0))
        end if
    end function component_restart_files

    subroutine component_write_output(self, file, t, stream, init)
        ! Default writer (core-design section 6): the exports listed in out_vars
        ! (unset: all), unpublished ones skipped. Stream 2D: the rank-2 and rank-3
        ! exports on the native grid (xc, yc [km], time; rank 3 also its axis),
        ! file <instance>.nc. Stream 1D: the 0-D exports, file <instance>_ts.nc.
        ! Values are written as float with the precision of the component
        ! kind's table input/kryos-variables-<kind>.md (lossless under
        ! [kryos] write_lossless; the 0-D rows of the tables are lossless: time
        ! series stay lossless); the units attribute is the export's own. A
        ! rank-2 or rank-3 export without a table row stops; a 0-D export without
        ! one (derived coefficients, named by the par file) is written lossless
        ! with its name as long_name. init = .true. creates the file; any other
        ! stream is a backend's, so the backend must have handled it.
        class(component_t), intent(inout) :: self
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init

        character(len=*), parameter :: proc = "kryos_component::write_output"
        type(var_io_type) :: row
        integer, allocatable :: hs(:)
        logical :: ts, found
        integer :: i, h, ncid, n
        real(sp), allocatable :: a2(:,:), a3(:,:,:)
        real(sp) :: a0

        select case (trim(stream))
        case ("2D")
            ts = .false.
        case ("1D")
            ts = .true.
        case default
            call fesm_error(proc, "stream `"//trim(stream)//"` is not served by the default writer", &
                            "instance "//trim(self%name)//" (backend "//trim(self%backend)// &
                            "): its write_output must handle every stream other than 2D and 1D")
            return
        end select

        hs = selected(self, ts)
        if (size(hs) == 0) then
            if (init) call log_msg(LOG_WARN, proc, "instance "//trim(self%name)//": no "// &
                                   merge("0-D    ", "2D, 3D ", ts)//"exports selected for stream "// &
                                   trim(stream)//"; "//trim(file)//" not written")
            return
        end if

        if (init) then
            if (ts) then
                call outfile_dims_ts(file, t)
            else
                if (.not. (allocated(self%xkm) .and. allocated(self%ykm))) then
                    call fesm_error(proc, "instance "//trim(self%name)//" has no grid axes for "//trim(file), &
                                    "component_identity sets them from the init context; a global instance has none")
                end if
                call outfile_dims_2d(file, self%xkm, self%ykm, t)
                do i = 1, size(hs)
                    if (self%exports%f(hs(i))%rank == 3) call axis_write(self, file, self%exports%f(hs(i))%zdim)
                end do
            end if
        end if

        call outfile_open(file, t, ncid, n)
        do i = 1, size(hs)
            h = hs(i)
            if (.not. self%exports%f(h)%published) cycle
            associate (f => self%exports%f(h))
                call outfile_row(self%kind, f%name, row, found)
                if (.not. found) then
                    if (f%rank > 0) then
                        call fesm_error(proc, "export `"//trim(f%name)//"` of instance "//trim(self%name)// &
                                        " has no row in "//outfile_table_file(self%kind), &
                                        "add it to the table (names and units: contracts-design section 2)")
                    end if
                    row = scalar_row(f%name, f%units)
                end if
                call row_check(self, row, f%rank, f%units, f%zdim, f%name)
                row%units = f%units
                select case (f%rank)
                case (0)
                    call fields_get(self%exports, h, a0)
                    call nc_write(file, trim(row%varname), a0, dim1="time", start=[n], count=[1], ncid=ncid, &
                                  units=trim(row%units), long_name=trim(row%long_name), &
                                  prec_abs=real(row%prec_abs, sp), prec_rel=real(row%prec_rel, sp))
                case (2)
                    allocate(a2(f%nx, f%ny))
                    call fields_get(self%exports, h, a2)
                    call var_io_write(file, row, a2, n, ncid=ncid)
                    deallocate(a2)
                case (3)
                    allocate(a3(f%nx, f%ny, f%nz))
                    call fields_get(self%exports, h, a3)
                    call var_io_write(file, row, a3, n, ncid=ncid)
                    deallocate(a3)
                end select
            end associate
        end do
        call nc_close(ncid)
    end subroutine component_write_output

    function selected(self, ts) result(hs)
        ! Handles of the exports a stream writes: the 0-D ones for 1D, ranks 2
        ! and 3 for 2D, restricted to out_vars when it is set.
        class(component_t), intent(in) :: self
        logical,            intent(in) :: ts
        integer, allocatable :: hs(:)
        integer :: h
        logical :: pick

        hs = [integer ::]
        do h = 1, self%exports%n
            pick = (self%exports%f(h)%rank == 0) .eqv. ts
            if (pick .and. allocated(self%out_vars)) then
                if (size(self%out_vars) > 0) pick = any(self%out_vars == self%exports%f(h)%name)
            end if
            if (pick) hs = [hs, h]
        end do
    end function selected

    function scalar_row(name, units) result(row)
        ! The row of a 0-D export the table does not list: lossless, no long name.
        character(len=*), intent(in) :: name, units
        type(var_io_type) :: row
        row%varname   = name
        row%long_name = name
        row%units     = units
        row%dimnames  = ""
        row%ndims     = 1
        allocate(row%dims(1))
        row%dims(1)   = ""
        row%prec_abs  = 0.0d0
        row%prec_rel  = 0.0d0
    end function scalar_row

    subroutine axis_write(self, file, name)
        ! The third axis `name` of a rank-3 export, from the instance's axes.
        class(component_t), intent(in) :: self
        character(len=*),   intent(in) :: file, name
        integer :: k
        k = 0
        if (allocated(self%axes)) k = axis_find(self%axes, name)
        if (k == 0) then
            call fesm_error("kryos_component::write_output", "axis `"//trim(name)//"` of instance "// &
                            trim(self%name)//" is not among its axes", "file "//trim(file))
            return
        end if
        call outfile_axis(file, trim(name), self%axes(k)%v, trim(self%axes(k)%units))
    end subroutine axis_write

    subroutine row_check(self, row, rank, units, zdim, name)
        ! The table row of an export agrees with the export: units, rank and,
        ! for rank 3, the third axis.
        class(component_t), intent(in) :: self
        type(var_io_type),  intent(in) :: row
        integer,            intent(in) :: rank
        character(len=*),   intent(in) :: units, zdim, name

        character(len=*), parameter :: proc = "kryos_component::write_output"
        character(len=:), allocatable :: where
        integer :: nd

        where = "table "//outfile_table_file(self%kind)//", variable `"//trim(name)//"`, instance "//trim(self%name)
        if (.not. units_equivalent(row%units, units)) then
            call fesm_error(proc, "units `"//trim(row%units)//"` of the table differ from the export's `"// &
                            trim(units)//"`", where)
        end if
        nd = 0
        if (len_trim(row%dims(1)) > 0) nd = row%ndims
        if (nd /= rank) then
            call fesm_error(proc, "dimensions `"//trim(row%dimnames)//"` of the table do not match the export (rank "// &
                            to_str(rank)//")", where)
        end if
        if (rank == 3) then
            if (trim(row%dims(3)) /= trim(zdim)) then
                call fesm_error(proc, "third dimension `"//trim(row%dims(3))//"` of the table, the export's is `"// &
                                trim(zdim)//"`", where)
            end if
        end if
    end subroutine row_check

    subroutine component_output(c, file, t, stream, init)
        ! The one call site of the write_output binding for a producer held as
        ! class(component_t) (the output module); the per-component wrappers
        ! (ice_write_output, ...) do the same for their slot.
        class(component_t), intent(inout) :: c
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init
        call component_enter(c)
        call c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine component_output

    subroutine component_prepare_call(c, st)
        ! The one call site of the prepare binding (domain_prepare, for any
        ! component held as class(component_t)); puts no export, so no
        ! window to close.
        class(component_t), intent(inout) :: c
        type(step_time_t),  intent(in)    :: st
        call component_enter(c)
        call c%prepare(st)
        call component_leave()
    end subroutine component_prepare_call

    subroutine component_end(self)
        ! Default: nothing to release.
        class(component_t), intent(inout) :: self
        associate (unused => self%nx); end associate
    end subroutine component_end

    ! =========================================================================
    ! Factories: identity and admission (decision 12)
    ! =========================================================================

    subroutine component_identity(c, ctx, kind, con)
        ! Copy the instance's identity from ctx and set the backend's
        ! constraints; the first step of every factory and of <comp>_adopt.
        class(component_t),  intent(inout) :: c
        type(init_ctx_t),    intent(in)    :: ctx
        character(len=*),    intent(in)    :: kind
        type(constraints_t), intent(in)    :: con

        c%name    = ctx%instance
        c%kind    = kind
        c%backend = ctx%backend
        c%grid    = ctx%grid_name
        c%domain  = ctx%domain
        c%scope   = ctx%scope
        c%con     = con
        if (allocated(ctx%grid%G%x) .and. allocated(ctx%grid%G%y)) then   ! a global instance has no grid
            c%xkm = ctx%grid%G%x
            c%ykm = ctx%grid%G%y
        end if
    end subroutine component_identity

    subroutine backend_admit(c, census)
        ! Enforce c%con with named stops, then count the instance. Called by
        ! the factory before init, so a stub stops here with the rebuild hint.
        class(component_t), intent(in)    :: c
        type(census_t),     intent(inout) :: census

        character(len=*), parameter :: proc = "kryos_component::backend_admit"
        character(len=:), allocatable :: who, key
        integer :: k

        who = "instance `"//trim(c%name)//"` ("//trim(c%kind)//", "//trim(c%scope)// &
              ", domain `"//trim(c%domain)//"`, grid `"//trim(c%grid)//"`)"

        if (.not. c%con%linked) then
            call fesm_error(proc, "instance `"//trim(c%name)//"` needs backend `"//trim(c%backend)// &
                "`: rebuild kryos with `"//trim(c%con%build_flag)//"`", who)
        end if
        if (trim(c%scope) /= "regional" .and. trim(c%scope) /= "global") then
            call fesm_error(proc, "scope must be `regional` or `global`", who)
        end if
        if (trim(c%scope) == "regional" .and. .not. c%con%regional) then
            call fesm_error(proc, "backend `"//trim(c%backend)//"` of "//trim(c%kind)// &
                " has no regional scope", who)
        end if
        if (trim(c%scope) == "global" .and. .not. c%con%global) then
            call fesm_error(proc, "backend `"//trim(c%backend)//"` of "//trim(c%kind)// &
                " has no global scope", who)
        end if
        if (len_trim(c%con%fixed_grid) > 0 .and. trim(c%grid) /= trim(c%con%fixed_grid)) then
            call fesm_error(proc, "backend `"//trim(c%backend)//"` runs only on grid `"// &
                trim(c%con%fixed_grid)//"`", who)
        end if
        if (len_trim(c%con%domain) > 0 .and. trim(c%domain) /= trim(c%con%domain)) then
            call fesm_error(proc, "backend `"//trim(c%backend)//"` runs only in domain `"// &
                trim(c%con%domain)//"`", who)
        end if

        key = trim(c%kind)//":"//trim(c%backend)//":"//trim(c%scope)
        k = census_find(census, key)
        if (c%con%max_instances > 0 .and. census_count(census, key) >= c%con%max_instances) then
            if (trim(c%scope) == "regional") then
                call fesm_error(proc, "backend `"//trim(c%backend)//"` allows "// &
                    count_word(c%con%max_instances)//"; another was requested by domain `"// &
                    trim(c%domain)//"`", who//new_line("a")//"first admitted: instance `"// &
                    trim(census%first(k))//"`")
            else
                call fesm_error(proc, "backend `"//trim(c%backend)//"` allows "// &
                    count_word(c%con%max_instances)//"; another was requested by instance `"// &
                    trim(c%name)//"`", who//new_line("a")//"first admitted: instance `"// &
                    trim(census%first(k))//"`")
            end if
        end if

        if (k == 0) then
            if (.not. allocated(census%key)) then
                allocate(census%key(0), census%n(0), census%first(0))
            end if
            census%key   = [character(len=3*NAME_LEN) :: census%key, key]
            census%n     = [census%n, 0]
            census%first = [character(len=NAME_LEN) :: census%first, c%name]
            k = size(census%key)
        end if
        census%n(k) = census%n(k) + 1

    contains

        function count_word(n) result(s)
            integer, intent(in) :: n
            character(len=:), allocatable :: s
            if (n == 1) then
                s = "one instance"
            else
                s = to_str(n)//" instances"
            end if
        end function count_word

    end subroutine backend_admit

    integer function census_count(census, key) result(n)
        ! Instances admitted so far under "kind:backend:scope".
        type(census_t),   intent(in) :: census
        character(len=*), intent(in) :: key
        integer :: k
        n = 0
        k = census_find(census, key)
        if (k > 0) n = census%n(k)
    end function census_count

    integer function census_find(census, key) result(k)
        type(census_t),   intent(in) :: census
        character(len=*), intent(in) :: key
        integer :: i
        k = 0
        if (.not. allocated(census%key)) return
        do i = 1, size(census%key)
            if (trim(census%key(i)) == trim(key)) then
                k = i
                return
            end if
        end do
    end function census_find

    ! =========================================================================
    ! Step time of an instance (core-design 3.4)
    ! =========================================================================

    function step_time_of(cadence, clk, clock_mode, init) result(st)
        ! The step_time_t of an instance with this cadence: its interval
        ! [last mark, now] (cadence_interval; init: t0 = t1 = now), t_eval = t1,
        ! or the forcing time of the clock (clock_time_rel) when its `clock`
        ! key is "forcing"; dt = the cadence step (dtt for every step, 0 off).
        type(cadence_t),  intent(in) :: cadence
        type(clock_t),    intent(in) :: clk
        character(len=*), intent(in) :: clock_mode      ! model | forcing
        logical,          intent(in), optional :: init
        type(step_time_t) :: st

        logical :: i0

        i0 = .false.
        if (present(init)) i0 = init
        if (i0) then
            st%t0 = clock_time(clk)
            st%t1 = st%t0
        else
            call cadence_interval(cadence, clk, st%t0, st%t1)
        end if
        st%t_eval = st%t1
        st%tick   = clk%n
        if (trim(clock_mode) == "forcing") st%t_eval = clock_time_rel(clk)
        st%dt  = cadence_dt_years(cadence, clk)
        st%dtt = clock_dtt(clk)
    end function step_time_of

    subroutine sec_year_check(ctx, instance)
        ! A backend that converts with the year length takes the run calendar's
        ! (ctx%sec_year, core-design 3.6 item 5); stops if the context has none.
        type(init_ctx_t), intent(in) :: ctx
        character(len=*), intent(in) :: instance
        if (ctx%sec_year > 0.0_dp) return
        call fesm_error("kryos_component::sec_year_check", "init_ctx_t%sec_year is not set", &
                        "instance "//trim(instance)//": the core fills it from the year of the run clock")
    end subroutine sec_year_check

    ! =========================================================================
    ! Wrapper helpers
    ! =========================================================================

    subroutine component_close(c, st)
        ! Publish the window means of the exports at st%t1: fields_finalize_
        ! interval when the backend put any MEAN or ACCUM export in this call.
        ! A call that put no window-field (an init pass that publishes nothing,
        ! a state without MEAN fields) leaves the state as it was; a partial
        ! put stops in fields_finalize_interval, naming the field.
        class(component_t), intent(inout) :: c
        type(step_time_t),  intent(in)    :: st
        integer :: h

        do h = 1, c%exports%n
            if (c%exports%f(h)%tsem == FIELD_INSTANT) cycle
            if (c%exports%f(h)%wsum > 0.0_dp) then
                call fields_finalize_interval(c%exports, st%t1)
                return
            end if
        end do
    end subroutine component_close

    subroutine component_enter(c)
        ! Log prefix = instance name around each backend call (core-design 8).
        class(component_t), intent(in) :: c
        call log_set_prefix(trim(c%name))
    end subroutine component_enter

    subroutine component_leave()
        call log_set_prefix("")
    end subroutine component_leave

    subroutine component_optional(c, name)
        ! Mark import `name` as optional: wire_resolve accepts it unwired.
        class(component_t), intent(inout) :: c
        character(len=*),   intent(in)    :: name
        if (.not. allocated(c%opt_imports)) allocate(c%opt_imports(0))
        c%opt_imports = [character(len=FIELD_NAME_LEN) :: c%opt_imports, name]
    end subroutine component_optional

    subroutine component_core_filled(c, name)
        ! Mark import `name` as filled by the core outside the wiring rows (the
        ! global geo's external ice and reference, phase 4): wire_resolve does
        ! not look for a producer of it.
        class(component_t), intent(inout) :: c
        character(len=*),   intent(in)    :: name
        if (.not. allocated(c%core_imports)) allocate(c%core_imports(0))
        c%core_imports = [character(len=FIELD_NAME_LEN) :: c%core_imports, name]
    end subroutine component_core_filled

    function component_axis(c, ctx, name) result(nz)
        ! Length of axis `name` for a rank-3 declaration: the instance's own
        ! (set by a producer in init), else a producer's from ctx%axes, else
        ! the canonical month axis; adopted into c%axes. Stops if unknown.
        class(component_t), intent(inout) :: c
        type(init_ctx_t),   intent(in)    :: ctx
        character(len=*),   intent(in)    :: name
        integer :: nz

        type(axis_t) :: ax
        integer :: k

        if (.not. allocated(c%axes)) allocate(c%axes(0))
        k = axis_find(c%axes, name)
        if (k > 0) then
            nz = size(c%axes(k)%v)
            return
        end if

        if (allocated(ctx%axes)) k = axis_find(ctx%axes, name)
        if (k > 0) then
            ax = ctx%axes(k)
        else if (trim(name) == "month") then
            ax = month_axis()
        else
            call fesm_error("kryos_component::component_axis", "axis `"//trim(name)// &
                "` is not declared by any producer initialised before instance `"//trim(c%name)//"`", &
                "the producer (climate for `depth`) must be initialised first and pass its axis in ctx%axes")
        end if
        c%axes = [c%axes, ax]
        nz = size(ax%v)
    end function component_axis

    function month_axis() result(ax)
        ! The canonical month axis: 12 values 1..12.
        type(axis_t) :: ax
        integer :: m
        ax%name  = "month"
        ax%units = "month"
        allocate(ax%v(12))
        ax%v = [(real(m, dp), m = 1, 12)]
    end function month_axis

    integer function axis_find(axes, name) result(k)
        type(axis_t),     intent(in) :: axes(:)
        character(len=*), intent(in) :: name
        integer :: i
        k = 0
        do i = 1, size(axes)
            if (trim(axes(i)%name) == trim(name)) then
                k = i
                return
            end if
        end do
    end function axis_find

    logical function axes_equal(a, b)
        ! Same name, units and values (bitwise).
        type(axis_t), intent(in) :: a, b
        axes_equal = .false.
        if (trim(a%name) /= trim(b%name) .or. trim(a%units) /= trim(b%units)) return
        if (.not. (allocated(a%v) .and. allocated(b%v))) return
        if (size(a%v) /= size(b%v)) return
        axes_equal = all(a%v == b%v)
    end function axes_equal

    ! =========================================================================
    ! Masks and fractions (contracts-design section 8)
    ! =========================================================================

    elemental integer function kryos_mask(f_ice, f_grnd) result(mask)
        ! Majority class of a cell from its ice and grounded fractions.
        real(dp), intent(in) :: f_ice, f_grnd
        if (f_ice >= 0.5_dp .and. f_grnd >= 0.5_dp) then
            mask = MASK_GRND
        else if (f_ice >= 0.5_dp) then
            mask = MASK_FLOAT
        else if (f_grnd >= 0.5_dp) then
            mask = MASK_LAND
        else
            mask = MASK_OCEAN
        end if
    end function kryos_mask

end module kryos_component
