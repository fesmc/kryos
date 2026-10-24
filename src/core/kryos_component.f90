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

    use precision,      only: sp, dp
    use coordinates,    only: grid_class
    use regions,        only: regions_class
    use phys_constants, only: phys_const_class
    use fesm_log,       only: LOG_INFO, log_msg, log_set_prefix, fesm_error, to_str
    use fesm_fields,    only: field_state_t, fields_finalize_interval, FIELD_INSTANT, &
                              FIELD_NAME_LEN, FIELD_UNITS_LEN
    use fesm_clock,     only: cadence_t

    implicit none

    private

    ! Initialisation passes of init_state (core-design section 3.2; PASS_RESTART
    ! after a restore, review outcome).
    integer, parameter, public :: PASS_GEOMETRY = 1, PASS_BOUNDARY = 2, PASS_STATE = 3, &
                                  PASS_FINAL = 4, PASS_RESTART = 5

    ! The canonical four-value mask (contracts-design section 8).
    integer, parameter, public :: MASK_OCEAN = 0, MASK_LAND = 1, MASK_GRND = 2, MASK_FLOAT = 3

    integer, parameter, public :: NAME_LEN = 64, PATH_LEN = 512, NLIB = 8

    type, public :: caps_t                      ! capability flags, set by the backend in init
        logical :: has_ocn = .false., has_ocn_shelf = .false., has_smb = .false., has_Qd = .false.
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
        type(grid_class)        :: grid         ! from maps/grid_<name>.txt
        type(regions_class)     :: reg          ! domain regions on this grid
        type(grid_ref_t)        :: ref          ! reference topography on this grid
        type(phys_const_class)  :: cnst
        real(dp)                :: time = 0.0_dp, time_rel = 0.0_dp, dtt = 0.0_dp   ! t0
        type(axis_t), allocatable :: axes(:)    ! axes declared by the producers initialised so far
    end type init_ctx_t

    type, public :: step_time_t
        real(dp) :: t0 = 0.0_dp, t1 = 0.0_dp    ! model-time interval; t0 = t1 = time_init during init
        real(dp) :: t_eval = 0.0_dp             ! t1, or the forcing time at t1 (clock = "forcing")
        real(dp) :: dt = 0.0_dp                 ! the instance's cadence step [yr]
        real(dp) :: dtt = 0.0_dp                ! the run's step [yr]
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
    contains
        procedure(c_init),   deferred :: init           ! (ctx): read group; set prec, nx, ny, caps, opt
        procedure                     :: init_state => component_init_state      ! (st, pass)
        procedure(c_update), deferred :: update         ! (st): integrate t0 -> t1; put every export
        procedure(c_rst),    deferred :: restart_write, restart_read             ! (dir, t)
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

    public :: backend_admit, component_identity, component_close, component_enter, component_leave
    public :: component_axis, component_optional, month_axis, axis_find, axes_equal
    public :: census_count, kryos_mask
    public :: component_init_state, component_write_output, component_end

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

    subroutine component_write_output(self, file, t, stream, init)
        ! Default writer: the exports on the native grid (core-design section 6).
        ! Logs only until the output module lands (core C8).
        class(component_t), intent(inout) :: self
        character(len=*),   intent(in)    :: file, stream
        real(dp),           intent(in)    :: t
        logical,            intent(in)    :: init
        call log_msg(LOG_INFO, "kryos_component", "write_output "//trim(self%name)//": stream "// &
                     trim(stream)//" -> "//trim(file)//" at t = "//to_str(t)//merge(" (init)", "       ", init)// &
                     " skipped: default writer pending (core C8)")
    end subroutine component_write_output

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
