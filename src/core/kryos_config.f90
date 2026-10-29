module kryos_config
    ! The par file (core-design sections 1 and 5): config_t holds the par file,
    ! the template file input/kryos_defaults.nml, the library list and the
    ! [kryos] run keys; the free procedures read a key with its template
    ! default, find a named group in the par file or the library, validate
    ! a group against its templates, and register an instance's cadence from
    ! its keys dt, phase and clock (config_cadence, globe and domains alike).
    !
    ! Templates: one group per kind in input/kryos_defaults.nml ([kryos],
    ! [clock], [schedule], [globe], [domain], [instance]) and one per backend,
    ! [<component>_<backend>] ([ice_yelmo], ...). An instance group is valid if
    ! a key is declared by [instance] or by its backend's template. Instance,
    ! domain and run groups must be in the par file; named library groups
    ! (snapshot_*, coef_*, data fields) are found in [par file, library(:)].
    !
    ! Style: free procedures on a plain type; no OOP. Errors stop through
    ! fesm_error (core-design section 8); the procedures that validate or look
    ! a group up take an optional stat, so a caller (and test_config) can
    ! observe the failure instead of the stop. config_read stops on a key its
    ! templates do not declare: that is a code error, not a par file error.

    use precision, only: dp
    use fesm_log,  only: LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR, fesm_error, to_str
    use nml,       only: nml_read, nml_validate, nml_has_param, nml_has_group, nml_find_group, &
                         nml_group_names, nml_group_len
    use fesm_clock, only: clock_t, cadence_t, cadence_init

    implicit none

    private

    ! Same values as NAME_LEN, PATH_LEN, NLIB in kryos_component (config uses no
    ! other kryos module); test_config checks they agree.
    integer, parameter, public :: CFG_NAME_LEN = 64, CFG_PATH_LEN = 512, CFG_NLIB = 8

    character(len=*), parameter, public :: CFG_DEFAULTS_FILE = "input/kryos_defaults.nml"

    integer, parameter :: MAX_DOMAINS = 32

    type, public :: config_t
        character(len=CFG_PATH_LEN) :: par_file = "", defaults_file = CFG_DEFAULTS_FILE
        character(len=CFG_PATH_LEN) :: library(CFG_NLIB) = ""     ! files of the named groups not in the par file
        ! the [kryos] group (template [kryos])
        character(len=CFG_NAME_LEN), allocatable :: domains(:)    ! domain groups, in order
        character(len=CFG_NAME_LEN) :: globe = "", clock = "", restart_out = "", phys_const = ""
        character(len=CFG_PATH_LEN) :: restart = "", out_dir = "", map_fldr = "", phys_const_file = ""
        integer :: log_level = LOG_INFO                           ! LOG_DEBUG .. LOG_ERROR
        logical :: write_lossless = .false., timing = .true.
    end type config_t

    interface config_read
        ! (cfg, group, key, value [, template [, extra]]): the value of `key` in
        ! `group`, found in the par file or the library. With a template, the
        ! key must be declared by it (or by `extra`, the backend template) and
        ! takes its default there when the group leaves it out; without one the
        ! key must be set in the group. Vectors: the value must fit.
        module procedure config_read_str, config_read_str1, config_read_dp, config_read_int, config_read_log
    end interface config_read

    public :: config_init, config_read, config_has, config_find
    public :: config_validate, config_validate_instance, config_validate_experiment, config_cadence

contains

    ! =========================================================================
    ! Init
    ! =========================================================================

    subroutine config_init(cfg, par_file, defaults_file, stat)
        ! Open the par file: check the files exist, validate [kryos] against its
        ! template, read the library list and the [kryos] keys (defaults filled).
        type(config_t),   intent(out)           :: cfg
        character(len=*), intent(in)            :: par_file
        character(len=*), intent(in),  optional :: defaults_file
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_config::config_init"
        character(len=CFG_PATH_LEN) :: lib(CFG_NLIB)
        character(len=CFG_NAME_LEN) :: dom(MAX_DOMAINS), level
        integer :: i

        if (present(stat)) stat = 0
        cfg%par_file = par_file
        if (present(defaults_file)) cfg%defaults_file = defaults_file
        allocate(cfg%domains(0))

        if (.not. file_exists(cfg%par_file)) then
            call fesm_error(proc, "par file not found: "//trim(cfg%par_file), stat=stat)
            return
        end if
        if (.not. file_exists(cfg%defaults_file)) then
            call fesm_error(proc, "template file not found: "//trim(cfg%defaults_file), &
                            "par file "//trim(cfg%par_file), stat)
            return
        end if

        call config_validate(cfg, "kryos", "kryos", stat=stat)
        if (failed(stat)) return

        lib = ""
        call config_read(cfg, "kryos", "library", lib, "kryos")
        do i = 1, CFG_NLIB
            if (len_trim(lib(i)) == 0) cycle
            if (.not. file_exists(lib(i))) then
                call fesm_error(proc, "library file not found: "//trim(lib(i)), &
                                "[kryos] library of "//trim(cfg%par_file), stat)
                return
            end if
        end do
        cfg%library = lib

        dom = ""
        call config_read(cfg, "kryos", "domains",         dom,                  "kryos")
        call config_read(cfg, "kryos", "globe",           cfg%globe,            "kryos")
        call config_read(cfg, "kryos", "clock",           cfg%clock,            "kryos")
        call config_read(cfg, "kryos", "restart",         cfg%restart,          "kryos")
        call config_read(cfg, "kryos", "restart_out",     cfg%restart_out,      "kryos")
        call config_read(cfg, "kryos", "out_dir",         cfg%out_dir,          "kryos")
        call config_read(cfg, "kryos", "map_fldr",        cfg%map_fldr,         "kryos")
        call config_read(cfg, "kryos", "phys_const",      cfg%phys_const,       "kryos")
        call config_read(cfg, "kryos", "phys_const_file", cfg%phys_const_file,  "kryos")
        call config_read(cfg, "kryos", "write_lossless",  cfg%write_lossless,   "kryos")
        call config_read(cfg, "kryos", "timing",          cfg%timing,           "kryos")
        call config_read(cfg, "kryos", "log_level",       level,                "kryos")
        cfg%domains = pack(dom, len_trim(dom) > 0)

        select case (trim(level))
        case ("debug"); cfg%log_level = LOG_DEBUG
        case ("info");  cfg%log_level = LOG_INFO
        case ("warn");  cfg%log_level = LOG_WARN
        case ("error"); cfg%log_level = LOG_ERROR
        case default
            call fesm_error(proc, "[kryos] log_level = '"//trim(level)//"' is not recognised", &
                            "use debug, info, warn or error; par file "//trim(cfg%par_file), stat)
        end select
    end subroutine config_init

    ! =========================================================================
    ! Lookup
    ! =========================================================================

    function config_find(cfg, group, what, stat) result(file)
        ! The first of [par file, library(:)] that holds `group` (nml_find_group:
        ! stops listing the files searched, or stat = 1 and ""). `what` names the
        ! reader in the message, e.g. "&snapshot_lgm tas".
        type(config_t),   intent(in)            :: cfg
        character(len=*), intent(in)            :: group
        character(len=*), intent(in),  optional :: what
        integer,          intent(out), optional :: stat
        character(len=CFG_PATH_LEN) :: file

        file = nml_find_group([character(len=CFG_PATH_LEN) :: cfg%par_file, cfg%library], &
                              trim(group), what, stat)
    end function config_find

    logical function config_has(cfg, group, key) result(found)
        ! Is `key` set in `group` (not defaulted by a template)? False when no
        ! file holds the group. For keys whose presence matters (time_const).
        type(config_t),   intent(in) :: cfg
        character(len=*), intent(in) :: group, key
        character(len=CFG_PATH_LEN) :: file

        file = locate(cfg, group)
        found = .false.
        if (len_trim(file) > 0) found = nml_has_param(trim(file), trim(group), trim(key))
    end function config_has

    function locate(cfg, group) result(file)
        ! As config_find, quietly: "" when no file holds the group.
        type(config_t),   intent(in) :: cfg
        character(len=*), intent(in) :: group
        character(len=CFG_PATH_LEN) :: file
        integer :: i

        file = ""
        if (nml_has_group(trim(cfg%par_file), trim(group))) then
            file = cfg%par_file
            return
        end if
        do i = 1, CFG_NLIB
            if (len_trim(cfg%library(i)) == 0) cycle
            if (nml_has_group(trim(cfg%library(i)), trim(group))) then
                file = cfg%library(i)
                return
            end if
        end do
    end function locate

    ! =========================================================================
    ! Validation
    ! =========================================================================

    subroutine config_validate(cfg, group, template, extra, required, stat)
        ! `group` must be in the par file and every key it sets declared by the
        ! template group `template` or, if given, `extra` (both groups of
        ! cfg%defaults_file); every key of `required` must be set. Stops naming
        ! the group, the keys and the templates, or sets stat = 1 and returns.
        type(config_t),   intent(in)            :: cfg
        character(len=*), intent(in)            :: group, template
        character(len=*), intent(in),  optional :: extra
        character(len=*), intent(in),  optional :: required(:)
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_config::config_validate"
        character(len=CFG_PATH_LEN) :: tfile(2)
        character(len=CFG_NAME_LEN) :: tgroup(2)
        integer :: i, nt, nbad

        if (present(stat)) stat = 0

        call require_group(cfg, group, template, proc, stat)
        if (failed(stat)) return

        nt = 1
        tgroup(1) = template
        if (present(extra)) then
            nt = 2
            tgroup(2) = extra
        end if
        tfile = cfg%defaults_file
        do i = 1, nt
            if (.not. nml_has_group(trim(cfg%defaults_file), trim(tgroup(i)))) then
                call fesm_error(proc, "no template &"//trim(tgroup(i))//" in "//trim(cfg%defaults_file), &
                                "validating group &"//trim(group)//" of "//trim(cfg%par_file), stat)
                return
            end if
        end do

        if (present(stat)) then
            call nml_validate(trim(cfg%par_file), trim(group), tfile(1:nt), tgroup(1:nt), &
                              halt_on_error=.false., nbad=nbad)
            if (nbad > 0) then
                stat = 1
                return
            end if
        else
            call nml_validate(trim(cfg%par_file), trim(group), tfile(1:nt), tgroup(1:nt))
        end if

        if (present(required)) then
            do i = 1, size(required)
                if (nml_has_param(trim(cfg%par_file), trim(group), trim(required(i)))) cycle
                call fesm_error(proc, "required key '"//trim(required(i))//"' is not set in group &"// &
                                trim(group), "par file "//trim(cfg%par_file), stat)
                return
            end do
        end if
    end subroutine config_validate

    subroutine config_validate_instance(cfg, group, component, backend, stat)
        ! Validate an instance group of `component` ("ice", "geo", ...): its
        ! `backend` key (required) names the template [<component>_<backend>],
        ! which with [instance] declares the valid keys. Returns the backend.
        type(config_t),   intent(in)            :: cfg
        character(len=*), intent(in)            :: group, component
        character(len=*), intent(out)           :: backend
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_config::config_validate_instance"
        character(len=CFG_NAME_LEN) :: extra

        backend = ""
        if (present(stat)) stat = 0

        ! the group and its backend key first: the template depends on it
        call require_group(cfg, group, "instance", proc, stat)
        if (failed(stat)) return
        if (.not. nml_has_param(trim(cfg%par_file), trim(group), "backend")) then
            call fesm_error(proc, "required key 'backend' is not set in instance group &"//trim(group), &
                            "component "//trim(component)//"; known backends: "// &
                            known_backends(cfg, component)//"; par file "//trim(cfg%par_file), stat)
            return
        end if

        call nml_read(trim(cfg%par_file), trim(group), "backend", backend, init=.true.)
        if (len_trim(backend) == 0) then
            call fesm_error(proc, "instance &"//trim(group)//" ("//trim(component)//") has an empty backend", &
                            "par file "//trim(cfg%par_file), stat)
            return
        end if

        extra = trim(component)//"_"//trim(backend)
        if (.not. nml_has_group(trim(cfg%defaults_file), trim(extra))) then
            call fesm_error(proc, "instance &"//trim(group)//": backend '"//trim(backend)// &
                            "' has no template &"//trim(extra)//" in "//trim(cfg%defaults_file), &
                            "known backends of "//trim(component)//": "//known_backends(cfg, component), stat)
            backend = ""
            return
        end if

        call config_validate(cfg, group, "instance", extra, stat=stat)
        if (failed(stat)) backend = ""
    end subroutine config_validate_instance

    subroutine config_validate_experiment(cfg, group, kind, stat)
        ! Validate an experiment group (core-design 3.9): its `kind` key
        ! (required) names the template [experiment_<kind>], which with
        ! [experiment] declares the valid keys. Returns the kind. The twin of
        ! config_validate_instance.
        type(config_t),   intent(in)            :: cfg
        character(len=*), intent(in)            :: group
        character(len=*), intent(out)           :: kind
        integer,          intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_config::config_validate_experiment"
        character(len=CFG_NAME_LEN) :: extra

        kind = ""
        if (present(stat)) stat = 0

        call require_group(cfg, group, "experiment", proc, stat)
        if (failed(stat)) return
        if (.not. nml_has_param(trim(cfg%par_file), trim(group), "kind")) then
            call fesm_error(proc, "required key 'kind' is not set in experiment group &"//trim(group), &
                            "known kinds: "//known_backends(cfg, "experiment")//"; par file "//trim(cfg%par_file), stat)
            return
        end if

        call nml_read(trim(cfg%par_file), trim(group), "kind", kind, init=.true.)
        if (len_trim(kind) == 0) then
            call fesm_error(proc, "experiment &"//trim(group)//" has an empty kind", &
                            "par file "//trim(cfg%par_file), stat)
            return
        end if

        extra = "experiment_"//trim(kind)
        if (.not. nml_has_group(trim(cfg%defaults_file), trim(extra))) then
            call fesm_error(proc, "experiment &"//trim(group)//": kind '"//trim(kind)// &
                            "' has no template &"//trim(extra)//" in "//trim(cfg%defaults_file), &
                            "known kinds: "//known_backends(cfg, "experiment"), stat)
            kind = ""
            return
        end if

        call config_validate(cfg, group, "experiment", extra, stat=stat)
        if (failed(stat)) kind = ""
    end subroutine config_validate_experiment

    ! =========================================================================
    ! Cadence of an instance
    ! =========================================================================

    subroutine config_cadence(cfg, group, component, backend, clk, cadence, clock_mode)
        ! The keys dt, phase and clock of an instance group (templates
        ! [instance] + [<component>_<backend>]): the instance's cadence on the
        ! clock (dt = 0: every step; dt > 0: the grid phase + k*dt; dt < 0: off)
        ! and its clock (model | forcing). Call after the instance's init and
        ! before clock_start (the cadence registers with the clock).
        type(config_t),   intent(in)    :: cfg
        character(len=*), intent(in)    :: group, component, backend
        type(clock_t),    intent(inout) :: clk
        type(cadence_t),  intent(inout) :: cadence
        character(len=*), intent(out)   :: clock_mode

        character(len=CFG_NAME_LEN) :: extra
        real(dp) :: dt, phase

        extra = trim(component)//"_"//trim(backend)
        dt = 0.0_dp;  phase = 0.0_dp;  clock_mode = ""
        call config_read(cfg, group, "dt",    dt,         "instance", extra)
        call config_read(cfg, group, "phase", phase,      "instance", extra)
        call config_read(cfg, group, "clock", clock_mode, "instance", extra)

        select case (trim(clock_mode))
        case ("model", "forcing")
        case default
            call fesm_error("kryos_config::config_cadence", "instance &"//trim(group)//": clock = '"// &
                            trim(clock_mode)//"' is not recognised", "use model or forcing; par file "// &
                            trim(cfg%par_file))
        end select

        if (dt == 0.0_dp) then
            call cadence_init(cadence, clk, trim(group), .true.)
        else
            call cadence_init(cadence, clk, trim(group), dt, phase)
        end if
    end subroutine config_cadence

    subroutine require_group(cfg, group, template, proc, stat)
        ! Instance, domain and run groups must be in the par file itself.
        type(config_t),   intent(in)            :: cfg
        character(len=*), intent(in)            :: group, template, proc
        integer,          intent(out), optional :: stat

        if (present(stat)) stat = 0
        if (nml_has_group(trim(cfg%par_file), trim(group))) return
        call fesm_error(proc, "group &"//trim(group)//" is missing from the par file "// &
                        trim(cfg%par_file), "template &"//trim(template)//" of "// &
                        trim(cfg%defaults_file), stat)
    end subroutine require_group

    function known_backends(cfg, component) result(list)
        ! The backends with a template [<component>_<backend>], blank separated
        ! (also the kinds of "experiment").
        type(config_t),   intent(in) :: cfg
        character(len=*), intent(in) :: component
        character(len=:), allocatable :: list
        character(len=nml_group_len), allocatable :: names(:)
        integer :: i, n

        list = ""
        names = nml_group_names(trim(cfg%defaults_file))
        n = len_trim(component) + 1
        do i = 1, size(names)
            if (len_trim(names(i)) <= n) cycle
            if (names(i)(1:n) /= trim(component)//"_") cycle
            list = list//" "//trim(names(i)(n+1:))
        end do
        list = trim(adjustl(list))
        if (len(list) == 0) list = "(none)"
    end function known_backends

    ! =========================================================================
    ! Reads
    ! =========================================================================

    subroutine resolve(cfg, group, key, template, extra, file, dgroup)
        ! Where `key` of `group` is read: `file` holds the group; `dgroup` is
        ! the template group that declares the key ("" without a template).
        type(config_t),   intent(in)  :: cfg
        character(len=*), intent(in)  :: group, key
        character(len=*), intent(in),  optional :: template, extra
        character(len=*), intent(out) :: file, dgroup

        character(len=*), parameter :: proc = "kryos_config::config_read"
        character(len=:), allocatable :: tmpl

        file = config_find(cfg, group, "&"//trim(group)//" "//trim(key))
        dgroup = ""

        if (.not. present(template)) then
            if (.not. nml_has_param(trim(file), trim(group), trim(key))) &
                call fesm_error(proc, "key '"//trim(key)//"' is not set in group &"//trim(group), &
                                "file "//trim(file)//"; no template given, so no default")
            return
        end if

        if (nml_has_param(trim(cfg%defaults_file), trim(template), trim(key))) then
            dgroup = template
        else if (present(extra)) then
            if (nml_has_param(trim(cfg%defaults_file), trim(extra), trim(key))) dgroup = extra
        end if
        if (len_trim(dgroup) == 0) then
            tmpl = "&"//trim(template)
            if (present(extra)) tmpl = tmpl//" &"//trim(extra)
            call fesm_error(proc, "key '"//trim(key)//"' of group &"//trim(group)// &
                            " is declared by no template", "templates "//tmpl//" in "// &
                            trim(cfg%defaults_file))
        end if
    end subroutine resolve

    subroutine config_read_str(cfg, group, key, value, template, extra)
        type(config_t),   intent(in)           :: cfg
        character(len=*), intent(in)           :: group, key
        character(len=*), intent(inout)        :: value
        character(len=*), intent(in), optional :: template, extra
        character(len=CFG_PATH_LEN) :: file
        character(len=CFG_NAME_LEN) :: dgroup

        call resolve(cfg, group, key, template, extra, file, dgroup)
        if (len_trim(dgroup) > 0) then
            call nml_read(trim(file), trim(group), trim(key), value, init=.true., &
                          defaults_file=trim(cfg%defaults_file), defaults_group=trim(dgroup))
        else
            call nml_read(trim(file), trim(group), trim(key), value, init=.true.)
        end if
    end subroutine config_read_str

    subroutine config_read_str1(cfg, group, key, value, template, extra)
        ! A list of strings; stops if the group sets more than size(value).
        type(config_t),   intent(in)           :: cfg
        character(len=*), intent(in)           :: group, key
        character(len=*), intent(inout)        :: value(:)
        character(len=*), intent(in), optional :: template, extra
        character(len=CFG_PATH_LEN) :: file
        character(len=CFG_NAME_LEN) :: dgroup
        character(len=len(value))   :: tmp(size(value)+1)

        call resolve(cfg, group, key, template, extra, file, dgroup)
        tmp = ""
        if (len_trim(dgroup) > 0) then
            call nml_read(trim(file), trim(group), trim(key), tmp, init=.true., &
                          defaults_file=trim(cfg%defaults_file), defaults_group=trim(dgroup))
        else
            call nml_read(trim(file), trim(group), trim(key), tmp, init=.true.)
        end if
        if (len_trim(tmp(size(value)+1)) > 0) &
            call fesm_error("kryos_config::config_read", "group &"//trim(group)//" key '"//trim(key)// &
                            "' lists more than "//to_str(size(value))//" values", "file "//trim(file))
        value = tmp(1:size(value))
    end subroutine config_read_str1

    subroutine config_read_dp(cfg, group, key, value, template, extra)
        type(config_t),   intent(in)           :: cfg
        character(len=*), intent(in)           :: group, key
        real(dp),         intent(inout)        :: value
        character(len=*), intent(in), optional :: template, extra
        character(len=CFG_PATH_LEN) :: file
        character(len=CFG_NAME_LEN) :: dgroup

        call resolve(cfg, group, key, template, extra, file, dgroup)
        if (len_trim(dgroup) > 0) then
            call nml_read(trim(file), trim(group), trim(key), value, init=.true., &
                          defaults_file=trim(cfg%defaults_file), defaults_group=trim(dgroup))
        else
            call nml_read(trim(file), trim(group), trim(key), value, init=.true.)
        end if
    end subroutine config_read_dp

    subroutine config_read_int(cfg, group, key, value, template, extra)
        type(config_t),   intent(in)           :: cfg
        character(len=*), intent(in)           :: group, key
        integer,          intent(inout)        :: value
        character(len=*), intent(in), optional :: template, extra
        character(len=CFG_PATH_LEN) :: file
        character(len=CFG_NAME_LEN) :: dgroup

        call resolve(cfg, group, key, template, extra, file, dgroup)
        if (len_trim(dgroup) > 0) then
            call nml_read(trim(file), trim(group), trim(key), value, init=.true., &
                          defaults_file=trim(cfg%defaults_file), defaults_group=trim(dgroup))
        else
            call nml_read(trim(file), trim(group), trim(key), value, init=.true.)
        end if
    end subroutine config_read_int

    subroutine config_read_log(cfg, group, key, value, template, extra)
        type(config_t),   intent(in)           :: cfg
        character(len=*), intent(in)           :: group, key
        logical,          intent(inout)        :: value
        character(len=*), intent(in), optional :: template, extra
        character(len=CFG_PATH_LEN) :: file
        character(len=CFG_NAME_LEN) :: dgroup

        call resolve(cfg, group, key, template, extra, file, dgroup)
        if (len_trim(dgroup) > 0) then
            call nml_read(trim(file), trim(group), trim(key), value, init=.true., &
                          defaults_file=trim(cfg%defaults_file), defaults_group=trim(dgroup))
        else
            call nml_read(trim(file), trim(group), trim(key), value, init=.true.)
        end if
    end subroutine config_read_log

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    logical function file_exists(file) result(found)
        character(len=*), intent(in) :: file
        inquire(file=trim(file), exist=found)
    end function file_exists

    logical function failed(stat)
        ! stat was set to a nonzero value by a callee (always false without stat).
        integer, intent(in), optional :: stat
        failed = .false.
        if (present(stat)) failed = stat /= 0
    end function failed

end module kryos_config
