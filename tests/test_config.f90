program test_config
    ! kryos_config (core-design section 5): defaults filled, library lookup
    ! across files, unknown key and missing group (and the other par file
    ! errors) rejected. Par files are written to tmp/test_config/ (run from the
    ! repo root, as make check does) and validated against the real
    ! input/kryos_defaults.nml.
    !
    ! Abort paths: fesm_error stops the program, so each failing call is made
    ! with `stat` (the procedures that validate or look a group up have it) and
    ! stat = 1 is checked; the stopping form (no stat) is run in a child
    ! process of this executable ("abort <case>") and must exit non-zero.

    use precision,       only: dp
    use fesm_log,        only: LOG_WARN, LOG_INFO
    use kryos_component, only: NAME_LEN, PATH_LEN, NLIB
    use kryos_config,    only: config_t, config_init, config_read, config_has, config_find, &
                               config_validate, config_validate_instance, config_validate_experiment, &
                               CFG_NAME_LEN, CFG_PATH_LEN, CFG_NLIB

    implicit none

    character(len=*), parameter :: dir = "tmp/test_config/"
    character(len=*), parameter :: par = dir//"par.nml", lib = dir//"lib.nml"
    integer :: nfail = 0

    character(len=512) :: arg

    call get_command_argument(1, arg)
    if (trim(arg) == "abort") then
        call get_command_argument(2, arg)
        call abort_case(trim(arg))          ! must stop; falling through is a failure
        stop 0
    end if

    call execute_command_line("mkdir -p "//dir)

    call check(CFG_NAME_LEN == NAME_LEN .and. CFG_PATH_LEN == PATH_LEN .and. CFG_NLIB == NLIB, &
               "lengths agree with kryos_component")

    call write_par_files()
    call test_defaults()
    call test_read()
    call test_library()
    call test_validate()
    call test_errors()
    call test_aborts()

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_config: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_config: OK"

contains

    ! -------------------------------------------------------------------------
    ! Tests
    ! -------------------------------------------------------------------------

    subroutine test_defaults()
        ! Keys the par file leaves out take their template value.
        type(config_t) :: cfg
        integer :: stat

        call config_init(cfg, par, stat=stat)
        call check(stat == 0, "config_init")
        call check(size(cfg%domains) == 2, "domains: two")
        call check(trim(cfg%domains(1)) == "Greenland" .and. trim(cfg%domains(2)) == "Antarctica", "domains: names")
        call check(trim(cfg%restart_out) == "out_rst", "set value: restart_out")
        call check(trim(cfg%library(1)) == lib .and. len_trim(cfg%library(2)) == 0, "library list")
        call check(cfg%log_level == LOG_WARN, "set value: log_level")
        call check(trim(cfg%restart) == "None", "default: restart")
        call check(trim(cfg%out_dir) == "./" .and. trim(cfg%map_fldr) == "maps", "default: out_dir, map_fldr")
        call check(trim(cfg%globe) == "globe" .and. trim(cfg%clock) == "clock", "default: globe, clock")
        call check(trim(cfg%phys_const) == "Earth" .and. trim(cfg%phys_const_file) == "input/phys_const.nml", &
                   "default: phys_const")
        call check(.not. cfg%write_lossless .and. cfg%timing, "default: write_lossless, timing")

        ! a minimal [kryos] group: everything defaulted, no domain, no library
        call config_init(cfg, dir//"min.nml", stat=stat)
        call check(stat == 0, "config_init: minimal [kryos]")
        call check(size(cfg%domains) == 0 .and. cfg%log_level == LOG_INFO .and. len_trim(cfg%library(1)) == 0, &
                   "minimal [kryos]: defaults only")
    end subroutine test_defaults

    subroutine test_read()
        ! config_read: template defaults, instance + backend template, set values.
        type(config_t) :: cfg
        character(len=64) :: s, names(4)
        real(dp) :: x
        integer :: i

        call config_init(cfg, par)

        ! instance group, core template and backend template
        x = -1.0_dp
        call config_read(cfg, "ice_grl", "dt", x, "instance", "ice_yelmo")
        call check(x == 0.0_dp, "default dp from [instance]")
        call config_read(cfg, "ice_grl", "init_equil_time", x, "instance", "ice_yelmo")
        call check(x == 10.0_dp, "set dp from the group, declared by the backend template")
        x = -1.0_dp
        call config_read(cfg, "ice_grl", "init_time_thrm", x, "instance", "ice_yelmo")
        call check(x == 0.0_dp, "init_time_thrm: default 0 from the backend template")
        call config_read(cfg, "ice_grl", "yelmo", s, "instance", "ice_yelmo")
        call check(trim(s) == "yelmo", "default string from the backend template")
        call config_read(cfg, "ice_grl", "out_1D", s, "instance", "ice_yelmo")
        call check(trim(s) == "none", "default string from [instance]")
        call config_read(cfg, "ice_grl", "out_2D", s, "instance", "ice_yelmo")
        call check(trim(s) == "out_2D", "set string")
        call config_read(cfg, "ice_grl", "out_vars", s, "instance", "ice_yelmo")
        call check(len_trim(s) == 0, "empty default string is reset")
        call config_read(cfg, "geo_glob", "dt", x, "instance", "geo_fastisostasy")
        call check(x == 10.0_dp, "set dp: dt")

        ! the domain's list of experiments: set, and empty by default
        names = "x"
        call config_read(cfg, "Greenland", "experiments", names, "domain")
        call check(trim(names(1)) == "exp_host" .and. trim(names(2)) == "exp_other" .and. len_trim(names(3)) == 0, &
                   "domain experiments: the list in order")
        names = "x"
        call config_read(cfg, "Bad_domain", "experiments", names, "domain")
        call check(all(len_trim(names) == 0), "domain experiments: default none")

        ! integer, logical, vector
        call config_read(cfg, "clock", "ticks_per_year", i, "clock")
        call check(i == 36000, "default integer")
        call config_read(cfg, "clock", "time_end", x, "clock")
        call check(x == 15.0e3_dp, "set dp with exponent")
        call config_read(cfg, "clock", "calendar", s, "clock")
        call check(trim(s) == "365_day", "default string: calendar")
        names = "x"
        call config_read(cfg, "kryos", "domains", names, "kryos")
        call check(trim(names(1)) == "Greenland" .and. trim(names(2)) == "Antarctica" .and. &
                   len_trim(names(3)) == 0, "string vector, shorter than the array")

        ! presence: set in the group, not defaulted
        call check(config_has(cfg, "clock", "time_const"), "config_has: set key")
        call check(.not. config_has(cfg, "clock", "time_const_cal"), "config_has: defaulted key")
        call check(.not. config_has(cfg, "no_such_group", "dt"), "config_has: no group")
        call check(config_has(cfg, "snapshot_pd", "tas_file"), "config_has: library group")
    end subroutine test_read

    subroutine test_library()
        ! Named groups are found in the par file first, then the library files.
        type(config_t) :: cfg
        character(len=CFG_PATH_LEN) :: f
        character(len=64) :: s
        integer :: stat

        call config_init(cfg, par)

        f = config_find(cfg, "Greenland")
        call check(trim(f) == par, "find: group in the par file")
        f = config_find(cfg, "snapshot_pd")
        call check(trim(f) == lib, "find: group only in the library")
        f = config_find(cfg, "ice_grl")
        call check(trim(f) == par, "find: the par file wins over the library")
        f = config_find(cfg, "coef_alpha")
        call check(trim(f) == lib, "find: second group of the library")

        call config_read(cfg, "snapshot_pd", "tas_file", s)
        call check(trim(s) == "tas_pd.nc", "read a library group without template")
        call config_read(cfg, "ice_grl", "backend", s, "instance")
        call check(trim(s) == "yelmo", "read from the par file where the library has the group too")

        f = config_find(cfg, "no_such_group", "test", stat)
        call check(stat == 1 .and. len_trim(f) == 0, "find: missing group sets stat = 1")
    end subroutine test_library

    subroutine test_validate()
        ! Valid groups pass; a key is valid if the core or the backend template declares it.
        type(config_t) :: cfg
        character(len=64) :: backend
        integer :: stat

        call config_init(cfg, par)

        call config_validate(cfg, "Greenland", "domain", required=["grid_hub"], stat=stat)
        call check(stat == 0, "validate: domain group")
        call config_validate(cfg, "clock", "clock", required=["time_init", "time_end ", "dtt      "], stat=stat)
        call check(stat == 0, "validate: clock group with required keys")
        call config_validate_instance(cfg, "ice_grl", "ice", backend, stat)
        call check(stat == 0 .and. trim(backend) == "yelmo", "validate: instance, core and backend keys")
        call config_validate_instance(cfg, "geo_glob", "geo", backend, stat)
        call check(stat == 0 .and. trim(backend) == "fastisostasy", "validate: global geo instance")
        call config_validate_experiment(cfg, "exp_host", backend, stat)
        call check(stat == 0 .and. trim(backend) == "host", "validate: experiment of kind host")
    end subroutine test_validate

    subroutine test_errors()
        ! Every par file error comes back as stat = 1.
        type(config_t) :: cfg, bad     ! bad: a config_init that fails leaves its cfg unusable
        character(len=64) :: backend
        integer :: stat

        call config_init(cfg, par)

        ! unknown key
        call config_validate_instance(cfg, "ice_bad", "ice", backend, stat)
        call check(stat == 1 .and. len_trim(backend) == 0, "unknown key in an instance group")
        call config_validate_instance(cfg, "ice_other", "ice", backend, stat)
        call check(stat == 1, "key of another backend's template (bsl in ice_yelmo)")
        call config_validate(cfg, "Bad_domain", "domain", stat=stat)
        call check(stat == 1, "unknown key in a domain group")
        call config_init(bad, dir//"bad_kryos.nml", stat=stat)
        call check(stat == 1, "unknown key in [kryos]")

        ! missing group
        call config_validate(cfg, "NoSuchDomain", "domain", stat=stat)
        call check(stat == 1, "missing domain group")
        call config_validate_instance(cfg, "no_such_inst", "ice", backend, stat)
        call check(stat == 1, "missing instance group")
        call config_validate(cfg, "snapshot_pd", "domain", stat=stat)
        call check(stat == 1, "a library group is not a par file group")
        call config_init(bad, dir//"no_kryos.nml", stat=stat)
        call check(stat == 1, "missing [kryos] group")

        ! required key, backend
        call config_validate(cfg, "Bad_domain", "domain", required=["grid_hub"], stat=stat)
        call check(stat == 1, "required key missing (also unknown key)")
        call config_validate(cfg, "Empty_domain", "domain", required=["grid_hub"], stat=stat)
        call check(stat == 1, "required key missing")
        call config_validate_instance(cfg, "ice_nobackend", "ice", backend, stat)
        call check(stat == 1, "instance without backend key")
        call config_validate_instance(cfg, "ice_unknown", "ice", backend, stat)
        call check(stat == 1 .and. len_trim(backend) == 0, "backend without template")
        call config_validate_instance(cfg, "geo_wrongkind", "ice", backend, stat)
        call check(stat == 1, "backend of another component")

        ! experiments
        call config_validate_experiment(cfg, "exp_badkey", backend, stat)
        call check(stat == 1 .and. len_trim(backend) == 0, "unknown key in an experiment group")
        call config_validate_experiment(cfg, "exp_nokind", backend, stat)
        call check(stat == 1, "experiment without kind key")
        call config_validate_experiment(cfg, "exp_emptykind", backend, stat)
        call check(stat == 1, "experiment with an empty kind")
        call config_validate_experiment(cfg, "exp_unknown", backend, stat)
        call check(stat == 1 .and. len_trim(backend) == 0, "kind without template")
        call config_validate_experiment(cfg, "no_such_exp", backend, stat)
        call check(stat == 1, "missing experiment group")

        ! files, values
        call config_init(bad, dir//"missing.nml", stat=stat)
        call check(stat == 1, "par file not found")
        call config_init(bad, par, defaults_file=dir//"missing_defaults.nml", stat=stat)
        call check(stat == 1, "template file not found")
        call config_init(bad, dir//"bad_lib.nml", stat=stat)
        call check(stat == 1, "library file not found")
        call config_init(bad, dir//"bad_level.nml", stat=stat)
        call check(stat == 1, "unknown log_level")
    end subroutine test_errors

    subroutine test_aborts()
        ! The stopping form: no stat, the child process must exit non-zero.
        call check(aborts("unknown-key"),   "abort: unknown key stops")
        call check(aborts("missing-group"), "abort: missing group stops")
        call check(aborts("missing-find"),  "abort: library lookup of a missing group stops")
        call check(aborts("undeclared"),    "abort: config_read of an undeclared key stops")
        call check(aborts("exp-unknown"),   "abort: an experiment of an unknown kind stops")
        call check(aborts("exp-badkey"),    "abort: an unknown key in an experiment group stops")
        call check(.not. aborts("valid"),   "abort: a valid group does not stop")
    end subroutine test_aborts

    subroutine abort_case(name)
        ! Runs in the child process.
        character(len=*), intent(in) :: name
        type(config_t) :: cfg
        character(len=64) :: backend, s

        call config_init(cfg, par)
        select case (name)
        case ("unknown-key");   call config_validate_instance(cfg, "ice_bad", "ice", backend)
        case ("missing-group"); call config_validate(cfg, "NoSuchDomain", "domain")
        case ("missing-find");  s = config_find(cfg, "no_such_group")
        case ("undeclared");    call config_read(cfg, "ice_grl", "bogus", s, "instance", "ice_yelmo")
        case ("exp-unknown");   call config_validate_experiment(cfg, "exp_unknown", backend)
        case ("exp-badkey");    call config_validate_experiment(cfg, "exp_badkey", backend)
        case ("valid");         call config_validate_instance(cfg, "ice_grl", "ice", backend)
        end select
    end subroutine abort_case

    logical function aborts(name)
        ! True if `<this program> abort <name>` exits non-zero.
        character(len=*), intent(in) :: name
        character(len=512) :: exe
        integer :: es

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" abort "//name//" > /dev/null 2>&1", exitstat=es)
        aborts = es /= 0
    end function aborts

    ! -------------------------------------------------------------------------
    ! Par files
    ! -------------------------------------------------------------------------

    subroutine write_par_files()
        call put(par, [character(len=100) :: &
            '&kryos', &
            '    domains     = "Greenland" "Antarctica"', &
            '    restart_out = "out_rst"', &
            '    library     = "'//lib//'"', &
            '    log_level   = "warn"', &
            '/', &
            '&clock', &
            '    method     = "const"', &
            '    time_const = 0.0', &
            '    time_init  = 0.0', &
            '    time_end   = 15e3', &
            '    dtt        = 10.0', &
            '/', &
            '&Greenland', &
            '    grid_hub    = "GRL-16KM"', &
            '    ice         = "ice_grl"', &
            '    experiments = "exp_host" "exp_other"', &
            '/', &
            '&exp_host', &
            '    kind = "host"', &
            '/', &
            '&exp_badkey', &
            '    kind  = "host"', &
            '    bogus = 1.0', &
            '/', &
            '&exp_nokind', &
            '/', &
            '&exp_emptykind', &
            '    kind = ""', &
            '/', &
            '&exp_unknown', &
            '    kind = "nonexistent"', &
            '/', &
            '&ice_grl', &
            '    backend         = "yelmo"', &
            '    init_equil_time = 10.0', &
            '    out_2D          = "out_2D"', &
            '/', &
            '&geo_glob', &
            '    backend = "fastisostasy"', &
            '    dt      = 10.0', &
            '    bsl     = "barysealevel"', &
            '/', &
            '&ice_bad', &
            '    backend = "yelmo"', &
            '    bogus   = 1.0', &
            '/', &
            '&ice_other', &
            '    backend = "yelmo"', &
            '    bsl     = "barysealevel"', &
            '/', &
            '&Bad_domain', &
            '    ice     = "ice_grl"', &
            '    bogus   = "x"', &
            '/', &
            '&Empty_domain', &
            '/', &
            '&ice_nobackend', &
            '    dt = 10.0', &
            '/', &
            '&ice_unknown', &
            '    backend = "nonexistent"', &
            '/', &
            '&geo_wrongkind', &
            '    backend = "fastisostasy"', &
            '/'])

        ! library: a group also in the par file (the par file wins), data groups
        call put(lib, [character(len=100) :: &
            '&snapshot_pd', &
            '    tas_file = "tas_pd.nc"', &
            '/', &
            '&coef_alpha', &
            '    method = "const"', &
            '/', &
            '&ice_grl', &
            '    backend = "library_loses"', &
            '/'])

        call put(dir//"min.nml",       [character(len=100) :: '&kryos', '/'])
        call put(dir//"no_kryos.nml",  [character(len=100) :: '&clock', '    dtt = 1.0', '/'])
        call put(dir//"bad_kryos.nml", [character(len=100) :: '&kryos', '    domains = "Greenland"', &
                                                              '    domain  = "typo"', '/'])
        call put(dir//"bad_lib.nml",   [character(len=100) :: '&kryos', '    library = "'//dir//'nothere.nml"', '/'])
        call put(dir//"bad_level.nml", [character(len=100) :: '&kryos', '    log_level = "loud"', '/'])
    end subroutine write_par_files

    subroutine put(file, lines)
        character(len=*), intent(in) :: file, lines(:)
        integer :: u, i

        open(newunit=u, file=file, status="replace", action="write")
        do i = 1, size(lines)
            write(u,'(a)') trim(lines(i))
        end do
        close(u)
    end subroutine put

    subroutine check(ok, msg)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) then
            nfail = nfail + 1
            write(*,'(a)') "test_config: FAIL: "//msg
        end if
    end subroutine check

end program test_config
