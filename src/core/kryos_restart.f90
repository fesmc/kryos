module kryos_restart
    ! Restart bundles (core-design sections 3.3, 7): the folder of one restart
    ! time, its files and its manifest. The run writes a bundle when its
    ! restart schedule is due (always at time_end) and restores from one when
    ! [kryos] restart names it; the order of the restore, with the exchanges
    ! between the instances, is the run sequence's (kryos_run_sequence).
    !
    !     restart-<time>/                     <time>: year_label, e.g. restart-500
    !       clock_restart.nc                  the clock and the cadence of every instance that steps
    !       <inst>_imports.nc  <inst>_exports.nc      every global instance: its field states
    !       <inst>_*restart*.nc                         a global instance's backend files (stateful backends)
    !       <domain>/
    !         <inst>_imports.nc  <inst>_exports.nc    every instance of the domain
    !         <inst>_*restart*.nc                       its backend files (stateful backends)
    !       manifest.nml                      [bundle] + one group per instance; written last
    !
    ! The field states (fesm_fields_restart_write/read, one file per state)
    ! are the core's: a published import or export survives the restart
    ! without the adapter (adapters-design A5); the backend files are the
    ! adapter's (restart_write/read through the component wrappers, nothing for
    ! a backend without state); a backend declares the files its restart_write
    ! creates (restart_files, <comp>_restart_files), and the writer checks that
    ! each exists. The hub is not stored: it is derived from the
    ! restored ice. Write order = restore order: the globe (forcing, geo),
    ! then each domain's instances in RESTORE_ORDER (ice, geo, climate,
    ! surface, shelf, geothermal, sediments); then the domain's stateful
    ! experiments (core-design 3.9), each a manifest group with component
    ! "experiment" and backend = its kind, its declared files in the domain's
    ! folder. A stateless experiment is not in the bundle: it may be added or
    ! dropped between chunks.
    !
    ! Restore of one instance (restart_read_instance, restart_read_global):
    ! <comp>_restart_read, then its imports and exports from the bundle, then
    ! <comp>_init_state(PASS_RESTART), which re-derives the backend's internal
    ! state from the restored state and puts no export (the exports are the
    ! ones the bundle holds, bit for bit).
    !
    ! manifest.nml (read back with nml): group [bundle] (time, ticks,
    ! ticks_per_year, kryos_git, fesm_utils_git, domains) and one group per
    ! instance, named as the instance (scope, domain, component, backend,
    ! grid, files, git; files: the two field files, then the backend's
    ! declared files). A restore checks the run's instances against it:
    ! the same instances, each with the same scope, domain, component,
    ! backend and grid, or it stops naming the instance (a changed cadence is
    ! allowed and logged by fesm_clock); then that every file the manifest
    ! lists exists in the bundle.
    !
    ! The folders are made with the shell (rm -rf, mkdir -p; one call, as yelmox).

    use precision,        only: dp
    use iso_fortran_env,  only: int64
    use fesm_log,         only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,      only: fields_restart_write, fields_restart_read
    use fesm_clock,       only: clock_t, clock_time, clock_restart_write, clock_restart_read, &
                                cadence_restart_write, cadence_restart_read, cadence_dt_years
    use nml,              only: nml_read, nml_has_group, nml_group_names, nml_group_len
    use kryos_component,  only: component_t, PASS_RESTART, NAME_LEN, PATH_LEN
    use kryos_globe,      only: globe_t, globe_component, globe_pass, globe_restart_write, globe_restart_read, &
                                globe_restart_files
    use kryos_domain,     only: domain_t, domain_component, domain_pass, domain_restart_write, &
                                domain_restart_read, domain_restart_files, RESTORE_ORDER, NCOMP
    use kryos_experiment, only: experiment_t, experiment_list_t, experiments_count, experiments_restart_write, &
                                experiments_restart_read, experiments_restart_files
    use kryos_build_info, only: GIT_KRYOS, LIB_NAME, LIB_GIT

    implicit none

    private
    public :: GLOBE_ORDER, year_label, restart_dir
    public :: restart_write, restart_read_clock, restart_check_manifest
    public :: restart_read_global, restart_read_instance, restart_read_experiments, restart_require_files

    character(len=*), parameter :: CLOCK_FILE    = "clock_restart.nc"
    character(len=*), parameter :: MANIFEST_FILE = "manifest.nml"
    character(len=*), parameter :: BUNDLE_GROUP  = "bundle"

    ! The global instances in write and restore order (core-design 3.3).
    character(len=8), parameter :: GLOBE_ORDER(2) = [character(len=8) :: "forcing", "geo"]

    ! The keys of an instance group the restore checks, in this order.
    character(len=9), parameter :: CHECK_KEYS(5) = [character(len=9) :: "scope", "domain", "component", &
                                                    "backend", "grid"]

    type :: entry_t                                 ! one instance of the run, as the manifest lists it
        character(len=NAME_LEN) :: name = "", scope = "", domain = "", component = "", backend = "", grid = ""
        character(len=PATH_LEN) :: dir = ""         ! its folder in the bundle (the bundle, or <domain>/)
        character(len=64)       :: git = ""
        character(len=PATH_LEN), allocatable :: files(:)
    end type entry_t

contains

    ! =========================================================================
    ! Names
    ! =========================================================================

    function year_label(t) result(s)
        ! `t` [yr] as the shortest decimal that reads back to the same value:
        ! 15000, 0.5, -21000 (restart-<time> folder names, log lines).
        real(dp), intent(in) :: t
        character(len=:), allocatable :: s

        character(len=64) :: buf, fmt
        real(dp) :: x
        integer  :: nd

        buf = ""
        do nd = 0, 17
            write(fmt, '(a,i0,a)') "(f0.", nd, ")"
            write(buf, fmt) t
            read(buf, *) x
            if (x == t) exit
        end do
        s = trim(adjustl(buf))
        if (s(1:1) == ".") s = "0"//s                      ! gfortran f0.d drops the leading zero
        if (len(s) > 1) then
            if (s(1:2) == "-.") s = "-0"//s(2:)
        end if
        if (s(len(s):len(s)) == ".") s = s(1:len(s)-1)
        if (s == "-0") s = "0"
    end function year_label

    function restart_dir(out_dir, t) result(dir)
        ! The bundle folder of time t under out_dir ("./" or ending in "/").
        character(len=*), intent(in) :: out_dir
        real(dp),         intent(in) :: t
        character(len=:), allocatable :: dir
        dir = trim(out_dir)//"restart-"//year_label(t)
    end function restart_dir

    ! =========================================================================
    ! Write
    ! =========================================================================

    subroutine restart_write(dir, clk, globe, dom, exps)
        ! The whole bundle at the clock's time into `dir` (replaced if it
        ! exists): the clock and the cadences, the globe's instances, each
        ! domain's instances in RESTORE_ORDER and its stateful experiments
        ! (exps(d), the list of dom(d)), the manifest last (a bundle
        ! without one is incomplete). The folders are made with one shell
        ! command; after each backend write, every file the backend declared
        ! (restart_files) must exist in its folder.
        character(len=*),       intent(in)    :: dir
        type(clock_t),          intent(in)    :: clk
        type(globe_t),  target, intent(inout) :: globe
        type(domain_t), target, intent(inout) :: dom(:)
        type(experiment_list_t), intent(inout) :: exps(:)

        type(entry_t), allocatable :: ent(:)
        class(component_t), pointer :: c
        character(len=:), allocatable :: ddir, file, cmd
        real(dp) :: t
        integer  :: d, j

        t = clock_time(clk)
        cmd = "rm -rf '"//dir//"' && mkdir -p '"//dir//"'"
        do d = 1, size(dom)
            cmd = cmd//" '"//dir//"/"//trim(dom(d)%name)//"'"
        end do
        call shell(cmd, "the bundle folder "//dir)

        ! clock and cadences
        file = dir//"/"//CLOCK_FILE
        call clock_restart_write(clk, file)
        do j = 1, size(GLOBE_ORDER)
            c => globe_component(globe, trim(GLOBE_ORDER(j)))
            if (cadenced(c, clk)) call cadence_restart_write(c%cadence, file)
        end do
        do d = 1, size(dom)
            do j = 1, NCOMP
                c => domain_component(dom(d), RESTORE_ORDER(j))
                if (cadenced(c, clk)) call cadence_restart_write(c%cadence, file)
            end do
        end do

        ! the instances
        allocate(ent(0))
        do j = 1, size(GLOBE_ORDER)
            c => globe_component(globe, trim(GLOBE_ORDER(j)))
            if (.not. associated(c)) cycle
            ent = [ent, entry_of(c, dir)]
            call backend_write(ent(size(ent)), globe=globe, which=trim(GLOBE_ORDER(j)), t=t)
            call fields_write(c, ent(size(ent)))
        end do
        do d = 1, size(dom)
            ddir = dir//"/"//trim(dom(d)%name)
            do j = 1, NCOMP
                c => domain_component(dom(d), RESTORE_ORDER(j))
                if (.not. associated(c)) cycle
                ent = [ent, entry_of(c, ddir)]
                call backend_write(ent(size(ent)), dom=dom(d), k=RESTORE_ORDER(j), t=t)
                call fields_write(c, ent(size(ent)))
            end do
            call experiments_write(exps(d), dom(d)%name, ddir, t, ent)
        end do

        call manifest_write(dir, clk, dom, ent)
        call log_msg(LOG_INFO, "kryos_restart", "restart bundle "//dir//": "//to_str(size(ent))//" instances")
    end subroutine restart_write

    subroutine experiments_write(list, domain, dir, t, ent)
        ! The stateful experiments of a domain: their restart_write into the
        ! domain's folder (which checks the files they declare), and one
        ! manifest entry each, in list order.
        type(experiment_list_t), intent(inout) :: list
        character(len=*),        intent(in)    :: domain, dir
        real(dp),                intent(in)    :: t
        type(entry_t), allocatable, intent(inout) :: ent(:)

        integer :: i

        call experiments_restart_write(list, dir, t)
        do i = 1, experiments_count(list)
            if (.not. list%s(i)%c%restart) cycle
            ent = [ent, entry_of_experiment(list%s(i)%c, domain, dir)]
            ent(size(ent))%files = experiments_restart_files(list, i)
        end do
    end subroutine experiments_write

    function entry_of_experiment(x, domain, dir) result(e)
        ! The manifest entry of a stateful experiment: component "experiment", backend = its kind.
        class(experiment_t), intent(in) :: x
        character(len=*),    intent(in) :: domain, dir
        type(entry_t) :: e

        e%name      = x%name
        e%scope     = "regional"
        e%domain    = domain
        e%component = "experiment"
        e%backend   = x%kind
        e%grid      = ""
        e%dir       = dir
        e%git       = GIT_KRYOS
        allocate(e%files(0))
    end function entry_of_experiment

    function entry_of(c, dir) result(e)
        ! The manifest entry of an instance written into the folder `dir`.
        class(component_t), intent(in) :: c
        character(len=*),   intent(in) :: dir
        type(entry_t) :: e

        e%name      = c%name
        e%scope     = c%scope
        e%domain    = c%domain
        e%component = c%kind
        e%backend   = c%backend
        e%grid      = c%grid
        e%dir       = dir
        e%git       = backend_git(c)
        allocate(e%files(0))
    end function entry_of

    subroutine backend_write(e, t, globe, which, dom, k)
        ! The backend's restart_write of one instance (globe: `which`; domain:
        ! component k), then the files it declares (restart_files): each must
        ! exist in the instance's folder, else stop naming instance, file and
        ! folder; they are recorded in e.
        type(entry_t),                        intent(inout)           :: e
        real(dp),                             intent(in)              :: t
        type(globe_t),                        intent(inout), optional :: globe
        character(len=*),                     intent(in),    optional :: which
        type(domain_t),                       intent(inout), optional :: dom
        integer,                              intent(in),    optional :: k

        character(len=PATH_LEN), allocatable :: declared(:)

        if (present(globe)) then
            call globe_restart_write(globe, which, trim(e%dir), t)
            declared = globe_restart_files(globe, which)
        else
            call domain_restart_write(dom, k, trim(e%dir), t)
            declared = domain_restart_files(dom, k)
        end if
        call restart_require_files(e%name, e%dir, declared, "its restart_files declares was not created by its "// &
                                   "restart_write (backend "//trim(e%backend)//")")
        e%files = [character(len=PATH_LEN) :: e%files, declared]
    end subroutine backend_write

    subroutine fields_write(c, e)
        ! The imports and exports of an instance: <inst>_imports.nc, <inst>_exports.nc
        ! (recorded in e before the backend's files).
        class(component_t), intent(in)    :: c
        type(entry_t),      intent(inout) :: e

        character(len=PATH_LEN) :: files(2)

        files = [character(len=PATH_LEN) :: trim(c%name)//"_imports.nc", trim(c%name)//"_exports.nc"]
        call fields_restart_write(c%imports, fields_file(e%dir, c%name, "imports"))
        call fields_restart_write(c%exports, fields_file(e%dir, c%name, "exports"))
        e%files = [character(len=PATH_LEN) :: files, e%files]
    end subroutine fields_write

    function fields_file(dir, name, which) result(file)
        character(len=*), intent(in) :: dir, name, which
        character(len=:), allocatable :: file
        file = trim(dir)//"/"//trim(name)//"_"//which//".nc"
    end function fields_file

    subroutine manifest_write(dir, clk, dom, ent)
        ! manifest.nml: [bundle], then one group per instance in write order.
        character(len=*), intent(in) :: dir
        type(clock_t),    intent(in) :: clk
        type(domain_t),   intent(in) :: dom(:)
        type(entry_t),    intent(in) :: ent(:)

        character(len=:), allocatable :: s
        integer :: u, i, d

        open(newunit=u, file=dir//"/"//MANIFEST_FILE, status="replace", action="write")
        write(u, '(a)') "! kryos restart bundle (core-design 7): written by kryos_restart, read back at a restore"
        write(u, '(a)') "&"//BUNDLE_GROUP
        call kv(u, "time",           year_label(clock_time(clk)), quote=.false.)
        call kv(u, "ticks",          int_str(clk%n),              quote=.false.)
        call kv(u, "ticks_per_year", int_str(clk%tpy),            quote=.false.)
        call kv(u, "kryos_git",      GIT_KRYOS)
        call kv(u, "fesm_utils_git", lib_version("fesm-utils"))
        s = ""
        do d = 1, size(dom)
            s = s//' "'//trim(dom(d)%name)//'"'
        end do
        if (size(dom) == 0) s = ' ""'
        call kv(u, "domains", s(2:), quote=.false.)
        write(u, '(a)') "/"
        do i = 1, size(ent)
            write(u, '(a)') ""
            write(u, '(a)') "&"//trim(ent(i)%name)
            call kv(u, "scope",     ent(i)%scope)
            call kv(u, "domain",    ent(i)%domain)
            call kv(u, "component", ent(i)%component)
            call kv(u, "backend",   ent(i)%backend)
            call kv(u, "grid",      ent(i)%grid)
            s = ""
            do d = 1, size(ent(i)%files)
                s = s//' "'//trim(ent(i)%files(d))//'"'
            end do
            call kv(u, "files",     s(2:), quote=.false.)
            call kv(u, "git",       ent(i)%git)
            write(u, '(a)') "/"
        end do
        close(u)

    contains

        subroutine kv(u, key, val, quote)
            integer,          intent(in) :: u
            character(len=*), intent(in) :: key, val
            logical,          intent(in), optional :: quote
            logical :: q
            q = .true.
            if (present(quote)) q = quote
            if (q) then
                write(u, '(a)') "    "//key//repeat(" ", max(0, 14 - len(key)))//' = "'//trim(val)//'"'
            else
                write(u, '(a)') "    "//key//repeat(" ", max(0, 14 - len(key)))//" = "//trim(val)
            end if
        end subroutine kv

    end subroutine manifest_write

    ! =========================================================================
    ! Restore
    ! =========================================================================

    subroutine restart_read_clock(dir, clk, globe, dom)
        ! The clock from the bundle (continuation if the bundle time is
        ! time_init, else a logged branch: fesm_clock), then the cadence of
        ! every instance that steps (continuation: each must have fired at the bundle
        ! time). Before clock_start.
        character(len=*),       intent(in)    :: dir
        type(clock_t),          intent(inout) :: clk
        type(globe_t),  target, intent(inout) :: globe
        type(domain_t), target, intent(inout) :: dom(:)

        class(component_t), pointer :: c
        character(len=:), allocatable :: file
        integer :: d, j

        file = dir//"/"//CLOCK_FILE
        if (.not. exists(file)) then
            call fesm_error("kryos_restart::restart_read_clock", "restart bundle has no "//CLOCK_FILE, &
                            "[kryos] restart = "//dir//"; looked for "//file)
        end if
        call clock_restart_read(clk, file)

        do j = 1, size(GLOBE_ORDER)
            c => globe_component(globe, trim(GLOBE_ORDER(j)))
            if (cadenced(c, clk)) call cadence_restart_read(c%cadence, clk, file)
        end do
        do d = 1, size(dom)
            do j = 1, NCOMP
                c => domain_component(dom(d), RESTORE_ORDER(j))
                if (cadenced(c, clk)) call cadence_restart_read(c%cadence, clk, file)
            end do
        end do
    end subroutine restart_read_clock

    subroutine restart_check_manifest(dir, globe, dom, exps)
        ! The run's instances against the bundle's manifest: every instance of
        ! the run (and every stateful experiment, component "experiment") is
        ! in the bundle with the same scope, domain, component,
        ! backend and grid, and the bundle holds no other instance; else stop
        ! naming the instance. Every file the manifest lists for an instance
        ! exists in its folder, else stop naming instance, file and folder
        ! (before any read). Logs the bundle's time and versions.
        character(len=*),       intent(in) :: dir
        type(globe_t),  target, intent(in) :: globe
        type(domain_t), target, intent(in) :: dom(:)
        type(experiment_list_t), intent(in) :: exps(:)

        character(len=*), parameter :: proc = "kryos_restart::restart_check_manifest"
        character(len=nml_group_len), allocatable :: groups(:)
        type(entry_t), allocatable :: ent(:)
        class(component_t), pointer :: c
        character(len=:), allocatable :: file, hint
        character(len=PATH_LEN) :: val, kgit, ftime
        integer :: d, j, i, q

        file = dir//"/"//MANIFEST_FILE
        if (.not. exists(file)) then
            call fesm_error(proc, "restart bundle has no "//MANIFEST_FILE//" (an incomplete bundle)", &
                            "[kryos] restart = "//dir//"; looked for "//file)
        end if
        if (.not. nml_has_group(file, BUNDLE_GROUP)) then
            call fesm_error(proc, MANIFEST_FILE//" has no group &"//BUNDLE_GROUP, "file "//file)
        end if
        hint = "the bundle "//dir//" was written by a run with other instances; restore it with the par file "// &
               "of that run (a changed cadence is allowed)"

        ! the run's instances
        allocate(ent(0))
        do j = 1, size(GLOBE_ORDER)
            c => globe_component(globe, trim(GLOBE_ORDER(j)))
            if (associated(c)) ent = [ent, entry_of(c, dir)]
        end do
        do d = 1, size(dom)
            do j = 1, NCOMP
                c => domain_component(dom(d), RESTORE_ORDER(j))
                if (associated(c)) ent = [ent, entry_of(c, dir//"/"//trim(dom(d)%name))]
            end do
            do j = 1, experiments_count(exps(d))
                if (.not. exps(d)%s(j)%c%restart) cycle
                ent = [ent, entry_of_experiment(exps(d)%s(j)%c, dom(d)%name, dir//"/"//trim(dom(d)%name))]
            end do
        end do

        do i = 1, size(ent)
            if (.not. nml_has_group(file, trim(ent(i)%name))) then
                call fesm_error(proc, "instance `"//trim(ent(i)%name)//"` ("//trim(ent(i)%component)//", "// &
                                trim(ent(i)%backend)//") is not in the restart bundle", hint)
            end if
            do q = 1, size(CHECK_KEYS)
                val = ""
                call nml_read(file, trim(ent(i)%name), trim(CHECK_KEYS(q)), val)
                if (trim(val) == trim(entry_value(ent(i), q))) cycle
                call fesm_error(proc, "instance `"//trim(ent(i)%name)//"`: "//trim(CHECK_KEYS(q))//" `"// &
                                trim(entry_value(ent(i), q))//"` in the run, `"//trim(val)//"` in the restart bundle", &
                                hint)
            end do
        end do

        ! the bundle's instances
        groups = nml_group_names(file)
        do i = 1, size(groups)
            if (trim(groups(i)) == BUNDLE_GROUP) cycle
            if (any(ent%name == groups(i))) cycle
            call fesm_error(proc, "instance `"//trim(groups(i))//"` of the restart bundle is not in the run", hint)
        end do

        ! the files the manifest lists
        do i = 1, size(ent)
            call check_listed_files(file, ent(i))
        end do

        ftime = "";  kgit = ""
        call nml_read(file, BUNDLE_GROUP, "time", ftime)
        call nml_read(file, BUNDLE_GROUP, "kryos_git", kgit)
        call log_msg(LOG_INFO, "kryos_restart", "restart bundle "//dir//": time "//trim(ftime)//" yr, "// &
                     to_str(size(ent))//" instances as in the run", "written by kryos "//trim(kgit)// &
                     "; this build "//GIT_KRYOS)
    end subroutine restart_check_manifest

    subroutine check_listed_files(manifest, e)
        ! Every file the manifest lists for the instance e exists in its folder.
        character(len=*), intent(in) :: manifest
        type(entry_t),    intent(in) :: e

        call restart_require_files(e%name, e%dir, manifest_files(manifest, trim(e%name)), &
                                   "the manifest lists is not in the restart bundle (an incomplete or damaged bundle)")
    end subroutine check_listed_files

    subroutine restart_require_files(instance, folder, files, what)
        ! Every file of `files` exists in `folder`, else stop naming the
        ! instance, the file and the folder; `what` completes "the file `f` ...".
        character(len=*),        intent(in) :: instance, folder, what
        character(len=PATH_LEN), intent(in) :: files(:)

        character(len=:), allocatable :: f
        integer :: i

        do i = 1, size(files)
            f = trim(files(i))
            if (exists(trim(folder)//"/"//f)) cycle
            call fesm_error("kryos_restart::restart_require_files", "instance `"//trim(instance)//"`: the file `"// &
                            f//"` "//what, "folder "//trim(folder))
        end do
    end subroutine restart_require_files

    function manifest_files(manifest, group) result(files)
        ! The `files` list of a manifest group, the quoted names of the line
        ! manifest_write wrote (not through nml: its string vector is 256
        ! characters long, a list of many backend files is longer).
        character(len=*), intent(in) :: manifest, group
        character(len=PATH_LEN), allocatable :: files(:)

        character(len=8192) :: line
        integer :: u, ios, p, i, i0
        logical :: in_group

        allocate(files(0))
        in_group = .false.
        open(newunit=u, file=manifest, status="old", action="read")
        do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (line(1:1) == "&") then
                in_group = trim(line(2:)) == trim(group)
            else if (line(1:1) == "/") then
                in_group = .false.
            else if (in_group) then
                p = index(line, "=")
                if (p == 0) cycle
                if (trim(line(1:p-1)) /= "files") cycle
                i0 = 0
                do i = p + 1, len_trim(line)
                    if (line(i:i) /= '"') cycle
                    if (i0 == 0) then
                        i0 = i + 1                                  ! opening quote
                    else
                        files = [character(len=PATH_LEN) :: files, line(i0:i-1)]
                        i0 = 0                                      ! closing quote
                    end if
                end do
                exit
            end if
        end do
        close(u)
    end function manifest_files

    function entry_value(e, q) result(s)
        ! The value of CHECK_KEYS(q) of an entry.
        type(entry_t), intent(in) :: e
        integer,       intent(in) :: q
        character(len=NAME_LEN) :: s
        select case (q)
        case (1); s = e%scope
        case (2); s = e%domain
        case (3); s = e%component
        case (4); s = e%backend
        case (5); s = e%grid
        end select
    end function entry_value

    subroutine restart_read_global(dir, globe, which, clk)
        ! Restore of the global instance `which` from the bundle folder:
        ! restart_read, its imports and exports, PASS_RESTART.
        character(len=*),      intent(in)    :: dir
        type(globe_t), target, intent(inout) :: globe
        character(len=*),      intent(in)    :: which
        type(clock_t),         intent(in)    :: clk

        class(component_t), pointer :: c

        c => globe_component(globe, which)
        if (.not. associated(c)) return
        call globe_restart_read(globe, which, dir, clock_time(clk))
        call fields_read(c, dir)
        call globe_pass(globe, which, clk, PASS_RESTART)
    end subroutine restart_read_global

    subroutine restart_read_instance(dir, dom, k, clk)
        ! Restore of component k of a domain from the bundle (its folder
        ! dir/<domain>): restart_read, its imports and exports, PASS_RESTART.
        character(len=*),       intent(in)    :: dir
        type(domain_t), target, intent(inout) :: dom
        integer,                intent(in)    :: k
        type(clock_t),          intent(in)    :: clk

        class(component_t), pointer :: c
        character(len=:), allocatable :: ddir

        c => domain_component(dom, k)
        if (.not. associated(c)) return
        ddir = dir//"/"//trim(dom%name)
        call domain_restart_read(dom, k, ddir, clock_time(clk))
        call fields_read(c, ddir)
        call domain_pass(dom, k, clk, PASS_RESTART)
    end subroutine restart_read_instance

    subroutine restart_read_experiments(dir, dom, list, clk)
        ! Restore of the stateful experiments of a domain from the bundle
        ! (their files in dir/<domain>); the stateless ones have none.
        character(len=*),        intent(in)    :: dir
        type(domain_t),          intent(in)    :: dom
        type(experiment_list_t), intent(inout) :: list
        type(clock_t),           intent(in)    :: clk

        call experiments_restart_read(list, dir//"/"//trim(dom%name), clock_time(clk))
    end subroutine restart_read_experiments

    subroutine fields_read(c, dir)
        ! The imports and exports of an instance as the bundle holds them.
        class(component_t), intent(inout) :: c
        character(len=*),   intent(in)    :: dir
        call fields_restart_read(c%imports, fields_file(dir, c%name, "imports"))
        call fields_restart_read(c%exports, fields_file(dir, c%name, "exports"))
    end subroutine fields_read

    ! =========================================================================
    ! Helpers
    ! =========================================================================

    logical function cadenced(c, clk)
        ! An instance whose cadence steps (not off): its last mark is state
        ! the bundle carries. An off cadence (a static instance) never marks;
        ! it is neither written nor read (its mark stays at time_init).
        class(component_t), pointer, intent(in) :: c
        type(clock_t),               intent(in) :: clk
        cadenced = associated(c)
        if (cadenced) cadenced = cadence_dt_years(c%cadence, clk) > 0.0_dp
    end function cadenced

    function backend_git(c) result(git)
        ! The version of the code a backend runs: the library of its build
        ! switch (constraints build_flag "yelmo=1" -> yelmo), else kryos (the
        ! adapters and backends compiled into kryos).
        class(component_t), intent(in) :: c
        character(len=64) :: git
        integer :: p
        git = GIT_KRYOS
        p = index(c%con%build_flag, "=")
        if (p > 1) git = lib_version(c%con%build_flag(1:p-1))
    end function backend_git

    function lib_version(name) result(v)
        ! The build's version of the library `name` (kryos_build_info).
        character(len=*), intent(in) :: name
        character(len=64) :: v
        integer :: i
        v = "unknown"
        do i = 1, size(LIB_NAME)
            if (trim(LIB_NAME(i)) == trim(name)) v = LIB_GIT(i)
        end do
    end function lib_version

    subroutine shell(cmd, what)
        ! A shell command that must succeed.
        character(len=*), intent(in) :: cmd, what
        integer :: es, cs
        call execute_command_line(cmd, exitstat=es, cmdstat=cs)
        if (cs /= 0 .or. es /= 0) then
            call fesm_error("kryos_restart", "shell command failed for "//what, cmd)
        end if
    end subroutine shell

    logical function exists(file)
        character(len=*), intent(in) :: file
        inquire(file=file, exist=exists)
    end function exists

    function int_str(n) result(s)
        integer(int64), intent(in) :: n
        character(len=:), allocatable :: s
        character(len=24) :: buf
        write(buf, '(i0)') n
        s = trim(buf)
    end function int_str

end module kryos_restart
