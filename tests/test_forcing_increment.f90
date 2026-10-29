program test_forcing_increment
    ! Forcing-increment restart bundles (phase 3, WP4; core-design section 7). A run
    ! of the forcing alone (no domain, as test_run_clock) with a ramp on alpha and
    ! restart_every_df = 0.25 writes a bundle restart-<n> each time alpha has moved
    ! by more than that since its last mark (yelmox coefs_restart_due, the counter of
    ! coefs_restart_fldr): at 30, 60 and 90 yr for alpha = t/100, besides the bundle
    ! of the restart schedule at time_end. A run restored from restart-1 (30 yr)
    ! continues bit for bit: the series of alpha, the bundle of the schedule at the
    ! end, and its own bundles are the original's, numbered on from the count in the
    ! bundle (restart-2 and restart-3). Restored in the out_dir of a run that stopped
    ! after restart-1 it writes restart-2 and restart-3 there, equal to those of the
    ! uninterrupted run; restored where restart-2 exists it stops (a bundle is never
    ! overwritten). Without restart_every_df no restart-<n> is written;
    ! restart_every_df on a derived coefficient stops (child processes of this
    ! executable).
    !
    ! Files go to tmp/test_forcing_increment/ (run from the repo root, as make check
    ! does).

    use precision,   only: sp, dp
    use fesm_log,    only: log_set_unit
    use ncio,        only: nc_read, nc_size
    use nml,         only: nml_read
    use kryos_run,   only: run_t, run_init, run_step, run_done, run_end

    implicit none

    character(len=*), parameter :: dir = "tmp/test_forcing_increment/"
    integer :: nfail = 0
    character(len=64) :: arg

    call get_command_argument(1, arg)
    if (len_trim(arg) > 0) then                 ! child process: a stop
        select case (trim(arg))
        case ("derived")
            call write_par(dir//"derived.nml", dir//"derived/", "None", 0.0_dp, 0.25_dp, derived=.true.)
            call run_all(dir//"derived.nml")
        case ("overwrite")                      ! restored from restart-1 of a, into a, where restart-2 exists
            call write_par(dir//"ow.nml", dir//"a/", dir//"a/restart-1", 30.0_dp, 0.25_dp)
            call run_all(dir//"ow.nml")
        end select
        write(*,'(a)') "SURVIVED"
        stop 0
    end if

    call execute_command_line("rm -rf "//dir//" && mkdir -p "//dir)

    call test_bundles()
    call test_restored()
    call test_same_dir()
    call test_off()
    call test_stop("derived", "restart_every_df > 0 on a derived coefficient", &
                   "derived: restart_every_df > 0 on a derived coefficient stops")
    call test_stop("overwrite", "forcing-increment bundle", &
                   "overwrite: a restored run never overwrites a bundle (stops where restart-2 exists)")

    if (nfail > 0) then
        write(*,'(a,i0,a)') "test_forcing_increment: FAILED (", nfail, " checks)"
        error stop 1
    end if
    write(*,'(a)') "test_forcing_increment: OK"

contains

    subroutine check(ok, label)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: label
        if (ok) then
            write(*,'(a)') "  ok    "//label
        else
            write(*,'(a)') "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    logical function exists(file)
        character(len=*), intent(in) :: file
        inquire(file=file, exist=exists)
    end function exists

    subroutine put(file, lines)
        character(len=*), intent(in) :: file, lines(:)
        integer :: u, i
        open(newunit=u, file=file, status="replace", action="write")
        do i = 1, size(lines)
            write(u,'(a)') trim(lines(i))
        end do
        close(u)
    end subroutine put

    function num(x) result(s)
        real(dp), intent(in) :: x
        character(len=:), allocatable :: s
        character(len=64) :: buf
        write(buf,'(f0.6)') x
        s = trim(adjustl(buf))
        if (s(1:1) == ".") s = "0"//s
    end function num

    subroutine write_par(file, out_dir, restart, time_init, every_df, derived, time_end)
        ! A zero-domain run, 0 .. 100 yr at dtt = 10: [globe] forcing, tsgen with alpha
        ! (ramp-time, alpha = t/100) and, with `derived`, a derived coefficient carrying
        ! restart_every_df.
        character(len=*), intent(in) :: file, out_dir, restart
        real(dp),         intent(in) :: time_init, every_df
        logical,          intent(in), optional :: derived
        real(dp),         intent(in), optional :: time_end

        logical :: drv
        real(dp) :: t_end
        drv = .false.
        if (present(derived)) drv = derived
        t_end = 100.0_dp
        if (present(time_end)) t_end = time_end
        call put(file, [character(len=200) :: &
            '&kryos', &
            '    out_dir        = "'//trim(out_dir)//'"', &
            '    restart        = "'//trim(restart)//'"', &
            '    restart_out    = "out_rst"', &
            '    write_lossless = True', &
            '    log_level      = "warn"', &
            '/', &
            '&clock', &
            '    method     = "rel"', &
            '    time_init  = '//num(time_init), &
            '    time_end   = '//num(t_end), &
            '    dtt        = 10.0', &
            '/', &
            '&out_1D',  '    method = "const"', '    dt = 10.0', '/', &
            '&out_rst', '    method = "const"', '    dt = 100.0', '/', &
            '&globe', '    forcing = "forcing"', '/', &
            '&forcing', '    backend = "tsgen"', '    clock = "forcing"', &
            merge('    coefs   = "alpha" "dTo"', '    coefs   = "alpha"      ', drv), '    out_1D  = "out_1D"', '/', &
            '&coef_alpha', '    method           = "ramp-time"', '    f_min            = 0.0', &
            '    f_max            = 1.0', '    dt_ramp          = 100.0', '    df_sign          = 1.0', &
            '    restart_every_df = '//merge(num(0.0_dp), num(every_df), drv), '/', &
            '&coef_dTo', '    method           = "derived"', '    from             = "air_anomaly"', &
            '    units            = "K"', '    restart_every_df = '//merge(num(every_df), num(0.0_dp), drv), '/'])
    end subroutine write_par

    subroutine run_all(par_file)
        character(len=*), intent(in) :: par_file
        type(run_t), target :: run
        call run_init(run, par_file)
        do while (.not. run_done(run))
            call run_step(run)
        end do
        call run_end(run)
    end subroutine run_all

    function series(file) result(f)
        character(len=*), intent(in) :: file
        real(sp), allocatable :: f(:)
        allocate(f(nc_size(file, "time")))
        call nc_read(file, "coef_alpha", f)
    end function series

    function bundle_time(folder) result(t)
        character(len=*), intent(in) :: folder
        real(dp) :: t
        call nml_read(folder//"/manifest.nml", "bundle", "time", t)
    end function bundle_time

    logical function same_bundle(a, b)
        ! Every file of the two bundles but the manifest (it holds the time of writing) equal.
        character(len=*), intent(in) :: a, b
        integer :: es
        call execute_command_line("diff -rq --exclude=manifest.nml "//a//" "//b//" > /dev/null", exitstat=es)
        same_bundle = es == 0
    end function same_bundle

    subroutine test_bundles()
        real(sp), allocatable :: f(:)
        character(len=*), parameter :: a = dir//"a/"

        call write_par(dir//"a.nml", a, "None", 0.0_dp, 0.25_dp)
        call run_all(dir//"a.nml")
        f = series(a//"forcing_ts.nc")
        call check(size(f) == 11 .and. abs(f(11) - 1.0_sp) < 1.0e-6_sp .and. abs(f(4) - 0.3_sp) < 1.0e-6_sp, &
                   "bundles: alpha ramps 0 .. 1 over 100 yr")
        call check(exists(a//"restart-1/manifest.nml") .and. exists(a//"restart-2/manifest.nml") .and. &
                   exists(a//"restart-3/manifest.nml") .and. .not. exists(a//"restart-4"), &
                   "bundles: restart-1, -2, -3 written (alpha moved by more than 0.25 at 30, 60, 90), no restart-4")
        call check(exists(a//"restart-100/manifest.nml"), "bundles: the bundle of the restart schedule at time_end besides")
        call check(bundle_time(a//"restart-1") == 30.0_dp .and. bundle_time(a//"restart-2") == 60.0_dp .and. &
                   bundle_time(a//"restart-3") == 90.0_dp, "bundles: restart-<n> hold the times 30, 60, 90")
        call check(exists(a//"restart-2/clock_restart.nc") .and. exists(a//"restart-2/forcing_alpha_restart.nc") .and. &
                   exists(a//"restart-2/forcing_exports.nc"), "bundles: a complete bundle (clock, forcing state and fields)")
    end subroutine test_bundles

    subroutine test_restored()
        ! Restored from restart-1 (30 yr), the run continues bit for bit.
        character(len=*), parameter :: a = dir//"a/", b = dir//"b/"
        real(sp), allocatable :: fa(:), fb(:)

        call write_par(dir//"b.nml", b, a//"restart-1", 30.0_dp, 0.25_dp)
        call run_all(dir//"b.nml")
        fa = series(a//"forcing_ts.nc")
        fb = series(b//"forcing_ts.nc")
        call check(size(fb) == 8 .and. all(fb == fa(4:11)), "restored: alpha of the restored run equals the original's from 30 yr on")
        call check(same_bundle(a//"restart-100", b//"restart-100"), "restored: the bundle at time_end equals the original's")
        call check(.not. exists(b//"restart-1") .and. exists(b//"restart-2/manifest.nml") .and. &
                   exists(b//"restart-3/manifest.nml") .and. .not. exists(b//"restart-4"), &
                   "restored: the restored run numbers on from the count in the bundle: restart-2 and restart-3 (60 and 90 yr)")
        call check(bundle_time(b//"restart-2") == 60.0_dp .and. bundle_time(b//"restart-3") == 90.0_dp, &
                   "restored: at the times of the original's restart-2 and restart-3")
        call check(same_bundle(a//"restart-2", b//"restart-2") .and. same_bundle(a//"restart-3", b//"restart-3"), &
                   "restored: those bundles equal the original's")
    end subroutine test_restored

    subroutine test_same_dir()
        ! A run that stopped at 40 yr (after restart-1), restored from its restart-1 in its own
        ! out_dir: it continues with restart-2 and restart-3, equal to the uninterrupted run's.
        character(len=*), parameter :: a = dir//"a/", s = dir//"s/"

        call write_par(dir//"s1.nml", s, "None", 0.0_dp, 0.25_dp, time_end=40.0_dp)
        call run_all(dir//"s1.nml")
        call check(exists(s//"restart-1/manifest.nml") .and. .not. exists(s//"restart-2"), "same out_dir: the first part wrote restart-1 only")
        call write_par(dir//"s2.nml", s, s//"restart-1", 30.0_dp, 0.25_dp)
        call run_all(dir//"s2.nml")
        call check(exists(s//"restart-1/manifest.nml") .and. exists(s//"restart-2/manifest.nml") .and. &
                   exists(s//"restart-3/manifest.nml") .and. .not. exists(s//"restart-4"), &
                   "same out_dir: the restored run continues with restart-2 and restart-3")
        call check(same_bundle(a//"restart-1", s//"restart-1") .and. same_bundle(a//"restart-2", s//"restart-2") .and. &
                   same_bundle(a//"restart-3", s//"restart-3"), "same out_dir: restart-1, -2, -3 equal the uninterrupted run's")
    end subroutine test_same_dir

    subroutine test_off()
        ! restart_every_df = 0 (the shipped configurations): only the schedule's bundle.
        character(len=*), parameter :: c = dir//"c/"
        integer :: es

        call write_par(dir//"c.nml", c, "None", 0.0_dp, 0.0_dp)
        call run_all(dir//"c.nml")
        call execute_command_line("ls "//c//" | grep -q '^restart-[0-9]*$' && test `ls "//c//" | grep -c '^restart-'` -eq 1", &
                                  exitstat=es)
        call check(exists(c//"restart-100/manifest.nml") .and. es == 0, &
                   "off: without restart_every_df only restart-100 is written")
    end subroutine test_off

    subroutine test_stop(case_name, text, label)
        ! A child process of this executable must print `text` on stderr, exit with 1 and not
        ! reach its last line.
        character(len=*), intent(in) :: case_name, text, label
        character(len=512) :: exe, line
        character(len=:), allocatable :: err, out
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > "//dir//"child.out 2> "//dir//"child.err", exitstat=es)
        err = ""; out = ""
        open(newunit=u, file=dir//"child.err", status="old", action="read", iostat=ios)
        do while (ios == 0)
            read(u, "(a)", iostat=ios) line
            if (ios == 0) err = err//trim(line)//" "
        end do
        close(u)
        open(newunit=u, file=dir//"child.out", status="old", action="read", iostat=ios)
        do while (ios == 0)
            read(u, "(a)", iostat=ios) line
            if (ios == 0) out = out//trim(line)//" "
        end do
        close(u)
        call check(es == 1 .and. index(err, text) > 0 .and. index(out, "SURVIVED") == 0, label)
    end subroutine test_stop

end program test_forcing_increment
