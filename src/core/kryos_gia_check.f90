module kryos_gia_check
    ! The init-time contract checks of the VILMA2 geo backend that need the whole run
    ! (phase 4 plan 3.2): which instances may be combined and what the ice instances
    ! must offer. Called by run_init after the instances exist (kryos_run), before
    ! anything is exchanged; every stop names the instance and the way out. Pure
    ! checks: nothing is changed.
    !
    ! Not here: the footprints of the domains must not overlap (checked where they are
    ! merged, kryos_globe_ice).

    use precision,       only: dp
    use fesm_log,        only: fesm_error, to_str
    use kryos_config,    only: config_t, config_has, config_read
    use kryos_geo,       only: geo_present
    use kryos_ice,       only: ice_present
    use kryos_globe,     only: globe_t
    use kryos_domain,    only: domain_t, C_GEO

    implicit none

    private
    public :: gia_check

    character(len=*), parameter :: VILMA2 = "vilma2", FASTISO = "fastisostasy"

contains

    subroutine gia_check(globe, dom, cfg)
        ! `vilma2` on the globe and in the domains go together, and not with FastIsostasy (a
        ! later package); with the global part every domain's ice offers `flotation`, `footprint`
        ! and `H_ice_pd`; a port's `dt`, when set, is the global part's.
        type(globe_t),  intent(in) :: globe
        type(domain_t), intent(in) :: dom(:)
        type(config_t), intent(in) :: cfg

        character(len=*), parameter :: proc = "kryos_gia_check::gia_check"
        logical :: glob_v, glob_f, any_port, any_f
        integer :: d

        glob_v = .false.; glob_f = .false.
        if (geo_present(globe%geo)) then
            glob_v = trim(globe%geo%c%backend) == VILMA2
            glob_f = trim(globe%geo%c%backend) == FASTISO
        end if
        any_port = .false.; any_f = .false.
        do d = 1, size(dom)
            any_port = any_port .or. trim(dom(d)%backend(C_GEO)) == VILMA2
            any_f    = any_f    .or. trim(dom(d)%backend(C_GEO)) == FASTISO
        end do
        if (globe%has_gi .and. .not. globe%host_gi) then       ! the observing hub: the ice offers what the gather reads
            do d = 1, size(dom)
                if (.not. ice_present(dom(d)%ice)) then
                    call fesm_error(proc, "domain `"//trim(dom(d)%name)//"` has no ice instance", &
                                    "[globe_ice] gathers the ice of every domain")
                end if
                call check_ice(dom(d))
            end do
        end if
        if (.not. (glob_v .or. any_port)) return

        if ((glob_v .and. any_f) .or. (any_port .and. glob_f)) then
            call fesm_error(proc, "FastIsostasy together with VILMA2 is not implemented", &
                            "phase 4, later package (VILMA2 change 5): use `geo = vilma2` or `geo = fastisostasy` "// &
                            "for the global part and every domain")
        end if

        if (any_port .and. .not. glob_v) then
            do d = 1, size(dom)
                if (trim(dom(d)%backend(C_GEO)) /= VILMA2) cycle
                call fesm_error(proc, "geo instance `"//trim(dom(d)%inst(C_GEO))//"` of domain `"//trim(dom(d)%name)// &
                                "` has backend vilma2, but the run has no global geo with vilma2", &
                                "[globe] geo names the global part (backend vilma2: the Earth the ports belong to)")
            end do
        end if

        if (glob_v) call check_spinup(globe, cfg)
        do d = 1, size(dom)
            call check_domain(globe, dom(d), cfg)
        end do
    end subroutine gia_check

    subroutine check_spinup(globe, cfg)
        ! The Earth is spun up under the t0 load before its seed (plan 7.4, F7: the dt = 0 seed solves with the
        ! previous load, so without the spin-up the t0 load, the external ice above all, enters with the first
        ! call): stops unless the instance sets `allow_no_spinup = True` (tests, present-day starts).
        type(globe_t),  intent(in) :: globe
        type(config_t), intent(in) :: cfg

        logical :: allow

        allow = .false.
        call config_read(cfg, trim(globe%geo%c%name), "allow_no_spinup", allow, "instance", "geo_vilma2")
        if (globe%geo%c%opt%host_spinup .or. allow) return
        call fesm_error("kryos_gia_check::gia_check", "global geo `"//trim(globe%geo%c%name)// &
                        "` has no VILMA2 spin-up before its seed", &
                        "set `pre_spinup_1d = .true.` or `equil_time_max > 0` in the &vilma group (plan 7.4): the seed "// &
                        "solves the previous load, so the t0 load would enter with the first call; "// &
                        "`allow_no_spinup = True` in the instance group accepts that (tests, present-day starts)")
    end subroutine check_spinup

    subroutine check_domain(globe, dom, cfg)
        type(globe_t),  intent(in) :: globe
        type(domain_t), intent(in) :: dom
        type(config_t), intent(in) :: cfg

        character(len=*), parameter :: proc = "kryos_gia_check::gia_check"
        character(len=:), allocatable :: who
        real(dp) :: dt_glob, dt_port

        who = "domain `"//trim(dom%name)//"`"
        if (len_trim(dom%inst(C_GEO)) == 0) then
            call fesm_error(proc, who//" has no geo instance", &
                            "the global geo is vilma2: every domain needs a geo instance with backend vilma2 (its port)")
        else if (trim(dom%backend(C_GEO)) /= VILMA2) then
            call fesm_error(proc, "geo instance `"//trim(dom%inst(C_GEO))//"` of "//who//" has backend "// &
                            trim(dom%backend(C_GEO)), "the global geo is vilma2: every domain's geo must be "// &
                            "vilma2 as well (the regional port of the Earth)")
        end if

        if (.not. ice_present(dom%ice)) then
            call fesm_error(proc, who//" has no ice instance", &
                            "the global geo is vilma2: its load is gathered from every domain's ice")
        end if
        call check_ice(dom)

        if (geo_present(globe%geo) .and. config_has(cfg, trim(dom%inst(C_GEO)), "dt")) then
            call config_read(cfg, trim(dom%inst(C_GEO)), "dt", dt_port, "instance")
            call config_read(cfg, trim(globe%geo%c%name), "dt", dt_glob, "instance")
            if (dt_port /= dt_glob) then
                call fesm_error(proc, "dt of geo instance `"//trim(dom%inst(C_GEO))//"` ("//to_str(dt_port)// &
                                " yr) differs from the global geo `"//trim(globe%geo%c%name)//"` ("//to_str(dt_glob)// &
                                " yr)", "a port steps with the Earth it belongs to: leave `dt` out of its group "// &
                                "or give it the global part's value")
            end if
        end if
    end subroutine check_domain

    subroutine check_ice(dom)
        ! What the gather reads of a domain's ice: the flotation binding, the static footprint and H_ice_pd.
        type(domain_t), intent(in) :: dom

        character(len=*), parameter :: proc = "kryos_gia_check::gia_check"
        character(len=:), allocatable :: who

        who = "domain `"//trim(dom%name)//"`"
        associate (ice => dom%ice%c)
            if (.not. ice%caps%has_flotation) then
                call fesm_error(proc, "ice instance `"//trim(ice%name)//"` ("//trim(ice%backend)//") of "//who// &
                                " has no flotation binding", "the host-resolved load needs the ice model's own "// &
                                "flotation rule outside its footprint (capability has_flotation)")
            end if
            if (.not. ice%opt%footprint) then
                call fesm_error(proc, "ice instance `"//trim(ice%name)//"` ("//trim(ice%backend)//") of "//who// &
                                " does not export `footprint`", "the host-resolved load needs the static footprint "// &
                                "of the ice model (where it holds ice)")
            end if
            if (.not. ice%opt%H_ice_pd) then
                call fesm_error(proc, "ice instance `"//trim(ice%name)//"` ("//trim(ice%backend)//") of "//who// &
                                " does not export `H_ice_pd`", "the reference of VILMA2 is built from the present-day "// &
                                "ice of every domain")
            end if
        end associate
    end subroutine check_ice

end module kryos_gia_check
