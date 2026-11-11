module kryos_orbit
    ! The orbit of the run ([orbit]): the orbital elements that apply at a
    ! time of the forcing time. The forcing resolves them (tsgen, key `orbit`)
    ! and exports them to every component that takes the orbit (the climate,
    ! the surface's insolation), so the run has one orbit; insol (fesmc/insol)
    ! holds the tables. Each consumer computes its own insolation, if any, from
    ! the elements, on its own grid; it recomputes only when
    ! the elements change. This module only resolves the elements:
    !
    !   method = "forcing_time": the elements at the forcing time time_rel (Laskar et al. 2004 tables)
    !            "fixed_time":   the elements at a fixed time [yr, model time: relative to 1950, past < 0]
    !            "fixed_orbit":  the given elements ecc, obl [deg], per [deg] (longitude of perihelion)
    !
    ! Read with input/kryos_defaults.nml (template group orbit): a par file lists only what differs.

    use nml,        only : nml_read, nml_validate
    use insolation, only : dp, insol_class, orbit_pars_class, orbit_pars_at, orbit_pars_fixed

    implicit none

    private

    ! Defaults (schema) of [orbit]: template group orbit of the kryos defaults file
    ! (yelmox: input/yelmox_defaults.nml)
    character(len=*), parameter :: def_file = "input/kryos_defaults.nml"

    type orbit_ctl_class
        character(len=16)  :: method    ! forcing_time | fixed_time | fixed_orbit
        real(dp)           :: time      ! [yr] model time of the orbit (fixed_time)
        real(dp)           :: ecc       ! [1]   eccentricity (fixed_orbit)
        real(dp)           :: obl       ! [deg] obliquity (fixed_orbit)
        real(dp)           :: per       ! [deg] longitude of perihelion (fixed_orbit)
        character(len=512) :: fldr      ! folder of the orbital tables
    end type orbit_ctl_class

    public :: orbit_ctl_class
    public :: orbit_pars_class
    public :: orbit_ctl_load
    public :: orbit_pars_now
    public :: orbit_pars_equal

contains

    subroutine orbit_ctl_load(ctl, filename, group)
        ! Read [orbit<sfx>] (optional group; the defaults without it).
        type(orbit_ctl_class), intent(out) :: ctl
        character(len=*),      intent(in)  :: filename, group

        type(orbit_pars_class) :: orbit

        call nml_validate(filename, def_file, group, defaults_group="orbit")
        call nml_read(filename, group, "method", ctl%method, defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "time",   ctl%time,   defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "ecc",    ctl%ecc,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "obl",    ctl%obl,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "per",    ctl%per,    defaults_file=def_file, defaults_group="orbit")
        call nml_read(filename, group, "fldr",   ctl%fldr,   defaults_file=def_file, defaults_group="orbit")

        select case(trim(ctl%method))
            case("forcing_time", "fixed_time")
                ! The tables are read by each host's insol_init.
            case("fixed_orbit")
                orbit = orbit_pars_fixed(ctl%ecc, ctl%obl, ctl%per)   ! stops on invalid elements
            case default
                write(*,*) "orbit_ctl_load:: error: ["//trim(group)//"] method must be &
                           &forcing_time, fixed_time or fixed_orbit; got '"//trim(ctl%method)//"'."
                error stop 1
        end select
    end subroutine orbit_ctl_load

    function orbit_pars_now(ctl, is, time_rel) result(orbit)
        ! The orbital elements at the forcing time time_rel [yr], from the
        ! tables of the host's insol_class (forcing_time, fixed_time).
        type(orbit_ctl_class), intent(in) :: ctl
        type(insol_class),     intent(in) :: is
        real(dp),              intent(in) :: time_rel
        type(orbit_pars_class) :: orbit

        select case(trim(ctl%method))
            case("forcing_time")
                orbit = orbit_pars_at(is, time_rel)
            case("fixed_time")
                orbit = orbit_pars_at(is, ctl%time)
            case default   ! fixed_orbit
                orbit = orbit_pars_fixed(ctl%ecc, ctl%obl, ctl%per)
        end select
    end function orbit_pars_now

    elemental function orbit_pars_equal(a, b) result(same)
        ! The same elements (the insolation need not be recomputed).
        type(orbit_pars_class), intent(in) :: a, b
        logical :: same

        same = (a%ecc == b%ecc .and. a%obl == b%obl .and. a%per == b%per)
    end function orbit_pars_equal

end module kryos_orbit
