module kryos_snap_coefs
    ! The forcing scalars of the snapshots backend without the registry
    ! (adapters-design 4.1): a private table of the coefficients its terms
    ! name, with the calls yelmox climate_terms and yelmox_climate make on the
    ! driver's coefs_class (coefs_has, coefs_is_derived, coefs_derived,
    ! coefs_value, coefs_set). Filled at init from the [coef_<name>] groups
    ! (as yelmox kryos_forcing coefs_init) and, at each update, from the
    ! climate's imports coef_<name>; `f` is sp as coef_class%f is in yelmox.

    use precision,       only: wp
    use fesm_log,        only: fesm_error
    use nml,             only: nml_read, nml_validate
    use kryos_snap_nml,  only: COEF_DEFAULTS, FILE_LEN, group_file

    implicit none

    private

    ! The definition of a derived coefficient ([coef_<name>] method = "derived").
    type, public :: derived_spec
        character(len=32) :: from    = ""        ! air_anomaly | obm_north | obm_south
        real(wp)          :: factor  = 1.0_wp    ! air_anomaly: value = factor * mean anomaly
        character(len=64) :: region  = "all"     ! air_anomaly: cells of the mean ([regions] selection)
        logical           :: anomaly = .true.    ! obm_*: the box temperature minus its initial value
    end type derived_spec

    type, public :: snap_coef_t
        character(len=64)  :: name  = ""
        character(len=32)  :: units = "1"
        logical            :: derived = .false.
        type(derived_spec) :: drv
        real(wp)           :: f = 0.0_wp
    end type snap_coef_t

    type, public :: snap_coef_table_t
        type(snap_coef_t), allocatable :: coef(:)
    end type snap_coef_table_t

    public :: coefs_add, coefs_has, coefs_is_derived, coefs_derived, coefs_value, coefs_set
    public :: coefs_position

contains

    subroutine coefs_add(table, par_file, library, name)
        ! Add coefficient `name` (once) from its group [coef_<name>]: the par
        ! file, else the first library file that has it; keys not set from the
        ! template `coef` of input/forcing_defaults.nml.
        type(snap_coef_table_t), intent(inout) :: table
        character(len=*),        intent(in)    :: par_file, library(:), name

        type(snap_coef_t) :: c
        character(len=FILE_LEN) :: gfile
        character(len=128) :: group
        character(len=32)  :: method

        if (.not. allocated(table%coef)) allocate(table%coef(0))
        if (coefs_has(table, name)) return

        c%name = trim(name)
        group  = "coef_"//trim(name)
        gfile  = group_file(par_file, library, trim(group), "a coefficient named by a term")
        call nml_validate(trim(gfile), COEF_DEFAULTS, trim(group), defaults_group="coef")
        call nml_read(trim(gfile), trim(group), "method", method, init=.true., &
                      defaults_file=COEF_DEFAULTS, defaults_group="coef")
        call nml_read(trim(gfile), trim(group), "units", c%units, init=.true., &
                      defaults_file=COEF_DEFAULTS, defaults_group="coef")
        if (trim(method) == "derived") then
            c%derived = .true.
            call nml_read(trim(gfile), trim(group), "from",    c%drv%from,    init=.true., &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "factor",  c%drv%factor, &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "region",  c%drv%region,  init=.true., &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            call nml_read(trim(gfile), trim(group), "anomaly", c%drv%anomaly, &
                          defaults_file=COEF_DEFAULTS, defaults_group="coef")
            select case (trim(c%drv%from))
            case ("air_anomaly", "obm_north", "obm_south")
            case default
                call fesm_error("kryos_snap_coefs::coefs_add", "["//trim(group)//"] from = '"// &
                                trim(c%drv%from)//"': derived coefficients come from air_anomaly, "// &
                                "obm_north or obm_south")
            end select
        end if
        table%coef = [table%coef, c]
    end subroutine coefs_add

    subroutine coefs_set(table, name, value)
        ! Set the value of the derived coefficient `name` (its producer, once
        ! per update); stops on an unknown or a tsgen coefficient.
        type(snap_coef_table_t), intent(inout) :: table
        character(len=*),        intent(in)    :: name
        real(wp),                intent(in)    :: value
        integer :: k
        k = coef_checked(table, name, "coefs_set")
        if (.not. table%coef(k)%derived) then
            call fesm_error("kryos_snap_coefs::coefs_set", "coefficient '"//trim(name)// &
                            "' is a tsgen coefficient, not derived")
        end if
        table%coef(k)%f = value
    end subroutine coefs_set

    logical function coefs_has(table, name) result(found)
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name
        found = coefs_position(table, name) > 0
    end function coefs_has

    logical function coefs_is_derived(table, name) result(derived)
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name
        derived = table%coef(coef_checked(table, name, "coefs_is_derived"))%derived
    end function coefs_is_derived

    function coefs_derived(table, name) result(drv)
        ! The definition (from, factor, region, anomaly) of the derived
        ! coefficient `name`; stops on an unknown or a tsgen coefficient.
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name
        type(derived_spec) :: drv
        integer :: k
        k = coef_checked(table, name, "coefs_derived")
        if (.not. table%coef(k)%derived) then
            call fesm_error("kryos_snap_coefs::coefs_derived", "coefficient '"//trim(name)// &
                            "' is a tsgen coefficient, not derived")
        end if
        drv = table%coef(k)%drv
    end function coefs_derived

    real(wp) function coefs_value(table, name) result(f)
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name
        f = table%coef(coef_checked(table, name, "coefs_value"))%f
    end function coefs_value

    integer function coef_checked(table, name, caller) result(k)
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name, caller
        k = coefs_position(table, name)
        if (k == 0) then
            call fesm_error("kryos_snap_coefs::"//trim(caller), "unknown coefficient '"//trim(name)//"'")
        end if
    end function coef_checked

    integer function coefs_position(table, name) result(k)
        ! Position of the coefficient `name` in the table, 0 if absent.
        type(snap_coef_table_t), intent(in) :: table
        character(len=*),        intent(in) :: name
        if (allocated(table%coef)) then
            do k = 1, size(table%coef)
                if (trim(table%coef(k)%name) == trim(name)) return
            end do
        end if
        k = 0
    end function coefs_position

end module kryos_snap_coefs
