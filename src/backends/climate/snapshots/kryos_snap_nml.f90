module kryos_snap_nml
    ! Namelist helpers shared by the private modules of the snapshots climate
    ! backend (coefficients, terms, snapshots, transforms): the schema files of
    ! their groups and the lookup of a named group (snapshot, term, coef,
    ! varslice, taper, bias) in the par file, then in the run's library files.
    ! Ports yelmox yelmox_nml.f90; the library list is [kryos] library
    ! (init_ctx_t%library), passed in by the adapter.

    use nml, only: nml_find_group

    implicit none

    private

    ! Schemas (input/): the groups of the backend and the forcing coefficients.
    character(len=*), parameter, public :: SNAP_DEFAULTS = "input/snapshots_defaults.nml"
    character(len=*), parameter, public :: COEF_DEFAULTS = "input/forcing_defaults.nml"

    integer, parameter, public :: FILE_LEN = 512

    public :: group_file

contains

    function group_file(par_file, library, group, what) result(filename)
        ! The par file if it holds `group` (searched first, so a par file never
        ! names itself and its groups override the library), else the first
        ! library file that does; stops naming the group and the files otherwise.
        ! `what` names the reader in the message (e.g. "&snapshot_lgm tas").
        character(len=*), intent(in) :: par_file, library(:), group, what
        character(len=FILE_LEN) :: filename

        character(len=FILE_LEN), allocatable :: files(:)

        files = [character(len=FILE_LEN) :: par_file, library]
        filename = nml_find_group(files, trim(group), what)
    end function group_file

end module kryos_snap_nml
