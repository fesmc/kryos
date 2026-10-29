module kryos_shelf_ocean
    ! The shelf component's own ocean-side state (phase 2 WP3, design Q2 = (b)): the
    ! thermal-forcing corrections of the shelf base, which a calibration (the
    ! experiment opt_tf) adjusts and every shelf backend adds to its thermal forcing,
    !
    !   tf_corr        [K] a field on the shelf grid: from a file (tf_correction), else 0,
    !                  then optimised by opt_tf (the control of yelmox's step_optimize)
    !   tf_corr_basin  [K] constants by region (corr_method = tf)
    !
    ! with the basins of the shelf grid (basin ids, also what the backend's own basin means
    ! use) and the grid spacing. They live in shelf_t%corr, not in a backend: the backend
    ! receives them as inputs before it computes the melt (marine_shelf adds
    ! (T_shlf - T_fp) + tf_corr + tf_corr_basin as it did when it held them).
    ! corr_method = bmb gives the backend its bmb correction by region (an initial field,
    ! bmb_corr, which stays the backend's state).
    !
    ! Keys: the group named by the shelf instance's `corr` key (default "shelf_corr"), all
    ! required (as the marine_shelf group they come from): basins, corr_method, corr_names,
    ! corr_values, corr_<name> (one selection expression per name), tf_correction, tf_path,
    ! tf_name. Restart: <instance>_corr_restart.nc (tf_corr, tf_corr_basin), written and
    ! read by the shelf wrappers (kryos_shelf).

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error, to_str
    use ncio,            only: nc_create, nc_write_attr, nc_write_dim, nc_open, nc_write, nc_close, nc_read
    use nml,             only: nml_read, nml_has_param, nml_has_group, nml_replace
    use regions,         only: regions_select, regions_basin_ids
    use kryos_component, only: init_ctx_t, NAME_LEN, PATH_LEN

    implicit none

    private
    public :: shelf_corr_init, shelf_corr_file, shelf_corr_restart_write, shelf_corr_restart_read

    character(len=*), parameter :: DEFAULT_GROUP = "shelf_corr"
    integer, parameter :: N_CORR_MAX = 50, LEN_EXPR = 1000

    type, public :: shelf_corr_t
        character(len=NAME_LEN) :: group = ""               ! the par-file group the keys came from
        character(len=NAME_LEN) :: domain = "", grid_name = ""
        integer  :: nx = 0, ny = 0
        real(sp) :: dx = 0.0_sp                             ! [m] grid spacing (as the backend's dx)
        real(sp), allocatable :: tf_corr(:,:)               ! [K] thermal-forcing correction, optimised (field)
        real(sp), allocatable :: tf_corr_basin(:,:)         ! [K] thermal-forcing correction by region (constants)
        real(sp), allocatable :: basins(:,:)                ! basin ids (0 = none), from the regions of the shelf grid
        real(sp), allocatable :: bmb_corr(:,:)              ! [m yr-1] bmb correction by region (initial field of the backend)
        real(sp), allocatable :: xc(:), yc(:)               ! [km] grid axes, for the restart file
    end type shelf_corr_t

contains

    subroutine shelf_corr_init(corr, ctx)
        ! The corrections of the shelf instance of ctx: the keys of its `corr` group, the
        ! basins and the region constants on the instance grid's regions, tf_corr from its
        ! file if tf_correction (yelmox marshelf_init).
        type(shelf_corr_t), intent(out) :: corr
        type(init_ctx_t),   intent(in)  :: ctx

        character(len=*), parameter :: proc = "kryos_shelf_ocean::shelf_corr_init"
        character(len=NAME_LEN) :: group, basins, method, names(N_CORR_MAX), tf_name
        character(len=LEN_EXPR) :: expr
        character(len=PATH_LEN) :: tf_path
        real(sp) :: values(N_CORR_MAX)
        logical  :: tf_correction
        integer  :: k, n, n_val

        group = DEFAULT_GROUP
        if (nml_has_param(ctx%par_file, ctx%instance, "corr")) call nml_read(ctx%par_file, ctx%instance, "corr", group)
        if (.not. nml_has_group(trim(ctx%par_file), trim(group))) then
            call fesm_error(proc, "shelf `"//trim(ctx%instance)//"`: the corrections group &"//trim(group)// &
                            " is not in the par file", "key `corr` of the instance; keys: basins, corr_method, "// &
                            "corr_names, corr_values, corr_<name>, tf_correction, tf_path, tf_name; par file "// &
                            trim(ctx%par_file))
        end if

        corr%group     = group
        corr%domain    = ctx%domain
        corr%grid_name = ctx%grid_name
        corr%nx        = ctx%grid%G%nx
        corr%ny        = ctx%grid%G%ny
        corr%dx        = real(ctx%grid%G%dx * ctx%grid%cs%xy_conv, sp)
        corr%xc        = real(ctx%grid%G%x, sp)
        corr%yc        = real(ctx%grid%G%y, sp)

        call nml_read(ctx%par_file, trim(group), "basins",      basins)
        call nml_read(ctx%par_file, trim(group), "corr_method", method)
        ! (values: a marker for the entries not given, to count them)
        names  = ""
        values = huge(1.0_sp)
        call nml_read(ctx%par_file, trim(group), "corr_names",  names)
        call nml_read(ctx%par_file, trim(group), "corr_values", values)
        n     = count(len_trim(names) > 0)
        n_val = count(values /= huge(1.0_sp))
        ! One value per name, in the same order. Without names corr_values is a
        ! placeholder and not checked.
        if (n > 0) then
            if (n_val /= n .or. any(len_trim(names(1:n)) == 0) .or. any(values(1:n) == huge(1.0_sp))) then
                call fesm_error(proc, "corr_names and corr_values must have the same number of entries", &
                                "group &"//trim(group)//": "//to_str(n)//" names, "//to_str(n_val)//" values")
            end if
        end if
        select case (trim(method))
        case ("none", "bmb", "tf")
        case default
            call fesm_error(proc, "corr_method must be none, bmb or tf", "group &"//trim(group)// &
                            ": corr_method = '"//trim(method)//"'")
        end select
        call nml_read(ctx%par_file, trim(group), "tf_correction", tf_correction)
        call nml_read(ctx%par_file, trim(group), "tf_path",       tf_path)
        call nml_read(ctx%par_file, trim(group), "tf_name",       tf_name)

        allocate(corr%tf_corr(corr%nx, corr%ny), corr%tf_corr_basin(corr%nx, corr%ny), &
                 corr%basins(corr%nx, corr%ny), corr%bmb_corr(corr%nx, corr%ny))
        corr%tf_corr       = 0.0_sp
        corr%tf_corr_basin = 0.0_sp
        corr%bmb_corr      = 0.0_sp
        corr%basins        = real(regions_basin_ids(ctx%reg, trim(basins)), sp)

        do k = 1, n
            call nml_read(ctx%par_file, trim(group), "corr_"//trim(names(k)), expr)
            select case (trim(method))
            case ("bmb")
                where (regions_select(ctx%reg, trim(expr))) corr%bmb_corr = values(k)
            case ("tf")
                where (regions_select(ctx%reg, trim(expr))) corr%tf_corr_basin = values(k)
            end select
        end do

        if (tf_correction) then                             ! start from a field of a file
            call nml_replace(tf_path, "{domain}",    trim(ctx%domain))
            call nml_replace(tf_path, "{grid_name}", trim(ctx%grid_name))
            call nc_read(trim(tf_path), trim(tf_name), corr%tf_corr)
        end if
    end subroutine shelf_corr_init

    function shelf_corr_file(instance) result(file)
        ! The base name of the corrections' restart file of shelf instance `instance`.
        character(len=*), intent(in) :: instance
        character(len=PATH_LEN) :: file
        file = trim(instance)//"_corr_restart.nc"
    end function shelf_corr_file

    subroutine shelf_corr_restart_write(corr, instance, dir, t)
        ! <instance>_corr_restart.nc: the variables of yelmox's marine_shelf restart
        ! (tf_corr, tf_corr_basin), same dimensions, units and attributes.
        type(shelf_corr_t), intent(in) :: corr
        character(len=*),   intent(in) :: instance, dir
        real(dp),           intent(in) :: t

        character(len=:), allocatable :: file
        integer :: ncid

        file = trim(dir)//"/"//trim(shelf_corr_file(instance))
        call nc_create(file)
        call nc_write_attr(file, "domain",    trim(corr%domain))
        call nc_write_attr(file, "grid_name", trim(corr%grid_name))
        call nc_write_dim(file, "xc", x=corr%xc, units="km")
        call nc_write_dim(file, "yc", x=corr%yc, units="km")
        call nc_write_dim(file, "month", x=1, dx=1, nx=12, units="month")
        call nc_write_dim(file, "time", x=t, dx=1.0_dp, nx=1, units="years", unlimited=.true.)
        call nc_open(file, ncid, writable=.true.)
        call nc_write(file, "tf_corr", corr%tf_corr, start=[1,1,1], units="K", &
                      dim1="xc", dim2="yc", dim3="time", ncid=ncid, grid_mapping="crs")
        call nc_write(file, "tf_corr_basin", corr%tf_corr_basin, start=[1,1,1], units="K", &
                      dim1="xc", dim2="yc", dim3="time", ncid=ncid, grid_mapping="crs")
        call nc_close(ncid)
    end subroutine shelf_corr_restart_write

    subroutine shelf_corr_restart_read(corr, instance, dir)
        type(shelf_corr_t), intent(inout) :: corr
        character(len=*),   intent(in)    :: instance, dir

        character(len=:), allocatable :: file

        file = trim(dir)//"/"//trim(shelf_corr_file(instance))
        call nc_read(file, "tf_corr",       corr%tf_corr,       start=[1,1,1], count=[corr%nx, corr%ny, 1])
        call nc_read(file, "tf_corr_basin", corr%tf_corr_basin, start=[1,1,1], count=[corr%nx, corr%ny, 1])
    end subroutine shelf_corr_restart_read

end module kryos_shelf_ocean
