program hub_areas
    ! Model areas of grids as the global geometry hub registers them (kryos_globe_ice: grid%area * xy_conv**2, true spherical
    ! cell areas), for the offline merge of gate G1b (4) on a real run (tools/levante/hub_offline.jl).
    !   hub_areas.x <map folder> <out.nc> <grid> [<grid> ...]
    ! Writes area_<grid>(x, y) [m2] with the dimensions x_<grid>, y_<grid>. Not part of `make check`.
    use precision,  only: dp
    use coords,     only: grid_class, grid_cdo_read_desc
    use ncio,       only: nc_create, nc_write_dim, nc_write
    implicit none
    character(len=512) :: fldr, out, name
    type(grid_class) :: g
    integer :: i, k
    real(dp), allocatable :: a(:,:)

    if (command_argument_count() < 3) error stop "usage: hub_areas.x <map folder> <out.nc> <grid> [<grid> ...]"
    call get_command_argument(1, fldr)
    call get_command_argument(2, out)
    call nc_create(trim(out))
    do i = 3, command_argument_count()
        call get_command_argument(i, name)
        call grid_cdo_read_desc(g, trim(name), trim(fldr))
        allocate(a(g%G%nx, g%G%ny))
        a = g%area*g%cs%xy_conv**2
        call nc_write_dim(trim(out), "x_"//trim(name), x=real([(k, k=1,g%G%nx)], dp), units="1")
        call nc_write_dim(trim(out), "y_"//trim(name), x=real([(k, k=1,g%G%ny)], dp), units="1")
        call nc_write(trim(out), "area_"//trim(name), a, dim1="x_"//trim(name), dim2="y_"//trim(name), units="m2")
        deallocate(a)
    end do
contains
end program hub_areas
