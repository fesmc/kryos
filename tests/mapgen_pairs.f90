program mapgen_pairs
    ! Map files for explicit grid pairs, from the grid descriptions of a map folder (coupler_prime; the folder is written):
    !   mapgen_pairs.x <map folder> <src>:<dst>:<con|bilin|nn> [...]
    ! The pairs of a kryos run are made by the run itself (the coupler registers the run's grids); this tool is for pairs the
    ! science case does not use (phase 4 WP7: GAUSS-L170 <-> ANT-16KM). Not part of `make check`.
    use coupler, only: coupler_class, coupler_init, coupler_prime
    implicit none
    character(len=512) :: fldr, arg
    character(len=64)  :: src, dst, meth
    type(coupler_class) :: cpl
    integer :: i, p1, p2

    if (command_argument_count() < 2) error stop "usage: mapgen_pairs.x <map folder> <src>:<dst>:<method> [...]"
    call get_command_argument(1, fldr)
    call coupler_init(cpl, trim(fldr))
    do i = 2, command_argument_count()
        call get_command_argument(i, arg)
        p1 = index(arg, ":"); p2 = index(arg(p1+1:), ":") + p1
        src = arg(1:p1-1); dst = arg(p1+1:p2-1); meth = trim(arg(p2+1:))
        write(*,'(a)') "map "//trim(src)//" -> "//trim(dst)//" ("//trim(meth)//")"
        call coupler_prime(cpl, trim(src), trim(dst), trim(meth))
    end do
end program mapgen_pairs
