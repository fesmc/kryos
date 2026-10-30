module g2a_mod
    ! Gate G2a (phase4-plan 13.2): the VILMA2 geo adapter (vilma_earth_t through the component wrappers, as the core drives
    ! it) against the standalone driver vilma.x on the same analytic load on the Gauss grid GAUSS-L<lmax>: a grounded land-ice
    ! cap inside a land box (the host region, w = 1; no coast, no external ice), growing by 25 % of its first size per step.
    ! Not part of `make check` (a Levante job: lmax 170, the 3D viscosity of the science case, minutes per call).
    !
    !   g2a_cap.x <nstep>          in a run folder with input/ (kryos), fesm-utils/, fesmdata/ and maps/ (writable)
    !
    ! Writes kryos_out.nc (dz_bed, dz_ss [= VILMA2's, i.e. the export plus dphi], bsl, bsl_ive, dphi, C_ocn; record 0 =
    ! the seed after the spin-up, records n = after step n) and the inputs of the standalone run: g2a_force.nc (h_ice, time),
    ! g2a_ref.nc (z_bed_eq, h_ice_eq), vilma_g2a.nml (&ctl, &vilma = the science case set with the run's densities).
    ! Orientation: everything here is on the host grid, latitude ascending (the exports, kryos_out.nc and the standalone
    ! inputs); vilma.x reads the inputs on the Gauss grid (remap_input = .false.), reverses their rows by their lat axis into
    ! VILMA2's layout (north first) and writes vilma_out.nc in that layout. The cap is northern and the viscosity (Bagge2021)
    ! hemispherically asymmetric, so a mirrored layout on either side shows; g2a_compare.jl aligns the files by lat.
    ! The load of step n is the same field in both runs: the host terms I_g = H, w = 1 in the box and h_ice = H in the
    ! standalone run; the reference has no ice, the start slice is the cap of size 1 (so both spin up under it).
    ! Timings of the adapter calls go to stdout ("g2a:" lines). Compared by tools/levante/g2a_compare.jl.

    use precision,       only: dp
    use coordinates,     only: grid_class
    use ncio,            only: nc_create, nc_write_dim, nc_write
    use constants,       only: pi
    use phys_constants,  only: phys_const_load
    use fesm_fields,     only: fields_get, fields_put
    use kryos_component, only: init_ctx_t, step_time_t, census_t, PASS_INIT, PASS_GEOMETRY, PASS_STATE
    use kryos_geo,       only: geo_slot_t, geo_create, geo_init, geo_init_state, geo_update, geo_end
    use kryos_geo_vilma2, only: vilma_earth_t

    implicit none

    integer, parameter :: LMAX = 170, NPHI = 680, NLAT = 342
    real(dp), parameter :: DTT = 10.0_dp
    character(len=*), parameter :: PAR = "g2a_par.nml"
    real(dp), parameter :: CLON = 90.0_dp, CLAT = 58.0_dp, RC = 12.0_dp, H0 = 2000.0_dp
    real(dp), parameter :: BOX(4) = [40.0_dp, 140.0_dp, 40.0_dp, 75.0_dp]      ! lon0 lon1 lat0 lat1 of the host region

    type(init_ctx_t), save :: ctx
    type(census_t),   save :: census
    type(geo_slot_t), save :: glob
    real(dp), allocatable, save :: lat(:,:), lon(:,:), z_ref(:,:), w(:,:), zero(:,:)
    real(dp), allocatable, save :: o_bed(:,:,:), o_ss(:,:,:), o_C(:,:,:), o_s(:,:), tt(:), forc(:,:,:)
    integer, save :: nstep
    integer(8), save :: c0, c1, rate

contains

    subroutine run()
        integer :: n
        character(len=64) :: arg
    call get_command_argument(1, arg)
    read(arg, *) nstep
    allocate(lat(NPHI,NLAT), lon(NPHI,NLAT), z_ref(NPHI,NLAT), w(NPHI,NLAT), zero(NPHI,NLAT))
    allocate(o_bed(NPHI,NLAT,0:nstep), o_ss(NPHI,NLAT,0:nstep), o_C(NPHI,NLAT,0:nstep), o_s(4,0:nstep), tt(0:nstep), &
             forc(NPHI,NLAT,0:nstep))
    call system_clock(count_rate=rate)
    zero = 0.0_dp

    call par_write()
    ctx%par_file  = PAR
    ctx%instance  = "geo_glob"
    ctx%component = "geo"
    ctx%backend   = "vilma2"
    ctx%scope     = "global"
    ctx%map_fldr  = "maps"
    ctx%time      = 0.0_dp
    ctx%time_rel  = 0.0_dp
    ctx%dtt       = DTT
    call phys_const_load(ctx%cnst, "fesm-utils/par/phys_const_earth.nml", "Earth")

    call system_clock(c0)
    call geo_create(glob, ctx, census)
    call geo_init(glob, ctx)
    call coordinates_of(glob)
    call terms()
    call system_clock(c1); write(*,'(a,f9.2,a)') "g2a: geo_init (grid, VILMA2 setup incl. 3D viscosity) ", real(c1-c0, dp)/rate, " s"

    call system_clock(c0)
    call fields_put(glob%c%imports, glob%c%hi_z_bed_ref,   z_ref, t1=0.0_dp)
    call fields_put(glob%c%imports, glob%c%hi_H_ice_ref,   zero,  t1=0.0_dp)
    call fields_put(glob%c%imports, glob%c%hi_host_w,      w,     t1=0.0_dp)
    call fields_put(glob%c%imports, glob%c%hi_host_Ig_ref, zero,  t1=0.0_dp)
    call fields_put(glob%c%imports, glob%c%hi_host_C_ref,  zero,  t1=0.0_dp)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_INIT)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_GEOMETRY)
    call put_load(0, 0.0_dp)
    call geo_init_state(glob, st(0.0_dp, 0.0_dp), PASS_STATE)       ! spin-up (pre_spinup_1d) and the seed
    call system_clock(c1); write(*,'(a,f9.2,a)') "g2a: reference, spin-up and seed ", real(c1-c0, dp)/rate, " s"
    tt(0) = 0.0_dp
    call record(0)
    do n = 1, nstep
        tt(n) = DTT*n
        call put_load(n, tt(n))
        call system_clock(c0)
        call geo_update(glob, st(tt(n) - DTT, tt(n)))
        call system_clock(c1)
        write(*,'(a,i0,a,f9.2,a)') "g2a: update ", n, " (dt = 10 yr) ", real(c1-c0, dp)/rate, " s"
        call record(n)
    end do
    call geo_end(glob)
    call write_outputs()
    call write_standalone_inputs()

    end subroutine run


    real(dp) function cap(k, i, j)
        integer, intent(in) :: k, i, j
        real(dp) :: a, b, dl, d
        a = lat(i,j)*pi/180.0_dp; b = CLAT*pi/180.0_dp; dl = (lon(i,j) - CLON)*pi/180.0_dp
        d = acos(max(-1.0_dp, min(1.0_dp, sin(a)*sin(b) + cos(a)*cos(b)*cos(dl))))*180.0_dp/pi
        cap = H0*(1.0_dp + 0.25_dp*real(k, dp))*max(0.0_dp, 1.0_dp - d/RC)
    end function cap

    subroutine terms()
        integer :: i, j, n
        w = 0.0_dp; z_ref = -3000.0_dp
        do j = 1, NLAT
            do i = 1, NPHI
                if (lon(i,j) >= BOX(1) .and. lon(i,j) <= BOX(2) .and. lat(i,j) >= BOX(3) .and. lat(i,j) <= BOX(4)) then
                    w(i,j) = 1.0_dp; z_ref(i,j) = 500.0_dp
                end if
            end do
        end do
        do n = 0, nstep
            do j = 1, NLAT
                do i = 1, NPHI
                    forc(i,j,n) = cap(n, i, j)             ! load n: the cap of size 1 + 0.25 n (zero outside the box: it ends inside)
                end do
            end do
        end do
    end subroutine terms

    subroutine put_load(k, t1)
        integer,  intent(in) :: k
        real(dp), intent(in) :: t1
        call fields_put(glob%c%imports, glob%c%hi_H_ice_ext, zero,       t1=t1)
        call fields_put(glob%c%imports, glob%c%hi_host_w,    w,          t1=t1)
        call fields_put(glob%c%imports, glob%c%hi_host_Ig,   forc(:,:,k), t1=t1)
        call fields_put(glob%c%imports, glob%c%hi_host_C,    zero,       t1=t1)
        call fields_put(glob%c%imports, glob%c%hi_host_F,    zero,       t1=t1)
    end subroutine put_load

    subroutine record(k)
        integer, intent(in) :: k
        real(dp) :: a(NPHI,NLAT), s
        call fields_get(glob%c%exports, glob%c%h_dz_bed, a);  o_bed(:,:,k) = a
        call fields_get(glob%c%exports, glob%c%h_dz_ss,  a);  o_ss(:,:,k) = a
        call fields_get(glob%c%exports, glob%c%h_C_ocn,  a);  o_C(:,:,k) = a
        call fields_get(glob%c%exports, glob%c%h_dphi,    s); o_s(1,k) = s
        call fields_get(glob%c%exports, glob%c%h_bsl,     s); o_s(2,k) = s
        call fields_get(glob%c%exports, glob%c%h_bsl_ive, s); o_s(3,k) = s
        call fields_get(glob%c%exports, glob%c%h_mass_resid, s); o_s(4,k) = s
        o_ss(:,:,k) = o_ss(:,:,k) + o_s(1,k)                   ! VILMA2's own dz_ss (the export is minus dphi)
    end subroutine record

    subroutine coordinates_of(g)
        type(geo_slot_t), intent(inout) :: g
        type(grid_class) :: gg
        integer :: i, j
        select type (c => g%c)
        type is (vilma_earth_t)
            gg = c%grid_obj
        end select
        do j = 1, NLAT
            do i = 1, NPHI
                lon(i,j) = gg%G%x(i)
                lat(i,j) = gg%G%y(j)
            end do
        end do
    end subroutine coordinates_of

    function st(t0, t1) result(s)
        real(dp), intent(in) :: t0, t1
        type(step_time_t) :: s
        s = step_time_t(t0=t0, t1=t1, t_eval=t1, dt=DTT, dtt=DTT)
    end function st

    subroutine vilma_group(u)
        integer, intent(in) :: u
        character(len=32) :: v
        write(u,'(a)') "&vilma", " lmax = 170", ' earth = "PREM"', " n_layer = 11", &
            " r_bot = 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3, 0.0", &
            " r_top = 6371.0e3, 6291.0e3, 6151.0e3, 5971.0e3, 5771.0e3, 5701.0e3, 5600.0e3, 4943.0e3, 4287.0e3, 3630.0e3, 3480.0e3", &
            " eta = 1.0e40, 4.0e20, 4.0e20, 4.0e20, 4.0e20, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 1.0e22, 0.0", &
            " rheology = 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2", " l_visc_3d = .true.", &
            ' visc_3d_file = "fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc"', " visc_3d_structure = 6", &
            " pre_spinup_1d = "//merge(".true. ", ".false.", spinup_on())
        call get_environment_variable("G2A_RATE_TOL", v)          ! equil_rate_tol [m/yr] of the spin-up in both runs
        if (len_trim(v) > 0) write(u,'(a)') " equil_rate_tol = "//trim(v)
    end subroutine vilma_group

    logical function spinup_on()
        ! G2A_SPINUP=0: no spin-up (the seed only), in both runs
        character(len=8) :: v
        call get_environment_variable("G2A_SPINUP", v)
        spinup_on = trim(v) /= "0"
    end function spinup_on

    subroutine par_write()
        integer :: u
        open(newunit=u, file=PAR, status="replace", action="write")
        write(u,'(a)') "&geo_glob", ' backend = "vilma2"', ' vilma = "vilma"', &
            ' allow_no_spinup = '//merge(".false.", ".true. ", spinup_on()), "/"
        call vilma_group(u)
        write(u,'(a)') "/"
        close(u)
    end subroutine par_write

    subroutine write_outputs()
        character(len=*), parameter :: f = "kryos_out.nc"
        call nc_create(f)
        call nc_write_dim(f, "lon", x=lon(:,1), units="degrees_east")
        call nc_write_dim(f, "lat", x=lat(1,:), units="degrees_north")
        call nc_write_dim(f, "time", x=tt, units="years", unlimited=.false.)
        call nc_write(f, "dz_bed", o_bed, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write(f, "dz_ss",  o_ss,  dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_write(f, "C_ocn",  o_C,   dim1="lon", dim2="lat", dim3="time", units="1")
        call nc_write(f, "dphi",    o_s(1,:), dim1="time", units="m")
        call nc_write(f, "bsl",     o_s(2,:), dim1="time", units="m")
        call nc_write(f, "bsl_ive", o_s(3,:), dim1="time", units="m")
        call nc_write(f, "mass_resid", o_s(4,:), dim1="time", units="1")
    end subroutine write_outputs

    subroutine write_standalone_inputs()
        character(len=*), parameter :: ff = "g2a_force.nc", fr = "g2a_ref.nc"
        integer :: u
        call nc_create(ff)
        call nc_write_dim(ff, "lon", x=lon(:,1), units="degrees_east")
        call nc_write_dim(ff, "lat", x=lat(1,:), units="degrees_north")
        call nc_write_dim(ff, "time", x=tt, units="years", unlimited=.false.)
        call nc_write(ff, "h_ice", forc, dim1="lon", dim2="lat", dim3="time", units="m")
        call nc_create(fr)
        call nc_write_dim(fr, "lon", x=lon(:,1), units="degrees_east")
        call nc_write_dim(fr, "lat", x=lat(1,:), units="degrees_north")
        call nc_write(fr, "bedrock_topography", z_ref, dim1="lon", dim2="lat", units="m")
        call nc_write(fr, "ice_thickness", zero, dim1="lon", dim2="lat", units="m")
        open(newunit=u, file="vilma_g2a.nml", status="replace", action="write")
        write(u,'(a)') "&ctl", ' file_forcing = "g2a_force.nc"', ' name_ice = "h_ice"', ' name_time = "time"', &
            ' file_ref = ""', ' name_zbed_eq = "z_bed_eq"', ' name_hice_ref = "h_ice_ref"', ' file_out = "vilma_out.nc"', &
            ' file_hor = ""', " write_lossless = .true.", " time_init = -1.0e30", " time_end = 1.0e30", &
            " remap_input = .false.", ' name_lon = "lon"', ' name_lat = "lat"', " i_eq = 2", ' z_bed_ref_file = ""', &
            ' h_ice_ref_file = ""', ' z_bed_eq_file = "g2a_ref.nc"', ' h_ice_eq_file = "g2a_ref.nc"', ' rsl_restart_file = ""', &
            ' name_z_bed_ref = "bedrock_topography"', ' name_h_ice_ref = "ice_thickness"', ' name_rsl = "rsl"', &
            ' restart_in_file = ""', "/"
        call vilma_group(u)
        write(u,'(a,es24.16)') " rho_ice = ", ctx%cnst%rho_ice
        write(u,'(a,es24.16)') " rho_water = ", ctx%cnst%rho_sw
        write(u,'(a)') "/"
        close(u)
    end subroutine write_standalone_inputs

end module g2a_mod

program g2a_cap
    use g2a_mod, only: run
    implicit none
    call run()
end program g2a_cap
