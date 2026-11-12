program test_forcing_tsgen
    ! forcing = tsgen (make check): coefficients from a par file, exported at
    ! t1 in sp widened; a derived coefficient set through forcing_set reaches
    ! the export one update later (the lag-1 copy of yelmox); restart write
    ! and read continue the series exactly. The orbit (key `orbit`): the
    ! elements exported at the forcing time are insol's (the tables, or the
    ! fixed elements of the template: 1950); the stream insol writes the
    ! monthly mean insolation on 1° latitudes.

    use precision,        only: sp, dp
    use fesm_fields,      only: fields_get
    use kryos_component,  only: init_ctx_t, step_time_t, census_t, PASS_INIT
    use kryos_forcing,    only: forcing_slot_t, forcing_create, forcing_init, forcing_init_state, &
                                forcing_update, forcing_set, forcing_restart_write, forcing_restart_read, &
                                forcing_kill, forcing_end, forcing_write_output
    use insolation,       only: insol_class, insol_init, orbit_pars_class, orbit_pars_at
    use ncio,             only: nc_read

    implicit none

    character(len=*), parameter :: par = "test_forcing_tsgen.nml", rdir = "test_forcing_tsgen_rst"

    type(forcing_slot_t) :: a, b, c, f
    type(init_ctx_t)     :: ctx
    type(census_t)       :: census, census_b, census_c, census_f
    type(insol_class)    :: tab
    type(orbit_pars_class) :: o
    real(dp) :: e(3), S(180, 12, 1)
    type(step_time_t)    :: st
    real(sp) :: v, w
    integer  :: u

    open(newunit=u, file=par, status="replace", action="write")
    write(u,'(a)') "&forcing"
    write(u,'(a)') "    backend = 'tsgen'"
    write(u,'(a)') "    coefs   = 'alpha' 'dTo'"
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_alpha"
    write(u,'(a)') "    method    = 'const'"
    write(u,'(a)') "    f_min     = 0.0"
    write(u,'(a)') "    f_max     = 1.0"
    write(u,'(a)') "    df_dt_max = 0.01"
    write(u,'(a)') "/"
    write(u,'(a)') "&coef_dTo"
    write(u,'(a)') "    method = 'derived'"
    write(u,'(a)') "    from   = 'air_anomaly'"
    write(u,'(a)') "    factor = 0.25"
    write(u,'(a)') "    units  = 'K'"
    write(u,'(a)') "/"
    close(u)

    ctx%par_file = par
    ctx%instance = "forcing"
    ctx%component = "forcing"
    ctx%backend  = "tsgen"
    ctx%scope    = "global"
    ctx%domain   = ""
    ctx%grid_name = ""

    call make(a, census)

    ! init: alpha at the start, dTo = 0 (derived)
    st = step_time_t(t0=0.0_dp, t1=0.0_dp, t_eval=0.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call forcing_init_state(a, st, PASS_INIT)
    call get(a, 1, v);  call check(v == 0.0_sp, "alpha at init")
    call get(a, 2, v);  call check(v == 0.0_sp, "dTo at init")

    ! step 1: alpha advances, dTo keeps 0 until a producer sets it
    st = step_time_t(t0=0.0_dp, t1=10.0_dp, t_eval=10.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call forcing_update(a, st)
    call get(a, 1, v);  call check(abs(v - 0.1_sp) < 1e-6_sp, "alpha after step 1")
    call get(a, 2, v);  call check(v == 0.0_sp, "dTo before it is set")

    ! the climate sets dTo at t = 10; the next update exports it
    call forcing_set(a, "dTo", real(2.5_sp, dp), 10.0_dp)
    st = step_time_t(t0=10.0_dp, t1=20.0_dp, t_eval=20.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call forcing_update(a, st)
    call get(a, 1, v);  call check(abs(v - 0.2_sp) < 1e-6_sp, "alpha after step 2")
    call get(a, 2, v);  call check(v == 2.5_sp, "dTo after forcing_set")
    call check(.not. forcing_kill(a), "kill switch off")

    ! restart: write at t = 20, restore into a fresh instance, continue both
    call forcing_restart_write(a, rdir, 20.0_dp)
    ctx%time_rel = 0.0_dp
    call make(b, census_b)
    call forcing_restart_read(b, rdir, 20.0_dp)
    st = step_time_t(t0=20.0_dp, t1=30.0_dp, t_eval=30.0_dp, dt=10.0_dp, dtt=10.0_dp)
    call forcing_update(a, st)
    call forcing_set(b, "dTo", real(2.5_sp, dp), 20.0_dp)
    call forcing_update(b, st)
    call get(a, 1, v);  call get(b, 1, w)
    call check(v == w, "alpha continues bit for bit after restart")
    call get(b, 2, w);  call check(w == 2.5_sp, "dTo after restart and set")

    call forcing_end(a)
    call forcing_end(b)

    ! the orbit at the forcing time: insol's tables
    open(newunit=u, file=par, status="replace", action="write")
    write(u,'(a)') "&forcing", "    backend = 'tsgen'", "    orbit = 'orbit'", "/", &
                   "&orbit", "    method = 'forcing_time'", "    fldr = 'input'", "/", &
                   "&forcing_fix", "    backend = 'tsgen'", "    orbit = 'orbit_fix'", "/", &
                   "&orbit_fix", "    method = 'fixed_orbit'", "/"
    close(u)
    call insol_init(tab, fldr="input")
    ctx%time_rel = -21000.0_dp
    ctx%sec_year = 365.0_dp*86400.0_dp
    call make(c, census_c)
    st = step_time_t(t0=-21000.0_dp, t1=-21000.0_dp, t_eval=-21000.0_dp, dt=100.0_dp, dtt=100.0_dp)
    call forcing_init_state(c, st, PASS_INIT)
    call orbit_get(c, e);  o = orbit_pars_at(tab, -21000.0_dp)
    call check(e(1) == o%ecc .and. e(2) == o%obl .and. e(3) == o%per, "orbit at init = insol's tables at -21000")
    st = step_time_t(t0=-21000.0_dp, t1=-20900.0_dp, t_eval=-20900.0_dp, dt=100.0_dp, dtt=100.0_dp)
    call forcing_update(c, st)
    call orbit_get(c, e);  o = orbit_pars_at(tab, -20900.0_dp)
    call check(e(1) == o%ecc .and. e(2) == o%obl .and. e(3) == o%per, "orbit after a step = the tables at -20900")
    call forcing_write_output(c, rdir//"/forcing_insol.nc", -20900.0_dp, "insol", .true.)
    call nc_read(rdir//"/forcing_insol.nc", "S", S)
    call check(all(S >= 0.0_dp) .and. abs(sum(S(90:91, :, 1))/24.0_dp - 418.0_dp) < 10.0_dp .and. &
               maxval(S(180, :, 1)) > 500.0_dp .and. minval(S(180, :, 1)) == 0.0_dp, &
               "stream insol: monthly S, equator about 418 W m-2, polar day and night")
    call forcing_end(c)

    ! a fixed orbit: the template's elements, those of 1950
    ctx%instance = "forcing_fix"
    call make(f, census_f)
    call forcing_init_state(f, st, PASS_INIT)
    call orbit_get(f, e)
    call check(e(1) == 0.016725_dp .and. e(2) == 23.4458_dp .and. e(3) == 102.063_dp, "fixed orbit: the 1950 elements")
    call forcing_end(f)

    call execute_command_line("rm -rf "//rdir//" "//par)
    write(*,'(a)') "test_forcing_tsgen: OK"

contains

    subroutine make(slot, cen)
        type(forcing_slot_t), intent(inout) :: slot
        type(census_t),       intent(inout) :: cen
        call forcing_create(slot, ctx, cen)
        call forcing_init(slot, ctx)
        call execute_command_line("mkdir -p "//rdir)
    end subroutine make

    subroutine get(slot, k, x)
        type(forcing_slot_t), intent(in)  :: slot
        integer,              intent(in)  :: k
        real(sp),             intent(out) :: x
        call fields_get(slot%c%exports, slot%c%h_coef(k), x)
    end subroutine get

    subroutine orbit_get(slot, x)
        type(forcing_slot_t), intent(in)  :: slot
        real(dp),             intent(out) :: x(3)
        integer :: k
        do k = 1, 3
            call fields_get(slot%c%exports, slot%c%h_orbit(k), x(k))
        end do
    end subroutine orbit_get

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            write(*,'(a)') "test_forcing_tsgen: FAILED ("//what//")"
            error stop 1
        end if
    end subroutine check

end program test_forcing_tsgen
