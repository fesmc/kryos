program test_forcing_tsgen
    ! forcing = tsgen (make check): coefficients from a par file, exported at
    ! t1 in sp widened; a derived coefficient set through forcing_set reaches
    ! the export one update later (the lag-1 copy of yelmox); restart write
    ! and read continue the series exactly.

    use precision,        only: sp, dp
    use fesm_fields,      only: fields_get
    use kryos_component,  only: init_ctx_t, step_time_t, census_t, PASS_INIT
    use kryos_forcing,    only: forcing_slot_t, forcing_create, forcing_init, forcing_init_state, &
                                forcing_update, forcing_set, forcing_restart_write, forcing_restart_read, &
                                forcing_kill, forcing_end

    implicit none

    character(len=*), parameter :: par = "test_forcing_tsgen.nml", rdir = "test_forcing_tsgen_rst"

    type(forcing_slot_t) :: a, b
    type(init_ctx_t)     :: ctx
    type(census_t)       :: census, census_b
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

    subroutine check(ok, what)
        logical,          intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            write(*,'(a)') "test_forcing_tsgen: FAILED ("//what//")"
            error stop 1
        end if
    end subroutine check

end program test_forcing_tsgen
