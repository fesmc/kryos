program test_wiring
    ! kryos_wiring on toy field states (no files, no real backends): the table,
    ! resolution (gaps, duplicates, contracts), identity and remapped rows, every
    ! transform kind in sp and dp, gather, LAG_INIT, masks, and every lag error.
    !
    ! Remapped rows use two aligned toy grids registered in the coupler (the maps
    ! are built in-package into a scratch folder, removed at the end). Abort
    ! paths (fesm_error stops the program) are covered two ways: the procedures
    ! take an optional stat, which the tests below use for every error case
    ! of resolution and of the lag check; and the program re-runs itself with an
    ! argument per case, the child aborting for real, the parent checking the exit
    ! status and the message on stderr.

    use precision,       only: sp, dp
    use coords,          only: grid_class, grid_init
    use coupler,         only: coupler_class, coupler_init, coupler_add_grid, remap
    use fesm_fields
    use kryos_component, only: axis_t, grid_ref_t, month_axis
    use kryos_wiring

    implicit none

    character(len=*), parameter :: fldr = "maps_wiring_test", errfile = "test_wiring_abort.err"
    character(len=*), parameter :: GF = "wt-fine", GC = "wt-coarse"   ! 12x12 at 50 km, 6x6 at 100 km
    integer, parameter :: NF = 12, NC = 6, NM = 12

    type(coupler_class) :: cpl
    type(grid_class)    :: gf_, gc_
    ! one producer and one consumer, rebuilt by every test
    type(field_state_t), target :: sx, si, dx, di
    type(axis_t),        target :: ax(1), ax2(1)
    integer :: nfail = 0
    ! the toy domain of tc_domain
    type(field_state_t), target :: ice_i, ice_e, geo_i, geo_e, gg_i, gg_e, hub_i, hub_e, cl_i, cl_e
    type(field_state_t), target :: sf_i, sf_e, sh_i, sh_e, gt_e, sd_e, fo_i, fo_e, e_i, e_e
    character(len=FIELD_NAME_LEN), target :: opt_ice(3), opt_sh(7)
    type(wire_node_t) :: nodes(10)
    integer, parameter :: I_ = FIELD_INSTANT, M_ = FIELD_MEAN, IT = FIELD_INTENSIVE, EX = FIELD_EXTENSIVE

    call grid_init(gf_, name=GF, mtype="cartesian", units="kilometers", x0=25.0_dp, dx=50.0_dp, nx=NF, &
                   y0=25.0_dp, dy=50.0_dp, ny=NF)
    call grid_init(gc_, name=GC, mtype="cartesian", units="kilometers", x0=50.0_dp, dx=100.0_dp, nx=NC, &
                   y0=50.0_dp, dy=100.0_dp, ny=NC)
    if (command_argument_count() == 0) call execute_command_line("rm -rf "//fldr)
    call coupler_init(cpl, map_fldr=fldr)
    call coupler_add_grid(cpl, GF, gf_)
    call coupler_add_grid(cpl, GC, gc_)
    ax(1) = month_axis()
    ax2(1) = month_axis()

    if (command_argument_count() > 0) then
        call abort_case()
        error stop 77                                   ! the case did not abort
    end if

    call tc_table()
    call tc_domain()
    call tc_resolve_errors()
    call tc_identity()
    call tc_remap()
    call tc_add_ref()
    call tc_add_scalar()
    call tc_sub_ref()
    call tc_lapse()
    call tc_gather()
    call tc_mean_to_instant()
    call tc_init_once()
    call tc_rank0()
    call tc_masked()
    call tc_lag()
    call tc_aborts()

    call execute_command_line("rm -rf "//fldr//" "//errfile)
    if (nfail > 0) then
        write(*,"(a,i0,a)") "test_wiring: ", nfail, " check(s) FAILED"
        error stop 1
    end if
    write(*,"(a)") "test_wiring: OK"

contains

    ! =========================================================================
    ! helpers
    ! =========================================================================

    subroutine check(label, ok)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: ok
        if (ok) then
            write(*,"(a)") "  ok    "//label
        else
            write(*,"(a)") "  FAIL  "//label
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine st_add(st, name, units, tsem, cons, grid, nz, masked)
        ! Declare a field: rank 0 without grid, rank 2 with grid, rank 3 with nz (month axis).
        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in) :: name, units
        integer,             intent(in) :: tsem, cons
        character(len=*),    intent(in), optional :: grid
        integer,             intent(in), optional :: nz
        logical,             intent(in), optional :: masked
        integer :: h, n

        if (.not. present(grid)) then
            call fields_add(st, name, units, tsem, cons, h)
            return
        end if
        n = merge(NF, NC, trim(grid) == GF)
        if (present(nz)) then
            call fields_add(st, name, units, tsem, cons, h, grid=grid, dims=[n, n, nz], zdim="month")
        else
            call fields_add(st, name, units, tsem, cons, h, grid=grid, dims=[n, n], masked=masked)
        end if
    end subroutine st_add

    function wave(n, nz, a, b) result(x)
        ! A smooth field whose values are not representable in sp.
        integer,  intent(in) :: n, nz
        real(dp), intent(in) :: a, b
        real(dp) :: x(n,n,nz)
        integer :: i, j, k
        do k = 1, nz
            do j = 1, n
                do i = 1, n
                    x(i,j,k) = a + b*sin(0.37_dp*real(i, dp) + 0.11_dp*real(k, dp))*cos(0.21_dp*real(j, dp)) &
                               + real(i, dp)/7.0_dp + 1.0_dp/3.0_dp
                end do
            end do
        end do
    end function wave

    function wave2(n, a, b) result(x)
        integer,  intent(in) :: n
        real(dp), intent(in) :: a, b
        real(dp) :: x(n,n), y(n,n,1)
        y = wave(n, 1, a, b)
        x = y(:,:,1)
    end function wave2

    elemental function rnd(a, p) result(b)
        real(dp), intent(in) :: a
        integer,  intent(in) :: p
        real(dp) :: b
        b = merge(real(real(a, sp), dp), a, p == sp)
    end function rnd

    subroutine put3(st, name, x, p, t1, t0)
        ! Put a rank-2 or rank-3 array as a producer of precision p would: narrowed to sp, or dp.
        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in) :: name
        real(dp),            intent(in) :: x(:,:,:)
        integer,             intent(in) :: p
        real(dp),            intent(in) :: t1
        real(dp),            intent(in), optional :: t0
        integer :: h
        h = fields_handle(st, name)
        if (st%f(h)%rank == 3) then
            if (p == sp) then
                call fields_put(st, h, real(x, sp), t1=t1, t0=t0)
            else
                call fields_put(st, h, x, t1=t1, t0=t0)
            end if
        else
            if (p == sp) then
                call fields_put(st, h, real(x(:,:,1), sp), t1=t1, t0=t0)
            else
                call fields_put(st, h, x(:,:,1), t1=t1, t0=t0)
            end if
        end if
    end subroutine put3

    subroutine put0(st, name, v, p, t1)
        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in) :: name
        real(dp),            intent(in) :: v, t1
        integer,             intent(in) :: p
        if (p == sp) then
            call fields_put(st, fields_handle(st, name), real(v, sp), t1=t1)
        else
            call fields_put(st, fields_handle(st, name), v, t1=t1)
        end if
    end subroutine put0

    function got(st, name) result(x)
        ! The published values of a field as rank 3.
        type(field_state_t), intent(in) :: st
        character(len=*),    intent(in) :: name
        real(dp), allocatable :: x(:,:,:)
        x = st%f(fields_handle(st, name))%v
    end function got

    function rmp(x, gs, gd, method) result(y)
        ! The reference remap, called directly.
        real(dp),         intent(in) :: x(:,:,:)
        character(len=*), intent(in) :: gs, gd, method
        real(dp), allocatable :: y(:,:,:)
        call remap(cpl, x, gs, y, gd, method=method)
    end function rmp

    subroutine empty(st, owner)
        type(field_state_t), intent(inout) :: st
        character(len=*),    intent(in) :: owner
        call fields_init(st, owner)
        call fields_lock(st)
    end subroutine empty

    subroutine resolve1(w, row, ps, pd, dts, dtd, stat)
        ! Resolve one row between the toy producer (sx) and consumer (di).
        type(wire_t), allocatable, intent(out) :: w(:)
        type(wire_t),              intent(in)  :: row
        integer,                   intent(in)  :: ps, pd
        real(dp),                  intent(in), optional :: dts, dtd
        integer,                   intent(out), optional :: stat
        type(wire_node_t) :: nodes(2)
        real(dp) :: a, b

        a = 10.0_dp;  b = 10.0_dp
        if (present(dts)) a = dts
        if (present(dtd)) b = dtd
        call empty(si, "src:imports")
        call empty(dx, "dst:exports")
        nodes(1) = wire_node_states(row%src_kind, "src", si, sx, ps, a, check=.false., axes=ax)
        nodes(2) = wire_node_states(row%dst_kind, "dst", di, dx, pd, b, check=.false., axes=ax)
        call wire_resolve(w, nodes, cpl, table=[row], stat=stat)
        if (present(stat)) then
            if (stat /= 0) return
        end if
        call wire_prime(w, cpl)
    end subroutine resolve1

    ! =========================================================================
    ! the table
    ! =========================================================================

    subroutine tc_table()
        type(wire_t), allocatable :: t(:)
        integer :: i, n_gather, n_init

        write(*,"(a)") "table"
        t = wiring_table()
        n_gather = 0;  n_init = 0
        do i = 1, size(t)
            if (t(i)%gather) n_gather = n_gather + 1
            if (t(i)%lag == LAG_INIT) n_init = n_init + 1
        end do
        call check("50 rows", size(t) == 50)
        call check("two gather rows (dV_bsl, M_ice)", n_gather == 2)
        call check("static rows: ref x2, z_srf_ref, Q_geo, H_sed", n_init == 5)
        call check("geo -> ice z_bed: ADD_REF z_bed, lag 0, auto", &
                   any([(wire_label(t(i)) == "geo.dz_bed -> ice.z_bed" .and. t(i)%transform == T_ADD_REF &
                         .and. t(i)%aux == "z_bed" .and. t(i)%lag == 0 .and. t(i)%method == "auto", i = 1, size(t))]))
        call check("geo -> ice z_sl: ADD_SCALAR bsl", &
                   any([(wire_label(t(i)) == "geo.dz_ss -> ice.z_sl" .and. t(i)%transform == T_ADD_SCALAR &
                         .and. t(i)%aux == "bsl", i = 1, size(t))]))
        call check("ice.z_bed_ref -> geo.z_bed_ref, lag init, auto", &
                   any([(wire_label(t(i)) == "ice.z_bed_ref -> geo.z_bed_ref" .and. t(i)%lag == LAG_INIT &
                         .and. t(i)%method == "auto", i = 1, size(t))]))
        call check("surface.smb -> ice.smb: con, lag 1", &
                   any([(wire_label(t(i)) == "surface.smb -> ice.smb" .and. t(i)%lag == 1 .and. t(i)%method == "con", &
                         i = 1, size(t))]))
        call check("no dz_bed_dt_corr row (dropped in review)", &
                   .not. any([(t(i)%src_field == "dz_bed_dt_corr", i = 1, size(t))]))
    end subroutine tc_table

    ! =========================================================================
    ! resolution of a whole toy domain
    ! =========================================================================

    subroutine tc_domain()
        ! The mirror domain: climate on the coarse grid, everything else on the fine one.
        type(wire_t), allocatable :: w(:)
        integer :: i, stat

        write(*,"(a)") "resolution of a toy domain"
        call build()
        call wire_resolve(w, nodes, cpl, stat=stat)
        call check("resolves", stat == 0)
        call check("43 of 50 rows active", count(w%active) == 43 .and. size(w) == 50)
        call check("rule rows expanded: coef_alpha and dTo", &
                   any([(wire_label(w(i)) == "forcing.alpha -> climate.coef_alpha" .and. w(i)%active, i = 1, size(w))]) &
                   .and. any([(wire_label(w(i)) == "climate.dTo -> forcing.dTo" .and. w(i)%active, i = 1, size(w))]))
        call check("anomaly rows inactive on the mirror hub", &
                   .not. any([(w(i)%dst_field == "dz_bed" .and. w(i)%active, i = 1, size(w))]))
        call check("shelf set A rows inactive", .not. any([(w(i)%dst_field == "to_ann" .and. w(i)%active, i = 1, size(w))]))
        call check("methods: hub -> climate H_ice (fine -> coarse) con", w(row_of(w, "hub.H_ice -> climate.H_ice"))%mres == "con")
        call check("methods: climate -> surface tas bilin", w(row_of(w, "climate.tas -> surface.tas"))%mres == "bilin")
        call check("methods: surface -> ice smb on one grid is an identity", w(row_of(w, "surface.smb -> ice.smb"))%mres == "")
        call check("lag and precision carried", w(row_of(w, "surface.smb -> ice.smb"))%lag == 1 .and. &
                   w(row_of(w, "ice.M_ice -> forcing.M_ice"))%prec_dst == dp)
        call check("derived and gather flags", count(w%gather .and. w%active) == 2)
        call wiring_describe(w)

        ! the same domain with a finer hub: anomaly rows instead of the mirror
        call fields_init(hub_i, "hub:imports")
        call st_add(hub_i, "dz_bed", "m", I_, IT, GF)
        call st_add(hub_i, "dH_ice", "m", I_, EX, GF)
        call st_add(hub_i, "z_sl", "m", I_, IT, GF)
        call fields_lock(hub_i)
        call wire_resolve(w, nodes, cpl, stat=stat)
        call check("anomaly hub resolves", stat == 0)
        call check("anomaly rows active, mirror rows inactive", &
                   count([(w(i)%dst_kind == "hub" .and. w(i)%active, i = 1, size(w))]) == 3)
        call check("dz_bed row: SUB_REF", w(row_of(w, "ice.z_bed -> hub.dz_bed"))%active .and. &
                   w(row_of(w, "ice.z_bed -> hub.dz_bed"))%transform == T_SUB_REF)

    end subroutine tc_domain

    subroutine build()
        opt_ice = [character(len=FIELD_NAME_LEN) :: "T_shlf", "Q_geo", "H_sed"]
        opt_sh  = [character(len=FIELD_NAME_LEN) :: "to_ann", "so_ann", "dto_ann", "to_shlf", "so_shlf", "dto_shlf", &
                                                   "dso_shlf"]
        call empty(e_i, "none:imports");  call empty(e_e, "none:exports")

        call fields_init(ice_i, "ice:imports")
        call st_add(ice_i, "z_bed", "m", I_, IT, GF);  call st_add(ice_i, "z_sl", "m", I_, IT, GF)
        call st_add(ice_i, "smb", "kg m-2 yr-1", I_, EX, GF);  call st_add(ice_i, "T_srf", "K", I_, IT, GF)
        call st_add(ice_i, "bmb_shlf", "kg m-2 yr-1", I_, EX, GF);  call st_add(ice_i, "T_shlf", "K", I_, IT, GF)
        call st_add(ice_i, "Q_geo", "W m-2", I_, IT, GF);  call st_add(ice_i, "H_sed", "m", I_, IT, GF)
        call fields_lock(ice_i)
        call fields_init(ice_e, "ice:exports")
        call st_add(ice_e, "H_ice", "m", I_, EX, GF);  call st_add(ice_e, "z_srf", "m", I_, IT, GF)
        call st_add(ice_e, "z_bed", "m", I_, IT, GF);  call st_add(ice_e, "z_sl", "m", I_, IT, GF)
        call st_add(ice_e, "f_ice", "1", I_, IT, GF);  call st_add(ice_e, "f_grnd", "1", I_, IT, GF)
        call st_add(ice_e, "z_bed_ref", "m", I_, IT, GF);  call st_add(ice_e, "H_ice_ref", "m", I_, EX, GF)
        call st_add(ice_e, "M_ice", "Gt", I_, EX)
        call fields_lock(ice_e)

        call fields_init(geo_i, "geo:imports")
        call st_add(geo_i, "H_ice", "m", I_, EX, GF);  call st_add(geo_i, "bsl", "m", I_, IT)
        call st_add(geo_i, "z_bed_ref", "m", I_, IT, GF);  call st_add(geo_i, "H_ice_ref", "m", I_, EX, GF)
        call fields_lock(geo_i)
        call fields_init(geo_e, "geo:exports")
        call st_add(geo_e, "dz_bed", "m", I_, IT, GF);  call st_add(geo_e, "dz_ss", "m", I_, IT, GF)
        call st_add(geo_e, "bsl", "m", I_, IT);  call st_add(geo_e, "dV_bsl", "m3", I_, EX)
        call fields_lock(geo_e)

        call fields_init(gg_i, "geo_glob:imports");  call st_add(gg_i, "dV_bsl", "m3", I_, EX);  call fields_lock(gg_i)
        call fields_init(gg_e, "geo_glob:exports")
        call st_add(gg_e, "bsl", "m", I_, IT);  call st_add(gg_e, "A_ocn", "m2", I_, EX);  call fields_lock(gg_e)

        call fields_init(hub_i, "hub:imports")
        call st_add(hub_i, "H_ice", "m", I_, EX, GF);  call st_add(hub_i, "z_bed", "m", I_, IT, GF)
        call st_add(hub_i, "z_srf", "m", I_, IT, GF);  call st_add(hub_i, "f_ice", "1", I_, IT, GF)
        call st_add(hub_i, "f_grnd", "1", I_, IT, GF);  call st_add(hub_i, "z_sl", "m", I_, IT, GF)
        call fields_lock(hub_i)
        call fields_init(hub_e, "hub:exports")
        call st_add(hub_e, "z_srf", "m", I_, IT, GF);  call st_add(hub_e, "H_ice", "m", I_, EX, GF)
        call st_add(hub_e, "z_bed", "m", I_, IT, GF);  call st_add(hub_e, "f_grnd", "1", I_, IT, GF)
        call st_add(hub_e, "z_sl", "m", I_, IT, GF);  call st_add(hub_e, "z_srf_ref", "m", I_, IT, GF)
        call fields_lock(hub_e)

        call fields_init(cl_i, "climate:imports")
        call st_add(cl_i, "z_srf", "m", I_, IT, GC);  call st_add(cl_i, "H_ice", "m", I_, EX, GC)
        call st_add(cl_i, "z_bed", "m", I_, IT, GC);  call st_add(cl_i, "f_grnd", "1", I_, IT, GC)
        call st_add(cl_i, "z_sl", "m", I_, IT, GC);  call st_add(cl_i, "z_srf_ref", "m", I_, IT, GC)
        call st_add(cl_i, "coef_alpha", "1", I_, IT)
        call fields_lock(cl_i)
        call fields_init(cl_e, "climate:exports")
        call st_add(cl_e, "tas", "K", I_, IT, GC, nz=NM);  call st_add(cl_e, "pr", "kg m-2 yr-1", I_, EX, GC, nz=NM)
        call st_add(cl_e, "z_srf_clim", "m", I_, IT, GC);  call st_add(cl_e, "lapse", "K m-1", I_, IT, GC, nz=NM)
        call st_add(cl_e, "beta_p", "K-1", I_, IT, GC)
        call st_add(cl_e, "to_shlf", "K", I_, IT, GC);  call st_add(cl_e, "so_shlf", "psu", I_, IT, GC)
        call st_add(cl_e, "dto_shlf", "K", I_, IT, GC);  call st_add(cl_e, "dso_shlf", "psu", I_, IT, GC)
        call st_add(cl_e, "dTo", "K", I_, IT)
        call fields_lock(cl_e)

        call fields_init(sf_i, "surface:imports")
        call st_add(sf_i, "z_srf", "m", I_, IT, GF);  call st_add(sf_i, "H_ice", "m", I_, EX, GF)
        call st_add(sf_i, "z_bed", "m", I_, IT, GF);  call st_add(sf_i, "z_sl", "m", I_, IT, GF)
        call st_add(sf_i, "tas", "K", I_, IT, GF, nz=NM);  call st_add(sf_i, "pr", "kg m-2 yr-1", I_, EX, GF, nz=NM)
        call fields_lock(sf_i)
        call fields_init(sf_e, "surface:exports")
        call st_add(sf_e, "smb", "kg m-2 yr-1", M_, EX, GF);  call st_add(sf_e, "T_srf", "K", M_, IT, GF)
        call fields_lock(sf_e)

        call fields_init(sh_i, "shelf:imports")
        call st_add(sh_i, "H_ice", "m", I_, EX, GF);  call st_add(sh_i, "z_bed", "m", I_, IT, GF)
        call st_add(sh_i, "f_grnd", "1", I_, IT, GF);  call st_add(sh_i, "z_sl", "m", I_, IT, GF)
        call st_add(sh_i, "z_srf", "m", I_, IT, GF)
        call st_add(sh_i, "to_shlf", "K", I_, IT, GF);  call st_add(sh_i, "so_shlf", "psu", I_, IT, GF)
        call st_add(sh_i, "dto_shlf", "K", I_, IT, GF);  call st_add(sh_i, "dso_shlf", "psu", I_, IT, GF)
        call fields_lock(sh_i)
        call fields_init(sh_e, "shelf:exports")
        call st_add(sh_e, "bmb_shlf", "kg m-2 yr-1", M_, EX, GF);  call st_add(sh_e, "T_shlf", "K", I_, IT, GF)
        call fields_lock(sh_e)

        call fields_init(gt_e, "ghf:exports");  call st_add(gt_e, "Q_geo", "W m-2", I_, IT, GF);  call fields_lock(gt_e)
        call fields_init(sd_e, "sed:exports");  call st_add(sd_e, "H_sed", "m", I_, IT, GF);  call fields_lock(sd_e)

        call fields_init(fo_i, "forcing:imports")
        call st_add(fo_i, "M_ice", "Gt", I_, EX);  call st_add(fo_i, "dTo", "K", I_, IT)
        call fields_lock(fo_i)
        call fields_init(fo_e, "forcing:exports");  call st_add(fo_e, "alpha", "1", I_, IT);  call fields_lock(fo_e)

        nodes(1)  = wire_node_states("ice", "ice_grl", ice_i, ice_e, sp, 10.0_dp, opt=opt_ice)
        nodes(2)  = wire_node_states("geo", "geo_grl", geo_i, geo_e, sp, 10.0_dp)
        nodes(3)  = wire_node_states("geo_g", "geo_glob", gg_i, gg_e, sp, 10.0_dp, check=.false.)
        nodes(4)  = wire_node_states("hub", "Greenland", hub_i, hub_e, sp, 10.0_dp)
        nodes(5)  = wire_node_states("climate", "clim_grl", cl_i, cl_e, sp, 10.0_dp, axes=ax)
        nodes(6)  = wire_node_states("surface", "surf_grl", sf_i, sf_e, sp, 10.0_dp, axes=ax)
        nodes(7)  = wire_node_states("shelf", "shelf_grl", sh_i, sh_e, sp, 10.0_dp, opt=opt_sh)
        nodes(8)  = wire_node_states("geothermal", "ghf_grl", e_i, gt_e, sp, 0.0_dp)
        nodes(9)  = wire_node_states("sediments", "sed_grl", e_i, sd_e, sp, 0.0_dp)
        nodes(10) = wire_node_states("forcing", "forcing", fo_i, fo_e, dp, 10.0_dp, check=.false.)
    end subroutine build


    integer function row_of(w, label)
        type(wire_t), intent(in) :: w(:)
        character(len=*), intent(in) :: label
        integer :: i
        row_of = 1
        do i = 1, size(w)
            if (wire_label(w(i)) == label) then
                row_of = i
                return
            end if
        end do
        write(*,"(a)") "  FAIL  no row "//label
        nfail = nfail + 1
    end function row_of

    subroutine tc_resolve_errors()
        ! Error cases of resolution, through stat; the messages are checked in the abort cases.
        integer :: stat, h
        type(wire_t), allocatable :: w(:)
        type(wire_node_t) :: nodes2(2)
        type(field_state_t), target :: ei, ee
        character(len=FIELD_NAME_LEN), target :: opt(1)

        write(*,"(a)") "resolution errors (stat)"
        call empty(ei, "e:imports");  call empty(ee, "e:exports")

        ! gap: a required import with no producer; optional ones may stay unwired
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", I_, IT, GF);  call fields_lock(sx)
        call fields_init(di, "dst:imports")
        call st_add(di, "a", "m", I_, IT, GF);  call st_add(di, "b", "m", I_, IT, GF);  call fields_lock(di)
        nodes2(1) = wire_node_states("ice", "src", ei, sx, sp, 10.0_dp, check=.false.)
        nodes2(2) = wire_node_states("geo", "dst", di, ee, sp, 10.0_dp)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a")], stat=stat)
        call check("gap: import b has no producer", stat /= 0)
        opt(1) = "b"
        nodes2(2) = wire_node_states("geo", "dst", di, ee, sp, 10.0_dp, opt=opt)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a")], stat=stat)
        call check("optional import may stay unwired", stat == 0)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a"), wire("ice", "a", "geo", "a")], stat=stat)
        call check("two active rows into one field", stat /= 0)
        nodes2(2) = wire_node_states("geo", "dst", di, ee, sp, 10.0_dp, check=.false.)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a"), wire("ice", "a", "geo", "a")], stat=stat)
        call check("duplicates stop even when the node is not checked for gaps", stat /= 0)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a", lag=3)], stat=stat)
        call check("bad lag", stat /= 0)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a", transform=T_ADD_REF, aux="nope")], stat=stat)
        call check("unknown reference field", stat /= 0)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a", transform=T_ADD_SCALAR, aux="zz")], stat=stat)
        call check("missing 0-D aux", stat /= 0)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a", transform=T_LAPSE_T)], stat=stat)
        call check("lapse needs rank 3", stat /= 0)

        ! one grid, two shapes
        call fields_init(di, "dst:imports")
        call fields_add(di, "a", "m", I_, IT, h, grid=GF, dims=[NF - 1, NF])
        call fields_lock(di)
        nodes2(2) = wire_node_states("geo", "dst", di, ee, sp, 10.0_dp, check=.false.)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a")], stat=stat)
        call check("one grid, two shapes", stat /= 0)

        ! MEAN -> INSTANT needs equal cadences
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", M_, IT, GF);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "m", I_, IT, GF);  call fields_lock(di)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp, dts=10.0_dp, dtd=10.0_dp, stat=stat)
        call check("MEAN -> INSTANT, equal cadences", stat == 0)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp, dts=10.0_dp, dtd=20.0_dp, stat=stat)
        call check("MEAN -> INSTANT, cadences differ", stat /= 0)

        ! masked across two grids; mask on one side only
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", I_, IT, GF, masked=.true.);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "m", I_, IT, GC, masked=.true.);  call fields_lock(di)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp, stat=stat)
        call check("masked field across two grids", stat /= 0)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "m", I_, IT, GF);  call fields_lock(di)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp, stat=stat)
        call check("mask on one side only", stat /= 0)

        ! rank 3 axes
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "K", I_, IT, GF, nz=NM);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "K", I_, IT, GF, nz=NM);  call fields_lock(di)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp, stat=stat)
        call check("rank 3, same month axis", stat == 0)
        nodes2(1) = wire_node_states("ice", "src", ei, sx, sp, 10.0_dp, check=.false., axes=ax)
        ax2(1)%v(3) = 3.5_dp
        nodes2(2) = wire_node_states("geo", "dst", di, ee, sp, 10.0_dp, check=.false., axes=ax2)
        call wire_resolve(w, nodes2, cpl, table=[wire("ice", "a", "geo", "a")], stat=stat)
        call check("rank 3, axes differ", stat /= 0)
        ax2(1) = month_axis()
    end subroutine tc_resolve_errors

    ! =========================================================================
    ! exchange
    ! =========================================================================

    subroutine tc_identity()
        ! One grid: a copy; the validity interval travels; the consumer's precision rounds.
        real(dp) :: x(NF,NF,1)
        type(wire_t), allocatable :: w(:)
        integer :: p

        write(*,"(a)") "identity row"
        x = wave(NF, 1, 100.0_dp, 20.0_dp)
        do p = sp, dp, dp - sp
            call fields_init(sx, "src:exports");  call st_add(sx, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, GF)
            call fields_lock(sx)
            call fields_init(di, "dst:imports");  call st_add(di, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, GF)
            call fields_lock(di)
            call resolve1(w, wire("ice", "H_ice", "geo", "H_ice", lag=1), dp, p)
            call put3(sx, "H_ice", x, dp, t1=10.0_dp, t0=0.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=10.0_dp)
            call check(prec_name(p)//": value is the source rounded to the consumer's precision", &
                       all(got(di, "H_ice") == rnd(x, p)))
            call check(prec_name(p)//": validity [0,10] forwarded", &
                       di%f(1)%t0 == 0.0_dp .and. di%f(1)%t1 == 10.0_dp .and. di%f(1)%published)
            call check(prec_name(p)//": identity row has no method", len_trim(w(1)%mres) == 0)
        end do
    end subroutine tc_identity

    function prec_name(p) result(s)
        integer, intent(in) :: p
        character(len=2) :: s
        s = merge("sp", "dp", p == sp)
    end function prec_name

    subroutine tc_remap()
        ! Remapped rows: con conserves, the result is the coupler's remap rounded to the consumer's precision.
        real(dp) :: x(NF,NF,1), yc(NC,NC,1), xi, yi
        type(wire_t), allocatable :: w(:)
        real(dp), allocatable :: y(:,:,:)
        integer :: p

        write(*,"(a)") "remapped rows"
        x = wave(NF, 1, 100.0_dp, 20.0_dp)
        do p = sp, dp, dp - sp
            call fields_init(sx, "src:exports");  call st_add(sx, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, GF)
            call fields_lock(sx)
            call fields_init(di, "dst:imports");  call st_add(di, "H_ice", "m", FIELD_INSTANT, FIELD_EXTENSIVE, GC)
            call fields_lock(di)
            call resolve1(w, wire("ice", "H_ice", "geo", "H_ice", method="auto"), dp, p)
            call put3(sx, "H_ice", x, dp, t1=0.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=0.0_dp, t_last=0.0_dp, due=.true.)
            y = rmp(sx%f(1)%v, GF, GC, "con")
            call check(prec_name(p)//": auto, fine -> coarse, resolves to con", trim(w(1)%mres) == "con")
            call check(prec_name(p)//": = coupler con remap, rounded", all(got(di, "H_ice") == rnd(y, p)))
            yc = got(di, "H_ice");  xi = sum(x)*50.0_dp**2;  yi = sum(yc)*100.0_dp**2
            call check(prec_name(p)//": con conserves the integral", abs(yi - xi) <= merge(1.0e-6_dp, 1.0e-12_dp, p == sp)*abs(xi))
        end do

        call fields_init(sx, "src:exports");  call st_add(sx, "T", "K", FIELD_INSTANT, FIELD_INTENSIVE, GC);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "T", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call fields_lock(di)
        call resolve1(w, wire("ice", "T", "geo", "T", method="auto"), sp, sp)
        call put3(sx, "T", wave(NC, 1, 250.0_dp, 10.0_dp), sp, t1=0.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=0.0_dp, t_last=0.0_dp, due=.true.)
        call check("auto, coarse -> fine, resolves to bilin", trim(w(1)%mres) == "bilin")
        call check("bilin = coupler remap, rounded to sp", all(got(di, "T") == rnd(rmp(sx%f(1)%v, GC, GF, "bilin"), sp)))
        call resolve1(w, wire("ice", "T", "geo", "T"), sp, sp)
        call check("method '' follows cons (intensive: bilin)", trim(w(1)%mres) == "bilin")
        call fields_init(sx, "src:exports");  call st_add(sx, "T", "K", FIELD_INSTANT, FIELD_EXTENSIVE, GC);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "T", "K", FIELD_INSTANT, FIELD_EXTENSIVE, GF);  call fields_lock(di)
        call resolve1(w, wire("ice", "T", "geo", "T"), sp, sp)
        call check("method '' follows cons (extensive: con)", trim(w(1)%mres) == "con")
    end subroutine tc_remap

    subroutine tc_add_ref()
        ! dst = ref_dst + remap(src), in the consumer's precision; geo (coarse) -> ice (fine).
        real(dp) :: dz(NC,NC,1), ref(NF,NF)
        type(grid_ref_t) :: rf, rbad
        type(wire_t), allocatable :: w(:)
        real(dp), allocatable :: r(:,:,:), e(:,:,:)
        integer :: p, stat

        write(*,"(a)") "T_ADD_REF"
        dz = wave(NC, 1, 0.0_dp, 5.0_dp)
        ref = rnd(wave2(NF, -500.0_dp, 300.0_dp), sp)
        rf%grid = GF;  allocate(rf%z_bed(NF,NF));  rf%z_bed = ref
        do p = sp, dp, dp - sp
            call fields_init(sx, "src:exports");  call st_add(sx, "dz_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, GC)
            call fields_lock(sx)
            call fields_init(di, "dst:imports");  call st_add(di, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(di)
            call resolve1(w, wire("geo", "dz_bed", "ice", "z_bed", "auto", transform=T_ADD_REF, aux="z_bed"), p, p)
            call put3(sx, "dz_bed", dz, p, t1=10.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, ref_dst=rf, due=.true.)
            r = rmp(sx%f(1)%v, GC, GF, "bilin");  r = rnd(r, p)
            if (p == sp) then
                e = reshape(real(real(ref, sp) + real(r(:,:,1), sp), dp), [NF, NF, 1])
            else
                e = reshape(ref + r(:,:,1), [NF, NF, 1])
            end if
            call check(prec_name(p)//": z_bed = ref + remap(dz_bed), in "//prec_name(p), all(got(di, "z_bed") == e))
        end do
        call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true., stat=stat)
        call check("ref_dst missing is an error", stat /= 0)
        rbad%grid = GC;  allocate(rbad%z_bed(NC,NC));  rbad%z_bed = 0.0_dp
        call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, ref_dst=rbad, due=.true., stat=stat)
        call check("reference on the wrong grid is an error", stat /= 0)
    end subroutine tc_add_ref

    subroutine tc_add_scalar()
        ! dst = bsl + remap(dz_ss), with bsl the producer's own 0-D export.
        real(dp) :: dz(NC,NC,1), bsl
        type(wire_t), allocatable :: w(:)
        real(dp), allocatable :: r(:,:,:), e(:,:,:)
        integer :: p, stat

        write(*,"(a)") "T_ADD_SCALAR"
        dz = wave(NC, 1, 0.0_dp, 5.0_dp)
        bsl = -123.456789012345_dp
        do p = sp, dp, dp - sp
            call fields_init(sx, "src:exports")
            call st_add(sx, "dz_ss", "m", FIELD_INSTANT, FIELD_INTENSIVE, GC);  call st_add(sx, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE)
            call fields_lock(sx)
            call fields_init(di, "dst:imports");  call st_add(di, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(di)
            call resolve1(w, wire("geo", "dz_ss", "ice", "z_sl", "auto", transform=T_ADD_SCALAR, aux="bsl"), p, p)
            call put3(sx, "dz_ss", dz, p, t1=10.0_dp)
            call put0(sx, "bsl", bsl, p, t1=10.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true.)
            r = rnd(rmp(sx%f(1)%v, GC, GF, "bilin"), p)
            if (p == sp) then
                e = real(real(sx%f(2)%v(1,1,1), sp) + real(r, sp), dp)
            else
                e = sx%f(2)%v(1,1,1) + r
            end if
            call check(prec_name(p)//": z_sl = bsl + remap(dz_ss), in "//prec_name(p), all(got(di, "z_sl") == e))
        end do
        call put0(sx, "bsl", bsl, sp, t1=0.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true., stat=stat)
        call check("bsl out of step with dz_ss is an error", stat /= 0)
    end subroutine tc_add_scalar

    subroutine tc_sub_ref()
        ! dst = remap(src - ref_src): the difference in the source's precision, then the remap
        ! (ice fine -> hub coarse), then the consumer's rounding. All four combinations.
        real(dp) :: x(NF,NF,1), ref(NF,NF)
        type(grid_ref_t) :: rf
        type(wire_t), allocatable :: w(:)
        real(dp), allocatable :: pre(:,:,:), e(:,:,:)
        integer :: ps, pd

        write(*,"(a)") "T_SUB_REF"
        rf%grid = GF;  allocate(rf%z_bed(NF,NF))
        do ps = sp, dp, dp - sp
            do pd = sp, dp, dp - sp
                x = wave(NF, 1, 100.0_dp, 20.0_dp);  ref = wave2(NF, 90.0_dp, 25.0_dp)
                if (ps == sp) then
                    x = rnd(x, sp);  ref = rnd(ref, sp)
                end if
                rf%z_bed = ref
                call fields_init(sx, "src:exports");  call st_add(sx, "z_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
                call fields_lock(sx)
                call fields_init(di, "dst:imports");  call st_add(di, "dz_bed", "m", FIELD_INSTANT, FIELD_INTENSIVE, GC)
                call fields_lock(di)
                call resolve1(w, wire("ice", "z_bed", "hub", "dz_bed", "bilin", transform=T_SUB_REF, aux="z_bed"), ps, pd)
                call put3(sx, "z_bed", x, ps, t1=0.0_dp)
                call wire_exchange(w(1), sx, di, cpl, t=0.0_dp, t_last=0.0_dp, ref_src=rf, due=.true.)
                pre = x
                if (ps == sp) then
                    pre(:,:,1) = real(real(x(:,:,1), sp) - real(ref, sp), dp)
                else
                    pre(:,:,1) = x(:,:,1) - ref
                end if
                e = rnd(rmp(pre, GF, GC, "bilin"), pd)
                call check("src "//prec_name(ps)//", dst "//prec_name(pd)//": remap(src - ref) with src-precision difference", &
                           all(got(di, "dz_bed") == e))
            end do
        end do
    end subroutine tc_sub_ref

    subroutine tc_lapse()
        ! tas/pr re-elevated to the consumer's z_srf, on two grids and on one (no identical-grid shortcut).
        character(len=16) :: gs, gd
        integer :: p, ng, m, tr, ns, nd
        type(wire_t), allocatable :: w(:)
        real(dp), allocatable :: tas(:,:,:), pr(:,:,:), zc(:,:,:), bp(:,:,:), lp(:,:,:), z(:,:,:)
        real(dp), allocatable :: tr_(:,:,:), zc_(:,:,:), bp_(:,:,:), lp_(:,:,:), e(:,:,:)
        real(dp), allocatable :: dz_dp(:,:)
        real(sp), allocatable :: dz_sp(:,:)
        character(len=16) :: name

        write(*,"(a)") "T_LAPSE_T, T_LAPSE_P"
        do ng = 1, 2
            if (ng == 1) then
                gs = GC;  gd = GF;  ns = NC;  nd = NF
            else
                gs = GF;  gd = GF;  ns = NF;  nd = NF
            end if
            tas = wave(ns, NM, 270.0_dp, 15.0_dp);  pr = abs(wave(ns, NM, 300.0_dp, 100.0_dp))
            zc = reshape(wave(ns, 1, 1000.0_dp, 200.0_dp), [ns, ns, 1])
            bp = reshape(wave(ns, 1, 0.05_dp, 0.01_dp)*1.0e-2_dp, [ns, ns, 1])
            allocate(lp(ns, ns, NM))
            do m = 1, NM
                lp(:,:,m) = -0.008_dp*(1.0_dp + 0.01_dp*real(m, dp))
            end do
            z = reshape(wave(nd, 1, 1300.0_dp, 150.0_dp), [nd, nd, 1])
            do p = sp, dp, dp - sp
                do tr = T_LAPSE_T, T_LAPSE_P
                    name = merge("tas", "pr ", tr == T_LAPSE_T)
                    call fields_init(sx, "src:exports")
                    call st_add(sx, "tas", "K", I_, IT, gs, nz=NM);  call st_add(sx, "pr", "kg m-2 yr-1", I_, EX, gs, nz=NM)
                    call st_add(sx, "z_srf_clim", "m", I_, IT, gs);  call st_add(sx, "beta_p", "K-1", I_, IT, gs)
                    call st_add(sx, "lapse", "K m-1", I_, IT, gs, nz=NM)
                    call fields_lock(sx)
                    call fields_init(di, "dst:imports")
                    call st_add(di, "z_srf", "m", I_, IT, gd);  call st_add(di, "tas", "K", I_, IT, gd, nz=NM)
                    call st_add(di, "pr", "kg m-2 yr-1", I_, EX, gd, nz=NM)
                    call fields_lock(di)
                    call resolve1(w, wire("climate", trim(name), "surface", trim(name), "bilin", transform=tr), p, p)
                    call put3(sx, "tas", tas, p, t1=10.0_dp);  call put3(sx, "pr", pr, p, t1=10.0_dp)
                    call put3(sx, "z_srf_clim", zc, p, t1=10.0_dp);  call put3(sx, "beta_p", bp, p, t1=10.0_dp)
                    call put3(sx, "lapse", lp, p, t1=10.0_dp)
                    call put3(di, "z_srf", z, p, t1=10.0_dp)
                    call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true.)

                    ! the reference: bilin remaps of the producer's fields as stored, rounded, then the formulas
                    tr_ = merge_arr(tr == T_LAPSE_T, sx%f(1)%v, sx%f(2)%v)
                    tr_ = rnd(rmp2(tr_, gs, gd), p)
                    zc_ = rnd(rmp2(sx%f(3)%v, gs, gd), p);  bp_ = rnd(rmp2(sx%f(4)%v, gs, gd), p)
                    lp_ = rnd(rmp2(sx%f(5)%v, gs, gd), p)
                    e = tr_
                    if (p == sp) then
                        allocate(dz_sp(nd, nd))
                        dz_sp = real(di%f(1)%v(:,:,1), sp) - real(zc_(:,:,1), sp)
                        do m = 1, NM
                            if (tr == T_LAPSE_T) then
                                e(:,:,m) = real(real(tr_(:,:,m), sp) - real(lp_(:,:,m), sp)*dz_sp, dp)
                            else
                                e(:,:,m) = real(real(tr_(:,:,m), sp)*exp(-real(bp_(:,:,1), sp)*real(lp_(:,:,m), sp)*dz_sp), dp)
                            end if
                        end do
                        deallocate(dz_sp)
                    else
                        allocate(dz_dp(nd, nd))
                        dz_dp = di%f(1)%v(:,:,1) - zc_(:,:,1)
                        do m = 1, NM
                            if (tr == T_LAPSE_T) then
                                e(:,:,m) = tr_(:,:,m) - lp_(:,:,m)*dz_dp
                            else
                                e(:,:,m) = tr_(:,:,m)*exp(-bp_(:,:,1)*lp_(:,:,m)*dz_dp)
                            end if
                        end do
                        deallocate(dz_dp)
                    end if
                    call check(trim(name)//", "//trim(gs)//" -> "//trim(gd)//", "//prec_name(p)//": formula in the consumer's precision", &
                               all(got(di, trim(name)) == e))
                    call check(trim(name)//", "//trim(gs)//" -> "//trim(gd)//", "//prec_name(p)//": elevation changes the field", &
                               any(got(di, trim(name)) /= tr_))
                    if (ng == 2 .and. p == dp) then          ! z = zc: the identity, exactly (dz = 0, no shortcut needed)
                        call put3(di, "z_srf", zc, p, t1=10.0_dp)
                        call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true.)
                        call check(trim(name)//", one grid, z = z_srf_clim: unchanged, bitwise", &
                                   all(got(di, trim(name)) == tr_))
                    end if
                end do
            end do
            deallocate(lp)
        end do
    end subroutine tc_lapse

    function merge_arr(c, a, b) result(r)
        logical,  intent(in) :: c
        real(dp), intent(in) :: a(:,:,:), b(:,:,:)
        real(dp), allocatable :: r(:,:,:)
        if (c) then
            r = a
        else
            r = b
        end if
    end function merge_arr

    function rmp2(x, gs, gd) result(y)
        ! bilin remap, or a copy on one grid
        real(dp),         intent(in) :: x(:,:,:)
        character(len=*), intent(in) :: gs, gd
        real(dp), allocatable :: y(:,:,:)
        if (trim(gs) == trim(gd)) then
            y = x
        else
            y = rmp(x, gs, gd, "bilin")
        end if
    end function rmp2

    subroutine tc_gather()
        ! Sum of N domains' 0-D values published at t, in order, in the consumer's precision; stale ones count 0.
        type(field_state_t), target :: ex(3), gi, ge
        type(wire_node_t) :: srcs(3), nodes(2), gnode
        type(wire_t), allocatable :: w(:)
        type(wire_t) :: ws(3)
        real(dp) :: v(3), expect
        integer :: p, i, stat

        write(*,"(a)") "gather"
        v = [1.0_dp/3.0_dp, 1.0e8_dp + 0.1_dp, 2.0_dp/7.0_dp]
        call empty(ge, "geo_glob:exports")
        call empty(si, "x:imports")
        do p = sp, dp, dp - sp
            do i = 1, 3
                call fields_init(ex(i), "geo:exports");  call st_add(ex(i), "dV_bsl", "m3", FIELD_INSTANT, FIELD_EXTENSIVE)
                call fields_lock(ex(i))
            end do
            call fields_init(gi, "geo_glob:imports");  call st_add(gi, "dV_bsl", "m3", FIELD_INSTANT, FIELD_EXTENSIVE)
            call fields_lock(gi)
            gnode = wire_node_states("geo_g", "geo_glob", gi, ge, p, 10.0_dp, check=.false.)
            do i = 1, 3
                nodes(1) = wire_node_states("geo", "geo_d", si, ex(i), p, 10.0_dp, check=.false.)
                nodes(2) = gnode
                call wire_resolve(w, nodes, cpl, table=[wire("geo", "dV_bsl", "geo_g", "dV_bsl", gather=.true.)])
                ws(i) = w(1)
                srcs(i) = nodes(1)
            end do
            call put0(ex(1), "dV_bsl", v(1), p, t1=10.0_dp)
            call put0(ex(2), "dV_bsl", v(2), p, t1=10.0_dp)
            call put0(ex(3), "dV_bsl", v(3), p, t1=0.0_dp)          ! stale
            call wire_gather(ws, srcs, gi, 10.0_dp)
            if (p == sp) then
                expect = real(real(v(1), sp) + real(v(2), sp), dp)
            else
                expect = v(1) + v(2)
            end if
            call check(prec_name(p)//": sum of the values published at t, stale one counts 0", &
                       gi%f(1)%v(1,1,1) == expect .and. gi%f(1)%t1 == 10.0_dp)
            call put0(ex(3), "dV_bsl", v(3), p, t1=10.0_dp)
            call wire_gather(ws, srcs, gi, 10.0_dp)
            if (p == sp) then
                expect = real(real(real(v(1), sp) + real(v(2), sp), sp) + real(v(3), sp), dp)
            else
                expect = (v(1) + v(2)) + v(3)
            end if
            call check(prec_name(p)//": all three, in order, in "//prec_name(p), gi%f(1)%v(1,1,1) == expect)
            call wire_exchange(ws(1), ex(1), gi, cpl, 10.0_dp, 10.0_dp, stat=stat)
            call check(prec_name(p)//": a gather row refuses wire_exchange", stat /= 0)
        end do
    end subroutine tc_gather

    subroutine tc_mean_to_instant()
        ! A MEAN export read as INSTANT, equal cadences: the window travels as the validity.
        real(dp) :: x(NF,NF,1)
        type(wire_t), allocatable :: w(:)

        write(*,"(a)") "MEAN -> INSTANT"
        x = rnd(wave(NF, 1, 300.0_dp, 50.0_dp), sp)
        call fields_init(sx, "surface:exports")
        call st_add(sx, "smb", "kg m-2 yr-1", FIELD_MEAN, FIELD_EXTENSIVE, GF);  call fields_lock(sx)
        call fields_init(di, "ice:imports")
        call st_add(di, "smb", "kg m-2 yr-1", FIELD_INSTANT, FIELD_EXTENSIVE, GF);  call fields_lock(di)
        call resolve1(w, wire("surface", "smb", "ice", "smb", "con", lag=1), sp, sp, dts=10.0_dp, dtd=10.0_dp)
        call put3(sx, "smb", x, sp, t1=10.0_dp, t0=0.0_dp)
        call put3(sx, "smb", x + 2.0_dp, sp, t1=20.0_dp, t0=10.0_dp)
        call fields_finalize_interval(sx, 20.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=30.0_dp, t_last=20.0_dp, due=.false.)
        call check("window mean forwarded", all(got(di, "smb") == rnd(sx%f(1)%v, sp)))
        call check("validity is the window [0,20]", di%f(1)%t0 == 0.0_dp .and. di%f(1)%t1 == 20.0_dp)
    end subroutine tc_mean_to_instant

    subroutine tc_init_once()
        ! A LAG_INIT row exchanges once: later calls do nothing, even with a changed source.
        real(dp) :: x(NF,NF,1)
        type(wire_t), allocatable :: w(:)

        write(*,"(a)") "LAG_INIT"
        x = rnd(wave(NF, 1, 0.05_dp, 0.01_dp), sp)
        call fields_init(sx, "ghf:exports");  call st_add(sx, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, GF)
        call fields_lock(sx)
        call fields_init(di, "ice:imports");  call st_add(di, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, GF)
        call fields_lock(di)
        call resolve1(w, wire("geothermal", "Q_geo", "ice", "Q_geo", "con", lag=LAG_INIT), sp, sp, dts=0.0_dp)
        call put3(sx, "Q_geo", x, sp, t1=0.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=0.0_dp, t_last=-5.0_dp, init=.true.)   ! static: t_last is not checked
        call check("first exchange copies", all(got(di, "Q_geo") == x) .and. w(1)%done)
        call put3(sx, "Q_geo", 2.0_dp*x, sp, t1=0.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=-5.0_dp)
        call check("second exchange is skipped", all(got(di, "Q_geo") == x))
    end subroutine tc_init_once

    subroutine tc_rank0()
        type(wire_t), allocatable :: w(:)

        write(*,"(a)") "0-D row"
        call fields_init(sx, "gg:exports");  call st_add(sx, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE);  call fields_lock(sx)
        call fields_init(di, "geo:imports");  call st_add(di, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE);  call fields_lock(di)
        call resolve1(w, wire("geo_g", "bsl", "geo", "bsl"), dp, sp)
        call fields_put(sx, 1, -120.0_dp + 1.0_dp/3.0_dp, t1=20.0_dp, t0=10.0_dp)
        call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=20.0_dp, due=.true.)
        call check("0-D value bit-exact, validity forwarded", di%f(1)%v(1,1,1) == -120.0_dp + 1.0_dp/3.0_dp &
                   .and. di%f(1)%t0 == 10.0_dp .and. di%f(1)%t1 == 20.0_dp)
    end subroutine tc_rank0

    subroutine tc_masked()
        type(wire_t), allocatable :: w(:)
        logical :: valid(NF,NF)
        integer :: i

        write(*,"(a)") "masked field, one grid"
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF, masked=.true.)
        call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF, masked=.true.)
        call fields_lock(di)
        call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
        valid = .true.
        do i = 1, NF
            valid(i,i) = .false.
        end do
        call fields_put(sx, 1, real(wave2(NF, 1.0_dp, 1.0_dp), sp), t1=0.0_dp, valid=valid)
        call wire_exchange(w(1), sx, di, cpl, t=0.0_dp, t_last=0.0_dp, due=.true.)
        call check("valid mask forwarded", all(di%f(1)%valid .eqv. valid))
    end subroutine tc_masked

    ! =========================================================================
    ! lag check: every case of wire_lag_status
    ! =========================================================================

    subroutine tc_lag()
        type(wire_t), allocatable :: w0(:), w1(:), wi(:), wa(:), wl(:)
        real(dp) :: x(NF,NF,1)
        real(dp), parameter :: t = 30.0_dp

        write(*,"(a)") "lag check"
        x = wave(NF, 1, 1.0_dp, 1.0_dp)
        call fields_init(sx, "src:exports")
        call st_add(sx, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
        call st_add(sx, "dz", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call st_add(sx, "bsl", "m", FIELD_INSTANT, FIELD_INTENSIVE)
        call fields_lock(sx)
        call fields_init(di, "dst:imports")
        call st_add(di, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call st_add(di, "dz", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
        call fields_lock(di)
        call resolve1(w0, wire("ice", "a", "geo", "a", lag=0), sp, sp)
        call resolve1(w1, wire("ice", "a", "geo", "a", lag=1), sp, sp)
        call resolve1(wi, wire("ice", "a", "geo", "a", lag=LAG_INIT), sp, sp)
        call resolve1(wa, wire("ice", "dz", "geo", "dz", transform=T_ADD_SCALAR, aux="bsl"), sp, sp)

        call check("unpublished source", len(wire_lag_status(w0(1), sx, di, t, t)) > 0)
        call put3(sx, "a", x, sp, t1=t)
        call check("lag 0, producer ran (t1 = t = t_last)", len(wire_lag_status(w0(1), sx, di, t, t, due=.true.)) == 0)
        call check("lag 0, not due and ran earlier", len(wire_lag_status(w0(1), sx, di, t + 10.0_dp, t)) == 0)
        call check("lag 0, due but has not run (swapped order)", &
                   len(wire_lag_status(w0(1), sx, di, t + 10.0_dp, t, due=.true.)) > 0)
        call check("t1 /= t_last: put without cadence_mark", len(wire_lag_status(w0(1), sx, di, t, t - 10.0_dp)) > 0)
        call check("t1 /= t_last: mark without put", len(wire_lag_status(w0(1), sx, di, t + 10.0_dp, t + 10.0_dp)) > 0)
        call check("t1 > t: producer ahead of the clock", len(wire_lag_status(w0(1), sx, di, t - 10.0_dp, t)) > 0)
        call check("lag 1, previous step's value", len(wire_lag_status(w1(1), sx, di, t + 10.0_dp, t)) == 0)
        call check("lag 1, producer already ran this step (consumer after producer)", &
                   len(wire_lag_status(w1(1), sx, di, t, t)) > 0)
        call check("lag 1 at the cold start counts as lag 0", len(wire_lag_status(w1(1), sx, di, t, t, init=.true.)) == 0)
        call check("lag 0 at the cold start unchanged", len(wire_lag_status(w0(1), sx, di, t + 10.0_dp, t, init=.true.)) == 0)
        call check("LAG_INIT: any t_last", len(wire_lag_status(wi(1), sx, di, t, -99.0_dp)) == 0)
        call check("LAG_INIT: t1 > t", len(wire_lag_status(wi(1), sx, di, t - 10.0_dp, -99.0_dp)) > 0)

        call put3(sx, "dz", x, sp, t1=t);  call put0(sx, "bsl", 1.0_dp, sp, t1=t)
        call check("aux in step", len(wire_lag_status(wa(1), sx, di, t, t, due=.true.)) == 0)
        call put0(sx, "bsl", 1.0_dp, sp, t1=t - 10.0_dp)
        call check("aux out of step", len(wire_lag_status(wa(1), sx, di, t, t, due=.true.)) > 0)

        ! lapse: the consumer's z_srf must have been exchanged
        call fields_init(sx, "src:exports")
        call st_add(sx, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
        call st_add(sx, "z_srf_clim", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call st_add(sx, "beta_p", "K-1", FIELD_INSTANT, FIELD_INTENSIVE, GF)
        call st_add(sx, "lapse", "K m-1", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
        call fields_lock(sx)
        call fields_init(di, "dst:imports")
        call st_add(di, "z_srf", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call st_add(di, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
        call fields_lock(di)
        call resolve1(wl, wire("climate", "tas", "surface", "tas", transform=T_LAPSE_T), sp, sp)
        call put3(sx, "tas", wave(NF, NM, 1.0_dp, 1.0_dp), sp, t1=t);  call put3(sx, "z_srf_clim", x, sp, t1=t)
        call put3(sx, "beta_p", x, sp, t1=t);  call put3(sx, "lapse", wave(NF, NM, 1.0_dp, 1.0_dp), sp, t1=t)
        call check("lapse: consumer z_srf not yet exchanged", len(wire_lag_status(wl(1), sx, di, t, t, due=.true.)) > 0)
        call put3(di, "z_srf", x, sp, t1=t)
        call check("lapse: ready", len(wire_lag_status(wl(1), sx, di, t, t, due=.true.)) == 0)
    end subroutine tc_lag

    ! =========================================================================
    ! the aborts, for real: the program re-runs itself, one case per child
    ! =========================================================================

    subroutine tc_aborts()
        write(*,"(a)") "aborts (child processes)"
        call expect_abort("lag_swapped", "runs before the producer")
        call expect_abort("lag_after", "runs after the producer")
        call expect_abort("lag_unmarked", "last update")
        call expect_abort("lag_unpublished", "never published")
        call expect_abort("lag_ahead", "ahead of the exchange time")
        call expect_abort("lapse_z", "not yet exchanged")
        call expect_abort("gap", "has no producer")
        call expect_abort("shelf_none", "needs ocean forcing")
        call expect_abort("units", "units differ")
        call expect_abort("no_ref", "needs ref_dst")
        call expect_abort("inactive", "is inactive")
    end subroutine tc_aborts

    subroutine expect_abort(case_name, text)
        character(len=*), intent(in) :: case_name, text
        character(len=512) :: exe, line
        character(len=:), allocatable :: all_
        integer :: es, u, ios

        call get_command_argument(0, exe)
        call execute_command_line(trim(exe)//" "//case_name//" > /dev/null 2> "//errfile, exitstat=es)
        all_ = ""
        open(newunit=u, file=errfile, status="old", action="read", iostat=ios)
        if (ios == 0) then
            do
                read(u, "(a)", iostat=ios) line
                if (ios /= 0) exit
                all_ = all_//trim(line)//" "
            end do
            close(u)
        end if
        call check(case_name//": aborts (exit 1) with '"//text//"'", es == 1 .and. index(all_, text) > 0)
    end subroutine expect_abort

    subroutine abort_case()
        ! Child: build the situation and call the procedure without stat.
        character(len=64) :: arg
        type(wire_t), allocatable :: w(:)
        type(wire_node_t) :: nodes(2)
        type(field_state_t), target :: ei, ee

        call get_command_argument(1, arg)
        call empty(ei, "e:imports");  call empty(ee, "e:exports")
        call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call fields_lock(sx)
        call fields_init(di, "dst:imports");  call st_add(di, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF);  call fields_lock(di)
        select case (trim(arg))
        case ("lag_swapped")
            call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
            call put3(sx, "a", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=10.0_dp, due=.true.)
        case ("lag_after")
            call resolve1(w, wire("ice", "a", "geo", "a", lag=1), sp, sp)
            call put3(sx, "a", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=20.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=20.0_dp)
        case ("lag_unmarked")
            call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
            call put3(sx, "a", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=20.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=10.0_dp)
        case ("lag_unpublished")
            call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=20.0_dp)
        case ("lag_ahead")
            call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
            call put3(sx, "a", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=30.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=20.0_dp, t_last=30.0_dp)
        case ("lapse_z")
            call fields_init(sx, "src:exports")
            call st_add(sx, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
            call st_add(sx, "z_srf_clim", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call st_add(sx, "beta_p", "K-1", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call st_add(sx, "lapse", "K m-1", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
            call fields_lock(sx)
            call fields_init(di, "dst:imports")
            call st_add(di, "z_srf", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call st_add(di, "tas", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF, nz=NM)
            call fields_lock(di)
            call resolve1(w, wire("climate", "tas", "surface", "tas", transform=T_LAPSE_T), sp, sp)
            call put3(sx, "tas", wave(NF, NM, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call put3(sx, "z_srf_clim", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call put3(sx, "beta_p", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call put3(sx, "lapse", wave(NF, NM, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true.)
        case ("gap")
            call fields_init(di, "dst:imports");  call st_add(di, "b", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(di)
            nodes(1) = wire_node_states("ice", "ice_grl", ei, sx, sp, 10.0_dp)
            nodes(2) = wire_node_states("geo", "geo_grl", di, ee, sp, 10.0_dp)
            call wire_resolve(w, nodes, cpl, table=[wire("ice", "a", "geo", "a")])
        case ("shelf_none")
            call fields_init(di, "shelf:imports");  call st_add(di, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(di)
            call fields_init(sx, "hub:exports");  call st_add(sx, "z_sl", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(sx)
            nodes(1) = wire_node_states("hub", "Greenland", ei, sx, sp, 10.0_dp)
            nodes(2) = wire_node_states("shelf", "shelf_grl", di, ee, sp, 10.0_dp)
            call wire_resolve(w, nodes, cpl, table=[wire("hub", "z_sl", "shelf", "z_sl")])
        case ("units")
            call fields_init(di, "dst:imports");  call st_add(di, "a", "K", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(di)
            call resolve1(w, wire("ice", "a", "geo", "a"), sp, sp)
        case ("no_ref")
            call fields_init(sx, "src:exports");  call st_add(sx, "a", "m", FIELD_INSTANT, FIELD_INTENSIVE, GF)
            call fields_lock(sx)
            call resolve1(w, wire("ice", "a", "geo", "a", transform=T_ADD_REF, aux="z_bed"), sp, sp)
            call put3(sx, "a", wave(NF, 1, 1.0_dp, 1.0_dp), sp, t1=10.0_dp)
            call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp, due=.true.)
        case ("inactive")
            call resolve1(w, wire("ice", "zz", "geo", "a"), sp, sp)
            call wire_exchange(w(1), sx, di, cpl, t=10.0_dp, t_last=10.0_dp)
        case default
            error stop 78
        end select
    end subroutine abort_case

end program test_wiring
