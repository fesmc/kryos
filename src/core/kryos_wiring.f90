module kryos_wiring
    ! The exchange engine and the fixed phase 1 wiring table (contracts-design
    ! section 4, core-design section 1): wire_t rows, the table, wire_resolve at
    ! init, wire_exchange on two field states, the lag check and wiring_describe.
    !
    ! Plain types and free procedures; no dependency on domain, globe or run.
    ! The engine works on field states and a coupler passed in: the producer's
    ! exports, the consumer's imports, and whatever the transform reads (a
    ! reference topography on a grid, grid_ref_t) is handed in by the caller. A
    ! wire_node_t describes one side for wire_resolve (an instance or the hub).
    !
    ! Exchanged values are `dp` in the field states. Transforms are evaluated in
    ! the consumer's declared precision with yelmox's operand order (section 9,
    ! decision 2026-10-25): a row never rounds (the consumer narrows once, after
    ! its own conversion); a T_SUB_REF difference is formed in the source's prec;
    ! every other transform runs in the destination's prec, an sp kernel
    ! narrowing its operands (the remap of widened sp values, narrowed, is the
    ! coupler's sp remap); T_LAPSE_P applies its sp factor to the value unnarrowed
    ! (pr carries a unit conversion the consumer reverses). The exchange never
    ! converts units. Failures stop through fesm_error naming the row; the procedures that
    ! can fail take an optional `stat` (fesm_error then returns instead).

    use precision,       only: sp, dp
    use fesm_log,        only: LOG_INFO, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t, field_t, fields_has, fields_handle, fields_put, &
                               fields_check_link, FIELD_INSTANT, FIELD_MEAN, FIELD_EXTENSIVE, &
                               FIELD_NAME_LEN
    use coupler,         only: coupler_class, coupler_method, coupler_prime, remap
    use kryos_component, only: component_t, axis_t, grid_ref_t, axis_find, axes_equal, NAME_LEN

    implicit none

    private

    integer, parameter, public :: T_NONE       = 0   ! dst = remap(src)
    integer, parameter, public :: T_ADD_REF    = 1   ! dst = ref_dst%<aux> + remap(src)             (decision 4)
    integer, parameter, public :: T_SUB_REF    = 2   ! dst = remap(src - ref_src%<aux>)            (hub anomalies)
    integer, parameter, public :: T_ADD_SCALAR = 3   ! dst = src_state%<aux> (0-D) + remap(src)     (z_sl = dphi + dz_ss)
    integer, parameter, public :: T_LAPSE_T    = 4   ! tas re-elevated to the consumer's z_srf
    integer, parameter, public :: T_LAPSE_P    = 5   ! pr  re-elevated to the consumer's z_srf
    integer, parameter, public :: T_GIA        = 6   ! the core's gather into a host-resolved global geo (phase 4)
    integer, parameter, public :: LAG_INIT     = -1  ! exchanged once after init (static fields)

    ! The other sources of the T_GIA rows: the step terms, and the reference terms.
    character(len=FIELD_NAME_LEN), parameter :: GIA_STEP(3) = [character(len=FIELD_NAME_LEN) :: "f_grnd", "z_bed", "z_sl"]
    character(len=FIELD_NAME_LEN), parameter :: GIA_REF(2)  = [character(len=FIELD_NAME_LEN) :: "z_bed_ref", "footprint"]

    ! Table rules: the destination's declared names select the rows (section 4.2).
    character(len=*), parameter, public :: RULE_COEF = "<coef>", RULE_DERIVED = "<derived>"

    type, public :: wire_t                           ! kinds: ice geo climate surface shelf geothermal sediments hub;
        character(len=16) :: src_kind = "", dst_kind = ""      ! globe: forcing, geo_g (global part of geo)
        character(len=FIELD_NAME_LEN) :: src_field = "", dst_field = "", aux = ""   ! aux: ref field or 0-D src field
        character(len=FIELD_NAME_LEN), allocatable :: srcs(:)   ! T_GIA: the other source fields the core reads
        character(len=8)  :: method = ""             ! "" (from cons) | con | bilin | auto
        integer           :: lag = 0, transform = T_NONE   ! lag: 0 | 1 | LAG_INIT
        logical           :: gather = .false.        ! N domains -> one global 0-D field, summed (wire_gather)
        ! resolved by wire_resolve
        logical           :: active = .false.
        character(len=160) :: reason = ""            ! why inactive
        character(len=NAME_LEN) :: src_name = "", dst_name = "", src_grid = "", dst_grid = ""
        character(len=32) :: mres = ""               ! method after resolution; "" for identity and 0-D rows
        integer           :: hs = 0, hd = 0, ha = 0, prec_src = 0, prec_dst = 0
        integer           :: hl(4) = 0               ! lapse rows: src z_srf_clim, beta_p, lapse; dst z_srf
        integer, allocatable :: hx(:)                   ! T_GIA rows: handles of the other source fields (srcs)
        logical           :: done = .false.          ! LAG_INIT row exchanged
        real(dp), allocatable :: work(:,:,:), pre(:,:,:)       ! persistent workspaces
        real(dp), allocatable :: zc(:,:,:), bp(:,:,:), lp(:,:,:) ! lapse: remapped z_srf_clim, beta_p, lapse
    end type wire_t

    type, public :: wire_node_t                      ! one side of the rows: an instance or the hub
        character(len=16)       :: kind = ""         ! table kind; geo of scope global is "geo_g"
        character(len=NAME_LEN) :: name = ""         ! instance (hub: domain) name, for messages
        integer                 :: prec = sp         ! kind of the instance's arrays
        real(dp)                :: dt = 0.0_dp       ! cadence step [yr] (dtt for every step, 0 off)
        logical                 :: check = .true.    ! required imports must be wired (.false.: globe node
                                                     ! seen from one domain)
        type(field_state_t), pointer :: imports => null(), exports => null()
        type(axis_t),        pointer :: axes(:) => null()
        character(len=FIELD_NAME_LEN), pointer :: opt(:) => null()   ! imports that may stay unwired
        character(len=FIELD_NAME_LEN), pointer :: core(:) => null()  ! imports the core fills outside the rows
    end type wire_node_t

    public :: wire, wiring_table, wire_node, wire_node_states
    public :: wire_resolve, wire_prime, wire_exchange, wire_gather, wire_lag_status
    public :: wiring_describe, wire_label, wire_runs

contains

    ! =========================================================================
    ! Rows and the table
    ! =========================================================================

    pure function wire(src_kind, src_field, dst_kind, dst_field, method, lag, transform, aux, gather, srcs) result(w)
        character(len=*), intent(in)           :: src_kind, src_field, dst_kind, dst_field
        character(len=*), intent(in), optional :: method, aux, srcs(:)
        integer,          intent(in), optional :: lag, transform
        logical,          intent(in), optional :: gather
        type(wire_t) :: w

        w%src_kind = src_kind;  w%src_field = src_field
        w%dst_kind = dst_kind;  w%dst_field = dst_field
        if (present(method))    w%method    = method
        if (present(aux))       w%aux       = aux
        if (present(srcs))      w%srcs      = srcs
        if (present(lag))       w%lag       = lag
        if (present(transform)) w%transform = transform
        if (present(gather))    w%gather    = gather
    end function wire

    function wiring_table() result(t)
        ! The table, contracts-design section 4.2, in table order (the order the
        ! core runs the rows of one src -> dst pair).
        type(wire_t), allocatable :: t(:)

        allocate(t(0))
        ! geo reads the ice (decision 4 refined: its reference is the ice instance's)
        call add(t, "ice", "H_ice", "geo", "H_ice", "auto", lag=1)
        call add_set(t, "ice", "z_bed_ref H_ice_ref", "geo", "auto", lag=LAG_INIT)
        call add(t, "geo_g", "bsl", "geo", "bsl")
        call add(t, "geo_g", "A_ocn", "geo", "A_ocn")
        call add(t, "geo", "dV_bsl", "geo_g", "dV_bsl", gather=.true.)
        ! geo = vilma2 (phase 4): the scatter into the regional ports and the gather of the loads. The rows
        ! exist for resolution, logging and lag checks; the arithmetic is the core's (gia_scatter, gia_gather).
        call add(t, "geo_g", "dz_bed", "geo", "dz_bed", "bilin")
        call add(t, "geo_g", "dz_ss", "geo", "dz_ss", "bilin")
        call add(t, "geo_g", "dphi", "geo", "dphi")
        ! (src_field is the principal source, srcs the others the gather reads)
        call add(t, "ice", "H_ice", "geo_g", "host_w", lag=1, transform=T_GIA, srcs=GIA_STEP)
        call add(t, "ice", "H_ice", "geo_g", "host_Ig", lag=1, transform=T_GIA, srcs=GIA_STEP)
        call add(t, "ice", "H_ice", "geo_g", "host_C", lag=1, transform=T_GIA, srcs=GIA_STEP)
        call add(t, "ice", "H_ice", "geo_g", "host_F", lag=1, transform=T_GIA, srcs=GIA_STEP)
        call add(t, "ice", "H_ice_pd", "geo_g", "host_Ig_ref", lag=LAG_INIT, transform=T_GIA, srcs=GIA_REF)
        call add(t, "ice", "H_ice_pd", "geo_g", "host_C_ref", lag=LAG_INIT, transform=T_GIA, srcs=GIA_REF)
        ! geo -> ice: bed and sea level in one datum
        call add(t, "geo", "dz_bed", "ice", "z_bed", "auto", transform=T_ADD_REF, aux="z_bed")
        call add(t, "geo", "dz_ss", "ice", "z_sl", "auto", transform=T_ADD_SCALAR, aux="dphi")
        ! ice -> hub: mirror set (same grid), anomaly set (finer hub); z_sl serves both (identity on one grid)
        call add_set(t, "ice", "H_ice z_bed z_srf f_ice f_grnd", "hub", "")
        call add(t, "ice", "z_sl", "hub", "z_sl", "bilin")
        call add(t, "ice", "z_bed", "hub", "dz_bed", "bilin", transform=T_SUB_REF, aux="z_bed")
        call add(t, "ice", "H_ice", "hub", "dH_ice", "bilin", transform=T_SUB_REF, aux="H_ice")
        ! hub -> climate, surface, shelf
        call add_set(t, "hub", "z_srf H_ice z_bed f_grnd z_sl", "climate", "auto")
        call add(t, "hub", "z_srf_ref", "climate", "z_srf_ref", "auto", lag=LAG_INIT)
        call add_set(t, "hub", "z_srf H_ice z_bed z_sl", "surface", "auto")
        call add_set(t, "hub", "H_ice z_bed f_grnd z_sl z_srf", "shelf", "bilin")
        ! climate -> surface (re-elevated), shelf
        call add(t, "climate", "tas", "surface", "tas", "bilin", transform=T_LAPSE_T)
        call add(t, "climate", "pr", "surface", "pr", "bilin", transform=T_LAPSE_P)
        ! surface = smb_simple: the sea-level air temperature, the reference ice of its mask
        call add(t, "climate", "tsl_ann", "surface", "tsl_ann", "bilin")
        call add(t, "ice", "H_ice_ref", "surface", "H_ice_ref", "bilin")
        call add_set(t, "climate", "to_ann so_ann dto_ann", "shelf", "bilin")
        call add_set(t, "climate", "to_shlf so_shlf dto_shlf dso_shlf", "shelf", "bilin")
        ! lagged boundary conditions into the ice
        call add_set(t, "surface", "smb T_srf", "ice", "con", lag=1)
        call add_set(t, "shelf", "bmb_shlf T_shlf", "ice", "con", lag=1)
        call add(t, "climate", "Qd", "ice", "Qd", "con", lag=1)
        call add(t, "geothermal", "Q_geo", "ice", "Q_geo", "con", lag=LAG_INIT)
        call add(t, "sediments", "H_sed", "ice", "H_sed", "con", lag=LAG_INIT)
        ! globe
        call add(t, "ice", "M_ice", "forcing", "M_ice", lag=1, gather=.true.)
        call add(t, "climate", RULE_DERIVED, "forcing", RULE_DERIVED, lag=1)
        call add(t, "forcing", RULE_COEF, "climate", RULE_COEF)
    end function wiring_table

    subroutine add(t, src_kind, src_field, dst_kind, dst_field, method, lag, transform, aux, gather, srcs)
        type(wire_t), allocatable, intent(inout) :: t(:)
        character(len=*), intent(in)           :: src_kind, src_field, dst_kind, dst_field
        character(len=*), intent(in), optional :: method, aux, srcs(:)
        integer,          intent(in), optional :: lag, transform
        logical,          intent(in), optional :: gather
        call append(t, wire(src_kind, src_field, dst_kind, dst_field, method, lag, transform, aux, gather, srcs))
    end subroutine add

    subroutine add_set(t, src_kind, fields, dst_kind, method, lag)
        ! One row per name of the blank-separated `fields`; same name on both sides.
        type(wire_t), allocatable, intent(inout) :: t(:)
        character(len=*), intent(in) :: src_kind, fields, dst_kind, method
        integer,          intent(in), optional :: lag
        integer :: p0, p1

        p0 = 1
        do while (p0 <= len_trim(fields))
            p1 = index(fields(p0:), " ")
            if (p1 == 0) p1 = len_trim(fields) - p0 + 2
            call append(t, wire(src_kind, fields(p0:p0+p1-2), dst_kind, fields(p0:p0+p1-2), method, lag))
            p0 = p0 + p1
        end do
    end subroutine add_set

    subroutine append(list, row)
        type(wire_t), allocatable, intent(inout) :: list(:)
        type(wire_t),              intent(in)    :: row
        type(wire_t), allocatable :: tmp(:)
        integer :: n

        if (.not. allocated(list)) allocate(list(0))
        n = size(list)
        allocate(tmp(n+1))
        tmp(1:n) = list
        tmp(n+1) = row
        call move_alloc(tmp, list)
    end subroutine append

    ! =========================================================================
    ! Nodes
    ! =========================================================================

    function wire_node(c, dt, check) result(n)
        ! The node of a component instance. `c` must have the target attribute
        ! (the node points into its field states); call after init (axes, opt).
        class(component_t), target, intent(in) :: c
        real(dp),                   intent(in) :: dt      ! cadence_dt_years
        logical,          intent(in), optional :: check
        type(wire_node_t) :: n

        n%kind = c%kind(1:16)
        if (trim(c%scope) == "global" .and. trim(c%kind) /= "forcing") n%kind = trim(c%kind)//"_g"
        n%name = c%name
        n%prec = c%prec
        n%dt   = dt
        if (present(check)) n%check = check
        n%imports => c%imports
        n%exports => c%exports
        if (allocated(c%axes))        n%axes => c%axes
        if (allocated(c%opt_imports)) n%opt  => c%opt_imports
        if (allocated(c%core_imports)) n%core => c%core_imports
    end function wire_node

    function wire_node_states(kind, name, imports, exports, prec, dt, check, axes, opt, core) result(n)
        ! A node from bare field states (the hub; toy states). The actual
        ! arguments must have the target attribute.
        character(len=*),            intent(in) :: kind, name
        type(field_state_t), target, intent(in) :: imports, exports
        integer,                     intent(in) :: prec
        real(dp),                    intent(in) :: dt
        logical,                     intent(in), optional :: check
        type(axis_t), target,        intent(in), optional :: axes(:)
        character(len=FIELD_NAME_LEN), target, intent(in), optional :: opt(:)
        character(len=FIELD_NAME_LEN), target, intent(in), optional :: core(:)
        type(wire_node_t) :: n

        n%kind = kind
        n%name = name
        n%prec = prec
        n%dt   = dt
        if (present(check)) n%check = check
        n%imports => imports
        n%exports => exports
        if (present(axes)) n%axes => axes
        if (present(opt))  n%opt  => opt
        if (present(core)) n%core => core
    end function wire_node_states

    integer function node_index(nodes, kind) result(k)
        type(wire_node_t), intent(in) :: nodes(:)
        character(len=*),  intent(in) :: kind
        integer :: i
        k = 0
        do i = 1, size(nodes)
            if (trim(nodes(i)%kind) == trim(kind)) then
                k = i
                return
            end if
        end do
    end function node_index

    ! =========================================================================
    ! Resolution at init (section 4.4)
    ! =========================================================================

    subroutine wire_resolve(w, nodes, cpl, table, stat)
        ! Bind the rows to the nodes of one domain (plus the globe nodes its rows
        ! reach). Each row is active iff both fields are declared, else inactive
        ! with a reason. Then, stopping with a named message on failure:
        !   1. every required import of a checked node has exactly one active
        !      row, an optional one at most one (gaps, duplicates); an import the
        !      node says the core fills (the global geo's external ice and
        !      background reference, filled by gia_gather / globe_ice) that no
        !      row of the table names is no gap, and is logged;
        !   2. fields_check_link per active row; rank 3: same zdim, identical axis;
        !   3. no masked field across two grids; MEAN -> INSTANT rows need equal
        !      cadences;
        !   4. the shelf's ocean set: exactly one of profiles and shelf-base
        !      values is complete and active;
        !   5. a diagnostic (restart = .false., not in the restart bundle; its
        !      backend puts it again on PASS_RESTART) is an INSTANT export that
        !      no active row reads (check_diagnostics).
        ! Then the methods are resolved and the workspaces allocated. coupler_prime
        ! of every remapped pair is the caller's (wire_prime).
        ! w: every row of the table, rules expanded, in table order.
        type(wire_t), allocatable, intent(out)   :: w(:)
        type(wire_node_t),         intent(in)    :: nodes(:)
        type(coupler_class),       intent(inout) :: cpl
        type(wire_t),              intent(in), optional :: table(:)
        integer,                   intent(out), optional :: stat

        type(wire_t), allocatable :: t(:)
        integer :: i

        if (present(stat)) stat = 0

        if (present(table)) then
            t = table
        else
            t = wiring_table()
        end if
        call expand_rules(w, t, nodes)

        do i = 1, size(w)
            if (bind_row(w(i), nodes, stat)) return
        end do
        if (check_gaps(w, nodes, stat)) return
        do i = 1, size(w)
            if (.not. w(i)%active) cycle
            if (check_row(w(i), nodes, stat)) return
        end do
        if (check_shelf(w, nodes, stat)) return
        if (check_diagnostics(w, nodes, stat)) return

        do i = 1, size(w)
            if (.not. w(i)%active) cycle
            call finish_row(w(i), nodes, cpl)
        end do
    end subroutine wire_resolve

    subroutine expand_rules(w, t, nodes)
        ! Copy the table, replacing each rule row by one row per matching name
        ! of the destination's declared imports (section 4.2).
        type(wire_t), allocatable, intent(out) :: w(:)
        type(wire_t),              intent(in)  :: t(:)
        type(wire_node_t),         intent(in)  :: nodes(:)

        type(wire_t) :: r
        integer :: i, k, id, j
        logical :: fed

        allocate(w(0))
        do i = 1, size(t)
            id = node_index(nodes, t(i)%dst_kind)
            if (trim(t(i)%src_field) == RULE_COEF .and. id > 0) then
                ! climate import coef_<name> reads forcing <name>
                do k = 1, nodes(id)%imports%n
                    associate (nm => nodes(id)%imports%f(k)%name)
                        if (nm(1:min(5, len(nm))) /= "coef_") cycle
                        r = t(i)
                        r%src_field = trim(nm(6:));  r%dst_field = nm
                        call append(w, r)
                    end associate
                end do
            else if (trim(t(i)%src_field) == RULE_DERIVED .and. id > 0) then
                ! forcing import <name> (not fed by another row) reads the climate export <name>
                do k = 1, nodes(id)%imports%n
                    associate (nm => nodes(id)%imports%f(k)%name)
                        fed = .false.
                        do j = 1, size(t)
                            if (trim(t(j)%dst_kind) == trim(t(i)%dst_kind) .and. trim(t(j)%dst_field) == trim(nm) &
                                .and. j /= i) fed = .true.
                        end do
                        if (fed) cycle
                        r = t(i)
                        r%src_field = nm;  r%dst_field = nm
                        call append(w, r)
                    end associate
                end do
            else
                call append(w, t(i))
            end if
        end do
    end subroutine expand_rules

    function bind_row(w, nodes, stat) result(bad)
        ! Activity, handles and the table's own consistency for one row.
        type(wire_t),      intent(inout) :: w
        type(wire_node_t), intent(in)    :: nodes(:)
        integer,           intent(inout), optional :: stat
        logical :: bad

        character(len=*), parameter :: proc = "kryos_wiring::wire_resolve"
        character(len=:), allocatable :: msg
        integer :: is, id, k

        bad = .false.
        w%active = .false.
        if (w%lag /= 0 .and. w%lag /= 1 .and. w%lag /= LAG_INIT) then
            bad = failed(proc, "row "//wire_label(w)//": lag must be 0, 1 or LAG_INIT", "lag = "//to_str(w%lag), stat)
            return
        end if

        is = node_index(nodes, w%src_kind)
        id = node_index(nodes, w%dst_kind)
        if (is > 0) w%src_name = nodes(is)%name
        if (id > 0) w%dst_name = nodes(id)%name
        if (is == 0) then
            w%reason = "domain has no "//trim(kind_word(w%src_kind))//" instance"
            return
        end if
        if (id == 0) then
            w%reason = "domain has no "//trim(kind_word(w%dst_kind))//" instance"
            return
        end if
        if (.not. (associated(nodes(is)%exports) .and. associated(nodes(id)%imports))) then
            bad = failed(proc, "row "//wire_label(w)//": node without field states", &
                         "nodes "//trim(w%src_name)//", "//trim(w%dst_name), stat)
            return
        end if
        if (.not. fields_has(nodes(is)%exports, w%src_field)) then
            w%reason = trim(w%src_name)//" does not export `"//trim(w%src_field)//"`"
            return
        end if
        if (.not. fields_has(nodes(id)%imports, w%dst_field)) then
            w%reason = trim(w%dst_name)//" has no import `"//trim(w%dst_field)//"`"
            return
        end if
        w%hs = fields_handle(nodes(is)%exports, w%src_field)
        w%hd = fields_handle(nodes(id)%imports, w%dst_field)
        w%src_grid = nodes(is)%exports%f(w%hs)%grid
        w%dst_grid = nodes(id)%imports%f(w%hd)%grid

        associate (fs => nodes(is)%exports%f(w%hs), fd => nodes(id)%imports%f(w%hd))
            if (w%gather .and. (fs%rank /= 0 .or. w%transform /= T_NONE)) then
                bad = failed(proc, "row "//wire_label(w)//": a gather row is 0-D without transform", "", stat)
                return
            end if
            if (fs%rank == 0 .and. w%transform /= T_NONE) then
                bad = failed(proc, "row "//wire_label(w)//": a 0-D row has no transform", "", stat)
                return
            end if

            select case (w%transform)
            case (T_NONE)
            case (T_ADD_REF, T_SUB_REF)
                select case (trim(w%aux))
                case ("z_bed", "H_ice", "z_srf", "z_bed_sd")
                case default
                    bad = failed(proc, "row "//wire_label(w)//": unknown reference field `"//trim(w%aux)//"`", &
                                 "grid_ref_t has z_bed, H_ice, z_srf, z_bed_sd", stat)
                    return
                end select
            case (T_ADD_SCALAR)
                if (.not. fields_has(nodes(is)%exports, w%aux)) then
                    bad = failed(proc, "row "//wire_label(w)//": "//trim(w%src_name)//" does not export `"// &
                                 trim(w%aux)//"`", "the 0-D field added to the remapped one (aux)", stat)
                    return
                end if
                w%ha = fields_handle(nodes(is)%exports, w%aux)
                if (nodes(is)%exports%f(w%ha)%rank /= 0) then
                    bad = failed(proc, "row "//wire_label(w)//": aux `"//trim(w%aux)//"` is not 0-D", "", stat)
                    return
                end if
            case (T_GIA)
                if (fs%rank /= 2 .or. fd%rank /= 2) then
                    bad = failed(proc, "row "//wire_label(w)//": a T_GIA row is rank 2 on both sides", "", stat)
                    return
                end if
                msg = gia_sources(w, nodes(is), fs)
                if (len(msg) > 0) then
                    if (msg(1:1) == "!") then           ! declared, but not on the principal source's grid
                        bad = failed(proc, "row "//wire_label(w)//": "//msg(2:), "", stat)
                    else
                        w%reason = msg                  ! the producer does not export one of them: inactive
                    end if
                    return
                end if
            case (T_LAPSE_T, T_LAPSE_P)
                if (fs%rank /= 3) then
                    bad = failed(proc, "row "//wire_label(w)//": a lapse row needs a rank-3 (monthly) field", "", stat)
                    return
                end if
                do k = 1, 3
                    associate (nm => lapse_name(k))
                        if (.not. fields_has(nodes(is)%exports, nm)) then
                            bad = failed(proc, "row "//wire_label(w)//": "//trim(w%src_name)//" does not export `"// &
                                         nm//"`", "the lapse transform needs z_srf_clim, beta_p and lapse", stat)
                            return
                        end if
                        w%hl(k) = fields_handle(nodes(is)%exports, nm)
                    end associate
                end do
                if (.not. fields_has(nodes(id)%imports, "z_srf")) then
                    bad = failed(proc, "row "//wire_label(w)//": "//trim(w%dst_name)//" has no import `z_srf`", &
                                 "the lapse transform re-elevates to the consumer's surface", stat)
                    return
                end if
                w%hl(4) = fields_handle(nodes(id)%imports, "z_srf")
            case default
                bad = failed(proc, "row "//wire_label(w)//": unknown transform "//to_str(w%transform), "", stat)
                return
            end select

            if (fs%rank >= 2 .and. trim(w%src_grid) == trim(w%dst_grid)) then
                if (fs%nx /= fd%nx .or. fs%ny /= fd%ny) then
                    bad = failed(proc, "row "//wire_label(w)//": one grid, two shapes", &
                                 "grid "//trim(w%src_grid)//": "//to_str(fs%nx)//"x"//to_str(fs%ny)//" -> "// &
                                 to_str(fd%nx)//"x"//to_str(fd%ny), stat)
                    return
                end if
            end if
        end associate

        w%active = .true.
        w%reason = ""
    end function bind_row

    function gia_sources(w, node, fs) result(msg)
        ! The other source fields of a T_GIA row (w%srcs): their handles into
        ! w%hx. msg = "" when the producer exports them all, on the grid of the
        ! principal source fs; else the reason the row is inactive (the first one
        ! it does not export), or, led by "!", a declaration error (one that is
        ! not rank 2 on that grid).
        type(wire_t),      intent(inout) :: w
        type(wire_node_t), intent(in)    :: node
        type(field_t),     intent(in)    :: fs
        character(len=:), allocatable :: msg

        integer :: k, h

        msg = ""
        if (allocated(w%hx)) deallocate(w%hx)
        allocate(w%hx(0))
        if (.not. allocated(w%srcs)) return
        do k = 1, size(w%srcs)
            if (.not. fields_has(node%exports, trim(w%srcs(k)))) then
                msg = trim(w%src_name)//" does not export `"//trim(w%srcs(k))//"`"
                return
            end if
            h = fields_handle(node%exports, trim(w%srcs(k)))
            associate (f => node%exports%f(h))
                if (f%rank /= 2 .or. trim(f%grid) /= trim(fs%grid) .or. f%nx /= fs%nx .or. f%ny /= fs%ny) then
                    msg = "!source `"//trim(w%srcs(k))//"` is not rank 2 on the grid of `"//trim(fs%name)//"` ("// &
                          trim(fs%grid)//")"
                    return
                end if
            end associate
            w%hx = [w%hx, h]
        end do
    end function gia_sources

    function check_gaps(w, nodes, stat) result(bad)
        ! Step 1: producers of every import.
        type(wire_t),      intent(in) :: w(:)
        type(wire_node_t), intent(in) :: nodes(:)
        integer,           intent(inout), optional :: stat
        logical :: bad

        character(len=*), parameter :: proc = "kryos_wiring::wire_resolve"
        character(len=:), allocatable :: who, rows, why
        integer :: n, h, i, nact, nany
        logical :: optional_, core_

        bad = .false.
        do n = 1, size(nodes)
            if (.not. associated(nodes(n)%imports)) cycle
            do h = 1, nodes(n)%imports%n
                associate (nm => nodes(n)%imports%f(h)%name)
                    optional_ = .false.
                    if (associated(nodes(n)%opt)) optional_ = any(trim(nm) == nodes(n)%opt)
                    core_ = .false.
                    if (associated(nodes(n)%core)) core_ = any(trim(nm) == nodes(n)%core)
                    nact = 0;  nany = 0
                    rows = "";  why = ""
                    do i = 1, size(w)
                        if (trim(w(i)%dst_kind) /= trim(nodes(n)%kind) .or. trim(w(i)%dst_field) /= trim(nm)) cycle
                        nany = nany + 1
                        if (w(i)%active) then
                            nact = nact + 1
                            rows = rows//" `"//wire_label(w(i))//"`"
                        else if (len(why) == 0) then
                            why = "row "//wire_label(w(i))//" inactive ("//trim(w(i)%reason)//")"
                        end if
                    end do
                    who = trim(nodes(n)%kind)//" `"//trim(nodes(n)%name)//"` import `"//trim(nm)//"`"
                    if (core_ .and. nany == 0) then     ! no row feeds it: the core does
                        call log_msg(LOG_INFO, proc, who//": no wiring row, filled by gia_gather / globe_ice")
                        cycle
                    end if
                    if (nact == 0 .and. .not. optional_ .and. nodes(n)%check) then
                        if (nany == 0) why = "no row of the wiring table feeds it"
                        if (nm(1:min(5, len(nm))) == "coef_") why = trim(why)//"; is the coefficient in `[forcing] coefs`?"
                        bad = failed(proc, who//" has no producer: "//why, "", stat)
                        return
                    else if (nact > 1) then
                        bad = failed(proc, who//" has "//to_str(nact)//" active rows:"//rows, "", stat)
                        return
                    end if
                end associate
            end do
        end do
    end function check_gaps

    function check_row(w, nodes, stat) result(bad)
        ! Steps 2 and 3 for one active row.
        type(wire_t),      intent(in) :: w
        type(wire_node_t), intent(in) :: nodes(:)
        integer,           intent(inout), optional :: stat
        logical :: bad

        character(len=*), parameter :: proc = "kryos_wiring::wire_resolve"
        character(len=:), allocatable :: msg
        integer :: is, id

        bad = .false.
        if (w%transform == T_GIA) return        ! the fields are related by the core's gather, not by value:
                                                ! units, nature and masks are not compared
        is = node_index(nodes, w%src_kind)
        id = node_index(nodes, w%dst_kind)
        call fields_check_link(nodes(is)%exports, w%hs, nodes(id)%imports, w%hd)   ! stops by itself
        associate (fs => nodes(is)%exports%f(w%hs), fd => nodes(id)%imports%f(w%hd))
            if (fs%rank == 3) then
                msg = axis_problem(fs, nodes(is), fd, nodes(id))
                if (len(msg) > 0) then
                    bad = failed(proc, "row "//wire_label(w)//": "//msg, "", stat)
                    return
                end if
            end if
            if (fs%masked .neqv. fd%masked) then
                bad = failed(proc, "row "//wire_label(w)//": mask on one side only", "", stat)
                return
            end if
            if (fs%masked .and. trim(fs%grid) /= trim(fd%grid)) then
                bad = failed(proc, "row "//wire_label(w)//": masked field across two grids", &
                             trim(fs%grid)//" -> "//trim(fd%grid)//" (not supported)", stat)
                return
            end if
            if (fs%tsem == FIELD_MEAN .and. fd%tsem == FIELD_INSTANT) then
                if (abs(nodes(is)%dt - nodes(id)%dt) > 0.0_dp) then
                    bad = failed(proc, "row "//wire_label(w)//": MEAN export read as INSTANT needs equal cadences", &
                                 "producer "//trim(nodes(is)%name)//": dt = "//to_str(nodes(is)%dt)//" yr; consumer "// &
                                 trim(nodes(id)%name)//": dt = "//to_str(nodes(id)%dt)//" yr", stat)
                    return
                end if
            end if
        end associate
    end function check_row

    function check_diagnostics(w, nodes, stat) result(bad)
        ! Step 5: the state of a restart is every import, every export an
        ! active row reads (its principal source, aux, lapse fields or T_GIA
        ! sources) and every export with a window (MEAN, ACCUM); a field
        ! declared a diagnostic (restart = .false.) must be none of them.
        type(wire_t),      intent(in) :: w(:)
        type(wire_node_t), intent(in) :: nodes(:)
        integer,           intent(inout), optional :: stat
        logical :: bad

        character(len=*), parameter :: proc = "kryos_wiring::wire_resolve"
        character(len=*), parameter :: why = "a diagnostic (restart = .false.) is not in the restart bundle"
        character(len=:), allocatable :: who
        integer :: n, h, i

        bad = .false.
        do n = 1, size(nodes)
            if (associated(nodes(n)%imports)) then
                do h = 1, nodes(n)%imports%n
                    if (nodes(n)%imports%f(h)%restart) cycle
                    bad = failed(proc, trim(nodes(n)%kind)//" `"//trim(nodes(n)%name)//"` import `"// &
                                 trim(nodes(n)%imports%f(h)%name)//"` is declared a diagnostic", why// &
                                 "; every import is state", stat)
                    return
                end do
            end if
            if (.not. associated(nodes(n)%exports)) cycle
            do h = 1, nodes(n)%exports%n
                associate (f => nodes(n)%exports%f(h))
                    if (f%restart) cycle
                    who = trim(nodes(n)%kind)//" `"//trim(nodes(n)%name)//"` export `"//trim(f%name)//"`"
                    if (f%tsem /= FIELD_INSTANT) then
                        bad = failed(proc, who//" is declared a diagnostic but is not INSTANT", why// &
                                     "; its open window is state", stat)
                        return
                    end if
                    do i = 1, size(w)
                        if (.not. w(i)%active) cycle
                        if (node_index(nodes, w(i)%src_kind) /= n) cycle
                        if (.not. reads(w(i), h)) cycle
                        bad = failed(proc, who//" is declared a diagnostic but row "//wire_label(w(i))// &
                                     " reads it", why//"; a row's source is state", stat)
                        return
                    end do
                end associate
            end do
        end do

    contains

        logical function reads(r, hh)
            ! Does the active row r read the producer's export hh?
            type(wire_t), intent(in) :: r
            integer,      intent(in) :: hh
            reads = r%hs == hh
            if (r%transform == T_ADD_SCALAR) reads = reads .or. r%ha == hh
            if (r%transform == T_LAPSE_T .or. r%transform == T_LAPSE_P) reads = reads .or. any(r%hl(1:3) == hh)
            if (allocated(r%hx)) reads = reads .or. any(r%hx == hh)
        end function reads

    end function check_diagnostics

    function axis_problem(fs, ns, fd, nd) result(msg)
        ! Rank 3: same zdim and identical axis on both sides.
        type(field_t),     intent(in) :: fs, fd
        type(wire_node_t), intent(in) :: ns, nd
        character(len=:), allocatable :: msg
        integer :: ks, kd

        msg = ""
        if (trim(fs%zdim) /= trim(fd%zdim)) then
            msg = "axis differs: `"//trim(fs%zdim)//"` -> `"//trim(fd%zdim)//"`"
            return
        end if
        ks = 0;  kd = 0
        if (associated(ns%axes)) ks = axis_find(ns%axes, fs%zdim)
        if (associated(nd%axes)) kd = axis_find(nd%axes, fd%zdim)
        if (ks == 0) msg = "axis `"//trim(fs%zdim)//"` is not declared by "//trim(ns%name)
        if (kd == 0) msg = "axis `"//trim(fd%zdim)//"` is not declared by "//trim(nd%name)
        if (ks == 0 .or. kd == 0) return
        if (.not. axes_equal(ns%axes(ks), nd%axes(kd))) then
            msg = "axis `"//trim(fs%zdim)//"` differs between "//trim(ns%name)//" and "//trim(nd%name)
        end if
    end function axis_problem

    function check_shelf(w, nodes, stat) result(bad)
        ! Step 4, contracts-design section 2.5: exactly one ocean set.
        type(wire_t),      intent(in) :: w(:)
        type(wire_node_t), intent(in) :: nodes(:)
        integer,           intent(inout), optional :: stat
        logical :: bad

        character(len=*), parameter :: proc = "kryos_wiring::wire_resolve"
        character(len=FIELD_NAME_LEN), parameter :: set_a(3) = [character(len=FIELD_NAME_LEN) :: &
            "to_ann", "so_ann", "dto_ann"]
        character(len=FIELD_NAME_LEN), parameter :: set_b(4) = [character(len=FIELD_NAME_LEN) :: &
            "to_shlf", "so_shlf", "dto_shlf", "dso_shlf"]
        character(len=:), allocatable :: clim
        integer :: n, ic, nset

        bad = .false.
        n = node_index(nodes, "shelf")
        if (n == 0) return
        nset = count([complete(set_a), complete(set_b)])
        if (nset == 1) return

        ic = node_index(nodes, "climate")
        if (ic > 0) then
            clim = "climate `"//trim(nodes(ic)%name)//"`"
        else
            clim = "the domain has no climate instance;"
        end if
        if (nset == 0) then
            bad = failed(proc, "shelf `"//trim(nodes(n)%name)//"` needs ocean forcing: "//clim// &
                         " declares neither profiles (has_ocn) nor shelf-base values (has_ocn_shelf)", "", stat)
        else
            bad = failed(proc, "shelf `"//trim(nodes(n)%name)//"`: both ocean sets (profiles and shelf-base values)"// &
                         " are wired; exactly one is allowed", "", stat)
        end if

    contains

        logical function complete(set)
            character(len=*), intent(in) :: set(:)
            integer :: i, j
            logical :: found
            complete = .true.
            do i = 1, size(set)
                found = .false.
                do j = 1, size(w)
                    if (trim(w(j)%dst_kind) == "shelf" .and. trim(w(j)%dst_field) == trim(set(i)) &
                        .and. w(j)%active) found = .true.
                end do
                if (.not. found) complete = .false.
            end do
        end function complete

    end function check_shelf

    subroutine finish_row(w, nodes, cpl)
        ! Precision, method and workspaces of an active row.
        type(wire_t),        intent(inout) :: w
        type(wire_node_t),   intent(in)    :: nodes(:)
        type(coupler_class), intent(inout) :: cpl

        integer :: is, id

        is = node_index(nodes, w%src_kind)
        id = node_index(nodes, w%dst_kind)
        w%prec_src = nodes(is)%prec
        w%prec_dst = nodes(id)%prec
        w%mres = ""
        if (w%transform == T_GIA) return        ! no remap, no workspace: the core's gather (gia_gather)
        associate (fs => nodes(is)%exports%f(w%hs), fd => nodes(id)%imports%f(w%hd))
            if (fs%rank == 0) return
            if (trim(w%src_grid) /= trim(w%dst_grid)) then
                select case (trim(w%method))
                case ("")
                    w%mres = merge("con  ", "bilin", fs%cons == FIELD_EXTENSIVE)
                case ("auto")
                    w%mres = coupler_method(cpl, trim(w%src_grid), trim(w%dst_grid), "auto")
                case default
                    w%mres = w%method
                end select
            end if
            ! Workspaces: the remapped or transformed value (none on one grid
            ! without a destination transform: the value is put as it is), the
            ! anomaly before the remap, the remapped lapse fields.
            if (len_trim(w%mres) > 0 .or. dst_transform(w)) allocate(w%work(fd%nx, fd%ny, fd%nz))
            if (w%transform == T_SUB_REF) allocate(w%pre(fs%nx, fs%ny, fs%nz))
            if ((w%transform == T_LAPSE_T .or. w%transform == T_LAPSE_P) .and. len_trim(w%mres) > 0) then
                allocate(w%zc(fd%nx, fd%ny, 1), w%bp(fd%nx, fd%ny, 1), w%lp(fd%nx, fd%ny, fd%nz))
            end if
        end associate
    end subroutine finish_row

    logical function dst_transform(w)
        ! The row transforms the value on the destination side (after the remap).
        type(wire_t), intent(in) :: w
        dst_transform = any(w%transform == [T_ADD_REF, T_ADD_SCALAR, T_LAPSE_T, T_LAPSE_P])
    end function dst_transform

    subroutine wire_prime(w, cpl)
        ! coupler_prime for every remapped pair of the active rows (the maps are
        ! built, or loaded, at init and not in the step).
        type(wire_t),        intent(in)    :: w(:)
        type(coupler_class), intent(inout) :: cpl
        integer :: i

        do i = 1, size(w)
            if (.not. w(i)%active .or. len_trim(w(i)%mres) == 0) cycle
            call coupler_prime(cpl, trim(w(i)%src_grid), trim(w(i)%dst_grid), method=trim(w(i)%mres))
        end do
    end subroutine wire_prime

    logical function failed(proc, msg, detail, stat)
        ! fesm_error: stops unless stat is present, then returns .true.
        character(len=*), intent(in)  :: proc, msg, detail
        integer,          intent(inout), optional :: stat
        if (len_trim(detail) > 0) then
            call fesm_error(proc, msg, detail, stat)
        else
            call fesm_error(proc, msg, stat=stat)
        end if
        failed = .true.
    end function failed

    function lapse_name(k) result(s)
        integer, intent(in) :: k
        character(len=FIELD_NAME_LEN) :: s
        select case (k)
        case (1);  s = "z_srf_clim"
        case (2);  s = "beta_p"
        case default;  s = "lapse"
        end select
    end function lapse_name

    function kind_word(kind) result(s)
        character(len=*), intent(in) :: kind
        character(len=:), allocatable :: s
        select case (trim(kind))
        case ("geo_g");  s = "global geo"
        case default;    s = trim(kind)
        end select
    end function kind_word

    function wire_label(w) result(s)
        type(wire_t), intent(in) :: w
        character(len=:), allocatable :: s
        s = trim(w%src_kind)//"."//trim(w%src_field)//" -> "//trim(w%dst_kind)//"."//trim(w%dst_field)
    end function wire_label

    logical function wire_runs(w, src_kind, dst_kind)
        ! Is `w` an active exchange row src_kind -> dst_kind (not a gather row,
        ! not a T_GIA row: the core runs those itself)?
        type(wire_t),     intent(in) :: w
        character(len=*), intent(in) :: src_kind, dst_kind
        wire_runs = w%active .and. .not. w%gather .and. w%transform /= T_GIA .and. trim(w%src_kind) == trim(src_kind) &
                    .and. trim(w%dst_kind) == trim(dst_kind)
    end function wire_runs

    ! =========================================================================
    ! Describe (logged once, after fields_describe of every state)
    ! =========================================================================

    subroutine wiring_describe(w, level)
        ! One line per row: src.field [grid] -> dst.field [grid]  method (resolved)  lag  transform
        ! active|inactive: reason.
        type(wire_t),  intent(in) :: w(:)
        integer,       intent(in), optional :: level
        character(len=*), parameter :: who = "kryos_wiring"
        character(len=:), allocatable :: line, meth, lag, tr
        integer :: i, lvl

        lvl = LOG_INFO
        if (present(level)) lvl = level
        call log_msg(lvl, who, "wiring: "//to_str(size(w))//" rows, "//to_str(count(w%active))//" active")
        do i = 1, size(w)
            associate (r => w(i))
                if (r%transform == T_GIA .or. (r%active .and. len_trim(r%src_grid) == 0)) then
                    meth = "-"
                else if (len_trim(r%mres) > 0) then
                    meth = method_word(r)//" ("//trim(r%mres)//")"
                else if (r%active) then
                    meth = "(identity)"
                else
                    meth = method_word(r)
                end if
                select case (r%lag)
                case (LAG_INIT);  lag = "init"
                case default;     lag = to_str(r%lag)
                end select
                tr = transform_name(r)
                line = trim(r%src_kind)//"."//trim(r%src_field)//" ["//trim(r%src_grid)//"] -> "// &
                       trim(r%dst_kind)//"."//trim(r%dst_field)//" ["//trim(r%dst_grid)//"]  "//meth// &
                       "  lag "//lag//"  "//tr
                if (r%active) then
                    line = line//"  active"
                else
                    line = line//"  inactive: "//trim(r%reason)
                end if
                call log_msg(lvl, who, line)
            end associate
        end do
    end subroutine wiring_describe

    function method_word(w) result(s)
        type(wire_t), intent(in) :: w
        character(len=:), allocatable :: s
        if (len_trim(w%method) == 0) then
            s = "cons"
        else
            s = trim(w%method)
        end if
    end function method_word

    function transform_name(w) result(s)
        type(wire_t), intent(in) :: w
        character(len=:), allocatable :: s
        integer :: k
        select case (w%transform)
        case (T_NONE);       s = merge("gather", "none  ", w%gather)
        case (T_ADD_REF);    s = "add_ref("//trim(w%aux)//")"
        case (T_SUB_REF);    s = "sub_ref("//trim(w%aux)//")"
        case (T_ADD_SCALAR); s = "add_scalar("//trim(w%aux)//")"
        case (T_LAPSE_T);    s = "lapse_t"
        case (T_LAPSE_P);    s = "lapse_p"
        case (T_GIA)
            s = "gia("
            if (allocated(w%srcs)) then
                do k = 1, size(w%srcs)
                    if (k > 1) s = s//" "
                    s = s//trim(w%srcs(k))
                end do
            end if
            s = s//")"
        case default;        s = "?"
        end select
        s = trim(s)
    end function transform_name

    ! =========================================================================
    ! Lag check
    ! =========================================================================

    function wire_lag_status(w, src, dst, t, t_last, due, init) result(msg)
        ! The checks of an exchange as a message: "" when the row may run. t is
        ! the clock time, t_last the producer's cadence_t_last (hub: hub_t).
        !   * the source was published, and its t1 is the producer's last update
        !     (the producer put what it produced, and marked it);
        !   * t1 <= t; lag 0 with the producer due at t: t1 = t (it has run);
        !     lag 1: t1 < t (the previous step's value);
        !   * cold start (init): every row counts as lag 0;
        !   * LAG_INIT rows: published, t1 <= t (static producers are not marked);
        !   * the aux 0-D field and the lapse fields were put in the same update,
        !     and the consumer's z_srf of a lapse row was exchanged;
        !   * a T_GIA row (run by the core's gather): the same for every source
        !     field it reads (the principal one and srcs).
        type(wire_t),        intent(in) :: w
        type(field_state_t), intent(in) :: src, dst
        real(dp),            intent(in) :: t, t_last
        logical,             intent(in), optional :: due, init
        character(len=:), allocatable :: msg

        logical :: d, i0
        integer :: lag, k

        msg = ""
        d = .false.;  i0 = .false.
        if (present(due))  d  = due
        if (present(init)) i0 = init
        lag = w%lag
        if (i0 .and. lag == 1) lag = 0

        associate (fs => src%f(w%hs))
            msg = source_problem(fs)
            if (len(msg) > 0) return
            if (w%transform == T_GIA .and. allocated(w%hx)) then
                do k = 1, size(w%hx)
                    msg = source_problem(src%f(w%hx(k)))
                    if (len(msg) > 0) then
                        if (index(msg, "never published") == 0) msg = "source `"//trim(src%f(w%hx(k))%name)//"`: "//msg
                        return
                    end if
                end do
            end if

            if (w%transform == T_ADD_SCALAR) then
                msg = aux_problem(src%f(w%ha), fs)
            else if (w%transform == T_LAPSE_T .or. w%transform == T_LAPSE_P) then
                do k = 1, 3
                    msg = aux_problem(src%f(w%hl(k)), fs)
                    if (len(msg) > 0) return
                end do
                if (.not. dst%f(w%hl(4))%published) then
                    msg = "consumer import `z_srf` not yet exchanged (hub -> "//trim(w%dst_kind)// &
                          " must run before the lapse row)"
                end if
            end if
        end associate

    contains

        function source_problem(f) result(m)
            ! Published, and at a time that fits the row's lag.
            type(field_t), intent(in) :: f
            character(len=:), allocatable :: m

            m = ""
            if (.not. f%published) then
                m = "source `"//trim(src%owner)//"."//trim(f%name)//"` was never published"
                return
            end if
            if (w%lag == LAG_INIT) then
                if (f%t1 > t) m = "published at t1 = "//to_str(f%t1)//", after the exchange time t = "//to_str(t)
            else if (.not. same_time(f%t1, t_last)) then
                m = "published at t1 = "//to_str(f%t1)//", but the producer's last update is at "//to_str(t_last)// &
                    " (a put without cadence_mark, or a mark without a put)"
            else if (f%t1 > t) then
                m = "published at t1 = "//to_str(f%t1)//", ahead of the exchange time t = "//to_str(t)
            else if (lag == 0 .and. d .and. f%t1 < t) then
                m = "published at t1 = "//to_str(f%t1)//", expected "//to_str(t)// &
                    " (lag 0, producer due): the row runs before the producer in the step"
            else if (lag == 1 .and. f%t1 >= t) then
                m = "published at t1 = "//to_str(f%t1)//", expected before "//to_str(t)// &
                    " (lag 1): the row runs after the producer in the step"
            end if
        end function source_problem

    end function wire_lag_status

    function aux_problem(fa, fs) result(msg)
        type(field_t), intent(in) :: fa, fs
        character(len=:), allocatable :: msg
        msg = ""
        if (.not. fa%published) then
            msg = "`"//trim(fa%name)//"` was never published"
        else if (.not. (same_time(fa%t1, fs%t1))) then
            msg = "`"//trim(fa%name)//"` is at t1 = "//to_str(fa%t1)//", `"//trim(fs%name)//"` at "//to_str(fs%t1)
        end if
    end function aux_problem

    logical function same_time(a, b)
        real(dp), intent(in) :: a, b
        same_time = .not. (a < b .or. a > b)
    end function same_time

    ! =========================================================================
    ! Exchange (section 4.3)
    ! =========================================================================

    subroutine wire_exchange(w, src, dst, cpl, t, t_last, ref_src, ref_dst, due, init, stat)
        ! One row: the published source value into the destination, remapped and
        ! transformed. t = clock time, t_last = the producer's cadence_t_last
        ! (hub: hub_t); ref_src, ref_dst: the references on the two grids (the
        ! domain's refs; needed by T_SUB_REF and T_ADD_REF); due: the producer is
        ! due at t; init: cold start (every row counts as lag 0). A LAG_INIT
        ! row that has run is skipped. Gather rows go through wire_gather.
        type(wire_t),        intent(inout) :: w
        type(field_state_t), intent(in)    :: src
        type(field_state_t), intent(inout) :: dst
        type(coupler_class), intent(inout) :: cpl
        real(dp),            intent(in)    :: t, t_last
        type(grid_ref_t),    intent(in), optional :: ref_src, ref_dst
        logical,             intent(in), optional :: due, init
        integer,             intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_wiring::wire_exchange"
        character(len=:), allocatable :: msg

        if (present(stat)) stat = 0
        if (.not. w%active) then
            call fesm_error(proc, "row "//wire_label(w)//" is inactive", trim(w%reason), stat)
            return
        end if
        if (w%gather) then
            call fesm_error(proc, "row "//wire_label(w)//" is a gather row", "use wire_gather", stat)
            return
        end if
        if (w%transform == T_GIA) then
            call fesm_error(proc, "row "//wire_label(w)//" is a T_GIA row", "the core's gia_gather runs it", stat)
            return
        end if
        if (w%lag == LAG_INIT .and. w%done) return

        msg = wire_lag_status(w, src, dst, t, t_last, due, init)
        if (len(msg) > 0) then
            call fesm_error(proc, "row "//wire_label(w)//": "//msg, &
                            "producer "//trim(w%src_name)//", consumer "//trim(w%dst_name)//", t = "//to_str(t), stat)
            return
        end if

        associate (fs => src%f(w%hs), fd => dst%f(w%hd))

            if (fs%rank == 0) then
                call fields_put(dst, w%hd, fs%v(1,1,1), t0=fs%t0, t1=fs%t1)
                w%done = .true.
                return
            end if

            ! references the transform needs
            msg = ""
            if (w%transform == T_SUB_REF) then
                if (.not. present(ref_src)) then
                    msg = "needs ref_src (the reference on grid "//trim(fs%grid)//")"
                else
                    msg = ref_problem(ref_src, fs, w%aux)
                end if
            else if (w%transform == T_ADD_REF) then
                if (.not. present(ref_dst)) then
                    msg = "needs ref_dst (the reference on grid "//trim(fd%grid)//")"
                else
                    msg = ref_problem(ref_dst, fd, w%aux)
                end if
            end if
            if (len(msg) > 0) then
                call fesm_error(proc, "row "//wire_label(w)//": "//msg, "", stat)
                return
            end if

            if (w%transform == T_SUB_REF) &                     ! anomaly in the source's precision
                call sub_ref(w%pre, fs%v, ref_src, w%aux, w%prec_src)

            ! One grid and no transform on the destination side: the value is
            ! published as it is, without a copy into the workspace.
            if (len_trim(w%mres) == 0 .and. .not. dst_transform(w)) then
                if (w%transform == T_SUB_REF) then
                    call put(w%pre)
                else
                    call put(fs%v)
                end if
                w%done = .true.
                return
            end if

            if (w%transform == T_SUB_REF) then
                call remap_row(w, cpl, w%pre, w%work)
            else
                call remap_row(w, cpl, fs%v, w%work)
            end if

            select case (w%transform)                           ! in the destination's precision (section 9)
            case (T_ADD_REF);    call add_ref(w%work, ref_dst, w%aux, w%prec_dst)
            case (T_ADD_SCALAR); call add_scalar(w%work, src%f(w%ha)%v(1,1,1), w%prec_dst)
            case (T_LAPSE_T, T_LAPSE_P)
                if (len_trim(w%mres) == 0) then                 ! one grid: the producer's fields
                    call lapse(w%work, src%f(w%hl(1))%v, src%f(w%hl(2))%v, src%f(w%hl(3))%v, &
                               dst%f(w%hl(4))%v, w%transform == T_LAPSE_P, w%prec_dst)
                else
                    call remap_row(w, cpl, src%f(w%hl(1))%v, w%zc)
                    call remap_row(w, cpl, src%f(w%hl(2))%v, w%bp)
                    call remap_row(w, cpl, src%f(w%hl(3))%v, w%lp)
                    call lapse(w%work, w%zc, w%bp, w%lp, dst%f(w%hl(4))%v, w%transform == T_LAPSE_P, w%prec_dst)
                end if
            end select

            call put(w%work)
            w%done = .true.
        end associate

    contains

        subroutine put(v)
            ! The row's value into the destination, with the source's interval
            ! and, for a masked field (one grid, checked at resolution), its mask.
            real(dp), intent(in) :: v(:,:,:)
            associate (fs => src%f(w%hs))
                if (dst%f(w%hd)%masked) then
                    call fields_put(dst, w%hd, v, t0=fs%t0, t1=fs%t1, valid=fs%valid)
                else
                    call fields_put(dst, w%hd, v, t0=fs%t0, t1=fs%t1)
                end if
            end associate
        end subroutine put

    end subroutine wire_exchange

    subroutine remap_row(w, cpl, a, b)
        ! b = a on the destination grid: a copy on one grid, else the row's map.
        type(wire_t),        intent(in)    :: w
        type(coupler_class), intent(inout) :: cpl
        real(dp),            intent(in)    :: a(:,:,:)
        real(dp), allocatable, intent(inout) :: b(:,:,:)

        if (len_trim(w%mres) == 0) then
            b = a
        else
            call remap(cpl, a, trim(w%src_grid), b, trim(w%dst_grid), method=trim(w%mres))
        end if
    end subroutine remap_row

    function ref_problem(ref, f, name) result(msg)
        ! The reference is on the field's grid and holds the named array.
        type(grid_ref_t), intent(in) :: ref
        type(field_t),    intent(in) :: f
        character(len=*), intent(in) :: name
        character(len=:), allocatable :: msg

        msg = ""
        if (trim(ref%grid) /= trim(f%grid)) then
            msg = "reference is on grid `"//trim(ref%grid)//"`, the field on `"//trim(f%grid)//"`"
            return
        end if
        select case (trim(name))
        case ("z_bed");    if (.not. allocated(ref%z_bed))    msg = "reference z_bed not allocated"
        case ("H_ice");    if (.not. allocated(ref%H_ice))    msg = "reference H_ice not allocated"
        case ("z_srf");    if (.not. allocated(ref%z_srf))    msg = "reference z_srf not allocated"
        case ("z_bed_sd"); if (.not. allocated(ref%z_bed_sd)) msg = "reference z_bed_sd not allocated"
        end select
        if (len(msg) > 0) return
        select case (trim(name))
        case ("z_bed");    if (any(shape(ref%z_bed)    /= [f%nx, f%ny])) msg = "reference z_bed has the wrong shape"
        case ("H_ice");    if (any(shape(ref%H_ice)    /= [f%nx, f%ny])) msg = "reference H_ice has the wrong shape"
        case ("z_srf");    if (any(shape(ref%z_srf)    /= [f%nx, f%ny])) msg = "reference z_srf has the wrong shape"
        case ("z_bed_sd"); if (any(shape(ref%z_bed_sd) /= [f%nx, f%ny])) msg = "reference z_bed_sd has the wrong shape"
        end select
    end function ref_problem

    ! =========================================================================
    ! Gather (0-D, N domains -> one global field)
    ! =========================================================================

    subroutine wire_gather(ws, srcs, dst, t, stat)
        ! Sum, in the order of srcs (the core passes domain-name order), the 0-D
        ! values published at exactly t, stale ones count 0, and put the sum at
        ! t (the caller chooses t: the previous step for a lag 1 row). ws(i) is
        ! the resolved gather row of srcs(i), all into one field; the sum is
        ! formed in the destination's precision.
        type(wire_t),        intent(in)    :: ws(:)
        type(wire_node_t),   intent(in)    :: srcs(:)
        type(field_state_t), intent(inout) :: dst
        real(dp),            intent(in)    :: t
        integer,             intent(out), optional :: stat

        character(len=*), parameter :: proc = "kryos_wiring::wire_gather"
        real(dp) :: s_dp
        real(sp) :: s_sp
        integer :: i, hd, prec

        if (present(stat)) stat = 0
        if (size(ws) /= size(srcs)) then
            call fesm_error(proc, "one row per source", "rows = "//to_str(size(ws))//", sources = "//to_str(size(srcs)), stat)
            return
        end if
        hd = 0;  prec = dp
        do i = 1, size(ws)
            if (.not. (ws(i)%active .and. ws(i)%gather)) then
                call fesm_error(proc, "row "//wire_label(ws(i))//" is not an active gather row", trim(ws(i)%reason), stat)
                return
            end if
            if (hd /= 0 .and. ws(i)%hd /= hd) then
                call fesm_error(proc, "gather rows feed different fields", wire_label(ws(i)), stat)
                return
            end if
            hd = ws(i)%hd;  prec = ws(i)%prec_dst
        end do
        if (hd == 0) then
            call fesm_error(proc, "no source", "a gather needs at least one row to name the field", stat)
            return
        end if

        s_dp = 0.0_dp;  s_sp = 0.0_sp
        do i = 1, size(ws)
            associate (f => srcs(i)%exports%f(ws(i)%hs))
                if (.not. (f%published .and. same_time(f%t1, t))) cycle
                if (prec == sp) then
                    s_sp = s_sp + real(f%v(1,1,1), sp)
                else
                    s_dp = s_dp + f%v(1,1,1)
                end if
            end associate
        end do
        if (prec == sp) then
            call fields_put(dst, hd, real(s_sp, dp), t1=t)
        else
            call fields_put(dst, hd, s_dp, t1=t)
        end if
    end subroutine wire_gather

    ! =========================================================================
    ! Transform kernels, in the stated precision: an sp kernel narrows its
    ! operands (sp values widened, or remaps of them: narrowed, the coupler's sp
    ! remap), does the arithmetic in sp in yelmox's order and widens the result
    ! exactly.
    ! =========================================================================

    subroutine sub_ref(pre, v, ref, name, prec)
        ! pre = v - ref%<name>, formed in prec (the source's).
        real(dp),         intent(inout) :: pre(:,:,:)
        real(dp),         intent(in)    :: v(:,:,:)
        type(grid_ref_t), intent(in)    :: ref
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: prec

        select case (trim(name))
        case ("z_bed");    call sub_ref_k(pre, v, ref%z_bed, prec)
        case ("H_ice");    call sub_ref_k(pre, v, ref%H_ice, prec)
        case ("z_srf");    call sub_ref_k(pre, v, ref%z_srf, prec)
        case ("z_bed_sd"); call sub_ref_k(pre, v, ref%z_bed_sd, prec)
        end select
    end subroutine sub_ref

    subroutine sub_ref_k(pre, v, r, prec)
        real(dp), intent(inout) :: pre(:,:,:)
        real(dp), intent(in)    :: v(:,:,:), r(:,:)
        integer,  intent(in)    :: prec
        integer :: k

        do k = 1, size(v, 3)
            if (prec == sp) then
                pre(:,:,k) = real(real(v(:,:,k), sp) - real(r, sp), dp)
            else
                pre(:,:,k) = v(:,:,k) - r
            end if
        end do
    end subroutine sub_ref_k

    subroutine add_ref(work, ref, name, prec)
        ! work = ref%<name> + work, in prec (the destination's).
        real(dp),         intent(inout) :: work(:,:,:)
        type(grid_ref_t), intent(in)    :: ref
        character(len=*), intent(in)    :: name
        integer,          intent(in)    :: prec

        select case (trim(name))
        case ("z_bed");    call add_ref_k(work, ref%z_bed, prec)
        case ("H_ice");    call add_ref_k(work, ref%H_ice, prec)
        case ("z_srf");    call add_ref_k(work, ref%z_srf, prec)
        case ("z_bed_sd"); call add_ref_k(work, ref%z_bed_sd, prec)
        end select
    end subroutine add_ref

    subroutine add_ref_k(work, r, prec)
        real(dp), intent(inout) :: work(:,:,:)
        real(dp), intent(in)    :: r(:,:)
        integer,  intent(in)    :: prec
        integer :: k

        do k = 1, size(work, 3)
            if (prec == sp) then
                work(:,:,k) = real(real(r, sp) + real(work(:,:,k), sp), dp)
            else
                work(:,:,k) = r + work(:,:,k)
            end if
        end do
    end subroutine add_ref_k

    subroutine add_scalar(work, s, prec)
        ! work = s + work, in prec (the destination's).
        real(dp), intent(inout) :: work(:,:,:)
        real(dp), intent(in)    :: s
        integer,  intent(in)    :: prec

        if (prec == sp) then
            work = real(real(s, sp) + real(work, sp), dp)
        else
            work = s + work
        end if
    end subroutine add_scalar

    subroutine lapse(work, zc, bp, lp, z, is_pr, prec)
        ! Re-elevate tas (or pr) from the surface it refers to (zc) to the
        ! consumer's z, per month m, in prec (the destination's):
        !   dz = z - zc;  tas = tas - lp(:,:,m)*dz;  pr = pr*exp(e),  e = -bp*lp(:,:,m)*dz
        ! A cell with e = 0 (+0 or -0) keeps pr without evaluating exp: exp(e)
        ! is exactly 1 there and pr*1 is pr, bit for bit (the common case of a
        ! consumer on the producer's surface, dz = 0).
        real(dp), intent(inout) :: work(:,:,:)            ! tas or pr, remapped (not rounded)
        real(dp), intent(in)    :: zc(:,:,:), bp(:,:,:), lp(:,:,:), z(:,:,:)
        logical,  intent(in)    :: is_pr
        integer,  intent(in)    :: prec

        real(sp) :: dz_sp, l_sp, e_sp
        real(dp) :: dz_dp, e_dp
        integer  :: i, j, m

        do m = 1, size(work, 3)
            do j = 1, size(work, 2)
                do i = 1, size(work, 1)
                    if (prec == sp) then
                        dz_sp = real(z(i,j,1), sp) - real(zc(i,j,1), sp)
                        l_sp  = real(lp(i,j,m), sp)
                        if (is_pr) then                 ! the sp factor on the dp value
                            e_sp = -real(bp(i,j,1), sp)*l_sp*dz_sp
                            if (e_sp /= 0.0_sp) work(i,j,m) = work(i,j,m)*real(exp(e_sp), dp)
                        else
                            work(i,j,m) = real(real(work(i,j,m), sp) - l_sp*dz_sp, dp)
                        end if
                    else
                        dz_dp = z(i,j,1) - zc(i,j,1)
                        if (is_pr) then
                            e_dp = -bp(i,j,1)*lp(i,j,m)*dz_dp
                            if (e_dp /= 0.0_dp) work(i,j,m) = work(i,j,m)*exp(e_dp)
                        else
                            work(i,j,m) = work(i,j,m) - lp(i,j,m)*dz_dp
                        end if
                    end if
                end do
            end do
        end do
    end subroutine lapse

end module kryos_wiring
