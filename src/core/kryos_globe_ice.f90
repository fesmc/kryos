module kryos_globe_ice
    ! The global ice field of the globe (phase4-plan sections 4, 5, 9.1): the
    ! gather of the domains' ice into the host-resolved terms of a global solid
    ! Earth, the reference they are measured against, and (later in the phase)
    ! the global geometry hub built on the same code. One code path for the
    ! reference and for a step, so that they cannot drift apart.
    !
    ! globe_ice_t knows no domain type: the caller registers each domain with
    ! its ice grid and its static exports (footprint, H_ice_pd, z_bed_ref) and
    ! passes arrays and the ice instance (for its flotation binding) at every
    ! call; the coupler is the caller's (grids by name, map cache). Three
    ! global grids take part: the background grid (the datasrc of the globe,
    ! `[globe_ice] grid`), the target grid of the Earth (an input: VILMA2's
    ! Gauss grid, or a small grid in a test) and each domain's ice grid.
    !
    ! Sequence (the run calls it as phase 4 WP4 wires it):
    !
    !     globe_ice_init(gi, cfg, cpl, cnst, grid_tgt)       group [globe_ice], datasrc, grids
    !     globe_ice_add_domain(gi, cpl, name, grid, ...)     per domain: statics, coverage, areas
    !     globe_ice_reference_begin / _add (per domain) / _end        section 5.1
    !     globe_ice_gather_begin(gi, cpl, t) / _add (per domain) / _end    section 4.1, every geo step
    !
    ! Outputs are components of globe_ice_t (the caller copies them into the
    ! imports of the global geo): on the target grid host_w, host_Ig, host_C,
    ! host_F, H_ice_ext (step), z_bed_ref, H_ice_ref, host_Ig_ref, host_C_ref
    ! (reference); on the background grid the spliced reference z_bed_spl,
    ! H_spl and its coverage cov_bg. The diagnostics of the last gather are in
    ! gi%diag.
    !
    ! Domains are merged and summed in the order of their ice grid names, so
    ! the result does not depend on the order of the add calls. Every error
    ! stops through fesm_error naming the domain, the grid or the key.
    !
    ! Metrics (decided with phase4-plan 4.4). The amounts (grounded ice, F, the
    ! merged ice) are per true area: each ice grid's model area is the true
    ! (ellipsoidal) cell area, so the merges conserve the ice model's own mass.
    ! The weights w and C are geometric on the sphere (the coupler's coverage
    ! fraction of the target cell, in its map metric). Near the poles an
    ! ellipsoidal cell is up to 0.9 % larger than the sphere's, so at the rim
    ! cells (0 < w < 1) the amount per area and w times the cell mean differ by
    ! up to that much (ANT-32KM, GRL-16KM); in the interior they agree to
    ! round-off, and the closure of the merges is round-off (4e-14).
    !
    ! One reference definition: in a host cell the ice is the domain's present
    ! day inside its footprint and the background X elsewhere (also in domain
    ! cells outside the footprint), the bed the domain's z_bed_ref. VILMA2's
    ! reference host terms and the spliced reference of the hub are both that.

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, LOG_WARN, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t, fields_init, fields_add, fields_lock, FIELD_INSTANT, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use fesm_datasrc,    only: datasrc_t, datasrc_init, datasrc_update, datasrc_get, datasrc_describe, datasrc_end
    use coords,          only: grid_class, mv_dp
    use coupler,         only: coupler_class, coupler_merge_t, coupler_add_grid, coupler_set_area, &
                               coupler_integral, coupler_merge_begin, coupler_merge_add, coupler_merge_end, &
                               coupler_owner_reduce, remap
    use nml,             only: nml_read
    use phys_constants,  only: phys_const_class, phys_const_isset
    use kryos_config,    only: config_t, config_read, config_validate, config_find, CFG_PATH_LEN
    use kryos_ice,       only: ice_slot_t, ice_flotation

    implicit none

    private

    integer, parameter, public :: GI_NAME_LEN = 64
    character(len=*), parameter, public :: GI_DATASRC_DEFAULTS = "input/datasrc_defaults.nml"

    real(dp), parameter :: KG_TO_GT = 1.0e-12_dp         ! kg -> Gt
    real(dp), parameter :: BOUND_TOL = 1.0e-9_dp         ! slack of the bounds of host_C, host_w
    real(dp), parameter :: W_TOL = 1.0e-12_dp            ! the coverage of a step against the reference

    type, public :: gi_domain_t                          ! one registered domain
        character(len=GI_NAME_LEN) :: name = ""          ! the domain
        character(len=GI_NAME_LEN) :: grid = ""          ! its ice grid, as the coupler knows it
        integer :: nx = 0, ny = 0
        real(dp), allocatable :: footprint(:,:)          ! 1 where the ice model holds ice (static)
        real(dp), allocatable :: H_pd(:,:)               ! present-day ice thickness [m] (static)
        real(dp), allocatable :: z_bed_ref(:,:)          ! reference bed [m] (static)
        real(dp), allocatable :: f_tgt(:,:), f_bg(:,:)   ! coverage of the domain (its grid) on the two global grids
        real(dp), allocatable :: fp_tgt(:,:), fp_bg(:,:) ! covered by the footprint
        real(dp), allocatable :: ext_ref(:,:)            ! background ice on the ice grid (reference)
        real(dp), allocatable :: fg_ref(:,:)             ! reference grounded fraction
        real(dp), allocatable :: I_ref(:,:), I_last(:,:) ! grounded ice H_tot f_g [m] of the reference, of the last gather
        real(dp), allocatable :: H_ref(:,:)              ! total ice of the reference m H_pd + (1 - m) X [m]
        logical  :: has_prio = .false.                   ! overlapping domains: priority given (owner weights)
        integer  :: priority = 0                         ! higher wins the cells two domains share
        real(dp), allocatable :: owner(:,:)              ! owner weights [0,1] on the ice grid (has_prio), else absent: 1
        real(dp) :: V_ref = 0.0_dp                       ! grounded ice volume of the reference [m3]
        real(dp) :: V_src = 0.0_dp                       ! of the last gather (source grid integral)
        logical  :: ref_added = .false., added = .false.
    end type gi_domain_t

    type, public :: gi_diag_t                            ! diagnostics of the last gather
        real(dp) :: t = 0.0_dp                           ! load time [yr]
        integer  :: nd = 0
        character(len=GI_NAME_LEN), allocatable :: name(:)   ! domains, in merge order
        real(dp), allocatable :: V_g(:)                  ! grounded ice volume [m3] (host terms of the domain)
        real(dp), allocatable :: dV_g(:)                 ! against the reference [m3]
        real(dp), allocatable :: dM(:)                   ! [Gt]
        real(dp), allocatable :: esl(:)                  ! sea-level equivalent [m] (negative for growth)
        real(dp) :: ext_dV = 0.0_dp, ext_dM = 0.0_dp, ext_esl = 0.0_dp   ! the external part (1 - w)
        real(dp) :: clip_Gt = 0.0_dp                     ! external mass clipped at zero thickness [Gt]
        real(dp) :: cons_I = 0.0_dp, cons_ext = 0.0_dp   ! relative closure of the two merges
        real(dp) :: a_ocn_ref = 0.0_dp                   ! ocean area of the reference [m2]
    end type gi_diag_t

    type, public :: globe_ice_t
        logical :: ready = .false., have_ref = .false.
        logical :: in_ref = .false., in_gather = .false., owners_done = .false.
        character(len=GI_NAME_LEN) :: group = "", bg = "", tgt = ""
        integer :: nxb = 0, nyb = 0, nxg = 0, nyg = 0
        character(len=8) :: ext = "ref"                  ! "ref" | "file": ice outside the footprints
        real(dp) :: ext_tref = 0.0_dp                    ! [yr] model time of the history's present day
        real(dp) :: overlap_tol = 1.0e-2_dp, cons_tol = 1.0e-9_dp
        real(dp) :: rho_ice = 0.0_dp, rho_sw = 0.0_dp, area_ocn = 0.0_dp      ! [kg m-3], [kg m-3], [m2]
        ! the datasrc of the globe (background grid): X at ext_tref, the history as an anomaly
        type(field_state_t) :: st_x, st_ext
        type(datasrc_t)     :: ds_x, ds_ext
        logical :: has_ds_ext = .false.
        ! background X on the background grid and on the target grid
        real(dp), allocatable :: z_bed_X(:,:), H_X(:,:)
        real(dp), allocatable :: z_bed_ref(:,:), H_ice_ref(:,:)                ! target grid: the 1 - w part
        ! reference host terms (target grid)
        real(dp), allocatable :: host_w(:,:), host_Ig_ref(:,:), host_C_ref(:,:)
        ! the spliced reference (background grid, section 5.1 item 3)
        real(dp), allocatable :: z_bed_spl(:,:), H_spl(:,:), cov_bg(:,:)
        ! the last gather
        real(dp), allocatable :: H_ext(:,:)                                    ! background grid
        real(dp), allocatable :: H_ice_ext(:,:), host_Ig(:,:), host_C(:,:), host_F(:,:)   ! target grid
        real(dp) :: t_load = 0.0_dp
        type(gi_diag_t) :: diag
        type(gi_domain_t), allocatable :: dom(:)
        integer, allocatable :: ord(:)                   ! dom indices in order of the ice grid names
        type(coupler_merge_t) :: acc_I, acc_C, acc_F, acc_zs, acc_hs
    end type globe_ice_t

    public :: globe_ice_init, globe_ice_add_domain, globe_ice_end
    public :: globe_ice_reference_begin, globe_ice_reference_add, globe_ice_reference_end
    public :: globe_ice_gather_begin, globe_ice_gather_add, globe_ice_gather_end

contains

    ! =========================================================================
    ! Init
    ! =========================================================================

    subroutine globe_ice_init(gi, cfg, cpl, cnst, grid_tgt, grid_bg, group, time_ref, time_span, datasrc_defaults)
        ! The group `group` ([globe_ice], template [globe_ice]): validated, read,
        ! the background grid (grid_bg, else grid_<grid>.txt of the map folder)
        ! and the target grid registered in the coupler, and the datasrc of the
        ! globe bound: z_bed_ref and H_ice_ref (X, read at ext_tref) and, for
        ! ext = "file", ext_field (the history as an anomaly against its own
        ! slice at ext_tref). time_ref = calendar year of model time zero
        ! (default 1950); time_span = the run [yr], checked against the history.
        type(globe_ice_t),      intent(out)   :: gi
        type(config_t),         intent(in)    :: cfg
        type(coupler_class),    intent(inout) :: cpl
        type(phys_const_class), intent(in)    :: cnst
        type(grid_class),       intent(in)    :: grid_tgt
        type(grid_class),       intent(in), optional :: grid_bg
        character(len=*),       intent(in), optional :: group
        real(dp),               intent(in), optional :: time_ref, time_span(2)
        character(len=*),       intent(in), optional :: datasrc_defaults

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_init"
        character(len=GI_NAME_LEN) :: key_grid, key_zb, key_hx, key_ext, ext_mode
        character(len=CFG_PATH_LEN) :: dsdef
        character(len=CFG_PATH_LEN), allocatable :: lib(:)
        type(grid_class) :: gbg
        real(dp) :: tref, span(2)
        integer :: h

        gi%group = "globe_ice"
        if (present(group)) gi%group = group
        dsdef = GI_DATASRC_DEFAULTS
        if (present(datasrc_defaults)) dsdef = datasrc_defaults
        tref = 1950.0_dp
        if (present(time_ref)) tref = time_ref

        call config_validate(cfg, trim(gi%group), "globe_ice")
        call config_read(cfg, trim(gi%group), "grid",        key_grid, "globe_ice")
        call config_read(cfg, trim(gi%group), "z_bed_ref",   key_zb,   "globe_ice")
        call config_read(cfg, trim(gi%group), "H_ice_ref",   key_hx,   "globe_ice")
        call config_read(cfg, trim(gi%group), "ext",         ext_mode, "globe_ice")
        call config_read(cfg, trim(gi%group), "ext_field",   key_ext,  "globe_ice")
        call config_read(cfg, trim(gi%group), "ext_tref",    gi%ext_tref,    "globe_ice")
        call config_read(cfg, trim(gi%group), "overlap_tol", gi%overlap_tol, "globe_ice")
        call config_read(cfg, trim(gi%group), "cons_tol",    gi%cons_tol,    "globe_ice")

        gi%ext = trim(ext_mode)
        if (trim(gi%ext) /= "ref" .and. trim(gi%ext) /= "file") then
            call fesm_error(proc, "[" //trim(gi%group)//"] ext = '"//trim(ext_mode)//"' is not recognised", &
                            "use ""ref"" or ""file""; par file "//trim(cfg%par_file))
        end if
        if (len_trim(key_grid) == 0) then
            call fesm_error(proc, "["//trim(gi%group)//"] grid is empty", "the background grid, e.g. GLOBAL-0.5DEG")
        end if
        if (len_trim(key_zb) == 0 .or. len_trim(key_hx) == 0) then
            call fesm_error(proc, "["//trim(gi%group)//"] z_bed_ref and H_ice_ref must name field groups (X)", &
                            "z_bed_ref = """"? par file "//trim(cfg%par_file))
        end if
        if (trim(gi%ext) == "file" .and. len_trim(key_ext) == 0) then
            call fesm_error(proc, "["//trim(gi%group)//"] ext = ""file"" needs ext_field (the history)", &
                            "par file "//trim(cfg%par_file))
        end if

        if (.not. phys_const_isset(cnst%rho_ice) .or. .not. phys_const_isset(cnst%rho_sw) .or. &
            .not. phys_const_isset(cnst%area_seasurf)) then
            call fesm_error(proc, "the run's constants lack rho_ice, rho_sw or area_seasurf", "phys_const group "// &
                            trim(cnst%label))
        end if
        gi%rho_ice  = cnst%rho_ice
        gi%rho_sw   = cnst%rho_sw
        gi%area_ocn = cnst%area_seasurf*1.0e6_dp         ! km2 -> m2

        ! --- the grids ---
        if (present(grid_bg)) then
            gbg = grid_bg
            if (trim(gbg%name) /= trim(key_grid)) then
                call fesm_error(proc, "the background grid is named '"//trim(gbg%name)//"', [" // &
                                trim(gi%group)//"] grid = '"//trim(key_grid)//"'", "they must agree")
            end if
        else
            call grid_from_disk(gbg, trim(key_grid), trim(cfg%map_fldr))
        end if
        gi%bg  = trim(key_grid)
        gi%tgt = trim(grid_tgt%name)
        if (trim(gi%bg) == trim(gi%tgt)) then
            call fesm_error(proc, "the background grid and the target grid have one name: '"//trim(gi%bg)//"'", &
                            "two distinct global grids are needed")
        end if
        call coupler_add_grid(cpl, trim(gi%bg),  gbg)
        call coupler_add_grid(cpl, trim(gi%tgt), grid_tgt)
        gi%nxb = gbg%G%nx;      gi%nyb = gbg%G%ny
        gi%nxg = grid_tgt%G%nx; gi%nyg = grid_tgt%G%ny

        ! --- the datasrc: X (read at ext_tref) and the history (anomaly) ---
        lib = pack(cfg%library, len_trim(cfg%library) > 0)
        call fields_init(gi%st_x, trim(gi%group)//":X")
        call fields_add(gi%st_x, "z_bed_ref", "m", FIELD_INSTANT, FIELD_INTENSIVE, h, grid=trim(gi%bg), dims=[gi%nxb, gi%nyb])
        call fields_add(gi%st_x, "H_ice_ref", "m", FIELD_INSTANT, FIELD_EXTENSIVE, h, grid=trim(gi%bg), dims=[gi%nxb, gi%nyb])
        call fields_lock(gi%st_x)
        call datasrc_init(gi%ds_x, trim(cfg%par_file), trim(gi%group), gi%st_x, &
                          [character(len=16) :: "z_bed_ref", "H_ice_ref"], gbg, trim(dsdef), var_defs=lib, &
                          time_ref=tref)
        call datasrc_describe(gi%ds_x)

        gi%has_ds_ext = (trim(gi%ext) == "file")
        if (gi%has_ds_ext) then
            call check_history_group(cfg, trim(key_ext), gi%ext_tref, trim(dsdef), trim(gi%group))
            call fields_init(gi%st_ext, trim(gi%group)//":ext")
            call fields_add(gi%st_ext, "ext_field", "m", FIELD_INSTANT, FIELD_EXTENSIVE, h, grid=trim(gi%bg), &
                            dims=[gi%nxb, gi%nyb])
            call fields_lock(gi%st_ext)
            if (present(time_span)) then
                span = [time_span(1), min(time_span(2), gi%ext_tref)]
                call datasrc_init(gi%ds_ext, trim(cfg%par_file), trim(gi%group), gi%st_ext, ["ext_field"], gbg, &
                                  trim(dsdef), var_defs=lib, time_ref=tref, time_span=span)
            else
                call datasrc_init(gi%ds_ext, trim(cfg%par_file), trim(gi%group), gi%st_ext, ["ext_field"], gbg, &
                                  trim(dsdef), var_defs=lib, time_ref=tref)
            end if
            call datasrc_describe(gi%ds_ext)
        end if

        allocate(gi%z_bed_X(gi%nxb, gi%nyb), gi%H_X(gi%nxb, gi%nyb), gi%H_ext(gi%nxb, gi%nyb))
        allocate(gi%z_bed_spl(gi%nxb, gi%nyb), gi%H_spl(gi%nxb, gi%nyb), gi%cov_bg(gi%nxb, gi%nyb))
        allocate(gi%z_bed_ref(gi%nxg, gi%nyg), gi%H_ice_ref(gi%nxg, gi%nyg), gi%host_w(gi%nxg, gi%nyg))
        allocate(gi%host_Ig_ref(gi%nxg, gi%nyg), gi%host_C_ref(gi%nxg, gi%nyg))
        allocate(gi%H_ice_ext(gi%nxg, gi%nyg), gi%host_Ig(gi%nxg, gi%nyg), gi%host_C(gi%nxg, gi%nyg), &
                 gi%host_F(gi%nxg, gi%nyg))
        allocate(gi%dom(0), gi%ord(0))
        gi%ready = .true.

        call log_msg(LOG_INFO, "kryos_globe_ice", "["//trim(gi%group)//"] background "//trim(gi%bg)//" ("// &
                     to_str(gi%nxb)//" x "//to_str(gi%nyb)//"), target "//trim(gi%tgt)//" ("//to_str(gi%nxg)// &
                     " x "//to_str(gi%nyg)//"), ice outside the footprints: "//trim(gi%ext)// &
                     hist_note(gi))
    end subroutine globe_ice_init

    function hist_note(gi) result(s)
        type(globe_ice_t), intent(in) :: gi
        character(len=:), allocatable :: s
        s = ""
        if (gi%has_ds_ext) s = ", anomaly against the history at t = "//to_str(gi%ext_tref)//" yr"
    end function hist_note

    subroutine grid_from_disk(g, name, fldr)
        ! The grid `name` from grid_<name>.txt of the map folder.
        use coords, only: grid_cdo_read_desc
        type(grid_class), intent(out) :: g
        character(len=*), intent(in)  :: name, fldr
        call grid_cdo_read_desc(g, name, fldr)
    end subroutine grid_from_disk

    subroutine check_history_group(cfg, key, ext_tref, dsdef, group)
        ! The field group of the history must be the anomaly against its own
        ! slice at ext_tref (anomaly_ref = slice, ref_time = ext_tref) and stop
        ! outside the file (out_of_range = error): the hold after ext_tref is
        ! the gather's (it evaluates at min(t, ext_tref)), so a "hold" in the
        ! group would also hide a run that starts before the file.
        type(config_t),   intent(in) :: cfg
        character(len=*), intent(in) :: key, dsdef, group
        real(dp),         intent(in) :: ext_tref

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_init"
        character(len=CFG_PATH_LEN) :: file
        character(len=32) :: s
        real(dp) :: rt

        file = config_find(cfg, key, "["//group//"] ext_field")
        s = ""
        call nml_read(trim(file), key, "mode", s, defaults_file=dsdef, defaults_group="datasrc_field")
        if (trim(s) /= "transient" .and. trim(s) /= "auto") then
            call fesm_error(proc, "the history group &"//trim(key)//" has mode = '"//trim(s)//"'", &
                            "the history is a transient field (mode = ""transient""); file "//trim(file))
        end if
        call nml_read(trim(file), key, "anomaly_ref", s, defaults_file=dsdef, defaults_group="datasrc_field")
        if (trim(s) /= "slice") then
            call fesm_error(proc, "the history group &"//trim(key)//" has anomaly_ref = '"//trim(s)//"'", &
                            "the history enters as an anomaly against its own present day: anomaly_ref = ""slice"", "// &
                            "ref_time = [globe_ice] ext_tref; file "//trim(file))
        end if
        call nml_read(trim(file), key, "ref_time", rt, defaults_file=dsdef, defaults_group="datasrc_field")
        if (rt /= ext_tref) then
            call fesm_error(proc, "the history group &"//trim(key)//" has ref_time = "//to_str(rt)// &
                            " but ["//group//"] ext_tref = "//to_str(ext_tref), &
                            "the two name one time (the history's present day); file "//trim(file))
        end if
        call nml_read(trim(file), key, "out_of_range", s, defaults_file=dsdef, defaults_group="datasrc_field")
        if (trim(s) /= "error") then
            call fesm_error(proc, "the history group &"//trim(key)//" has out_of_range = '"//trim(s)//"'", &
                            "leave ""error"": the gather holds the history after ext_tref itself and a run before "// &
                            "the first record must stop; file "//trim(file))
        end if
    end subroutine check_history_group

    ! =========================================================================
    ! Domains
    ! =========================================================================

    subroutine globe_ice_add_domain(gi, cpl, name, grid, footprint, H_ice_pd, z_bed_ref, priority)
        ! Register a domain: its ice grid in the coupler with the true spherical
        ! cell areas as model areas (section 4.4: the load per true area, global
        ! integrals spherical; the plane dx*dy is Yelmo's), its static exports
        ! (copied), and its coverage and footprint fraction on the two global
        ! grids. Call before the reference; any order. priority (optional): for
        ! domains that overlap; the domain with the higher priority owns the cells
        ! two share (owner weights, coupler_owner_reduce, section 5.4), so the
        ! result does not depend on the order of the calls either. Domains without
        ! a priority, or with equal ones, must not overlap (the merge stops).
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(grid_class),    intent(in)    :: grid
        real(dp),            intent(in)    :: footprint(:,:), H_ice_pd(:,:), z_bed_ref(:,:)
        integer,             intent(in), optional :: priority

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_add_domain"
        type(gi_domain_t) :: dm
        type(gi_domain_t), allocatable :: tmp(:)
        real(dp), allocatable :: ones(:,:), v(:,:)
        integer :: k, n

        call require_ready(gi, proc)
        if (gi%have_ref .or. gi%in_ref .or. gi%in_gather) then
            call fesm_error(proc, "domain `"//trim(name)//"` registered after the reference was started", &
                            "register every domain first")
        end if
        do k = 1, size(gi%dom)
            if (trim(gi%dom(k)%name) == trim(name)) then
                call fesm_error(proc, "domain `"//trim(name)//"` is registered twice", "")
            end if
            if (trim(gi%dom(k)%grid) == trim(grid%name)) then
                call fesm_error(proc, "domains `"//trim(gi%dom(k)%name)//"` and `"//trim(name)// &
                                "` are on one ice grid ("//trim(grid%name)//")", "the merge names its sources by grid")
            end if
        end do
        if (trim(grid%name) == trim(gi%bg) .or. trim(grid%name) == trim(gi%tgt)) then
            call fesm_error(proc, "domain `"//trim(name)//"`: the ice grid has the name of a global grid ("// &
                            trim(grid%name)//")", "")
        end if

        dm%name = trim(name)
        dm%grid = trim(grid%name)
        dm%nx = grid%G%nx
        dm%ny = grid%G%ny
        call check_shape(footprint, dm, "footprint", proc)
        call check_shape(H_ice_pd,  dm, "H_ice_pd",  proc)
        call check_shape(z_bed_ref, dm, "z_bed_ref", proc)
        call check_finite(footprint, dm, "footprint", proc)
        call check_finite(H_ice_pd,  dm, "H_ice_pd",  proc)
        call check_finite(z_bed_ref, dm, "z_bed_ref", proc)
        if (any(footprint /= 0.0_dp .and. footprint /= 1.0_dp)) then
            call fesm_error(proc, "domain `"//trim(name)//"`: the footprint holds values other than 0 and 1", "")
        end if
        if (any(H_ice_pd < 0.0_dp)) then
            call fesm_error(proc, "domain `"//trim(name)//"`: H_ice_pd holds negative thickness", "")
        end if
        dm%footprint = footprint
        dm%H_pd      = H_ice_pd
        dm%z_bed_ref = z_bed_ref
        if (present(priority)) then
            dm%has_prio = .true.
            dm%priority = priority
        end if

        call coupler_add_grid(cpl, trim(dm%grid), grid)
        call coupler_set_area(cpl, trim(dm%grid), grid%area*grid%cs%xy_conv**2)

        ! coverage of the grid and of the footprint on the global grids (map metric, as the merge's)
        allocate(ones(dm%nx, dm%ny), v(dm%nx, dm%ny))
        ones = 1.0_dp
        v = merge(1.0_dp, mv_dp, footprint > 0.5_dp)
        allocate(dm%f_tgt(gi%nxg, gi%nyg), dm%fp_tgt(gi%nxg, gi%nyg), dm%f_bg(gi%nxb, gi%nyb), dm%fp_bg(gi%nxb, gi%nyb))
        call coverage(cpl, ones, dm%grid, gi%tgt, dm%f_tgt, mv_dp)
        call coverage(cpl, ones, dm%grid, gi%bg,  dm%f_bg,  mv_dp)
        call coverage(cpl, v,    dm%grid, gi%tgt, dm%fp_tgt, mv_dp)
        call coverage(cpl, v,    dm%grid, gi%bg,  dm%fp_bg,  mv_dp)

        n = size(gi%dom)
        allocate(tmp(n + 1))
        tmp(1:n) = gi%dom
        tmp(n + 1) = dm
        call move_alloc(tmp, gi%dom)
        call sort_order(gi)
    end subroutine globe_ice_add_domain

    subroutine merge_add(acc, cpl, dm, v)
        ! coupler_merge_add of a domain, with its owner weights when it has them.
        type(coupler_merge_t), intent(inout) :: acc
        type(coupler_class),   intent(inout) :: cpl
        type(gi_domain_t),     intent(in)    :: dm
        real(dp),              intent(in)    :: v(:,:)
        if (allocated(dm%owner)) then
            call coupler_merge_add(acc, cpl, trim(dm%grid), v, owner=dm%owner)
        else
            call coupler_merge_add(acc, cpl, trim(dm%grid), v)
        end if
    end subroutine merge_add

    function owned_integral(cpl, dm, v) result(tot)
        ! Integral of v over the part of the domain it owns (true areas).
        type(coupler_class), intent(in) :: cpl
        type(gi_domain_t),   intent(in) :: dm
        real(dp),            intent(in) :: v(:,:)
        real(dp) :: tot
        if (allocated(dm%owner)) then
            tot = coupler_integral(cpl, trim(dm%grid), dm%owner*v)
        else
            tot = coupler_integral(cpl, trim(dm%grid), v)
        end if
    end function owned_integral

    subroutine finalize_owners(gi, cpl)
        ! Owner weights of the domains that have a priority (once, before the
        ! first reference): in the order of descending priority (ties by ice
        ! grid name), a domain's owner is 1 on its grid minus the part already
        ! owned by each domain of higher priority that shares a global cell with
        ! it (coupler_owner_reduce, which needs the ice-to-ice map: pairs that
        ! share no cell are skipped). The coverage and footprint fractions are
        ! then the owned ones.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl

        integer, allocatable :: ord(:)
        real(dp), allocatable :: ones(:,:), v(:,:), tmp(:,:)
        integer :: i, j, k, a, b, n
        logical :: share

        gi%owners_done = .true.
        n = size(gi%dom)
        if (.not. any(gi%dom(:)%has_prio)) return

        ord = pack(gi%ord, gi%dom(gi%ord)%has_prio)          ! name order of the domains with a priority
        do i = 2, size(ord)                                   ! stable insertion sort, priority descending
            k = ord(i); j = i - 1
            do while (j >= 1)
                if (gi%dom(ord(j))%priority >= gi%dom(k)%priority) exit
                ord(j + 1) = ord(j); j = j - 1
            end do
            ord(j + 1) = k
        end do

        do i = 1, size(ord)
            associate (lo => gi%dom(ord(i)))
                allocate(lo%owner(lo%nx, lo%ny))
                lo%owner = 1.0_dp
                do j = 1, i - 1
                    associate (hi => gi%dom(ord(j)))
                        if (hi%priority == lo%priority) cycle
                        share = any(hi%f_tgt > 0.0_dp .and. lo%f_tgt > 0.0_dp) .or. any(hi%f_bg > 0.0_dp .and. lo%f_bg > 0.0_dp)
                        if (.not. share) cycle
                        call coupler_owner_reduce(cpl, trim(hi%grid), hi%owner, trim(lo%grid), lo%owner)
                    end associate
                end do
            end associate
        end do
        a = 0; b = 0
        do i = 1, size(ord)
            associate (dm => gi%dom(ord(i)))
                allocate(ones(dm%nx, dm%ny), v(dm%nx, dm%ny))
                ones = dm%owner
                v = merge(dm%owner, mv_dp, dm%footprint > 0.5_dp)
                call owned_coverage(cpl, ones, dm%grid, gi%tgt, dm%f_tgt)
                call owned_coverage(cpl, ones, dm%grid, gi%bg,  dm%f_bg)
                call owned_coverage(cpl, v,    dm%grid, gi%tgt, dm%fp_tgt)
                call owned_coverage(cpl, v,    dm%grid, gi%bg,  dm%fp_bg)
                deallocate(ones, v)
                if (any(dm%owner < 1.0_dp)) a = a + 1
                b = b + 1
            end associate
        end do
        if (allocated(tmp)) deallocate(tmp)
        call log_msg(LOG_INFO, "kryos_globe_ice", "owner weights for "//to_str(b)//" domains with a priority, "// &
                     to_str(a)//" of them share cells with a domain of higher priority")
    end subroutine finalize_owners

    subroutine owned_coverage(cpl, v, src, dst, f)
        ! Coverage of a cell of `dst` by the valid cells of `src` (v /= mv_dp), weighted by v: the
        ! coverage times the mean owner weight, as the merge's.
        type(coupler_class), intent(inout) :: cpl
        real(dp),            intent(in)    :: v(:,:)
        character(len=*),    intent(in)    :: src, dst
        real(dp),            intent(out)   :: f(:,:)
        real(dp), allocatable :: mean(:,:), cov(:,:)
        allocate(cov(size(f,1), size(f,2)))
        call remap(cpl, v, trim(src), mean, trim(dst), method="con", missing_value=mv_dp, frac2=cov)
        f = cov*merge(mean, 0.0_dp, mean /= mv_dp)
    end subroutine owned_coverage

    subroutine coverage(cpl, v, src, dst, f, miss)
        ! Fraction of each cell of `dst` covered by the valid cells (v /= miss) of `src`.
        type(coupler_class), intent(inout) :: cpl
        real(dp),            intent(in)    :: v(:,:), miss
        character(len=*),    intent(in)    :: src, dst
        real(dp),            intent(out)   :: f(:,:)
        real(dp), allocatable :: tmp(:,:)
        call remap(cpl, v, trim(src), tmp, trim(dst), method="con", missing_value=miss, frac2=f)
    end subroutine coverage

    subroutine sort_order(gi)
        ! gi%ord: the domains in the order of their ice grid names (the merge's order).
        type(globe_ice_t), intent(inout) :: gi
        integer :: i, j, t, n
        n = size(gi%dom)
        if (allocated(gi%ord)) deallocate(gi%ord)
        allocate(gi%ord(n))
        gi%ord = [(i, i = 1, n)]
        do i = 2, n
            t = gi%ord(i); j = i - 1
            do while (j >= 1)
                if (.not. llt(gi%dom(t)%grid, gi%dom(gi%ord(j))%grid)) exit
                gi%ord(j + 1) = gi%ord(j); j = j - 1
            end do
            gi%ord(j + 1) = t
        end do
    end subroutine sort_order

    integer function find_domain(gi, name, proc) result(d)
        type(globe_ice_t), intent(in) :: gi
        character(len=*),  intent(in) :: name, proc
        integer :: k
        d = 0
        do k = 1, size(gi%dom)
            if (trim(gi%dom(k)%name) == trim(name)) then
                d = k
                return
            end if
        end do
        call fesm_error(proc, "domain `"//trim(name)//"` is not registered", "globe_ice_add_domain first")
    end function find_domain

    ! =========================================================================
    ! Reference (section 5.1)
    ! =========================================================================

    subroutine globe_ice_reference_begin(gi, cpl)
        ! Background X on the background grid (datasrc at ext_tref) and on the
        ! target grid (z_bed_ref, H_ice_ref: the 1 - w part of VILMA2's
        ! reference), then the accumulators of the reference host terms and of
        ! the spliced reference. May be called again: the reference is rebuilt
        ! and is bit identical.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_reference_begin"
        real(dp), allocatable :: tmp(:,:)
        integer :: k

        call require_ready(gi, proc)
        if (gi%in_ref .or. gi%in_gather) then
            call fesm_error(proc, "a reference or a gather is open", "close it with its _end first")
        end if

        if (.not. gi%owners_done) call finalize_owners(gi, cpl)
        call datasrc_update(gi%ds_x, gi%ext_tref)
        call datasrc_get(gi%ds_x, "z_bed_ref", gi%z_bed_X)
        call datasrc_get(gi%ds_x, "H_ice_ref", gi%H_X)
        call check_bg(gi%z_bed_X, "z_bed (X)", proc)
        call check_bg(gi%H_X,     "H_ice (X)", proc)
        if (any(gi%H_X < 0.0_dp)) then
            call fesm_error(proc, "the background ice X holds negative thickness", "minval = "//to_str(minval(gi%H_X)))
        end if

        call remap(cpl, gi%z_bed_X, trim(gi%bg), tmp, trim(gi%tgt), method="con")
        gi%z_bed_ref = tmp
        call ext_to_tgt(cpl, trim(gi%bg), trim(gi%tgt), gi%H_X, gi%H_ice_ref)

        call coupler_merge_begin(gi%acc_I,  cpl, trim(gi%tgt), "extensive")
        call coupler_merge_begin(gi%acc_C,  cpl, trim(gi%tgt), "intensive")
        call coupler_merge_begin(gi%acc_zs, cpl, trim(gi%bg),  "intensive")
        call coupler_merge_begin(gi%acc_hs, cpl, trim(gi%bg),  "extensive")
        do k = 1, size(gi%dom)
            gi%dom(k)%ref_added = .false.
        end do
        gi%in_ref = .true.
    end subroutine globe_ice_reference_begin

    subroutine globe_ice_reference_add(gi, cpl, name, ice)
        ! Domain `name` in the reference: its present day with the background ice
        ! outside the footprint, the reference flotation of the ice instance
        ! (sea level 0, the bed z_bed_ref in every cell), the host terms of the
        ! same procedure as a step, and the same ice and bed for the spliced reference.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(ice_slot_t),    intent(in)    :: ice

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_reference_add"
        real(dp), allocatable :: H_tot(:,:), f_g(:,:), I(:,:), Cc(:,:), zero(:,:)
        integer :: d

        if (.not. gi%in_ref) call fesm_error(proc, "no reference is open", "globe_ice_reference_begin first")
        d = find_domain(gi, name, proc)
        associate (dm => gi%dom(d))
            if (dm%ref_added) call fesm_error(proc, "domain `"//trim(name)//"` was already added", "")
            call ext_to_ice(cpl, trim(gi%bg), trim(dm%grid), gi%H_X, dm%ext_ref)
            allocate(zero(dm%nx, dm%ny))
            zero = 0.0_dp
            call host_arrays(dm, ice, dm%ext_ref, dm%H_pd, dm%z_bed_ref, zero, H_tot, f_g)
            dm%fg_ref = f_g
            I  = H_tot*f_g
            Cc = 1.0_dp - f_g
            call merge_add(gi%acc_I, cpl, dm, I)
            call merge_add(gi%acc_C, cpl, dm, Cc)
            call merge_add(gi%acc_zs, cpl, dm, dm%z_bed_ref)
            call merge_add(gi%acc_hs, cpl, dm, H_tot)
            dm%I_ref = I
            dm%H_ref = H_tot
            dm%V_ref = owned_integral(cpl, dm, I)
            dm%ref_added = .true.
        end associate
    end subroutine globe_ice_reference_add

    subroutine globe_ice_reference_end(gi, cpl)
        ! Close the reference: host_w, host_Ig_ref, host_C_ref on the target grid
        ! (F_ref = 0), the spliced reference on the background grid
        ! ((1 - sum f_d) X + sum f_d H_ref,d; the bed likewise), the ocean area of the reference, the log.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_reference_end"
        real(dp), allocatable :: zero(:,:), fg_bg(:,:), c_bg(:,:), c_tgt(:,:), tmp(:,:), wtmp(:,:)
        real(dp) :: a_ocn
        integer :: k, d
        character(len=:), allocatable :: line

        if (.not. gi%in_ref) call fesm_error(proc, "no reference is open", "globe_ice_reference_begin first")
        call all_added(gi, gi%dom(:)%ref_added, proc)

        allocate(zero(gi%nxg, gi%nyg), wtmp(gi%nxg, gi%nyg))
        zero = 0.0_dp
        call coupler_merge_end(gi%acc_I, cpl, zero, tmp, frac_tot=gi%host_w, overlap_tol=gi%overlap_tol)
        gi%host_Ig_ref = tmp
        call coupler_merge_end(gi%acc_C, cpl, zero, tmp, frac_tot=wtmp, overlap_tol=gi%overlap_tol)
        gi%host_C_ref = tmp
        call check_host(gi, gi%host_Ig_ref, gi%host_C_ref, gi%host_w, "the reference", proc)

        call coupler_merge_end(gi%acc_zs, cpl, gi%z_bed_X, tmp, frac_tot=gi%cov_bg, overlap_tol=gi%overlap_tol)
        gi%z_bed_spl = tmp
        call coupler_merge_end(gi%acc_hs, cpl, gi%H_X, tmp, frac_tot=gi%cov_bg, overlap_tol=gi%overlap_tol)
        gi%H_spl = tmp

        ! ocean area of the reference: the host ocean function plus the background's (1 - w) part
        ! (the background flotation is the hydrostatic rule, for this log only)
        allocate(fg_bg(gi%nxb, gi%nyb))
        fg_bg = merge(1.0_dp, 0.0_dp, gi%rho_ice*gi%H_X/gi%rho_sw + gi%z_bed_X >= 0.0_dp)
        c_bg = 1.0_dp - fg_bg
        call remap(cpl, c_bg, trim(gi%bg), tmp, trim(gi%tgt), method="con", norm="dst")
        c_tgt = gi%host_C_ref + max(0.0_dp, 1.0_dp - gi%host_w)*tmp
        a_ocn = coupler_integral(cpl, trim(gi%tgt), c_tgt)
        gi%diag%a_ocn_ref = a_ocn

        line = "reference: coverage of the host region up to "//to_str(maxval(gi%host_w))
        do k = 1, size(gi%ord)
            d = gi%ord(k)
            line = line//new_line("a")//"    domain "//trim(gi%dom(d)%name)//" ("//trim(gi%dom(d)%grid)// &
                   "): grounded ice "//to_str(gi%dom(d)%V_ref*gi%rho_ice*KG_TO_GT)//" Gt, footprint "// &
                   to_str(count(gi%dom(d)%footprint > 0.5_dp))//" of "//to_str(gi%dom(d)%nx*gi%dom(d)%ny)//" cells"
        end do
        line = line//new_line("a")//"    ocean area of the reference "//to_str(a_ocn)//" m2 ("// &
               to_str(100.0_dp*(a_ocn - gi%area_ocn)/gi%area_ocn)//" % of the sea surface "//to_str(gi%area_ocn)//" m2)"
        call log_msg(LOG_INFO, "kryos_globe_ice", line)
        if (abs(a_ocn - gi%area_ocn) > 0.01_dp*gi%area_ocn) then
            call log_msg(LOG_WARN, "kryos_globe_ice", "the ocean area of the reference differs from the sea surface by "// &
                         "more than 1 %: "//to_str(a_ocn)//" m2 against "//to_str(gi%area_ocn)//" m2")
        end if

        gi%in_ref = .false.
        gi%have_ref = .true.
    end subroutine globe_ice_reference_end

    ! =========================================================================
    ! A step (section 4.1)
    ! =========================================================================

    subroutine globe_ice_gather_begin(gi, cpl, t_load)
        ! The external ice at the load time on the background grid,
        ! H_ext = max(0, X + dH_hist) ("ref": X), and on the target grid
        ! (H_ice_ext); opens the three accumulators. The history is evaluated at
        ! min(t_load, ext_tref): held at its present day after it, a stop before
        ! the first record.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl
        real(dp),            intent(in)    :: t_load

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_gather_begin"
        real(dp), allocatable :: dH(:,:), raw(:,:)
        integer :: k

        call require_ready(gi, proc)
        if (.not. gi%have_ref) call fesm_error(proc, "no reference yet", "run the reference first")
        if (gi%in_ref .or. gi%in_gather) then
            call fesm_error(proc, "a reference or a gather is open", "close it with its _end first")
        end if
        gi%t_load = t_load

        gi%diag%clip_Gt = 0.0_dp
        if (gi%has_ds_ext) then
            allocate(dH(gi%nxb, gi%nyb), raw(gi%nxb, gi%nyb))
            call datasrc_update(gi%ds_ext, min(t_load, gi%ext_tref))
            call datasrc_get(gi%ds_ext, "ext_field", dH)
            call check_bg(dH, "the history anomaly", proc)
            raw = gi%H_X + dH
            gi%H_ext = max(0.0_dp, raw)
            gi%diag%clip_Gt = coupler_integral(cpl, trim(gi%bg), gi%H_ext - raw)*gi%rho_ice*KG_TO_GT
        else
            gi%H_ext = gi%H_X
        end if
        call ext_to_tgt(cpl, trim(gi%bg), trim(gi%tgt), gi%H_ext, gi%H_ice_ext)

        call coupler_merge_begin(gi%acc_I, cpl, trim(gi%tgt), "extensive")
        call coupler_merge_begin(gi%acc_C, cpl, trim(gi%tgt), "intensive")
        call coupler_merge_begin(gi%acc_F, cpl, trim(gi%tgt), "extensive")
        do k = 1, size(gi%dom)
            gi%dom(k)%added = .false.
        end do
        gi%in_gather = .true.
    end subroutine globe_ice_gather_begin

    subroutine globe_ice_gather_add(gi, cpl, name, ice, H_ice, f_grnd, z_bed, z_sl, footprint)
        ! Domain `name` at the load time: the ice model's thickness, grounded
        ! fraction, bed and sea level (its exports of the previous step), and its
        ! footprint (compared with the registered one: a change stops).
        !   H_tot = m H + (1 - m) ext_i;  f_g = m f_Y + (1 - m) flot(H_tot, z_bed, z_sl)
        !   I = H_tot f_g,  C = 1 - f_g,  F = zeta0 (f_g,ref - f_g),  zeta0 = z_bed_ref
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(ice_slot_t),    intent(in)    :: ice
        real(dp),            intent(in)    :: H_ice(:,:), f_grnd(:,:), z_bed(:,:), z_sl(:,:), footprint(:,:)

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_gather_add"
        real(dp), allocatable :: ext_i(:,:), H_tot(:,:), f_g(:,:), I(:,:), Cc(:,:), F(:,:)
        integer :: d

        if (.not. gi%in_gather) call fesm_error(proc, "no gather is open", "globe_ice_gather_begin first")
        d = find_domain(gi, name, proc)
        associate (dm => gi%dom(d))
            if (dm%added) call fesm_error(proc, "domain `"//trim(name)//"` was already added", "")
            call check_shape(H_ice,  dm, "H_ice",  proc)
            call check_shape(f_grnd, dm, "f_grnd", proc)
            call check_shape(z_bed,  dm, "z_bed",  proc)
            call check_shape(z_sl,   dm, "z_sl",   proc)
            call check_shape(footprint, dm, "footprint", proc)
            call check_finite(H_ice,  dm, "H_ice",  proc)
            call check_finite(f_grnd, dm, "f_grnd", proc)
            call check_finite(z_bed,  dm, "z_bed",  proc)
            call check_finite(z_sl,   dm, "z_sl",   proc)
            if (any(footprint /= dm%footprint)) then
                call fesm_error(proc, "the footprint of domain `"//trim(name)//"` changed ("// &
                                to_str(count(footprint /= dm%footprint))//" cells differ from the one registered)", &
                                "the footprint is fixed (phase4-plan 4.2)")
            end if
            if (any(H_ice < 0.0_dp)) then
                call fesm_error(proc, "domain `"//trim(name)//"`: negative ice thickness (min "//to_str(minval(H_ice))//")", "")
            end if

            call ext_to_ice(cpl, trim(gi%bg), trim(dm%grid), gi%H_ext, ext_i)
            call host_arrays(dm, ice, ext_i, H_ice, z_bed, z_sl, H_tot, f_g, f_grnd)
            I  = H_tot*f_g
            Cc = 1.0_dp - f_g
            F  = dm%z_bed_ref*(dm%fg_ref - f_g)
            call merge_add(gi%acc_I, cpl, dm, I)
            call merge_add(gi%acc_C, cpl, dm, Cc)
            call merge_add(gi%acc_F, cpl, dm, F)
            dm%I_last = I
            dm%V_src = owned_integral(cpl, dm, I)
            dm%added = .true.
        end associate
    end subroutine globe_ice_gather_add

    subroutine globe_ice_gather_end(gi, cpl)
        ! The host terms on the target grid (background 0: the weighted fields
        ! per unit target area), the checks and the log of section 4.5, the
        ! diagnostics in gi%diag.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl

        character(len=*), parameter :: proc = "kryos_globe_ice::globe_ice_gather_end"
        real(dp), allocatable :: zero(:,:), tmp(:,:), w(:,:), wc(:,:), wf(:,:), ext_d(:,:)
        real(dp) :: v_src, v_tgt, v_ext_bg, v_ext_tgt, area_ocn
        integer :: k, d, n
        character(len=:), allocatable :: line

        if (.not. gi%in_gather) call fesm_error(proc, "no gather is open", "globe_ice_gather_begin first")
        call all_added(gi, gi%dom(:)%added, proc)

        allocate(zero(gi%nxg, gi%nyg), w(gi%nxg, gi%nyg), wc(gi%nxg, gi%nyg), wf(gi%nxg, gi%nyg))
        zero = 0.0_dp
        call coupler_merge_end(gi%acc_I, cpl, zero, tmp, frac_tot=w,  overlap_tol=gi%overlap_tol)
        gi%host_Ig = tmp
        call coupler_merge_end(gi%acc_C, cpl, zero, tmp, frac_tot=wc, overlap_tol=gi%overlap_tol)
        gi%host_C = tmp
        call coupler_merge_end(gi%acc_F, cpl, zero, tmp, frac_tot=wf, overlap_tol=gi%overlap_tol)
        gi%host_F = tmp

        ! the coverage of the host region is fixed with the footprint (section 4.2)
        if (maxval(abs(w - gi%host_w)) > W_TOL) then
            call fesm_error(proc, "the coverage of the host region differs from the reference (max "// &
                            to_str(maxval(abs(w - gi%host_w)))//")", "the footprints and grids are fixed")
        end if
        call check_host(gi, gi%host_Ig, gi%host_C, gi%host_w, "the gather at t = "//to_str(gi%t_load), proc)
        if (.not. all(ieee_is_finite(gi%host_F))) then
            call fesm_error(proc, "host_F is not finite at t = "//to_str(gi%t_load), "")
        end if
        if (.not. all(ieee_is_finite(gi%H_ice_ext)) .or. any(gi%H_ice_ext < 0.0_dp)) then
            call fesm_error(proc, "H_ice_ext is not finite or negative at t = "//to_str(gi%t_load), "")
        end if

        ! conservation: the merged grounded ice against the sources, the external ice against its source grid
        n = size(gi%ord)
        v_src = 0.0_dp
        do k = 1, n
            v_src = v_src + gi%dom(gi%ord(k))%V_src
        end do
        v_tgt = coupler_integral(cpl, trim(gi%tgt), gi%host_Ig)
        gi%diag%cons_I = abs(v_tgt - v_src)/max(abs(v_src), tiny(1.0_dp))
        v_ext_bg  = coupler_integral(cpl, trim(gi%bg),  gi%H_ext)
        v_ext_tgt = coupler_integral(cpl, trim(gi%tgt), gi%H_ice_ext)
        gi%diag%cons_ext = abs(v_ext_tgt - v_ext_bg)/max(abs(v_ext_bg), tiny(1.0_dp))
        if (gi%diag%cons_I > gi%cons_tol .or. gi%diag%cons_ext > gi%cons_tol) then
            call fesm_error(proc, "the merge does not conserve the ice at t = "//to_str(gi%t_load), &
                            "grounded ice: relative difference "//to_str(gi%diag%cons_I)//", external ice: "// &
                            to_str(gi%diag%cons_ext)//", cons_tol = "//to_str(gi%cons_tol))
        end if

        ! per domain and external: grounded volume, mass and sea-level equivalent against the reference
        if (allocated(gi%diag%name)) deallocate(gi%diag%name, gi%diag%V_g, gi%diag%dV_g, gi%diag%dM, gi%diag%esl)
        allocate(gi%diag%name(n), gi%diag%V_g(n), gi%diag%dV_g(n), gi%diag%dM(n), gi%diag%esl(n))
        gi%diag%nd = n
        gi%diag%t  = gi%t_load
        area_ocn = gi%area_ocn
        line = "t = "//to_str(gi%t_load)//" yr, sea-level equivalent of the ice change [m]:"
        do k = 1, n
            d = gi%ord(k)
            gi%diag%name(k) = gi%dom(d)%name
            gi%diag%V_g(k)  = gi%dom(d)%V_src
            gi%diag%dV_g(k) = gi%dom(d)%V_src - gi%dom(d)%V_ref
            gi%diag%dM(k)   = gi%diag%dV_g(k)*gi%rho_ice*KG_TO_GT
            gi%diag%esl(k)  = -(gi%rho_ice/gi%rho_sw)*gi%diag%dV_g(k)/area_ocn
            line = line//" "//trim(gi%dom(d)%name)//" "//to_str(gi%diag%esl(k))
        end do
        ext_d = max(0.0_dp, 1.0_dp - gi%host_w)*(gi%H_ice_ext - gi%H_ice_ref)
        gi%diag%ext_dV  = coupler_integral(cpl, trim(gi%tgt), ext_d)
        gi%diag%ext_dM  = gi%diag%ext_dV*gi%rho_ice*KG_TO_GT
        gi%diag%ext_esl = -(gi%rho_ice/gi%rho_sw)*gi%diag%ext_dV/area_ocn
        line = line//", external "//to_str(gi%diag%ext_esl)//"; clipped external ice "//to_str(gi%diag%clip_Gt)// &
               " Gt; merge closure "//to_str(max(gi%diag%cons_I, gi%diag%cons_ext))
        call log_msg(LOG_INFO, "kryos_globe_ice", line)

        gi%in_gather = .false.
    end subroutine globe_ice_gather_end

    subroutine globe_ice_end(gi)
        type(globe_ice_t), intent(inout) :: gi
        if (gi%ready) then
            call datasrc_end(gi%ds_x)
            if (gi%has_ds_ext) call datasrc_end(gi%ds_ext)
        end if
        gi%ready = .false.
    end subroutine globe_ice_end

    ! =========================================================================
    ! The one procedure of the reference and the step
    ! =========================================================================

    subroutine host_arrays(dm, ice, ext_i, H, z_bed, z_sl, H_tot, f_g, f_in)
        ! H_tot = m H + (1 - m) ext_i; f_g = m f_in + (1 - m) flot(H_tot, z_bed, z_sl)
        ! (the reference has no f_in: flot everywhere, as the same rule on the footprint).
        type(gi_domain_t), intent(in)  :: dm
        type(ice_slot_t),  intent(in)  :: ice
        real(dp),          intent(in)  :: ext_i(:,:), H(:,:), z_bed(:,:), z_sl(:,:)
        real(dp), allocatable, intent(out) :: H_tot(:,:), f_g(:,:)
        real(dp), intent(in), optional :: f_in(:,:)
        real(dp), allocatable :: flot(:,:)

        allocate(flot(dm%nx, dm%ny))
        H_tot = dm%footprint*H + (1.0_dp - dm%footprint)*ext_i
        call ice_flotation(ice, H_tot, z_bed, z_sl, flot)
        if (present(f_in)) then
            f_g = dm%footprint*f_in + (1.0_dp - dm%footprint)*flot
        else
            f_g = flot
        end if
    end subroutine host_arrays

    ! =========================================================================
    ! Remaps
    ! =========================================================================

    subroutine ext_to_tgt(cpl, bg, tgt, H_bg, H_tgt)
        ! Ice on the background grid to the target grid: an amount per target
        ! area (conserves the volume in the model areas), the whole-cell mean
        ! that the 1 - w part takes (D17).
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: bg, tgt
        real(dp),            intent(in)    :: H_bg(:,:)
        real(dp),            intent(out)   :: H_tgt(:,:)
        real(dp), allocatable :: tmp(:,:)
        call remap(cpl, H_bg, bg, tmp, tgt, method="con", norm="dst")
        H_tgt = tmp
    end subroutine ext_to_tgt

    subroutine ext_to_ice(cpl, bg, grid, H_bg, H_ice)
        ! Ice on the background grid to the ice grid of a domain (amount per true area).
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: bg, grid
        real(dp),            intent(in)    :: H_bg(:,:)
        real(dp), allocatable, intent(out) :: H_ice(:,:)
        real(dp), allocatable :: tmp(:,:)
        call remap(cpl, H_bg, bg, tmp, grid, method="con", norm="dst")
        H_ice = tmp
    end subroutine ext_to_ice

    ! =========================================================================
    ! Checks
    ! =========================================================================

    subroutine require_ready(gi, proc)
        type(globe_ice_t), intent(in) :: gi
        character(len=*),  intent(in) :: proc
        if (.not. gi%ready) call fesm_error(proc, "globe_ice_t is not initialised", "globe_ice_init first")
    end subroutine require_ready

    subroutine all_added(gi, flags, proc)
        type(globe_ice_t), intent(in) :: gi
        logical,           intent(in) :: flags(:)
        character(len=*),  intent(in) :: proc
        integer :: k
        do k = 1, size(gi%ord)
            if (.not. flags(gi%ord(k))) then
                call fesm_error(proc, "domain `"//trim(gi%dom(gi%ord(k))%name)//"` was not added", &
                                "every registered domain takes part in the reference and in every gather")
            end if
        end do
    end subroutine all_added

    subroutine check_shape(a, dm, what, proc)
        real(dp),          intent(in) :: a(:,:)
        type(gi_domain_t), intent(in) :: dm
        character(len=*),  intent(in) :: what, proc
        if (size(a,1) /= dm%nx .or. size(a,2) /= dm%ny) then
            call fesm_error(proc, "domain `"//trim(dm%name)//"`: "//trim(what)//" is "//to_str(size(a,1))//" x "// &
                            to_str(size(a,2))//", the ice grid "//trim(dm%grid)//" is "//to_str(dm%nx)//" x "//to_str(dm%ny), "")
        end if
    end subroutine check_shape

    subroutine check_finite(a, dm, what, proc)
        real(dp),          intent(in) :: a(:,:)
        type(gi_domain_t), intent(in) :: dm
        character(len=*),  intent(in) :: what, proc
        if (.not. all(ieee_is_finite(a))) then
            call fesm_error(proc, "domain `"//trim(dm%name)//"`: "//trim(what)//" is not finite", &
                            to_str(count(.not. ieee_is_finite(a)))//" cells")
        end if
    end subroutine check_finite

    subroutine check_bg(a, what, proc)
        real(dp),         intent(in) :: a(:,:)
        character(len=*), intent(in) :: what, proc
        if (.not. all(ieee_is_finite(a)) .or. any(a == mv_dp)) then
            call fesm_error(proc, trim(what)//" on the background grid is not finite or holds missing values", &
                            to_str(count(.not. ieee_is_finite(a) .or. a == mv_dp))//" cells")
        end if
    end subroutine check_bg

    subroutine check_host(gi, Ig, C, w, what, proc)
        ! finite; 0 <= C <= w; Ig >= 0; w <= 1 + tolerance (the restart of VILMA2 repeats it)
        type(globe_ice_t), intent(in) :: gi
        real(dp),          intent(in) :: Ig(:,:), C(:,:), w(:,:)
        character(len=*),  intent(in) :: what, proc
        if (.not. all(ieee_is_finite(Ig)) .or. .not. all(ieee_is_finite(C)) .or. .not. all(ieee_is_finite(w))) then
            call fesm_error(proc, "host_Ig, host_C or host_w is not finite in "//trim(what), "group ["//trim(gi%group)//"]")
        end if
        if (any(Ig < 0.0_dp)) then
            call fesm_error(proc, "host_Ig is negative in "//trim(what), "min "//to_str(minval(Ig)))
        end if
        if (any(C < -BOUND_TOL) .or. any(C > w + BOUND_TOL)) then
            call fesm_error(proc, "host_C is outside [0, host_w] in "//trim(what), &
                            "min "//to_str(minval(C))//", max (C - w) "//to_str(maxval(C - w)))
        end if
        if (any(w > 1.0_dp + gi%overlap_tol)) then
            call fesm_error(proc, "host_w exceeds 1 in "//trim(what), "max "//to_str(maxval(w)))
        end if
    end subroutine check_host

end module kryos_globe_ice
