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

    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use precision,       only: dp
    use fesm_log,        only: LOG_INFO, LOG_WARN, log_msg, fesm_error, to_str
    use fesm_fields,     only: field_state_t, fields_init, fields_add, fields_lock, FIELD_INSTANT, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE
    use fesm_datasrc,    only: datasrc_t, datasrc_init, datasrc_update, datasrc_get, datasrc_describe, datasrc_end
    use coords,          only: grid_class, mv_dp
    use coupler,         only: coupler_class, coupler_merge_t, coupler_add_grid, coupler_set_area, &
                               coupler_integral, coupler_merge_begin, coupler_merge_add, coupler_merge_end, remap
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
        logical :: in_ref = .false., in_gather = .false.
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

    subroutine globe_ice_add_domain(gi, cpl, name, grid, footprint, H_ice_pd, z_bed_ref)
        ! Register a domain: its ice grid in the coupler with the true spherical
        ! cell areas as model areas (section 4.4: the load per true area, global
        ! integrals spherical; the plane dx*dy is Yelmo's), its static exports
        ! (copied), and its coverage and footprint fraction on the two global
        ! grids. Call before the reference; any order.
        type(globe_ice_t),   intent(inout) :: gi
        type(coupler_class), intent(inout) :: cpl
        character(len=*),    intent(in)    :: name
        type(grid_class),    intent(in)    :: grid
        real(dp),            intent(in)    :: footprint(:,:), H_ice_pd(:,:), z_bed_ref(:,:)

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

    subroutine globe_ice_end(gi)
        type(globe_ice_t), intent(inout) :: gi
        if (gi%ready) then
            call datasrc_end(gi%ds_x)
            if (gi%has_ds_ext) call datasrc_end(gi%ds_ext)
        end if
        gi%ready = .false.
    end subroutine globe_ice_end

    ! =========================================================================
    ! Checks
    ! =========================================================================

    subroutine require_ready(gi, proc)
        type(globe_ice_t), intent(in) :: gi
        character(len=*),  intent(in) :: proc
        if (.not. gi%ready) call fesm_error(proc, "globe_ice_t is not initialised", "globe_ice_init first")
    end subroutine require_ready

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

end module kryos_globe_ice
