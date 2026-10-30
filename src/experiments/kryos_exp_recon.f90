module kryos_exp_recon
    ! Experiment `recon` (phase 2 WP4): the cold start from an ice reconstruction (yelmox
    ! domain_init_recon, kryos_startup.f90:289-349; domain_ctl_load, kryos.f90:629-645). Built
    ! with yelmo=1; otherwise kryos_exp_recon_stub.f90 provides this module.
    !
    ! Hook COLD_POST_STATE (after Yelmo's state is initialised, before PASS_SPINUP):
    !   1. the reconstruction's slice nearest recon_time of recon_var is Yelmo's reference ice
    !      thickness (bnd%H_ice_ref);
    !   2. the initial ice is that field where the bed is above recon_zb_min (default -500 m) and the
    !      cell is in the recon_regions ("all": every cell, no regions file read), smoothed (Gaussian, sigma =
    !      recon_smooth_sigma dx, default 2; <= 0: not smoothed), the topography initialised from it
    !      (yelmo_init_topo, load_topo = .false.);
    !   3. a one-year fixed-topography equilibration, thin floating ice (mask_bed = 5, H < recon_H_float_min,
    !      default 50 m; <= 0: kept) removed, another year;
    !   4. the reference ice thickness is the clean ice thickness (it drives the smb masks);
    !   5. the ice's INSTANT exports are published (yelmo_ice_publish), the domain's boundary
    !      refreshed (domain_refresh_boundary: hub, climate, surface, surface -> ice, prepare);
    !   6. 200 yr of equilibration with beta_min raised to 100 (stabilises the dynamic fields).
    ! recon_ref (the reconstruction as the reference ice only) is not ported. Stateless: the state
    ! is Yelmo's. Reaches Yelmo by `select type` on the domain's ice (yelmo_ice_t).

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error, to_str
    use ncio,            only: nc_read, nc_size, nc_read_attr
    use nml,             only: nml_replace
    use regions,         only: regions_class, regions_select
    use coords,          only: grid_class, grid_cdo_read_desc
    use yelmo,           only: yelmo_update_equil, yelmo_init_topo
    use yelmo_tools,     only: smooth_gauss_2D
    use kryos_component, only: NAME_LEN, PATH_LEN
    use kryos_config,    only: config_read
    use kryos_domain,    only: domain_t, C_ICE, domain_regions_load, domain_refresh_boundary
    use kryos_experiment, only: experiment_t, exp_ctx_t, hook_time_t, HOOK_COLD_POST_STATE
    use kryos_ice_yelmo, only: yelmo_ice_t, yelmo_ice_publish

    implicit none

    private
    public :: recon_exp_t, recon_set_initial_ice, recon_remove_thin_float

    integer, parameter :: LEN_EXPR = 1000

    type, extends(experiment_t) :: recon_exp_t
        character(len=PATH_LEN) :: path = ""            ! the reconstruction file ({domain}, {grid_name} replaced)
        character(len=NAME_LEN) :: var = ""             ! its ice-thickness variable
        real(dp) :: time = 0.0_dp                       ! [yr] model time of the state (nearest slice)
        character(len=LEN_EXPR) :: regions = ""         ! where its ice is imposed (a regions selection, or "all")
        real(dp) :: zb_min = -500.0_dp                  ! [m] its ice is imposed where z_bed > zb_min
        real(dp) :: smooth_sigma = 2.0_dp               ! [grid cells] Gaussian smoothing of the initial ice (<= 0: none)
        real(dp) :: H_float_min = 50.0_dp               ! [m] floating ice thinner than this is removed after year 1 (<= 0: kept)
        logical, allocatable :: mask_recon(:,:)         ! the same on the ice grid
        character(len=PATH_LEN) :: par_file = ""        ! (yelmo_init_topo reads its group from the par file)
    contains
        procedure :: init => recon_init
        procedure :: hook => recon_hook
    end type recon_exp_t

contains

    subroutine recon_init(self, ctx, dom)
        class(recon_exp_t), intent(inout) :: self
        type(exp_ctx_t),    intent(in)    :: ctx
        type(domain_t),     intent(inout) :: dom

        character(len=*), parameter :: tmpl = "experiment", extra = "experiment_recon"
        character(len=*), parameter :: proc = "kryos_exp_recon::init"
        type(regions_class) :: reg
        type(grid_class) :: grid
        logical :: found

        call config_read(ctx%cfg, ctx%instance, "recon_path",    self%path,    tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_var",     self%var,     tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_time",    self%time,    tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_regions", self%regions, tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_zb_min",  self%zb_min,  tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_smooth_sigma", self%smooth_sigma, tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "recon_H_float_min",  self%H_float_min,  tmpl, extra)
        if (len_trim(self%path) == 0 .or. len_trim(self%var) == 0 .or. len_trim(self%regions) == 0) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs recon_path, recon_var and recon_regions", &
                            "group &"//trim(ctx%instance)//" of "//trim(ctx%par_file))
        end if
        if (self%smooth_sigma > 0.0_dp .and. self%smooth_sigma < 1.0_dp) then       ! (smooth_gauss_2D needs f_sigma >= 1)
            call fesm_error(proc, "experiment `"//trim(self%name)//"`: recon_smooth_sigma must be 0 (no smoothing) or >= 1 grid cell", &
                            "recon_smooth_sigma = "//to_str(self%smooth_sigma)//", group &"//trim(ctx%instance)//" of "//trim(ctx%par_file))
        end if
        if (.not. allocated(dom%ice%c)) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs an ice instance", &
                            "domain &"//trim(dom%name)//" has none")
        end if
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            continue
        class default
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs ice backend yelmo", &
                            "domain &"//trim(dom%name)//": the ice instance `"//trim(dom%ice%c%name)//"` is backend `"// &
                            trim(dom%ice%c%backend)//"`")
        end select

        call nml_replace(self%path, "{domain}",    trim(dom%name))
        call nml_replace(self%path, "{grid_name}", trim(dom%grid(C_ICE)))
        inquire(file=trim(self%path), exist=found)
        if (.not. found) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"`: the reconstruction file does not exist", &
                            "recon_path = "//trim(self%path))
        end if

        call grid_cdo_read_desc(grid, trim(dom%grid(C_ICE)), trim(ctx%cfg%map_fldr))
        if (trim(adjustl(self%regions)) == "all") then
            allocate(self%mask_recon(grid%G%nx, grid%G%ny))     ! the whole domain: no regions file
            self%mask_recon = .true.
        else
            call domain_regions_load(reg, dom, ctx%cfg, grid)
            self%mask_recon = regions_select(reg, trim(self%regions))
        end if
        self%par_file   = ctx%par_file

        self%wants(HOOK_COLD_POST_STATE) = .true.
    end subroutine recon_init

    subroutine recon_hook(self, dom, id, st)
        class(recon_exp_t), intent(inout) :: self
        type(domain_t),     intent(inout) :: dom
        integer,            intent(in)    :: id
        type(hook_time_t),  intent(in)    :: st

        integer  :: nx, ny, k
        real(sp) :: beta_min_save

        associate (unused => id); end associate
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            nx = c%ylmo%tpo%par%nx
            ny = c%ylmo%tpo%par%ny

            ! the reconstruction at recon_time is the reference ice thickness
            k = recon_time_index(self%path, self%time)
            call nc_read(trim(self%path), trim(self%var), c%ylmo%bnd%H_ice_ref, start=[1,1,k], count=[nx,ny,1])

            ! initial ice: the reconstruction on the regions, smoothed; the topography from it
            call recon_set_initial_ice(c%ylmo%tpo%now%H_ice, c%ylmo%bnd%H_ice_ref, c%ylmo%bnd%z_bed, self%mask_recon, &
                                       self%zb_min, self%smooth_sigma, real(c%ylmo%grd%G%dx, sp))
            call yelmo_init_topo(c%ylmo, trim(self%par_file), c%ylmo%par%nml_init_topo, st%t1, load_topo=.false.)

            ! update the surface topography (fixed H), remove thin floating ice, update again
            call yelmo_update_equil(c%ylmo, st%t1, time_tot=1.0_sp, dt=1.0_sp, topo_fixed=.true.)
            call recon_remove_thin_float(c%ylmo%tpo%now%H_ice, c%ylmo%tpo%now%mask_bed, self%H_float_min)
            call yelmo_update_equil(c%ylmo, st%t1, time_tot=1.0_sp, dt=1.0_sp, topo_fixed=.true.)

            ! the clean thickness is the reference state (drives the smb masks)
            c%ylmo%bnd%H_ice_ref = c%ylmo%tpo%now%H_ice

            ! the new geometry and reference to the exports, the boundary fields on it
            call yelmo_ice_publish(c, st%t1)
            call domain_refresh_boundary(dom, st%t1, st%clk)

            ! stabilise the dynamic fields with a raised beta_min
            beta_min_save = c%ylmo%dyn%par%beta_min
            c%ylmo%dyn%par%beta_min = 100.0_sp
            call yelmo_update_equil(c%ylmo, st%t1, time_tot=2.0e2_sp, dt=5.0_sp, topo_fixed=.false.)
            c%ylmo%dyn%par%beta_min = beta_min_save
        end select
    end subroutine recon_hook

    subroutine recon_set_initial_ice(H_ice, H_ref, z_bed, mask, zb_min, sigma, dx)
        ! The initial ice: the reconstruction `H_ref` where the bed is above zb_min [m] and `mask` holds (H_ice
        ! elsewhere as it is), then smoothed with a Gaussian of `sigma` grid cells of size dx [m]; sigma <= 0: no
        ! smoothing (the field is the reconstruction's exactly in the selected cells).
        real(sp), intent(inout) :: H_ice(:,:)
        real(sp), intent(in)    :: H_ref(:,:), z_bed(:,:), dx
        logical,  intent(in)    :: mask(:,:)
        real(dp), intent(in)    :: zb_min, sigma

        where (z_bed > real(zb_min, sp) .and. mask) H_ice = H_ref
        if (sigma > 0.0_dp) call smooth_gauss_2D(H_ice, dx=dx, f_sigma=real(sigma, sp))
    end subroutine recon_set_initial_ice

    subroutine recon_remove_thin_float(H_ice, mask_bed, H_min)
        ! Floating ice (mask_bed = 5) thinner than H_min [m] is removed; H_min <= 0: nothing is.
        real(sp), intent(inout) :: H_ice(:,:)
        integer,  intent(in)    :: mask_bed(:,:)
        real(dp), intent(in)    :: H_min

        if (H_min > 0.0_dp) then
            where (mask_bed == 5 .and. H_ice < real(H_min, sp)) H_ice = 0.0_sp
        end if
    end subroutine recon_remove_thin_float

    function recon_time_index(filename, time) result(k)
        ! Index of the slice of a reconstruction file nearest to the model time `time` [yr]. The
        ! time axis is converted to model time from its units: "...ago"/"...BP" count positive
        ! into the past; a leading "k" (kiloyears, kyr, ka) means thousands of years. Stops if
        ! `time` lies outside the axis (yelmox recon_time_index).
        character(len=*), intent(in) :: filename
        real(dp),         intent(in) :: time
        integer :: k

        real(dp), allocatable :: t(:)
        character(len=256) :: units
        integer  :: nt
        real(dp) :: scale

        nt = nc_size(filename, "time")
        allocate(t(nt))
        call nc_read(filename, "time", t)
        units = ""
        call nc_read_attr(filename, "time", "units", units)
        units = adjustl(units)

        scale = 1.0_dp
        if (units(1:1) == "k") scale = 1e3_dp
        if (index(units, "ago") > 0 .or. index(units, "BP") > 0) scale = -scale
        t = t*scale

        if (time < minval(t) .or. time > maxval(t)) then
            call fesm_error("kryos_exp_recon::recon_time_index", "recon_time = "//to_str(time)// &
                            " is outside the time axis of "//trim(filename), "["//to_str(minval(t))//", "// &
                            to_str(maxval(t))//"] yr (units: "//trim(units)//")")
        end if

        k = minloc(abs(t - time), dim=1)
    end function recon_time_index

end module kryos_exp_recon
