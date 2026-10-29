module kryos_exp_opt_tf
    ! Experiment `opt_tf` (phase 2 WP3): optimisation of the shelf component's
    ! thermal-forcing correction tf_corr towards the present-day ice (Lipscomb et
    ! al. 2021; yelmox step_optimize, kryos_coupling.f90:88-115). Built with
    ! yelmo=1 (the calibration routines are Yelmo's ice_optimization); otherwise
    ! kryos_exp_opt_tf_stub.f90 provides this module.
    !
    ! Hook STEP_BEGIN, in the window [tf_time_init, tf_time_end] of the elapsed
    ! time: the ice's exports of the previous step (H_ice, H_grnd, dHidt and the
    ! static H_ice_pd, and for method L21-points H_grnd_pd) are remapped ice grid
    ! -> shelf grid by the domain coupler (bilinear where the shelf grid is finer,
    ! else conservative; the map is primed at init), then optimize_tf_corr_basin
    ! (L21, one correction per basin) or optimize_tf_corr (L21-points) updates
    ! shelf_t%corr%tf_corr in place, on the basins and the grid spacing of the
    ! shelf component. No `select type`: only the contract fields of the ice and
    ! the shelf component are used. Stateless: tf_corr is in the shelf component's
    ! restart (<shelf>_corr_restart.nc).

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: field_state_t, fields_has, fields_get, fields_handle
    use nml,             only: nml_read
    use coupler,         only: coupler_prime, remap
    use ice_optimization, only: optimize_tf_corr, optimize_tf_corr_basin
    use kryos_config,    only: config_read
    use kryos_domain,    only: domain_t, C_ICE, C_SHELF
    use kryos_experiment, only: experiment_t, exp_ctx_t, hook_time_t, HOOK_STEP_BEGIN, &
                                experiment_require_ice_every_step

    implicit none

    private
    public :: opt_tf_exp_t

    integer, parameter :: N_BASINS = 100             ! length of Yelmo's tf_basins list

    type, extends(experiment_t) :: opt_tf_exp_t
        character(len=56) :: method = "L21-points"    ! L21 | L21-points
        real(sp) :: tf_time_init = 0.0_sp, tf_time_end = 0.0_sp   ! [yr] window of the elapsed time
        real(sp) :: H_grnd_lim = 0.0_sp, tf_sigma = 0.0_sp, tau_m = 0.0_sp, m_temp = 0.0_sp
        real(sp) :: tf_min = 0.0_sp, tf_max = 0.0_sp
        integer  :: tf_basins(N_BASINS) = 0          ! L21: the basins to optimise (first < 0: all)
        logical  :: basin_fill = .false.             ! L21-points
        character(len=32) :: remap_method = "con"    ! ice grid -> shelf grid
        character(len=64) :: grid_ice = "", grid_shelf = ""
    contains
        procedure :: init => opttf_init
        procedure :: hook => opttf_hook
    end type opt_tf_exp_t

contains

    subroutine opttf_init(self, ctx, dom)
        class(opt_tf_exp_t), intent(inout) :: self
        type(exp_ctx_t),     intent(in)    :: ctx
        type(domain_t),      intent(inout) :: dom

        character(len=*), parameter :: tmpl = "experiment", extra = "experiment_opt_tf"
        character(len=*), parameter :: proc = "kryos_exp_opt_tf::init"
        real(dp) :: x
        real(sp) :: dx_ice, dx_shelf

        call experiment_require_ice_every_step(ctx, "opt_tf")

        call config_read(ctx%cfg, ctx%instance, "method", self%method, tmpl, extra)
        select case (trim(self%method))
        case ("L21", "L21-points")
        case default
            call fesm_error(proc, "experiment `"//trim(self%name)//"`: method = '"//trim(self%method)// &
                            "' is not recognised", "use L21 or L21-points")
        end select
        call config_read(ctx%cfg, ctx%instance, "tf_time_init", x, tmpl, extra);  self%tf_time_init = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tf_time_end",  x, tmpl, extra);  self%tf_time_end  = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "H_grnd_lim",   x, tmpl, extra);  self%H_grnd_lim   = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tf_sigma",     x, tmpl, extra);  self%tf_sigma     = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tau_m",        x, tmpl, extra);  self%tau_m        = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "m_temp",       x, tmpl, extra);  self%m_temp       = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tf_min",       x, tmpl, extra);  self%tf_min       = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tf_max",       x, tmpl, extra);  self%tf_max       = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "basin_fill", self%basin_fill, tmpl, extra)
        self%tf_basins = 0                                 ! (the entries the list does not give)
        call nml_read(ctx%par_file, trim(ctx%instance), "tf_basins", self%tf_basins, init=.true., &
                      defaults_file=trim(ctx%cfg%defaults_file), defaults_group=extra)

        ! the instances and the fields the calibration reads
        if (.not. allocated(dom%ice%c)) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs an ice instance", &
                            "domain &"//trim(dom%name)//" has none")
        end if
        if (.not. allocated(dom%shelf%c)) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs a shelf instance", &
                            "domain &"//trim(dom%name)//" has none")
        end if
        call need(dom%ice%c%exports, "H_ice");  call need(dom%ice%c%exports, "H_grnd")
        call need(dom%ice%c%exports, "dHidt");  call need(dom%ice%c%exports, "H_ice_pd")
        if (trim(self%method) == "L21-points") call need(dom%ice%c%exports, "H_grnd_pd")

        ! ice grid -> shelf grid: bilinear to refine (shelf finer), conservative to coarsen
        self%grid_ice   = dom%grid(C_ICE)
        self%grid_shelf = dom%grid(C_SHELF)
        dx_ice   = real(dom%ice%c%xkm(2) - dom%ice%c%xkm(1), sp)
        dx_shelf = dom%shelf%c%corr%xc(2) - dom%shelf%c%corr%xc(1)
        self%remap_method = "con"
        if (dx_shelf < dx_ice) self%remap_method = "bilin"
        call coupler_prime(dom%cpl, trim(self%grid_ice), trim(self%grid_shelf), method=trim(self%remap_method))

        self%wants(HOOK_STEP_BEGIN) = .true.

    contains

        subroutine need(exports, name)
            type(field_state_t), intent(in) :: exports
            character(len=*),    intent(in) :: name
            if (.not. fields_has(exports, name)) then
                call fesm_error(proc, "experiment `"//trim(self%name)//"` needs the ice export `"//name//"`", &
                                "ice instance `"//trim(dom%ice%c%name)//"` (backend "//trim(dom%ice%c%backend)//")")
            end if
        end subroutine need

    end subroutine opttf_init

    subroutine opttf_hook(self, dom, id, st)
        class(opt_tf_exp_t), intent(inout) :: self
        type(domain_t),      intent(inout) :: dom
        integer,             intent(in)    :: id
        type(hook_time_t),   intent(in)    :: st

        real(sp), allocatable :: H_ice_m(:,:), H_grnd_m(:,:), dHidt_m(:,:), H_ice_pd_m(:,:), H_grnd_pd_m(:,:)

        associate (unused => id); end associate
        if (st%t_elapsed < self%tf_time_init .or. st%t_elapsed > self%tf_time_end) return

        call ice_field(self, dom, "H_ice",    H_ice_m)
        call ice_field(self, dom, "H_grnd",   H_grnd_m)
        call ice_field(self, dom, "dHidt",    dHidt_m)
        call ice_field(self, dom, "H_ice_pd", H_ice_pd_m)

        associate (corr => dom%shelf%c%corr)
            select case (trim(self%method))
            case ("L21")
                call optimize_tf_corr_basin(corr%tf_corr, H_ice_m, H_grnd_m, dHidt_m, H_ice_pd_m, corr%basins, &
                                            self%H_grnd_lim, self%tau_m, self%m_temp, self%tf_min, self%tf_max, &
                                            self%tf_basins, dt=st%dtt)
            case ("L21-points")
                call ice_field(self, dom, "H_grnd_pd", H_grnd_pd_m)
                call optimize_tf_corr(corr%tf_corr, H_ice_m, H_grnd_m, dHidt_m, H_ice_pd_m, H_grnd_pd_m, &
                                      self%H_grnd_lim, corr%basins, self%basin_fill, self%tau_m, self%m_temp, &
                                      self%tf_min, self%tf_max, corr%dx, sigma=self%tf_sigma, dt=st%dtt)
            end select
        end associate
    end subroutine opttf_hook

    subroutine ice_field(self, dom, name, v_m)
        ! The ice's export `name` (narrowed to sp by fields_get, as the ice adapter widened an sp
        ! field) remapped to the shelf grid; a copy where the grids are the same.
        class(opt_tf_exp_t),   intent(in)    :: self
        type(domain_t),        intent(inout) :: dom
        character(len=*),      intent(in)    :: name
        real(sp), allocatable, intent(inout) :: v_m(:,:)

        real(sp), allocatable :: v(:,:)

        allocate(v(dom%ice%c%nx, dom%ice%c%ny))
        call fields_get(dom%ice%c%exports, fields_handle(dom%ice%c%exports, name), v)
        call remap(dom%cpl, v, trim(self%grid_ice), v_m, trim(self%grid_shelf), method=trim(self%remap_method))
    end subroutine ice_field

end module kryos_exp_opt_tf
