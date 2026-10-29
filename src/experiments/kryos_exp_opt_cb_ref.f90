module kryos_exp_opt_cb_ref
    ! Experiment `opt_cb_ref` (phase 2 WP2): optimisation of Yelmo's basal
    ! friction cb_ref towards present-day observations (Lipscomb et al. 2021;
    ! yelmox step_optimize, kryos_coupling.f90:50-86; domain_opt_init and
    ! domain_opt_init_cb_ref, kryos.f90:516-542, kryos_startup.f90:197-216).
    ! Built with yelmo=1; otherwise kryos_exp_opt_cb_ref_stub.f90 provides it.
    !
    ! init: Yelmo takes cb_ref from outside (ydyn till_method = -1) and the
    !   friction bounds cf_min, cf_max are seeded from ytill (also on a restore,
    !   where cb_ref comes from Yelmo's restart).
    ! HOOK_COLD_PRE_STATE: the cold-start guess: cf_init > 0 a constant, else
    !   the till friction of the bed (calc_cb_ref).
    ! HOOK_STEP_BEGIN, in the window [cf_time_init, cf_time_end] of the elapsed
    !   time: optimize_cb_ref on the Yelmo state of the previous step.
    ! method = "none" keeps the first two (an external, constant cb_ref) and
    ! never optimises. Stateless (cb_ref is in Yelmo's restart). Reaches Yelmo
    ! by `select type` on the domain's ice (yelmo_ice_t).

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error
    use basal_dragging,  only: calc_cb_ref
    use ice_optimization, only: optimize_cb_ref
    use kryos_config,    only: config_read
    use kryos_domain,    only: domain_t
    use kryos_experiment, only: experiment_t, exp_ctx_t, hook_time_t, HOOK_STEP_BEGIN, HOOK_COLD_PRE_STATE, &
                                experiment_require_ice_every_step
    use kryos_ice_yelmo, only: yelmo_ice_t

    implicit none

    private
    public :: opt_cb_ref_exp_t

    type, extends(experiment_t) :: opt_cb_ref_exp_t
        character(len=56) :: method      = "L21"      ! none | L21
        real(sp) :: cf_time_init = 0.0_sp, cf_time_end = 0.0_sp   ! [yr] window of the elapsed time
        real(sp) :: cf_init = 0.0_sp                  ! cold-start cb_ref everywhere (<= 0: the till friction)
        real(sp) :: tau_c = 0.0_sp, H0 = 0.0_sp, sigma_err = 0.0_sp, sigma_vel = 0.0_sp
        character(len=56) :: fill_method = "cf_min"
        real(sp), allocatable :: cf_min(:,:), cf_max(:,:)   ! bounds, from ytill (init)
    contains
        procedure :: init => optcf_init
        procedure :: hook => optcf_hook
    end type opt_cb_ref_exp_t

contains

    subroutine optcf_init(self, ctx, dom)
        class(opt_cb_ref_exp_t), intent(inout) :: self
        type(exp_ctx_t),         intent(in)    :: ctx
        type(domain_t),          intent(inout) :: dom

        character(len=*), parameter :: tmpl = "experiment", extra = "experiment_opt_cb_ref"
        character(len=*), parameter :: proc = "kryos_exp_opt_cb_ref::init"
        real(dp) :: x

        call experiment_require_ice_every_step(ctx, "opt_cb_ref")

        call config_read(ctx%cfg, ctx%instance, "method", self%method, tmpl, extra)
        select case (trim(self%method))
        case ("none", "L21")
        case default
            call fesm_error(proc, "experiment `"//trim(self%name)//"`: method = '"//trim(self%method)// &
                            "' is not recognised", "use none or L21")
        end select
        call config_read(ctx%cfg, ctx%instance, "cf_time_init", x, tmpl, extra);  self%cf_time_init = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "cf_time_end",  x, tmpl, extra);  self%cf_time_end  = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "cf_init",      x, tmpl, extra);  self%cf_init      = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tau_c",        x, tmpl, extra);  self%tau_c        = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "H0",           x, tmpl, extra);  self%H0           = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "sigma_err",    x, tmpl, extra);  self%sigma_err    = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "sigma_vel",    x, tmpl, extra);  self%sigma_vel    = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "fill_method", self%fill_method, tmpl, extra)

        if (.not. allocated(dom%ice%c)) then
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs an ice instance", &
                            "domain &"//trim(dom%name)//" has none")
        end if
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            if (trim(self%method) /= "none" .and. .not. c%ylmo%dta%par%pd_vel_load) then
                call fesm_error(proc, "experiment `"//trim(self%name)//"` (L21) needs the present-day velocity", &
                                "ice instance `"//trim(c%name)//"`: set pd_vel_load = True in Yelmo's group &yelmo_data")
            end if
            allocate(self%cf_min(c%ylmo%grd%G%nx, c%ylmo%grd%G%ny), self%cf_max(c%ylmo%grd%G%nx, c%ylmo%grd%G%ny))
            self%cf_min = c%ylmo%dyn%par%till_cf_min
            self%cf_max = c%ylmo%dyn%par%till_cf_ref
            c%ylmo%dyn%par%till_method = -1               ! cb_ref is set from outside from here on
        class default
            call fesm_error(proc, "experiment `"//trim(self%name)//"` needs ice backend yelmo", &
                            "domain &"//trim(dom%name)//": the ice instance `"//trim(dom%ice%c%name)//"` is backend `"// &
                            trim(dom%ice%c%backend)//"`")
        end select

        self%wants(HOOK_COLD_PRE_STATE) = .true.
        if (trim(self%method) /= "none") self%wants(HOOK_STEP_BEGIN) = .true.
    end subroutine optcf_init

    subroutine optcf_hook(self, dom, id, st)
        class(opt_cb_ref_exp_t), intent(inout) :: self
        type(domain_t),          intent(inout) :: dom
        integer,                 intent(in)    :: id
        type(hook_time_t),       intent(in)    :: st

        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            select case (id)
            case (HOOK_COLD_PRE_STATE)
                if (self%cf_init > 0.0_sp) then
                    c%ylmo%dyn%now%cb_ref = self%cf_init
                else
                    associate (par => c%ylmo%dyn%par, bnd => c%ylmo%bnd)
                        call calc_cb_ref(c%ylmo%dyn%now%cb_ref, bnd%z_bed, bnd%z_bed_sd, bnd%z_sl, bnd%H_sed, &
                                         par%till_f_sed, par%till_sed_min, par%till_sed_max, par%till_cf_ref, &
                                         par%till_cf_min, par%till_z0, par%till_z1, par%till_n_sd, &
                                         par%till_scale_zb, par%till_scale_sed)
                    end associate
                end if
            case (HOOK_STEP_BEGIN)
                if (st%t_elapsed < self%cf_time_init .or. st%t_elapsed > self%cf_time_end) return
                call optimize_cb_ref(c%ylmo%dyn%now%cb_ref, c%ylmo%tpo%now%H_ice, c%ylmo%tpo%now%dHidt, &
                                     c%ylmo%bnd%z_bed, c%ylmo%bnd%z_sl, c%ylmo%dyn%now%ux_s, c%ylmo%dyn%now%uy_s, &
                                     c%ylmo%dta%pd%H_ice, c%ylmo%dta%pd%uxy_s, c%ylmo%dta%pd%H_grnd, &
                                     self%cf_min, self%cf_max, c%ylmo%tpo%par%dx, &
                                     self%sigma_err, self%sigma_vel, self%tau_c, self%H0, &
                                     dt=st%dtt, fill_method=self%fill_method, fill_dist=self%sigma_err, &
                                     cb_tgt=c%ylmo%dyn%now%cb_tgt)
            end select
        end select
    end subroutine optcf_hook

end module kryos_exp_opt_cb_ref
