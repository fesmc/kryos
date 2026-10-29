module kryos_exp_relax
    ! Experiment `relax` (phase 2 WP2): relaxation of Yelmo's topography
    ! towards the reference, ytopo.topo_rel set from the time elapsed since
    ! the start of the run (yelmox step_relax, kryos_coupling.f90:36-47;
    ! domain_relax_init, kryos.f90:501-514). Built with yelmo=1; otherwise
    ! kryos_exp_relax_stub.f90 provides this module.
    !
    ! Hook STEP_BEGIN: relax_update (yelmo's ice_optimization) on the ytopo
    ! parameters, window time st%t_elapsed. Stateless: a function of time.
    ! Reaches Yelmo by `select type` on the domain's ice (yelmo_ice_t).

    use precision,       only: sp, dp
    use fesm_log,        only: fesm_error
    use ice_optimization, only: relax_params, relax_update
    use kryos_config,    only: config_read
    use kryos_domain,    only: domain_t
    use kryos_experiment, only: experiment_t, exp_ctx_t, hook_time_t, HOOK_STEP_BEGIN, &
                                experiment_require_ice_every_step
    use kryos_ice_yelmo, only: yelmo_ice_t

    implicit none

    private
    public :: relax_exp_t

    type, extends(experiment_t) :: relax_exp_t
        type(relax_params) :: par                   ! topo_rel, tau1, tau2, time1, time2, m
        integer  :: topo_rel0     = 0               ! ytopo.topo_rel and topo_rel_tau of the par file,
        real(sp) :: topo_rel_tau0 = 0.0_sp          ! which apply again after time2
    contains
        procedure :: init => relax_init
        procedure :: hook => relax_hook
    end type relax_exp_t

contains

    subroutine relax_init(self, ctx, dom)
        class(relax_exp_t), intent(inout) :: self
        type(exp_ctx_t),    intent(in)    :: ctx
        type(domain_t),     intent(inout) :: dom

        character(len=*), parameter :: tmpl = "experiment", extra = "experiment_relax"
        real(dp) :: x

        call experiment_require_ice_every_step(ctx, "relax")

        call config_read(ctx%cfg, ctx%instance, "topo_rel", self%par%topo_rel, tmpl, extra)
        call config_read(ctx%cfg, ctx%instance, "tau1",  x, tmpl, extra);  self%par%tau1  = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "tau2",  x, tmpl, extra);  self%par%tau2  = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "time1", x, tmpl, extra);  self%par%time1 = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "time2", x, tmpl, extra);  self%par%time2 = real(x, sp)
        call config_read(ctx%cfg, ctx%instance, "m",     x, tmpl, extra);  self%par%m     = real(x, sp)

        if (.not. allocated(dom%ice%c)) then
            call fesm_error("kryos_exp_relax::init", "experiment `"//trim(self%name)//"` needs an ice instance", &
                            "domain &"//trim(dom%name)//" has none")
        end if
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            self%topo_rel0     = c%ylmo%tpo%par%topo_rel
            self%topo_rel_tau0 = c%ylmo%tpo%par%topo_rel_tau
        class default
            call fesm_error("kryos_exp_relax::init", "experiment `"//trim(self%name)//"` needs ice backend yelmo", &
                            "domain &"//trim(dom%name)//": the ice instance `"//trim(dom%ice%c%name)//"` is backend `"// &
                            trim(dom%ice%c%backend)//"`")
        end select

        self%wants(HOOK_STEP_BEGIN) = .true.
    end subroutine relax_init

    subroutine relax_hook(self, dom, id, st)
        class(relax_exp_t), intent(inout) :: self
        type(domain_t),     intent(inout) :: dom
        integer,            intent(in)    :: id
        type(hook_time_t),  intent(in)    :: st

        associate (unused => id); end associate
        select type (c => dom%ice%c)
        type is (yelmo_ice_t)
            call relax_update(self%par, c%ylmo%tpo%par, st%t_elapsed, self%topo_rel0, self%topo_rel_tau0)
        end select
    end subroutine relax_hook

end module kryos_exp_relax
