module kryos_exp_relax
    ! Stub of the experiment `relax`, compiled when kryos is built with
    ! yelmo=0: same module and type as kryos_exp_relax.f90, no libyelmo. Its
    ! init stops with the rebuild hint.

    use fesm_log,        only: fesm_error
    use kryos_domain,    only: domain_t
    use kryos_experiment, only: experiment_t, exp_ctx_t, hook_time_t

    implicit none

    private
    public :: relax_exp_t

    type, extends(experiment_t) :: relax_exp_t
    contains
        procedure :: init => stub_init
        procedure :: hook => stub_hook
    end type relax_exp_t

contains

    subroutine stub_init(self, ctx, dom)
        class(relax_exp_t), intent(inout) :: self
        type(exp_ctx_t),   intent(in)    :: ctx
        type(domain_t),    intent(inout) :: dom
        associate (unused => ctx%time + len(dom%name)); end associate
        call fesm_error("kryos_exp_relax", "experiment `"//trim(self%name)// &
                        "` needs ice backend yelmo: rebuild kryos with `yelmo=1`")
    end subroutine stub_init

    subroutine stub_hook(self, dom, id, st)
        class(relax_exp_t), intent(inout) :: self
        type(domain_t),    intent(inout) :: dom
        integer,           intent(in)    :: id
        type(hook_time_t), intent(in)    :: st
        associate (unused => len(dom%name) + id + st%t1); end associate
        call fesm_error("kryos_exp_relax", "experiment `"//trim(self%name)// &
                        "` needs ice backend yelmo: rebuild kryos with `yelmo=1`")
    end subroutine stub_hook

end module kryos_exp_relax
