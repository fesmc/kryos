submodule (kryos_experiment) kryos_experiment_factory
    ! The experiment factory: the only place that names the experiment kinds:
    ! "host" (an experiment the host program passes in) and the plugins of
    ! phase 2, one `case` each (the module of a plugin is its yelmo=1 module or
    ! its stub, config/common.mk).

    use kryos_exp_relax,      only: relax_exp_t
    use kryos_exp_opt_cb_ref, only: opt_cb_ref_exp_t
    use kryos_exp_opt_tf,     only: opt_tf_exp_t
    use kryos_exp_recon,      only: recon_exp_t

    implicit none

contains

    module subroutine experiment_create(slot, ctx, protos)
        type(experiment_slot_t),  intent(out) :: slot
        type(exp_ctx_t),          intent(in)  :: ctx
        type(experiment_proto_t), intent(in)  :: protos(:)

        integer :: i, k

        select case (trim(ctx%kind))
        case ("host")
            k = 0
            do i = 1, size(protos)
                if (trim(protos(i)%name) == trim(ctx%instance)) k = i
            end do
            if (k == 0) then
                call fesm_error("kryos_experiment::experiment_create", "no host experiment named `"// &
                                trim(ctx%instance)//"`", "experiment group &"//trim(ctx%instance)// &
                                " has kind = host: the host passes an experiment_proto_t of that name to run_init")
            end if
            call experiment_adopt(slot, protos(k)%c, ctx)
        case ("relax")
            allocate(relax_exp_t :: slot%c)
            call experiment_identity(slot, ctx)
        case ("opt_cb_ref")
            allocate(opt_cb_ref_exp_t :: slot%c)
            call experiment_identity(slot, ctx)
        case ("opt_tf")
            allocate(opt_tf_exp_t :: slot%c)
            call experiment_identity(slot, ctx)
        case ("recon")
            allocate(recon_exp_t :: slot%c)
            call experiment_identity(slot, ctx)
        case default
            call fesm_error("kryos_experiment::experiment_create", "unknown experiment kind '"//trim(ctx%kind)//"'", &
                            "experiment "//trim(ctx%instance)//"; known: host, relax, opt_cb_ref, opt_tf, recon")
        end select
    end subroutine experiment_create

end submodule kryos_experiment_factory
