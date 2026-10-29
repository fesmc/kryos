submodule (kryos_experiment) kryos_experiment_factory
    ! The experiment factory: the only place that names the experiment kinds.
    ! WP1 has the kind "host" (an experiment the host program passes in); the
    ! plugins of WP2 and later add a `case` each.

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
        case default
            call fesm_error("kryos_experiment::experiment_create", "unknown experiment kind '"//trim(ctx%kind)//"'", &
                            "experiment "//trim(ctx%instance)//"; known: host")
        end select
    end subroutine experiment_create

end submodule kryos_experiment_factory
