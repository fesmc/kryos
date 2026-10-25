module kryos_forcing
    ! The forcing component, run-level scalars (contracts-design sections 2.7,
    ! 3.3, 7): abstract forcing_t (one 0-D export per coefficient), its slot,
    ! the wrappers incl. forcing_kill, forcing_restart_due and forcing_set, and
    ! the factory kryos_forcing_factory. Global only.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE, fields_put
    use kryos_component, only: component_t, PATH_LEN, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis

    implicit none

    private

    type, public :: forcing_opt_t                   ! optional imports (set in init)
        logical :: M_ice = .false.                  ! a coefficient uses the ice-mass feedback
    end type forcing_opt_t

    type, abstract, extends(component_t), public :: forcing_t
        type(forcing_opt_t) :: opt
        ! Set by the backend in init: one coefficient per name of [forcing] coefs,
        ! its units, and whether it is derived (set by a producer, section 7).
        character(len=32), allocatable :: coef_names(:), coef_units(:)
        logical,           allocatable :: coef_derived(:)
        integer, allocatable :: h_coef(:)           ! exports, one per coefficient
        integer, allocatable :: hi_derived(:)       ! imports, per coefficient (0 unless derived)
        real(dp), allocatable :: t_set(:)           ! forcing_set: time of the last set
        integer :: hi_M_ice = 0
    contains
        procedure :: kill        => forcing_kill_default          ! () -> logical: tsgen kill switch
        procedure :: restart_due => forcing_restart_due_default   ! () -> logical: restart_every_df
    end type forcing_t

    type, public :: forcing_slot_t                      ! what domain_t (or globe_t) holds
        class(forcing_t), allocatable :: c
    end type forcing_slot_t

    interface                                       ! factory: submodule kryos_forcing_factory
        module subroutine forcing_create(forcing, ctx, census)
            type(forcing_slot_t), intent(out)   :: forcing
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine forcing_create
    end interface

    public :: forcing_create, forcing_adopt, forcing_present
    public :: forcing_init, forcing_init_state, forcing_update, forcing_restart_write, forcing_restart_read
    public :: forcing_restart_files
    public :: forcing_write_output, forcing_end, forcing_kill, forcing_restart_due, forcing_set

contains

    subroutine forcing_adopt(forcing, proto, con, ctx, census)
        ! A host-provided extension of forcing_t (design decision 1), admitted as a factory would.
        type(forcing_slot_t),     intent(out)   :: forcing
        class(forcing_t), allocatable, intent(in) :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        ! proto is allocatable, so the actual is polymorphic: ifx 2023.2 builds the
        ! descriptor of a non-polymorphic actual on the caller's stack, and an object
        ! allocated with mold= from it keeps a pointer into that frame (its
        ! deallocation later reads a dead stack).
        allocate(forcing%c, mold=proto)
        call component_identity(forcing%c, ctx, "forcing", con)
        call backend_admit(forcing%c, census)
    end subroutine forcing_adopt

    logical function forcing_present(forcing)
        type(forcing_slot_t), intent(in) :: forcing
        forcing_present = allocated(forcing%c)
    end function forcing_present

    subroutine forcing_init(forcing, ctx)
        type(forcing_slot_t), intent(inout) :: forcing
        type(init_ctx_t),     intent(in)    :: ctx
        call component_enter(forcing%c)
        call forcing%c%init(ctx)
        call forcing_declare(forcing%c, ctx)
        call component_leave()
    end subroutine forcing_init

    subroutine forcing_init_state(forcing, st, pass)
        type(forcing_slot_t),  intent(inout) :: forcing
        type(step_time_t),     intent(in)    :: st
        integer,               intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(forcing%c)
        call forcing%c%init_state(s, pass)
        call component_close(forcing%c, s)
        call component_leave()
    end subroutine forcing_init_state

    subroutine forcing_update(forcing, st)
        type(forcing_slot_t),  intent(inout) :: forcing
        type(step_time_t),     intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(forcing%c)
        call forcing%c%update(st)
        call component_close(forcing%c, st)
        call component_leave()
    end subroutine forcing_update

    subroutine forcing_restart_write(forcing, dir, t)
        type(forcing_slot_t), intent(inout) :: forcing
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. forcing%c%con%restart) return             ! backend keeps no restart file
        call component_enter(forcing%c)
        call forcing%c%restart_write(dir, t)
        call component_leave()
    end subroutine forcing_restart_write

    subroutine forcing_restart_read(forcing, dir, t)
        type(forcing_slot_t), intent(inout) :: forcing
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        if (.not. forcing%c%con%restart) return
        call component_enter(forcing%c)
        call forcing%c%restart_read(dir, t)
        call component_leave()
    end subroutine forcing_restart_read

    function forcing_restart_files(forcing) result(files)
        ! The base names the backend's restart_write creates in the bundle folder.
        type(forcing_slot_t), intent(in) :: forcing
        character(len=PATH_LEN), allocatable :: files(:)
        associate (c => forcing%c)                          ! (gfortran: ICE on a call through forcing%c itself)
            files = c%restart_files()
        end associate
    end function forcing_restart_files

    subroutine forcing_write_output(forcing, file, t, stream, init)
        type(forcing_slot_t), intent(inout) :: forcing
        character(len=*),     intent(in)    :: file, stream
        real(dp),             intent(in)    :: t
        logical,              intent(in)    :: init
        call component_enter(forcing%c)
        call forcing%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine forcing_write_output

    subroutine forcing_end(forcing)
        type(forcing_slot_t), intent(inout) :: forcing
        call component_enter(forcing%c)
        call forcing%c%end()
        call component_leave()
    end subroutine forcing_end

    logical function forcing_kill_default(self)
        class(forcing_t), intent(in) :: self
        forcing_kill_default = .false.
        associate (unused => self%nx); end associate
    end function forcing_kill_default

    logical function forcing_restart_due_default(self)
        class(forcing_t), intent(in) :: self
        forcing_restart_due_default = .false.
        associate (unused => self%nx); end associate
    end function forcing_restart_due_default

    logical function forcing_kill(forcing)
        ! The run stops at the end of this step (tsgen kill switch).
        type(forcing_slot_t), intent(in) :: forcing
        forcing_kill = forcing%c%kill()
    end function forcing_kill

    logical function forcing_restart_due(forcing)
        ! A forcing-increment restart bundle is due (restart_every_df).
        type(forcing_slot_t), intent(in) :: forcing
        forcing_restart_due = forcing%c%restart_due()
    end function forcing_restart_due

    subroutine forcing_set(forcing, name, value, t)
        ! Set the derived coefficient `name` at t (its import; the backend
        ! copies it in update). Stops on an unknown or a tsgen name and on a
        ! second set of one name at one t (two producers, study section 4.3).
        type(forcing_slot_t), intent(inout) :: forcing
        character(len=*),     intent(in)    :: name
        real(dp),             intent(in)    :: value, t
        integer :: k
        character(len=*), parameter :: proc = "kryos_forcing::forcing_set"

        associate (c => forcing%c)
            k = findloc([(trim(c%coef_names(k)) == trim(name), k = 1, size(c%coef_names))], .true., dim=1)
            if (k == 0) then
                call fesm_error(proc, "unknown coefficient `"//trim(name)//"`", "instance "//trim(c%name))
            end if
            if (.not. c%coef_derived(k)) then
                call fesm_error(proc, "coefficient `"//trim(name)//"` is not derived (tsgen)", &
                                "instance "//trim(c%name))
            end if
            if (.not. (t > c%t_set(k)) .and. .not. (t < c%t_set(k))) then
                call fesm_error(proc, "coefficient `"//trim(name)//"` set twice at one time", &
                                "instance "//trim(c%name)//"; two producers of one derived coefficient")
            end if
            call fields_put(c%imports, c%hi_derived(k), value, t1=t)
            c%t_set(k) = t
        end associate
    end subroutine forcing_set

    subroutine forcing_declare(c, ctx)
        ! Contracts-design section 2.7: 0-D only; exports every coefficient,
        ! imports M_ice (feedback) and every derived coefficient.
        class(forcing_t), intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: k, n

        associate (unused => ctx%time); end associate
        if (.not. allocated(c%coef_names)) allocate(c%coef_names(0), c%coef_units(0), c%coef_derived(0))
        n = size(c%coef_names)
        allocate(c%h_coef(n), c%hi_derived(n), c%t_set(n))
        c%hi_derived = 0
        c%t_set = -huge(1.0_dp)

        call fields_init(c%imports, trim(c%name)//":imports")
        if (c%opt%M_ice) call fields_add(c%imports, "M_ice", "Gt", FIELD_INSTANT, FIELD_EXTENSIVE, c%hi_M_ice)
        do k = 1, n
            if (c%coef_derived(k)) call fields_add(c%imports, c%coef_names(k), c%coef_units(k), &
                                                   FIELD_INSTANT, FIELD_INTENSIVE, c%hi_derived(k))
        end do
        call fields_lock(c%imports)

        call fields_init(c%exports, trim(c%name)//":exports")
        do k = 1, n
            call fields_add(c%exports, c%coef_names(k), c%coef_units(k), FIELD_INSTANT, FIELD_INTENSIVE, c%h_coef(k))
        end do
        call fields_lock(c%exports)
    end subroutine forcing_declare

end module kryos_forcing
