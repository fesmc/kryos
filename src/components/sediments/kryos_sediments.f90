module kryos_sediments
    ! The sediments component (contracts-design sections 2.6, 3.3): abstract
    ! sediments_t (export H_sed [m]), its slot, wrappers, factory, and the
    ! backend "data" (adapters-design 7), which lives here.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE, fields_get, fields_put
    use fesm_datasrc,    only: datasrc_t, datasrc_init, datasrc_update, datasrc_fill, datasrc_describe
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis, PASS_INIT

    implicit none

    private

    ! Schema of the field groups (synchronised copy of fesm-utils par/datasrc_defaults.nml)
    character(len=*), parameter, public :: DATASRC_DEFAULTS = "input/datasrc_defaults.nml"

    type, abstract, extends(component_t), public :: sediments_t
        integer :: h_H_sed = 0
    end type sediments_t

    type, public :: sediments_slot_t                      ! what domain_t (or globe_t) holds
        class(sediments_t), allocatable :: c
    end type sediments_slot_t

    type, extends(sediments_t), public :: data_sediments_t    ! backend "data": component_t + datasrc_t (adapters-design 7)
        type(init_ctx_t) :: ctx                     ! kept from init until the field is served
        type(datasrc_t)  :: ds
    contains
        procedure :: init          => data_init_
        procedure :: init_state    => data_init_state_
        procedure :: update        => data_update_
        procedure :: restart_write => data_restart_
        procedure :: restart_read  => data_restart_
    end type data_sediments_t

    interface                                       ! factory: submodule kryos_sediments_factory
        module subroutine sediments_create(sediments, ctx, census)
            type(sediments_slot_t), intent(out)   :: sediments
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine sediments_create
    end interface

    public :: sediments_create, sediments_adopt, sediments_present
    public :: sediments_init, sediments_init_state, sediments_update, sediments_restart_write, sediments_restart_read
    public :: sediments_write_output, sediments_end, data_sediments_constraints

contains

    subroutine sediments_adopt(sediments, proto, con, ctx, census)
        ! A host-provided extension of sediments_t (design decision 1), admitted as a factory would.
        type(sediments_slot_t),       intent(out)   :: sediments
        class(sediments_t),           intent(in)    :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        allocate(sediments%c, mold=proto)
        call component_identity(sediments%c, ctx, "sediments", con)
        call backend_admit(sediments%c, census)
    end subroutine sediments_adopt

    logical function sediments_present(sediments)
        type(sediments_slot_t), intent(in) :: sediments
        sediments_present = allocated(sediments%c)
    end function sediments_present

    subroutine sediments_init(sediments, ctx)
        type(sediments_slot_t), intent(inout) :: sediments
        type(init_ctx_t),       intent(in)    :: ctx
        call component_enter(sediments%c)
        call sediments%c%init(ctx)
        call sediments_declare(sediments%c, ctx)
        call component_leave()
    end subroutine sediments_init

    subroutine sediments_init_state(sediments, st, pass)
        type(sediments_slot_t),  intent(inout) :: sediments
        type(step_time_t),       intent(in)    :: st
        integer,                 intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(sediments%c)
        call sediments%c%init_state(s, pass)
        call component_close(sediments%c, s)
        call component_leave()
    end subroutine sediments_init_state

    subroutine sediments_update(sediments, st)
        type(sediments_slot_t),  intent(inout) :: sediments
        type(step_time_t),       intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(sediments%c)
        call sediments%c%update(st)
        call component_close(sediments%c, st)
        call component_leave()
    end subroutine sediments_update

    subroutine sediments_restart_write(sediments, dir, t)
        type(sediments_slot_t), intent(inout) :: sediments
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t
        if (.not. sediments%c%con%restart) return             ! backend keeps no restart file
        call component_enter(sediments%c)
        call sediments%c%restart_write(dir, t)
        call component_leave()
    end subroutine sediments_restart_write

    subroutine sediments_restart_read(sediments, dir, t)
        type(sediments_slot_t), intent(inout) :: sediments
        character(len=*),       intent(in)    :: dir
        real(dp),               intent(in)    :: t
        if (.not. sediments%c%con%restart) return
        call component_enter(sediments%c)
        call sediments%c%restart_read(dir, t)
        call component_leave()
    end subroutine sediments_restart_read

    subroutine sediments_write_output(sediments, file, t, stream, init)
        type(sediments_slot_t), intent(inout) :: sediments
        character(len=*),       intent(in)    :: file, stream
        real(dp),               intent(in)    :: t
        logical,                intent(in)    :: init
        call component_enter(sediments%c)
        call sediments%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine sediments_write_output

    subroutine sediments_end(sediments)
        type(sediments_slot_t), intent(inout) :: sediments
        call component_enter(sediments%c)
        call sediments%c%end()
        call component_leave()
    end subroutine sediments_end

    ! ========================================================================
    ! Backend "data" (H_sed from a field group through fesm_datasrc)
    ! ========================================================================

    function data_sediments_constraints() result(con)
        ! Regional or global, any grid; no backend restart file (re-served at init).
        type(constraints_t) :: con
        con = constraints_t(global=.true., restart=.false.)
    end function data_sediments_constraints

    subroutine data_init_(self, ctx)
        ! Keep the context; the field is bound and served in init_state
        ! (PASS_INIT), after the wrapper has declared the exports.
        class(data_sediments_t), intent(inout) :: self
        type(init_ctx_t),     intent(in)    :: ctx
        self%ctx = ctx
        self%nx  = ctx%grid%G%nx
        self%ny  = ctx%grid%G%ny
        self%prec = dp
    end subroutine data_init_

    subroutine data_init_state_(self, st, pass)
        ! PASS_INIT: datasrc_init, datasrc_update(t1) and datasrc_fill,
        ! once; the export is INSTANT, valid at [t1, t1]. Other passes: nothing.
        class(data_sediments_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        integer,             intent(in)    :: pass
        real(dp) :: v(self%nx, self%ny)
        integer  :: n
        if (pass /= PASS_INIT) return
        n = count(len_trim(self%ctx%library) > 0)
        call datasrc_init(self%ds, self%ctx%par_file, self%ctx%instance, self%exports, ["H_sed"], &
                          self%ctx%grid, DATASRC_DEFAULTS, var_defs=pack(self%ctx%library, &
                          len_trim(self%ctx%library) > 0), domain=self%ctx%domain)
        call datasrc_describe(self%ds)
        call datasrc_update(self%ds, st%t1)
        call datasrc_fill(self%ds, self%exports)
        associate (unused => n + v(1,1)); end associate
    end subroutine data_init_state_

    subroutine data_update_(self, st)
        ! Static: nothing to integrate (cadence off).
        class(data_sediments_t), intent(inout) :: self
        type(step_time_t),    intent(in)    :: st
        associate (unused => st%t1 + self%nx); end associate
    end subroutine data_update_

    subroutine data_restart_(self, dir, t)
        ! Never called (constraints restart = .false.): the export state is the core's.
        class(data_sediments_t), intent(inout) :: self
        character(len=*),     intent(in)    :: dir
        real(dp),             intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine data_restart_

    subroutine sediments_declare(c, ctx)
        ! Contracts-design section 2.6: export H_sed only.
        class(sediments_t), intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        associate (unused => ctx%time); end associate
        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_lock(c%imports)
        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "H_sed", "m", FIELD_INSTANT, FIELD_INTENSIVE, c%h_H_sed, grid=g, dims=d2)
        call fields_lock(c%exports)
    end subroutine sediments_declare

end module kryos_sediments
