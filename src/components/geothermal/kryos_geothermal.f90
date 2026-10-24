module kryos_geothermal
    ! The geothermal component (contracts-design sections 2.6, 3.3): abstract
    ! geothermal_t (export Q_geo [W m-2, + into the ice]), its slot, wrappers,
    ! factory, and the backend "data" (adapters-design 7), which lives here.

    use precision,       only: dp
    use fesm_log,        only: fesm_error
    use fesm_fields,     only: fields_init, fields_add, fields_lock, FIELD_INSTANT, FIELD_MEAN, &
                               FIELD_INTENSIVE, FIELD_EXTENSIVE, fields_get, fields_put
    use fesm_datasrc,    only: datasrc_t, datasrc_init, datasrc_update, datasrc_fill, datasrc_describe
    use nml,             only: nml_read, nml_has_param
    use kryos_component, only: component_t, init_ctx_t, step_time_t, census_t, constraints_t, &
                               component_identity, backend_admit, component_close, &
                               component_enter, component_leave, component_optional, component_axis, PASS_GEOMETRY

    implicit none

    private

    ! Schema of the field groups (synchronised copy of fesm-utils par/datasrc_defaults.nml)
    character(len=*), parameter, public :: DATASRC_DEFAULTS = "input/datasrc_defaults.nml"

    type, abstract, extends(component_t), public :: geothermal_t
        integer :: h_Q_geo = 0
    end type geothermal_t

    type, public :: geothermal_slot_t                      ! what domain_t (or globe_t) holds
        class(geothermal_t), allocatable :: c
    end type geothermal_slot_t

    type, extends(geothermal_t), public :: data_geothermal_t    ! backend "data": component_t + datasrc_t (adapters-design 7)
        type(init_ctx_t) :: ctx                     ! kept from init until the field is served
        type(datasrc_t)  :: ds
        real(dp) :: Q_geo_min = -huge(1.0_dp)     ! floor after the fill ([instance] Q_geo_min)
    contains
        procedure :: init          => data_init_
        procedure :: init_state    => data_init_state_
        procedure :: update        => data_update_
        procedure :: restart_write => data_restart_
        procedure :: restart_read  => data_restart_
    end type data_geothermal_t

    interface                                       ! factory: submodule kryos_geothermal_factory
        module subroutine geothermal_create(geothermal, ctx, census)
            type(geothermal_slot_t), intent(out)   :: geothermal
            type(init_ctx_t), intent(in)    :: ctx
            type(census_t),   intent(inout) :: census
        end subroutine geothermal_create
    end interface

    public :: geothermal_create, geothermal_adopt, geothermal_present
    public :: geothermal_init, geothermal_init_state, geothermal_update, geothermal_restart_write, geothermal_restart_read
    public :: geothermal_write_output, geothermal_end, data_geothermal_constraints

contains

    subroutine geothermal_adopt(geothermal, proto, con, ctx, census)
        ! A host-provided extension of geothermal_t (design decision 1), admitted as a factory would.
        type(geothermal_slot_t),        intent(out)   :: geothermal
        class(geothermal_t),            intent(in)    :: proto
        type(constraints_t), intent(in)    :: con
        type(init_ctx_t),    intent(in)    :: ctx
        type(census_t),      intent(inout) :: census
        allocate(geothermal%c, mold=proto)
        call component_identity(geothermal%c, ctx, "geothermal", con)
        call backend_admit(geothermal%c, census)
    end subroutine geothermal_adopt

    logical function geothermal_present(geothermal)
        type(geothermal_slot_t), intent(in) :: geothermal
        geothermal_present = allocated(geothermal%c)
    end function geothermal_present

    subroutine geothermal_init(geothermal, ctx)
        type(geothermal_slot_t), intent(inout) :: geothermal
        type(init_ctx_t),        intent(in)    :: ctx
        call component_enter(geothermal%c)
        call geothermal%c%init(ctx)
        call geothermal_declare(geothermal%c, ctx)
        call component_leave()
    end subroutine geothermal_init

    subroutine geothermal_init_state(geothermal, st, pass)
        type(geothermal_slot_t),  intent(inout) :: geothermal
        type(step_time_t),        intent(in)    :: st
        integer,                  intent(in)    :: pass
        type(step_time_t) :: s
        s = st
        s%t0 = st%t1 - st%dt                            ! virtual previous step for MEAN puts
        call component_enter(geothermal%c)
        call geothermal%c%init_state(s, pass)
        call component_close(geothermal%c, s)
        call component_leave()
    end subroutine geothermal_init_state

    subroutine geothermal_update(geothermal, st)
        type(geothermal_slot_t),  intent(inout) :: geothermal
        type(step_time_t),        intent(in)    :: st
        if (.not. (st%t1 > st%t0)) return               ! idempotent at dt = 0
        call component_enter(geothermal%c)
        call geothermal%c%update(st)
        call component_close(geothermal%c, st)
        call component_leave()
    end subroutine geothermal_update

    subroutine geothermal_restart_write(geothermal, dir, t)
        type(geothermal_slot_t), intent(inout) :: geothermal
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t
        if (.not. geothermal%c%con%restart) return             ! backend keeps no restart file
        call component_enter(geothermal%c)
        call geothermal%c%restart_write(dir, t)
        call component_leave()
    end subroutine geothermal_restart_write

    subroutine geothermal_restart_read(geothermal, dir, t)
        type(geothermal_slot_t), intent(inout) :: geothermal
        character(len=*),        intent(in)    :: dir
        real(dp),                intent(in)    :: t
        if (.not. geothermal%c%con%restart) return
        call component_enter(geothermal%c)
        call geothermal%c%restart_read(dir, t)
        call component_leave()
    end subroutine geothermal_restart_read

    subroutine geothermal_write_output(geothermal, file, t, stream, init)
        type(geothermal_slot_t), intent(inout) :: geothermal
        character(len=*),        intent(in)    :: file, stream
        real(dp),                intent(in)    :: t
        logical,                 intent(in)    :: init
        call component_enter(geothermal%c)
        call geothermal%c%write_output(file, t, stream, init)
        call component_leave()
    end subroutine geothermal_write_output

    subroutine geothermal_end(geothermal)
        type(geothermal_slot_t), intent(inout) :: geothermal
        call component_enter(geothermal%c)
        call geothermal%c%end()
        call component_leave()
    end subroutine geothermal_end

    ! ========================================================================
    ! Backend "data" (Q_geo from a field group through fesm_datasrc)
    ! ========================================================================

    function data_geothermal_constraints() result(con)
        ! Regional or global, any grid; no backend restart file (re-served at init).
        type(constraints_t) :: con
        con = constraints_t(global=.true., restart=.false.)
    end function data_geothermal_constraints

    subroutine data_init_(self, ctx)
        ! Keep the context; the field is bound and served in init_state
        ! (PASS_GEOMETRY), after the wrapper has declared the exports.
        class(data_geothermal_t), intent(inout) :: self
        type(init_ctx_t),      intent(in)    :: ctx
        self%ctx = ctx
        self%nx  = ctx%grid%G%nx
        self%ny  = ctx%grid%G%ny
        self%prec = dp
        if (nml_has_param(ctx%par_file, ctx%instance, "Q_geo_min")) &
            call nml_read(ctx%par_file, ctx%instance, "Q_geo_min", self%Q_geo_min)
    end subroutine data_init_

    subroutine data_init_state_(self, st, pass)
        ! PASS_GEOMETRY: datasrc_init, datasrc_update(t1) and datasrc_fill,
        ! once; the export is INSTANT, valid at [t1, t1]. Other passes: nothing.
        class(data_geothermal_t), intent(inout) :: self
        type(step_time_t),     intent(in)    :: st
        integer,              intent(in)    :: pass
        real(dp) :: v(self%nx, self%ny)
        integer  :: n
        if (pass /= PASS_GEOMETRY) return
        n = count(len_trim(self%ctx%library) > 0)
        call datasrc_init(self%ds, self%ctx%par_file, self%ctx%instance, self%exports, ["Q_geo"], &
                          self%ctx%grid, DATASRC_DEFAULTS, var_defs=pack(self%ctx%library, &
                          len_trim(self%ctx%library) > 0), domain=self%ctx%domain)
        call datasrc_describe(self%ds)
        call datasrc_update(self%ds, st%t1)
        call datasrc_fill(self%ds, self%exports)
        if (self%Q_geo_min > -huge(1.0_dp)) then            ! the floor (today's ghf_min)
            call fields_get(self%exports, self%h_Q_geo, v)
            call fields_put(self%exports, self%h_Q_geo, max(v, self%Q_geo_min), t1=st%t1)
        end if
        associate (unused => n + v(1,1)); end associate
    end subroutine data_init_state_

    subroutine data_update_(self, st)
        ! Static: nothing to integrate (cadence off).
        class(data_geothermal_t), intent(inout) :: self
        type(step_time_t),     intent(in)    :: st
        associate (unused => st%t1 + self%nx); end associate
    end subroutine data_update_

    subroutine data_restart_(self, dir, t)
        ! Never called (constraints restart = .false.): the export state is the core's.
        class(data_geothermal_t), intent(inout) :: self
        character(len=*),      intent(in)    :: dir
        real(dp),              intent(in)    :: t
        associate (unused => len(dir) + t + self%nx); end associate
    end subroutine data_restart_

    subroutine geothermal_declare(c, ctx)
        ! Contracts-design section 2.6: export Q_geo only.
        class(geothermal_t), intent(inout) :: c
        type(init_ctx_t), intent(in)    :: ctx
        integer :: d2(2)
        character(len=:), allocatable :: g

        d2 = [c%nx, c%ny]
        g  = trim(c%grid)

        associate (unused => ctx%time); end associate
        call fields_init(c%imports, trim(c%name)//":imports")
        call fields_lock(c%imports)
        call fields_init(c%exports, trim(c%name)//":exports")
        call fields_add(c%exports, "Q_geo", "W m-2", FIELD_INSTANT, FIELD_INTENSIVE, c%h_Q_geo, grid=g, dims=d2)
        call fields_lock(c%exports)
    end subroutine geothermal_declare

end module kryos_geothermal
