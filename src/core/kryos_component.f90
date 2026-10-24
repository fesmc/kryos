module kryos_component
    ! The component contract (design section 4): every component instance
    ! extends component_t; one abstract type per component adds the canonical
    ! field handles; each backend is one extension of that type.
    !
    ! placeholder: phase 1 fills this in (exchange, cadence and capabilities
    ! are stand-ins until fesm_fields/fesm_clock exist).

    use precision,          only: dp
    use coordinates,        only: grid_class
    use phys_constants,     only: phys_const_class
    use kryos_fesm_pending, only: field_state_t, cadence_t

    implicit none

    private
    public :: component_t, caps_t

    ! Capability flags, declared at init and checked against the wiring.
    type :: caps_t
        logical :: has_smb       = .false.
        logical :: has_ocn_shelf = .false.
    end type caps_t

    type, abstract :: component_t
        character(len=64) :: name, kind, backend, grid     ! instance, role, backend, native grid tag
        character(len=16) :: scope                          ! "regional" | "global"
        type(field_state_t) :: imports, exports             ! declared at init (section 5)
        type(cadence_t)     :: cadence                      ! dt, phase, t_last
        type(caps_t)        :: caps                         ! capability flags (has_smb, has_ocn_shelf, ...)
    contains
        procedure(c_init),    deferred :: init          ! (par_file, group, grid, time, cnst)
        procedure(c_update),  deferred :: update        ! (t_target): integrate from t_last to t_target; idempotent at dt=0
        procedure(c_rst_w),   deferred :: restart_write ! (dir, time)
        procedure(c_rst_r),   deferred :: restart_read  ! (dir, time)
        procedure                      :: write_output => component_write_output
        procedure                      :: end => component_end
    end type component_t

    abstract interface

        subroutine c_init(self, par_file, group, grid, time, cnst)
            import :: component_t, grid_class, phys_const_class, dp
            class(component_t),     intent(inout) :: self
            character(len=*),       intent(in)    :: par_file
            character(len=*),       intent(in)    :: group
            type(grid_class),       intent(in)    :: grid
            real(dp),               intent(in)    :: time
            type(phys_const_class), intent(in)    :: cnst
        end subroutine c_init

        subroutine c_update(self, t_target)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            real(dp),           intent(in)    :: t_target
        end subroutine c_update

        subroutine c_rst_w(self, dir, time)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: dir
            real(dp),           intent(in)    :: time
        end subroutine c_rst_w

        subroutine c_rst_r(self, dir, time)
            import :: component_t, dp
            class(component_t), intent(inout) :: self
            character(len=*),   intent(in)    :: dir
            real(dp),           intent(in)    :: time
        end subroutine c_rst_r

    end interface

contains

    subroutine component_write_output(self, time)
        ! Default: write the exports on the native grid from a variable table.
        ! placeholder: phase 1 fills this in.
        class(component_t), intent(inout) :: self
        real(dp),           intent(in)    :: time
        write(*,'(a,f10.2)') "kryos: write_output (placeholder) "//trim(self%name)//" t=", time
    end subroutine component_write_output

    subroutine component_end(self)
        class(component_t), intent(inout) :: self
        write(*,'(a)') "kryos: end (placeholder) "//trim(self%name)
    end subroutine component_end

end module kryos_component
