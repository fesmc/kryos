module kryos_fesm_pending
    ! placeholder: phase 1 fills this in.
    ! Types kryos needs from fesm-utils additions that do not exist yet
    ! (design section 12): field_state_t from fesm_fields (row 2), cadence_t
    ! from fesm_clock (row 4). Delete this module when they are available.

    use precision, only: dp

    implicit none

    private
    public :: field_state_t, cadence_t

    type :: field_state_t
        integer :: n = 0                    ! number of declared fields
    end type field_state_t

    type :: cadence_t
        real(dp) :: dt     = 1.0_dp         ! [yr]
        real(dp) :: phase  = 0.0_dp         ! [yr]
        real(dp) :: t_last = 0.0_dp         ! [yr]
    end type cadence_t

end module kryos_fesm_pending
