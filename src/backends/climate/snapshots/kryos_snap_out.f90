module kryos_snap_out
    ! The adapter's internal now/ref record of the climate (yelmox climate_out,
    ! snapshots branch): the snapshots backend fills `now` and `ref` after each
    ! update, the terms and transforms work on it, and the adapter exports from
    ! it. Not a framework type. Holding both `now` and `ref` lets the
    ! transforms and the exports form anomalies (now - ref).
    !
    ! The atmosphere is valid at the surface elevation z_srf of the climate
    ! grid; the monthly lapse rate and the precipitation sensitivity beta_p
    ! move it to another surface (the wiring, contracts-design 4.1):
    !   tas(z) = tas - lapse*(z - z_srf),  pr(z) = pr*exp(-beta_p*lapse*(z - z_srf)).

    use precision, only: wp

    implicit none

    private

    type, public :: clim_state_class
        ! atmosphere
        real(wp), allocatable :: tas(:,:,:)     ! near-surface air temperature, monthly [K]
        real(wp), allocatable :: pr(:,:,:)      ! precipitation, monthly [mm/d]
        real(wp), allocatable :: tsl_ann(:,:)   ! sea-level air temperature, annual [K]
        real(wp), allocatable :: ta_ann(:,:)    ! near-surface air temperature, annual [K]
        real(wp), allocatable :: pr_ann(:,:)    ! precipitation, annual [mm/a]
        real(wp), allocatable :: z_srf(:,:)     ! surface elevation of tas/pr [m]
        real(wp), allocatable :: beta_p(:,:)    ! precipitation sensitivity to temperature [1/K]
        ! ocean as depth profiles
        real(wp), allocatable :: to_ann(:,:,:)  ! ocean temperature over depth, annual [K]
        real(wp), allocatable :: so_ann(:,:,:)  ! ocean salinity over depth, annual [psu]
        real(wp), allocatable :: depth(:)       ! ocean depth axis [m]
    end type clim_state_class

    type, public :: climate_out_class
        type(clim_state_class) :: now           ! current climate
        type(clim_state_class) :: ref           ! reference climate (anomaly baseline)
        real(wp) :: lapse(12) = 0.0_wp          ! monthly lapse rate of tas [K/m]
    end type climate_out_class

    public :: summer_mean

contains

    function summer_mean(var, south) result(var_sum)
        ! Mean of a monthly field over the summer months: DJF in the southern
        ! hemisphere, JJA in the northern.
        real(wp), intent(in) :: var(:,:,:)
        logical,  intent(in) :: south
        real(wp) :: var_sum(size(var,1),size(var,2))

        if (south) then
            var_sum = (var(:,:,1)+var(:,:,2)+var(:,:,12)) / 3.0_wp
        else
            var_sum = (var(:,:,6)+var(:,:,7)+var(:,:,8)) / 3.0_wp
        end if

    end function summer_mean

end module kryos_snap_out
