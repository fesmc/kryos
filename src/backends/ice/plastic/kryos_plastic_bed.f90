module kryos_plastic_bed
    ! The bed factor lambda of the ice backend `plastic` (docs/backend-plastic.md), shared by the sliding law "coulomb_n"
    ! (kryos_plastic_friction) and the Coulomb yield law of tau_b (taub_method = "coulomb", kryos_plastic_taub). No
    ! library: compiled whatever `plastic=`.
    !
    !   lambda = lam_min + (1 - lam_min) clamp((z_bed - z_sl - zb0) / (zb1 - zb0), 0, 1)
    !
    ! linear in the bed height above the sea level of the time, from lam_min at zb0 or below to 1 at zb1 or above (weak
    ! marine beds; Yelmo's form, Robinson et al. 2020); 1 everywhere when lam_min >= 1.
    !
    ! relaxed_bed: the present-day bed relaxed by local isostasy, z_bed + (rho_i/RHO_M) H_ice (the bed without its
    ! ice load): with taub_lambda_bed = "relaxed" the Coulomb law's lambda is a fixed (geological) field of it,
    ! relative to present sea level, so that isostatically depressed interior beds do not read as marine.

    use precision, only: dp

    implicit none

    private
    public :: bed_lambda, relaxed_bed, RHO_M

    real(dp), parameter :: RHO_M = 3300.0_dp     ! [kg m-3] mantle density of relaxed_bed

contains

    elemental function bed_lambda(z_bed, z_sl, lam_min, zb0, zb1) result(lam)
        real(dp), intent(in) :: z_bed, z_sl, lam_min, zb0, zb1
        real(dp) :: lam
        lam = 1.0_dp
        if (lam_min < 1.0_dp) lam = lam_min + (1.0_dp - lam_min)*min(max((z_bed - z_sl - zb0)/(zb1 - zb0), 0.0_dp), 1.0_dp)
    end function bed_lambda

    elemental function relaxed_bed(z_bed, H_ice, rho_i) result(z)
        ! [m] z_bed + (rho_i/RHO_M) max(H_ice, 0): the bed relaxed by local isostasy.
        real(dp), intent(in) :: z_bed, H_ice, rho_i
        real(dp) :: z
        z = z_bed + (rho_i/RHO_M)*max(H_ice, 0.0_dp)
    end function relaxed_bed

end module kryos_plastic_bed
