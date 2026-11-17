module kryos_plastic_bed
    ! The bed factor lambda of the ice backend `plastic` (docs/backend-plastic.md), shared by the sliding law "coulomb_n"
    ! (kryos_plastic_friction) and the Coulomb yield law of tau_b (taub_method = "coulomb", kryos_plastic_taub). No
    ! library: compiled whatever `plastic=`.
    !
    !   lambda = lam_min + (1 - lam_min) clamp((z_bed - z_sl - zb0) / (zb1 - zb0), 0, 1)
    !
    ! linear in the bed height above the sea level of the time, from lam_min at zb0 or below to 1 at zb1 or above (weak
    ! marine beds; Yelmo's form, Robinson et al. 2020); 1 everywhere when lam_min >= 1.

    use precision, only: dp

    implicit none

    private
    public :: bed_lambda

contains

    elemental function bed_lambda(z_bed, z_sl, lam_min, zb0, zb1) result(lam)
        real(dp), intent(in) :: z_bed, z_sl, lam_min, zb0, zb1
        real(dp) :: lam
        lam = 1.0_dp
        if (lam_min < 1.0_dp) lam = lam_min + (1.0_dp - lam_min)*min(max((z_bed - z_sl - zb0)/(zb1 - zb0), 0.0_dp), 1.0_dp)
    end function bed_lambda

end module kryos_plastic_bed
