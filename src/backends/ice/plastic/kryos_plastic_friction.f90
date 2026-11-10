module kryos_plastic_friction
    ! The sliding laws of the ice backend `plastic` whose coefficient follows the imposed tau_b, per cell
    ! (docs/dev/plastic-backend.md, "Implied mass balance"; drive = "taub"). No library: compiled whatever `plastic=`.
    !
    !   "yield"      the library's Weertman law, beta = tau_b / u_y^(1/m), u_y = uy (tau_b / 100 kPa)^uy_b: the ice
    !                slides at u_y where the stress is at its yield value tau_b (uy_b < 0: the weaker the bed, the
    !                faster);
    !   "coulomb_n"  the library's regularised Coulomb law with C = tau_b lambda(z_bed - z_sl) N / (rho_i g Hc):
    !                N = rho_i g H (1 - H_f/H)^p the effective pressure (Leguy et al. 2014; H_f = rho_w/rho_i
    !                max(0, z_sl - z_bed) the flotation thickness; 0 at or below it), lambda linear from lam_min at
    !                z_bed - z_sl <= zb0 to 1 at >= zb1. At tau = tau_b the speed u_0 x/(1-x), x = (Hc/(lambda N/
    !                (rho_i g)))^m, depends on N and the bed alone: slow under thick ice, fast where the ice thins
    !                below Hc/lambda or nears flotation; u_max where C <= tau.
    !
    ! Cells off the ice or without a positive tau_b get no sliding (huge coefficient; a margin face takes the
    ! coefficient of its ice cell).

    use precision, only: dp

    implicit none

    private

    type, public :: friction_t
        character(len=16) :: kind = ""                   ! "" (the law's uniform coefficient) | "yield" | "coulomb_n"
        real(dp) :: uy = 0.0_dp, uy_b = 0.0_dp           ! yield: [m yr-1], exponent
        real(dp) :: Hc = 0.0_dp                          ! coulomb_n: [m] the thickness scale of the capacity
        real(dp) :: lam_min = 1.0_dp                     ! coulomb_n: lambda on low beds
        real(dp) :: zb0 = -500.0_dp, zb1 = 0.0_dp        ! coulomb_n: [m] the bed above sea level of lam_min, of 1
        real(dp) :: p = 0.5_dp                           ! coulomb_n: the exponent of N
    end type friction_t

    public :: friction_law, friction_map

contains

    function friction_law(fr) result(law)
        ! The library law of fr%kind: "weertman" (yield), "coulomb" (coulomb_n); "" for any other kind.
        type(friction_t), intent(in) :: fr
        character(len=8) :: law
        select case (trim(fr%kind))
        case ("yield");      law = "weertman"
        case ("coulomb_n");  law = "coulomb"
        case default;        law = ""
        end select
    end function friction_law

    subroutine friction_map(fr, m, tau_b, H, z_bed, z_sl, mask, rho_i, rho_w, g, fric)
        ! The coefficient per cell of fr%kind (the library's slide%fric_map) on the geometry H, z_bed, z_sl, mask
        ! and tau_b [Pa]; m the law's stress exponent.
        type(friction_t), intent(in)  :: fr
        real(dp),         intent(in)  :: m, rho_i, rho_w, g
        real(dp),         intent(in)  :: tau_b(:,:), H(:,:), z_bed(:,:), z_sl(:,:)
        logical,          intent(in)  :: mask(:,:)
        real(dp),         intent(out) :: fric(:,:)

        real(dp), allocatable :: Hf(:,:), N(:,:), lam(:,:)

        select case (trim(fr%kind))
        case ("yield")
            where (tau_b > 0.0_dp)
                fric = tau_b/(fr%uy*(tau_b/1.0e5_dp)**fr%uy_b)**(1.0_dp/m)
            elsewhere
                fric = huge(1.0_dp)
            end where
        case ("coulomb_n")
            Hf = rho_w/rho_i*max(z_sl - z_bed, 0.0_dp)
            allocate(N(size(H, 1), size(H, 2)), lam(size(H, 1), size(H, 2)))
            N = 0.0_dp
            where (H > Hf) N = rho_i*g*H*(1.0_dp - Hf/H)**fr%p
            lam = 1.0_dp
            if (fr%lam_min < 1.0_dp) lam = fr%lam_min + (1.0_dp - fr%lam_min)* &
                                           min(max((z_bed - z_sl - fr%zb0)/(fr%zb1 - fr%zb0), 0.0_dp), 1.0_dp)
            where (mask .and. tau_b > 0.0_dp)
                fric = tau_b*lam*N/(rho_i*g*fr%Hc)
            elsewhere
                fric = huge(1.0_dp)
            end where
        case default
            fric = huge(1.0_dp)
        end select
    end subroutine friction_map

end module kryos_plastic_friction
