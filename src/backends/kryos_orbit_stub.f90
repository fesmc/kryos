module kryos_orbit
    ! Stub of the orbit control (insol=0): kryos_orbit.f90 needs libinsol. Its
    ! users (the chion host, the climate paleoemu) are not built without insol
    ! (config/common.mk stops chion=1 or paleoemu=1 with insol=0), so the stub
    ! holds nothing.
    implicit none
    private
end module kryos_orbit
