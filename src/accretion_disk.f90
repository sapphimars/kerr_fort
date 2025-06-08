module accretion_disk
    use kerr_geometry, only: mass_bh, spin_bh
    implicit none
    private

    public :: isco_radius

contains

    ! Calculates the ISCO radius
    ! for a given black hole mass and spin - prograde orbit
    function isco_radius() result(r_isco)
        real(kind=8) :: r_isco
        real(kind=8) :: a_norm, z1, z2

        ! Spin parameter normalized by mass
        a_norm = spin_bh / mass_bh

        ! Auxiliary quantities from Bardeen, Press, & Teukolsky (1972)
        z1 = 1.0d0 + (1.0d0 - a_norm**2)**(1.0d0/3.0d0) * &
             ((1.0d0 + a_norm)**(1.0d0/3.0d0) + (1.0d0 - a_norm)**(1.0d0/3.0d0))
        z2 = sqrt(3.0d0 * a_norm**2 + z1**2)

        ! ISCO radius for a prograde orbit (co-rotating with the spin)
        r_isco = mass_bh * (3.0d0 + z2 - sqrt((3.0d0 - z1)*(3.0d0 + z1 + 2.0d0*z2)))

    end function

end module