module kerr_geometry
   implicit none
   private

   ! Lots of help from here: https://arxiv.org/abs/0904.4184

   ! Physical constants
   real(kind = 8), parameter :: c = 1.0d0  ! c=G=1 

   ! Black hole parameters
   real(kind = 8) :: mass_bh = 1.0d0      ! M 
   real(kind = 8) :: spin_bh = 0.5d0      ! a 

   public :: set_bh_params, kerr_metric, kerr_metric_contravariant, rho_squared, delta_kerr, &
   kerr_metric_determinant, kerr_metric_det_closed, kerr_christoffel, kerr_geodesic_rhs, &
   kerr_constants_of_motion, set_circular_orbit_ic, debug_circular_orbit, mass_bh, spin_bh

   contains

   subroutine set_bh_params(m, a)
      real(kind = 8), intent(in) :: m, a
      mass_bh = m
      spin_bh = a
   end subroutine

   ! Helper functions
   function rho_squared(r, theta)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8) :: rho_squared
      rho_squared = r**2 + spin_bh**2 * cos(theta)**2
   end function

   function delta_kerr(r)
      real(kind = 8), intent(in) :: r
      real(kind = 8) :: delta_kerr
      delta_kerr = r**2 - 2.0d0 * mass_bh * r + spin_bh**2
   end function

   ! Kerr metric tensor g_munu in covariant form
   subroutine kerr_metric(r, theta, g_covar)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8), intent(out) :: g_covar(0:3, 0:3)

      real(kind = 8) :: rho2, delta, sin2_theta, cos2_theta

      ! Initialize to zero
      g_covar = 0.0d0

      ! Compute helper quantities
      rho2 = rho_squared(r, theta)
      delta = delta_kerr(r)
      sin2_theta = sin(theta)**2
      cos2_theta = cos(theta)**2

      ! Diagonal components
      g_covar(0, 0) = -(1.0d0 - 2.0d0 * mass_bh * r / rho2)  ! g_tt
      g_covar(1, 1) = rho2 / delta                         ! g_rr  
      g_covar(2, 2) = rho2                               ! g_thetatheta
      g_covar(3, 3) = sin2_theta * (r**2 + spin_bh**2 + &
      2.0d0 * mass_bh * spin_bh**2 * r * sin2_theta / rho2)  ! g_phiphi

      ! Off-diagonal components (only non-zero is t-φ coupling)
      g_covar(0, 3) = -2.0d0 * mass_bh * spin_bh * r * sin2_theta / rho2  ! g_tphi
      g_covar(3, 0) = g_covar(0, 3)                               ! g_phit

   end subroutine

   ! Kerr metric tensor g^munu in contravariant form
   subroutine kerr_metric_contravariant(r, theta, g_contravar)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8), intent(out) :: g_contravar(0:3, 0:3)

      real(kind = 8) :: rho2, delta, sin2_theta
      real(kind = 8) :: g_tt, g_rr, g_thth, g_phph, g_tph
      real(kind = 8) :: det_tph_block

      ! Initialize to zero
      g_contravar = 0.0d0

      ! Compute helper quantities
      rho2 = rho_squared(r, theta)
      delta = delta_kerr(r)
      sin2_theta = sin(theta)**2

      ! Covariant components (for inversion)
      g_tt = -(1.0d0 - 2.0d0 * mass_bh * r / rho2)
      g_rr = rho2 / delta
      g_thth = rho2
      g_phph = sin2_theta * (r**2 + spin_bh**2 + &
      2.0d0 * mass_bh * spin_bh**2 * r * sin2_theta / rho2)
      g_tph = -2.0d0 * mass_bh * spin_bh * r * sin2_theta / rho2

      ! Diagonal components (simple inverse)
      g_contravar(1, 1) = 1.0d0 / g_rr     ! g^rr = 1/g_rr
      g_contravar(2, 2) = 1.0d0 / g_thth   ! g^θθ = 1/g_thetatheta

      ! (t,phi) block inversion
      ! For 2x2 block: det = g_tt * g_phiphi - g_tphi²
      det_tph_block = g_tt * g_phph - g_tph**2

      g_contravar(0, 0) = g_phph / det_tph_block    ! g^tt
      g_contravar(3, 3) = g_tt / det_tph_block      ! g^φφ
      g_contravar(0, 3) = -g_tph / det_tph_block    ! g^tφ
      g_contravar(3, 0) = g_contravar(0, 3)          ! g^φt

   end subroutine

   ! Kerr metric determinant (using block structure for edge cases)
   function kerr_metric_determinant(r, theta)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8) :: kerr_metric_determinant

      real(kind = 8) :: rho2, delta, sin2_theta
      real(kind = 8) :: g_tt, g_phph, g_tph, det_tph_block

      ! Compute helper quantities
      rho2 = rho_squared(r, theta)
      delta = delta_kerr(r)
      sin2_theta = sin(theta)**2

      ! Method 1: Use block structure
      ! det(g) = g_rr * g_thth * det(t-ph block)

      ! (t,ph) block determinant
      g_tt = -(1.0d0 - 2.0d0 * mass_bh * r / rho2)
      g_phph = sin2_theta * (r**2 + spin_bh**2 + &
      2.0d0 * mass_bh * spin_bh**2 * r * sin2_theta / rho2)
      g_tph = -2.0d0 * mass_bh * spin_bh * r * sin2_theta / rho2

      det_tph_block = g_tt * g_phph - g_tph**2

      ! Full determinant
      kerr_metric_determinant = (rho2 / delta) * rho2 * det_tph_block

      ! Note: This simplifies to the known result

   end function

   ! Direct closed form (more efficient for repeated calls)
   function kerr_metric_det_closed(r, theta)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8) :: kerr_metric_det_closed

      real(kind = 8) :: rho2, sin2_theta, sigma_squared

      rho2 = rho_squared(r, theta)
      sin2_theta = sin(theta)**2

      ! Sig² = (r2 + a2)2 - a2*Delta sin2theta
      sigma_squared = (r**2 + spin_bh**2)**2 - &
      spin_bh**2 * delta_kerr(r) * sin2_theta

      kerr_metric_det_closed = -rho2 * sin2_theta * sigma_squared

   end function

   ! Christoffels for Kerr metric
   subroutine kerr_christoffel(r, theta, christoffel)
      real(kind = 8), intent(in) :: r, theta
      real(kind = 8), intent(out) :: christoffel(0:3, 0:3, 0:3)

      real(kind = 8) :: rho2, delta, sin_theta, cos_theta, sin2_theta, &
      cos2_theta
      real(kind = 8) :: m, a, a2, r_s
      real(kind = 8) :: a_term ! As defined in reference notes for Γ^th_phph

      ! Initialize to zero
      christoffel = 0.0d0

      ! Shortcuts and standard notation from reference
      m = mass_bh
      a = spin_bh
      a2 = a**2
      r_s = 2.0d0 * m ! Schwarzschild radius
      sin_theta = sin(theta)
      cos_theta = cos(theta)
      sin2_theta = sin_theta**2
      cos2_theta = cos_theta**2

      rho2 = r**2 + a2 * cos2_theta ! This is Sigma in the reference
      delta = r**2 - r_s * r + a2

      ! --- ^t components ---
      ! ^t_tt (2.14.5a, right)
      christoffel(0, 0, 0) = (r_s * a2 * r * sin_theta * cos_theta) / rho2**3
      ! ^t_tr (2.14.5b, left)
      christoffel(0, 0, 1) = (r_s * (r**2 + a2) * (r**2 - a2 * cos2_theta)) / (2.0d0 * rho2**2 * delta)
      christoffel(0, 1, 0) = christoffel(0, 0, 1)
      ! ^t_tth (2.14.5c, left)
      christoffel(0, 0, 2) = -(r_s * a2 * r * sin_theta * cos_theta) / rho2**2
      christoffel(0, 2, 0) = christoffel(0, 0, 2)
      ! ^t_rph (2.14.5i, left)
      christoffel(0, 1, 3) = (r_s * a * sin2_theta * (a2 * cos2_theta * (a2 - r**2) - r**2 * (a2 + 3.0d0 * r**2)))&
       / (2.0d0 * rho2**2 * delta)
      christoffel(0, 3, 1) = christoffel(0, 1, 3)
      ! ^t_thph (2.14.5h, right)
      christoffel(0, 2, 3) = (r_s * a**3 * r * sin_theta**3 * cos_theta) / rho2**2
      christoffel(0, 3, 2) = christoffel(0, 2, 3)

      ! --- ^r components ---
      ! ^r_tt (2.14.5a, left)
      christoffel(1, 0, 0) = (r_s * delta * (r**2 - a2 * cos2_theta)) / (2.0d0 * rho2**3)
      ! ^r_rr (2.14.5e, left)
      christoffel(1, 1, 1) = (2.0d0 * r * a2 * sin2_theta - r_s * (r**2 - a2 * cos2_theta)) / (2.0d0 * rho2 * delta)
      ! ^r_rth (2.14.5f, left)
      christoffel(1, 1, 2) = -(a2 * sin_theta * cos_theta) / rho2
      christoffel(1, 2, 1) = christoffel(1, 1, 2)
      ! Γ^r_thth (2.14.5g, left)
      christoffel(1, 2, 2) = -(r * delta) / rho2
      ! Γ^r_tph (2.14.5d, left)
      christoffel(1, 0, 3) = -(delta * r_s * a * sin2_theta * (r**2 - a2 * cos2_theta)) / (2.0d0 * rho2**3)
      christoffel(1, 3, 0) = christoffel(1, 0, 3)
      ! Γ^r_phph (2.14.5k, left)
      christoffel(1, 3, 3) = (delta * sin2_theta / (2.0d0 * rho2**3)) * (-2.0d0 * r * rho2**2 + r_s * a2 *&
      sin2_theta * (r**2 - a2 * cos2_theta))

      ! --- ^th components ---
      ! ^th_tt (2.14.5a, right)
      christoffel(2, 0, 0) = -(r_s * a2 * r * sin_theta * cos_theta) / rho2**3
      ! ^th_rr (2.14.5e, right)
      christoffel(2, 1, 1) = (a2 * sin_theta * cos_theta) / (rho2 * delta)
      ! ^th_rth (2.14.5f, right)
      christoffel(2, 1, 2) = r / rho2
      christoffel(2, 2, 1) = christoffel(2, 1, 2)
      ! ^th_thth (2.14.5g, right)
      christoffel(2, 2, 2) = -(a2 * sin_theta * cos_theta) / rho2
      ! ^th_tph (2.14.5d, right)
      christoffel(2, 0, 3) = (r_s * a * r * (r**2 + a2) * sin_theta * cos_theta) / rho2**3
      christoffel(2, 3, 0) = christoffel(2, 0, 3)
      ! ^th_phph (2.14.5l, left)
      a_term = (r**2 + a2)**2 - delta * a2 * sin2_theta
      christoffel(2, 3, 3) = -(sin_theta * cos_theta / rho2**3) * (a_term * delta + (r**2 + a2) * r_s * a2 * r * sin2_theta)

      ! --- ^ph components ---
      ! ^ph_tr (2.14.5b, right)
      christoffel(3, 0, 1) = (r_s * a * (r**2 - a2 * cos2_theta)) / (2.0d0 * rho2**2 * delta)
      christoffel(3, 1, 0) = christoffel(3, 0, 1)
      ! ^ph_tth (2.14.5c, right)
      christoffel(3, 0, 2) = -(r_s * a * r * cos_theta) / (rho2**2 * sin_theta)
      christoffel(3, 2, 0) = christoffel(3, 0, 2)
      ! ^ph_rph (2.14.5j, left)
      christoffel(3, 1, 3) = (2.0d0 * r * rho2**2 + r_s * (a**4 * sin2_theta * cos2_theta - r**2 * (rho2 + r**2 + a2)))&
       / (2.0d0 * rho2**2 * delta)
      christoffel(3, 3, 1) = christoffel(3, 1, 3)
      ! ^ph_thph (2.14.5h, left)
      christoffel(3, 2, 3) = (cos_theta / (rho2**2 * sin_theta)) * (rho2**2 + r_s * a2 * r * sin2_theta)
      christoffel(3, 3, 2) = christoffel(3, 2, 3)

   end subroutine

   ! Geodesic equations: d^2 x^mu/dtau^2 + Γ^mu_nu,lambda (dx^nu/dtau)(dx^lambda/dtau) = 0
   subroutine kerr_geodesic_rhs(tau, y, dydt)
      real(kind = 8), intent(in) :: tau
      real(kind = 8), intent(in) :: y(:)      ! (t,r,th,ph,dt/dtau,dr/dτ\tau,dth/dtau,dph/dtau)
      real(kind = 8), intent(out) :: dydt(:)

      real(kind = 8) :: pos(0:3), vel(0:3)
      real(kind = 8) :: christoffel(0:3, 0:3, 0:3)
      integer :: mu, nu, lambda

      ! pos and vel
      pos(0:3) = y(1:4)   ! (t, r, th, ph)
      vel(0:3) = y(5:8)   ! (dt/dtau, dr/dtau, dth/dtau, dph/dtau)

   
      ! ! DEBUG: Print first call
      ! if (tau < 1.0d-6) then
      !    write(*, *) 'First geodesic call:'
      !    write(*, *) 'pos =', pos
      !    write(*, *) 'vel =', vel
      ! end if

      ! Get Christoffels at this pos
      call kerr_christoffel(pos(1), pos(2), christoffel)

      !
      ! ! DEBUG: Check specific Christoffel components for first call
      ! if (tau < 1.0d-6) then
      !    write(*, *) 'Key Christoffel symbols:'
      !    write(*, *) 'Γ^r_tt =', christoffel(1, 0, 0)
      !    write(*, *) 'Γ^r_φφ =', christoffel(1, 3, 3)
      !    write(*, *) 'Γ^r_tφ =', christoffel(1, 0, 3)
      !    write(*, *) 'Terms in d²r/dτ²:'
      !    write(*, *) '  -Γ^r_tt (u^t)² =', -christoffel(1, 0, 0) * vel(0)**2
      !    write(*, *) '  -Γ^r_φφ (u^φ)² =', -christoffel(1, 3, 3) * vel(3)**2
      !    write(*, *) '  -2Γ^r_tφ u^t u^φ =', -2.0d0 * christoffel(1, 0, 3) * vel(0) * vel(3)
      !    write(*, *) '  Total =', dydt(6)
      ! end if
      ! ! DEBUG: Check more Christoffel components for first call
      ! if (tau < 1.0d-6) then
      !    write(*, *) 'More Christoffel symbols:'
      !    write(*, *) 'Γ^t_tt =', christoffel(0, 0, 0)
      !    write(*, *) 'Γ^t_tr =', christoffel(0, 0, 1)
      !    write(*, *) 'Γ^t_tφ =', christoffel(0, 0, 3)
      !    write(*, *) 'Γ^φ_tt =', christoffel(3, 0, 0)
      !    write(*, *) 'Γ^φ_tr =', christoffel(3, 0, 1)
      !    write(*, *) 'Γ^φ_tφ =', christoffel(3, 0, 3)
      ! end if

      ! First half: dx^mu/dtau = v^mu
      dydt(1:4) = vel(0:3)

      ! Second half: d^2 x^mu/dtau^2 = -Γ^mu_nu,lambda (dx^nu/dtau)(dx^lamb/dtau)
      do mu = 0, 3
         dydt(5 + mu) = 0.0d0
         do nu = 0, 3
            do lambda = 0, 3
               dydt(5 + mu) = dydt(5 + mu) - christoffel(mu, nu, lambda) * &
               vel(nu) * vel(lambda)
            end do
         end do
      end do

      ! ! DEBUG: Print accelerations for first call
      ! ! In kerr_geodesic_rhs, add this after computing all accelerations:
      ! if (tau < 1.0d-6) then
      !    write(*, *) 'All accelerations:'
      !    write(*, *) 'd²t/dτ² =', dydt(5)
      !    write(*, *) 'd²r/dτ² =', dydt(6)
      !    write(*, *) 'd²θ/dτ² =', dydt(7)
      !    write(*, *) 'd²φ/dτ² =', dydt(8)
      !    write(*, *) 'Position after step should be:', y + dydt * 0.01d0
      ! end if
      ! 

   end subroutine

   ! Constants of motion for Kerr geodesics (for accuracy checks)
   subroutine kerr_constants_of_motion(r, theta, vel, energy, ang_mom, carter_const)
      real(kind = 8), intent(in) :: r, theta, vel(0:3)
      real(kind = 8), intent(out) :: energy, ang_mom, carter_const

      real(kind = 8) :: g_covar(0:3, 0:3)
      real(kind = 8) :: rho2, delta, a2

      call kerr_metric(r, theta, g_covar)

      ! Energy: E = -g_tt * dt/dtau - g_tph * dph/dtau
      energy = -(g_covar(0, 0) * vel(0) + g_covar(0, 3) * vel(3))

      ! Angular momentum: L = g_tph * dt/dtau + g_phph * dph/dtau
      ang_mom = g_covar(0, 3) * vel(0) + g_covar(3, 3) * vel(3)

      ! Carter constant (simplified? cant tell)
      rho2 = rho_squared(r, theta)
      delta = delta_kerr(r)
      a2 = spin_bh**2

      carter_const = rho2**2 * vel(2)**2 + cos(theta)**2 * &
      (ang_mom**2 / sin(theta)**2 - a2 * energy**2)

   end subroutine

   subroutine debug_circular_orbit(r_orbit)
      real(kind = 8), intent(in) :: r_orbit
      real(kind = 8) :: omega, energy_kerr, l_kerr
      real(kind = 8) :: sqrt_mr, sqrt_m, r32

      write(*, *) 'Debug Kerr circular orbit at r =', r_orbit

      sqrt_m = sqrt(mass_bh)
      sqrt_mr = sqrt(mass_bh * r_orbit)
      r32 = r_orbit**(1.5d0)

      ! Kerr circular orbit frequency (prograde)
      omega = sqrt_m / (r32 + spin_bh * sqrt_m)

      ! Kerr energy for circular orbit
      energy_kerr = (r_orbit**2 - 2.0d0 * mass_bh * r_orbit + spin_bh * sqrt_mr) / &
      (r_orbit * sqrt(r_orbit**2 - 3.0d0 * mass_bh * r_orbit + &
      2.0d0 * spin_bh * sqrt_mr))

      ! Kerr angular momentum for circular orbit
      l_kerr = sqrt_mr * (r_orbit**2 - 2.0d0 * spin_bh * sqrt_mr + spin_bh**2) / &
      (r32 - 2.0d0 * spin_bh * sqrt_mr + spin_bh * sqrt_m)

      write(*, *) 'Kerr Omega =', omega
      write(*, *) 'Kerr E =', energy_kerr
      write(*, *) 'Kerr L =', l_kerr
      write(*, *) 'Spin a =', spin_bh

   end subroutine

   subroutine set_circular_orbit_ic(r_orbit, y_init)
      real(kind = 8), intent(in) :: r_orbit
      real(kind = 8), intent(out) :: y_init(8)

      real(kind = 8) :: omega, ut, uphi
      real(kind = 8), parameter :: theta_eq = 1.570796326794897d0  ! π/2
      real(kind = 8) :: sqrt_m, r32, g_tt, g_tphi, g_phiphi
      real(kind = 8) :: det_metric_2x2, rho2

      sqrt_m = sqrt(mass_bh)
      r32 = r_orbit**(1.5d0)
      rho2 = r_orbit**2 + spin_bh**2  ! At θ=π/2, cos^2 θ=0

      ! Kerr circular orbit frequency (prograde)
      omega = sqrt_m / (r32 + spin_bh * sqrt_m)

      ! Kerr metric components at (r, θ=π/2)
      g_tt = -(1.0d0 - 2.0d0 * mass_bh * r_orbit / rho2)
      g_tphi = -2.0d0 * mass_bh * spin_bh * r_orbit / rho2
      g_phiphi = (r_orbit**2 + spin_bh**2 + &
      2.0d0 * mass_bh * spin_bh**2 * r_orbit / rho2)


      det_metric_2x2 = g_tt + 2.0d0 * g_tphi * omega + g_phiphi * omega**2
      ut = sqrt(-1.0d0 / det_metric_2x2)
      uphi = omega * ut

      ! Initial conditions: (t,r,θ,φ,dt/dτ,dr/dτ,dθ/dτ,dφ/dτ)
      y_init = [0.0d0, r_orbit, theta_eq, 0.0d0, ut, 0.0d0, 0.0d0, uphi]

      write(*, *) 'Set ut =', ut, ', uphi =', uphi
      write(*, *) 'Omega =', omega

   end subroutine

end module