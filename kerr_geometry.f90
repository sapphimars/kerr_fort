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
   kerr_metric_determinant, kerr_metric_det_closed, kerr_christoffel

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

      real(kind = 8) :: rho2, delta, sin_theta, cos_theta, sin2_theta
      real(kind = 8) :: m, a, a2
      real(kind = 8) :: rho2_inv, delta_inv

      ! Initialize to zero
      christoffel = 0.0d0

      ! Shortcuts
      m = mass_bh
      a = spin_bh
      a2 = a**2
      sin_theta = sin(theta)
      cos_theta = cos(theta)
      sin2_theta = sin_theta**2

      rho2 = r**2 + a2 * cos_theta**2
      delta = r**2 - 2.0d0 * m * r + a2
      rho2_inv = 1.0d0 / rho2
      delta_inv = 1.0d0 / delta

      ! Non-zero Christoffel symbols (only the essential ones)

      ! ^t components
      christoffel(0, 0, 1) = m * (r**2 - a2 * cos_theta**2) * delta_inv / rho2**2  ! ^t_tr
      christoffel(0, 1, 0) = christoffel(0, 0, 1)                            ! ^t_rt
      christoffel(0, 0, 2) = -2.0d0 * m * a2 * r * sin_theta * cos_theta / rho2**2     ! ^t_tth  
      christoffel(0, 2, 0) = christoffel(0, 0, 2)                            ! ^t_tht
      christoffel(0, 1, 3) = m * a * sin2_theta * (r**2 - a2 * cos_theta**2) * delta_inv / rho2**2  ! ^t_rph
      christoffel(0, 3, 1) = christoffel(0, 1, 3)                            ! ^t_phr
      christoffel(0, 2, 3) = -2.0d0 * m * a**3 * r * sin_theta**3 * cos_theta / rho2**2  ! ^t_thph
      christoffel(0, 3, 2) = christoffel(0, 2, 3)                            ! ^t_phth

      ! ^r components
      christoffel(1, 0, 0) = m * delta * (r**2 - a2 * cos_theta**2) / rho2**3      ! ^r_tt
      christoffel(1, 1, 1) = (r - m) * delta_inv + m / rho2                    ! ^r_rr
      christoffel(1, 2, 2) = -r * delta_inv                                  ! ^r_thth
      christoffel(1, 3, 3) = -r * sin2_theta * delta_inv + &
      m * a2 * sin2_theta * (2.0d0 * r**2 - a2 * cos_theta**2) / rho2**3  ! ^r_phph
      christoffel(1, 0, 3) = -m * a * sin2_theta * delta * (r**2 - a2 * cos_theta**2) / rho2**3  ! ^r_tph
      christoffel(1, 3, 0) = christoffel(1, 0, 3)                            ! ^r_pht

      ! ^th components
      christoffel(2, 1, 2) = r * rho2_inv                                    ! ^th_rth
      christoffel(2, 2, 1) = christoffel(2, 1, 2)                           ! ^th_thr  
      christoffel(2, 0, 0) = -2.0d0 * m * a2 * r * sin_theta * cos_theta / rho2**3    ! ^th_tt
      christoffel(2, 3, 3) = -sin_theta * cos_theta * (1.0d0 + &
      2.0d0 * m * r * (r**2 + a2) / rho2**3)               ! ^th_phph
      christoffel(2, 0, 3) = 2.0d0 * m * a**3 * r * sin_theta**3 * cos_theta / rho2**3  ! ^th_tph
      christoffel(2, 3, 0) = christoffel(2, 0, 3)                           ! ^th_pht

      ! ^ph components
      christoffel(3, 0, 1) = m * a * (r**2 - a2 * cos_theta**2) / (rho2**2 * delta)  ! ^ph_tr
      christoffel(3, 1, 0) = christoffel(3, 0, 1)                            ! ^ph_rt
      christoffel(3, 2, 3) = cos_theta / (sin_theta * rho2) * &
      (r**2 + a2 + 2.0d0 * m * a2 * r * sin2_theta / rho2)    ! ^ph_thph
      christoffel(3, 3, 2) = christoffel(3, 2, 3)                            ! ^ph_phth
      christoffel(3, 0, 0) = -2.0d0 * m * a * r * (r**2 - a2 * cos_theta**2) / rho2**3  ! ^ph_tt
      christoffel(3, 1, 3) = (r * delta_inv - m * a2 * sin2_theta / rho2) / rho2     ! ^ph_rph
      christoffel(3, 3, 1) = christoffel(3, 1, 3)                            ! ^ph_phr

   end subroutine

end module