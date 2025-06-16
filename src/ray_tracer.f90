module ray_tracer
   use kerr_geometry, only: mass_bh, spin_bh
   implicit none
   private

   public :: setup_camera, get_photon_ic, calculate_photon_spheres

   ! cam params
   real(kind = 8) :: r_cam, theta_cam, phi_cam
   real(kind = 8) :: fov_rad

   contains

   subroutine setup_camera(r_in, theta_in, phi_in, fov_deg)
      real(kind = 8), intent(in) :: r_in, theta_in, phi_in, fov_deg
      r_cam = r_in
      theta_cam = theta_in
      phi_cam = phi_in
      fov_rad = fov_deg * (4.0d0 * atan(1.0d0)) / 180.0d0
   end subroutine

   ! calc init conds using a ZAMO tetrad
   subroutine get_photon_ic(x_screen, y_screen, y_init)
      real(kind = 8), intent(in) :: x_screen, y_screen
      real(kind = 8), intent(out) :: y_init(8)

      ! local (ZAMO) frame quantities are hatted
      real(kind = 8) :: k_hat(0:3)
      ! camera setup
      real(kind = 8) :: cam_right(3), cam_up(3), cam_forward(3)
      real(kind = 8) :: sin_theta, cos_theta, sin_phi, cos_phi
      real(kind = 8) :: dir_camera(3), dir_global(3)
      real(kind = 8) :: r_comp, theta_comp, phi_comp, r_mag, xy_mag
      ! glob (Boyer-Lindquist) frame quantities are unhatted
      real(kind = 8) :: k(0:3)
      ! ZAMO tetrad transformation matrix components e^mu_(alpha)
      real(kind = 8) :: e_t_0, e_phi_0, e_r_1, e_theta_2, e_phi_3

      ! helper quantities from reference (Eq. 2.14.8)
      real(kind = 8) :: r, theta, a, m, r_s
      real(kind = 8) :: rho2, delta, a_term, omega
      real(kind = 8) :: alpha, beta ! Local sky angles
      real(kind = 8) :: g_tt, g_rr, g_thth, g_pp, g_tp, norm_check
      real(kind = 8) :: impact_param
      real(kind = 8) :: f_schw
      real(kind = 8) :: k_spatial


      k_hat = 0.0d0
      k = 0.0d0

      ! --- init Position ---
      y_init(1) = 0.0d0
      y_init(2) = r_cam
      y_init(3) = theta_cam
      y_init(4) = phi_cam

      ! --- def Photon Momentum in the Local ZAMO Frame ---
      ! convert screen coordinates to local sky angles.
      alpha = -x_screen * fov_rad / 2.0d0
      beta = y_screen * fov_rad / 2.0d0


      ! define mom in local frame using these angles
      ! ensures null vector
      k_hat(0) = 1.0d0 ! Photon energy in the local frame is 1
      k_hat(1) = cos(beta) * cos(alpha) ! k_r_hat (radial, points to origin)
      k_hat(2) = sin(beta)              ! k_theta_hat (polar)
      k_hat(3) = cos(beta) * sin(alpha) ! k_phi_hat (azimuthal)

      ! --- ZAMO Tetrad Transformation Matrix ---
      ! transforms vectors from the ZAMO frame to the global BL frame.
      ! These formulas are a direct transcription of Eq. 2.14.8. (i think)

      k_spatial = sqrt(k_hat(1)**2 + k_hat(2)**2 + k_hat(3)**2)
      k_hat(1:3) = k_hat(1:3) / k_spatial ! Normalize spatial part
      k_hat(0) = 1.0d0                    ! Energy = 1 for unit spatial vector

      r = r_cam
      theta = theta_cam
      a = spin_bh
      m = mass_bh
      r_s = 2.0d0 * m

      rho2 = r**2 + a**2 * cos(theta)**2
      delta = r**2 - r_s * r + a**2
      a_term = (r**2 + a**2)**2 - delta * a**2 * sin(theta)**2
      omega = (r_s * a * r) / a_term

      ! e_t_0 = sqrt(a_term / (rho2 * delta))
      ! e_phi_0 = omega * e_t_0
      ! e_r_1 = sqrt(delta / rho2)
      ! e_theta_2 = 1.0d0 / sqrt(rho2)
      ! e_phi_3 = sqrt(rho2 / a_term) / sin(theta)
      ! ! --- trans the local mom to glob coords ---
      ! ! k^mu = e^mu_(alpha) * k^(alpha)
      ! k(0) = e_t_0 * k_hat(0)
      ! k(1) = e_r_1 * k_hat(1)
      ! k(2) = e_theta_2 * k_hat(2)
      ! k(3) = e_phi_0 * k_hat(0) + e_phi_3 * k_hat(3)


      f_schw = 1.0d0 - 2.0d0 * mass_bh / r_cam

      ! Static tetrad (much simpler for a=0)
      e_t_0 = sqrt(f_schw)                      ! √(1-rs/r)
      e_phi_0 = 0.0d0                           ! No rotation
      e_r_1 = 1.0d0 / sqrt(f_schw)              ! 1/√(1-rs/r)  
      e_theta_2 = 1.0d0 / r_cam                 ! 1/r
      e_phi_3 = 1.0d0 / (r_cam * sin(theta_cam)) ! 1/(r sin θ)

      ! Transformation
      k(0) = e_t_0 * k_hat(0)
      k(1) = e_r_1 * k_hat(1)
      k(2) = e_theta_2 * k_hat(2)
      k(3) = e_phi_3 * k_hat(3)


      ! DEBUG
      impact_param = sqrt(k(2)**2 * rho2**2 + (k(3) / sin(theta))**2 * rho2**2)
      if (x_screen == 0.0d0 .and. y_screen == 0.0d0) then
         write(*, *) 'Central ray impact parameter:', impact_param
      end if



      g_tt = -(1.0d0 - r_s * r / rho2)
      g_rr = rho2 / delta
      g_thth = rho2
      g_pp = sin(theta)**2 * a_term / rho2
      g_tp = -r_s * a * r * sin(theta)**2 / rho2

      norm_check = g_tt * k(0) * k(0) + g_rr * k(1) * k(1) + g_thth * k(2) * k(2) + &
      g_pp * k(3) * k(3) + 2.0d0 * g_tp * k(0) * k(3)

      if (abs(norm_check) > 1.0d-1) then
         write(*, *) 'Warning: photon not null, norm =', norm_check
      end if

      ! set final initial 4-velocity for the integrator
      y_init(5:8) = k(0:3)

   end subroutine

   subroutine calculate_photon_spheres(a, r_pro, r_retro)
      real(kind = 8), intent(in) :: a
      real(kind = 8), intent(out) :: r_pro, r_retro
      real(kind = 8) :: rs, z_pro, z_retro

      rs = 2.0d0 * mass_bh

      ! For prograde (+a)
      z_pro = cos(2.0d0 / 3.0d0 * acos(a / mass_bh))
      r_pro = rs * (1.0d0 + z_pro)

      ! For retrograde (-a)
      z_retro = cos(2.0d0 / 3.0d0 * acos(-a / mass_bh))
      r_retro = rs * (1.0d0 + z_retro)
   end subroutine

end module