program kerr_simulator
   use kerr_geometry
   use accretion_disk
   use ray_tracer
   use integrators
   implicit none

   ! --- Parameters ---
   integer, parameter :: nx = 768
   integer, parameter :: ny = 768
   real(kind = 8), parameter :: m = 1.0d0
   real(kind = 8), parameter :: a = 0.0d0
   real(kind = 8), parameter :: r_cam = 85.0d0
   real(kind = 8), parameter :: theta_cam_deg = 90.0d0
   real(kind = 8), parameter :: fov_deg = 15.0d0
   integer, parameter :: max_steps = 10000000


   ! --- Variables ---
   real(kind = 8) :: image(nx, ny)
   integer :: i, j, k
   real(kind = 8) :: x_screen, y_screen, y(8), y_old(8), y_new(8), t, h
   real(kind = 8) :: r_isco, r_horizon
   real(kind = 8), parameter :: pi = 4.0d0 * atan(1.0d0)
   real(kind = 8) :: theta_cam_rad, phi_cam_rad
   real(kind = 8) :: alpha, r_cross, intensity_factor, doppler_boost, g_factor
   real(kind = 8) :: u_disk(0:3), u_photon(0:3), g_covar(0:3, 0:3)

   ! --- Initialization ---
   write(*, *) '=== Kerr Ray Tracer (Thin Disk + Doppler - Corrected) ==='
   call set_bh_params(m, a)
   theta_cam_rad = theta_cam_deg * pi / 180.0d0
   phi_cam_rad = 0.0d0
   call setup_camera(r_cam, theta_cam_rad, phi_cam_rad, fov_deg)
   r_isco = isco_radius()
   r_horizon = mass_bh + sqrt(mass_bh**2 - spin_bh**2)
   write(*, *) 'Event Horizon:', r_horizon, ' ISCO:', r_isco
   image = 0.0d0
   write(*, *) 'r_isco =', r_isco
   write(*, *) 'r_horizon =', r_horizon
   write(*, *) 'Disk range: [', r_isco, ', 30.0]'


   ! --- Main Ray-Tracing Loop ---
   do j = 1, ny
      if (mod(j, 10) == 0) then
         write(*, '("Row ", I4, " / ", I4)') j, ny
      end if
      do i = 1, nx


         call get_photon_ic((2.d0 * i - nx - 1.d0) / real(nx, 8), &
         (2.d0 * j - ny - 1.d0) / real(ny, 8), y)
         t = 0.0d0
         do k = 1, max_steps
            h = get_adaptive_step(y, r_horizon)
            call rk4_step(8, t, y, h, kerr_geodesic_rhs, y_new)
            y_old = y
            y = y_new
            t = t + h

            if (y(2) < r_horizon) then
               image(i, j) = 0.0d0
               exit
            end if

            if ((abs(y(3) - pi / 2.d0) < 1.0d-6) .and. (y(2) >= r_isco .and. y(2) <= 30.0d0)) then
               r_cross = y(2)
               call get_disk_velocity(r_cross, a, u_disk)
               u_photon(0:3) = y(5:8)
               call kerr_metric(r_cross, pi / 2.d0, g_covar)
               g_factor = g_covar(0, 0) * u_disk(0) * u_photon(0) + &
               g_covar(3, 3) * u_disk(3) * u_photon(3) + &
               g_covar(0, 3) * (u_disk(0) * u_photon(3) + u_disk(3) * u_photon(0))
               doppler_boost = g_factor**(-4)
               intensity_factor = (r_isco / r_cross)**(0.75d0)
               image(i, j) = intensity_factor * doppler_boost
               exit
            end if

            if ((y_old(3) - pi / 2.d0) * (y(3) - pi / 2.d0) < 0.0d0) then
               alpha = abs(y_old(3) - pi / 2.d0) / abs(y(3) - y_old(3))
               r_cross = y_old(2) + alpha * (y(2) - y_old(2))

               if (r_cross >= r_isco .and. r_cross <= 30.0d0) then
                  call get_disk_velocity(r_cross, a, u_disk)
                  u_photon(0:3) = y(5:8)
                  call kerr_metric(r_cross, pi / 2.d0, g_covar)
                  g_factor = g_covar(0, 0) * u_disk(0) * u_photon(0) + &
                  g_covar(3, 3) * u_disk(3) * u_photon(3) + &
                  g_covar(0, 3) * (u_disk(0) * u_photon(3) + u_disk(3) * u_photon(0))
                  doppler_boost = g_factor**(-4)
                  intensity_factor = (r_isco / r_cross)**(0.75d0)
                  image(i, j) = intensity_factor * doppler_boost
               else
                  image(i, j) = 0.0d0
               end if
               exit
            end if

            if (y(2) > r_cam * 5.0d0) then
               image(i, j) = 0.0d0
               exit
            end if
         end do
      end do
   end do

   write(*, *) 'Min/Max intensity:', minval(image), maxval(image)
   call write_image_file('image.dat', image, nx, ny)
   write(*, *) 'Ray tracing complete.'

   contains

   ! CORRECTED: Calculates the 4-velocity for a circular prograde orbit
   subroutine get_disk_velocity(r, a, u_disk)
      real(kind = 8), intent(in) :: r, a
      real(kind = 8), intent(out) :: u_disk(0:3)
      real(kind = 8) :: omega, ut, g_tt, g_tphi, g_phiphi, rho2

      ! Angular velocity of a circular orbit in Kerr metric
      omega = sqrt(mass_bh) / (r**(1.5d0) + a * sqrt(mass_bh))

      ! Metric components at the equator (theta=pi/2)
      rho2 = r**2 ! cos(theta) is 0
      g_tt = -(1.0d0 - 2.0d0 * mass_bh * r / rho2)
      g_tphi = -2.0d0 * mass_bh * a * r / rho2
      g_phiphi = (r**2 + a**2 + 2.d0 * mass_bh * a**2 * r / rho2)

      ! Normalize from u.u = -1
      ut = 1.0d0 / sqrt(-(g_tt + 2.d0 * g_tphi * omega + g_phiphi * omega**2))

      u_disk(0) = ut
      u_disk(1) = 0.0d0
      u_disk(2) = 0.0d0
      u_disk(3) = omega * ut
   end subroutine

   ! CORRECTED: Uses robust ES format specifier
   subroutine write_image_file(filename, img_data, w, h)
      character(len = *), intent(in) :: filename
      integer, intent(in) :: w, h
      real(kind = 8), intent(in) :: img_data(w, h)
      integer :: unit_num, i, j
      open(newunit = unit_num, file = filename, status = 'replace')
      write(unit_num, '(A,2I6)') '# nx ny', w, h
      do j = 1, h
         do i = 1, w
            write(unit_num, '(G0.15, 1X)', advance = 'no') img_data(i, j) ! '(ES23.15,1X)'
         end do
         write(unit_num, *)
      end do
      close(unit_num)
   end subroutine

end program