program kerr_simulator
   use kerr_geometry
   use accretion_disk
   use ray_tracer
   use integrators
   implicit none

   ! --- sim params ---
   integer, parameter :: nx = 256
   integer, parameter :: ny = 256
   real(kind = 8), parameter :: m = 1.0d0
   real(kind = 8), parameter :: a = 0.5d0
   real(kind = 8), parameter :: r_cam = 75.0d0
   real(kind = 8), parameter :: theta_cam_deg = 75.0d0
   real(kind = 8), parameter :: phi_cam_deg = 0.0d0
   real(kind = 8), parameter :: fov_deg = 25.0d0

   ! Integration settings
   integer, parameter :: max_steps = 1000000
   real(kind = 8), parameter :: t_start = 0.0d0
   real(kind = 8), parameter :: h_fixed = -0.005d0 ! Fixed step size (negative for reverse extrapolation)

   ! --- Variable Declarations ---
   real(kind = 8) :: image(nx, ny)
   integer :: image_quantized(nx, ny)
   integer :: i, j, k
   real(kind = 8) :: x_screen, y_screen
   real(kind = 8) :: y(8), y_old(8), y_new(8)
   real(kind = 8) :: t, r_isco, r_horizon
   real(kind = 8), parameter :: pi = 4.0d0 * atan(1.0d0)
   real(kind = 8) :: theta_cam_rad, phi_cam_rad, alpha, r_cross, intensity_factor
   real(kind = 8), parameter :: disk_thickness = 0.1d0  ! Small thickness
   integer :: pixel_counts(0:9)

   ! real(kind = 8), parameter :: r_photon_pro = 2.6d0   ! Approximate prograde photon sphere
   ! real(kind = 8), parameter :: r_photon_retro = 4.2d0 ! Approximate retrograde photon sphere
   real(kind = 8) :: r_photon_exact_pro, r_photon_exact_retro
   real(kind = 8) :: r_avg, dr_dt_avg, orbital_motion

   ! --- init ---
   write(*, *) '=== Kerr Ray Tracer (Fixed Step) ==='
   write(*, *) 'Image size: ', nx, 'x', ny
   write(*, *) 'Spin a = ', a

   call set_bh_params(m, a)
   theta_cam_rad = theta_cam_deg * pi / 180.0d0
   phi_cam_rad = phi_cam_deg * pi / 180.0d0
   call setup_camera(r_cam, theta_cam_rad, phi_cam_rad, fov_deg)

   r_isco = isco_radius()
   r_horizon = mass_bh + sqrt(mass_bh**2 - spin_bh**2)

   call calculate_photon_spheres(a, r_photon_exact_pro, r_photon_exact_retro)

   write(*, *) 'Event Horizon at r = ', r_horizon
   write(*, *) 'ISCO at r = ', r_isco
   write(*, *) 'Prograde Photon Sphere at r = ', r_photon_exact_pro
   write(*, *) 'Retrograde Photon Sphere at r = ', r_photon_exact_retro
   write(*, *)

   image = 0
   image_quantized = 0

   ! --- Main Ray-Tracing Loop ---
   do j = 1, ny
      if (mod(j, 10) == 0) then
         write(*, '("Row ", I4, " / ", I4)') j, ny
      end if
      do i = 1, nx
         x_screen = (2.0d0 * i - nx - 1.0d0) / real(nx, 8)
         y_screen = (2.0d0 * j - ny - 1.0d0) / real(ny, 8)

         call get_photon_ic(x_screen, y_screen, y)

         ! integrate this ray using fixed-step RK4
         t = t_start




         ! Replace your integration loop:
         do k = 1, max_steps
            y_old = y
            call rk4_step(8, t, y, h_fixed, kerr_geodesic_rhs, y_new)
            y = y_new
            t = t + h_fixed

            ! Horizon check
            if (y(2) <= r_horizon) then
               image(i, j) = 0.0d0
               exit
            end if

            ! REFINED photon sphere detection

            if (k > 1000) then  ! Let orbit stabilize
               ! Check for near-circular orbit at photon sphere radii

               ! Average radius over last few steps
               r_avg = (y(2) + y_old(2)) * 0.5d0
               dr_dt_avg = abs(y(6))  ! |dr/dtau|
               orbital_motion = abs(y(4))  ! Total phi accumulated

               ! Criteria for photon sphere orbit:
               ! 1. Near photon sphere radius
               ! 2. Small radial velocity (nearly circular)
               ! 3. Significant orbital motion
               ! 4. Near equatorial plane

               if (orbital_motion > 4.0d0 * pi .and. dr_dt_avg < 0.02d0 .and. &
               abs(y(3) - pi / 2.0d0) < 0.2d0) then

                  if (abs(r_avg - r_photon_exact_pro) < 0.3d0) then
                     image(i, j) = 4.0d0  ! Prograde photon sphere
                     image_quantized(i, j) = 9
                     exit
                  else if (abs(r_avg - r_photon_exact_retro) < 0.3d0) then
                     image(i, j) = 3.5d0  ! Retrograde photon sphere  
                     image_quantized(i, j) = 8
                     exit
                  end if
               end if
            end if

            ! Disk intersection (unchanged)
            if ((y_old(3) - pi / 2.0d0) * (y(3) - pi / 2.0d0) < 0.0d0) then
               alpha = abs(y_old(3) - pi / 2.0d0) / abs(y(3) - y_old(3))
               r_cross = y_old(2) + alpha * (y(2) - y_old(2))

               if (r_cross >= r_isco .and. r_cross <= 30.0d0) then
                  intensity_factor = (r_isco / r_cross)**(0.75d0)
                  image_quantized(i, j) = max(1, min(7, int(intensity_factor * 5.0d0)))
                  image(i, j) = intensity_factor * 5.0d0
               else
                  image(i, j) = 0.0d0
               end if
               exit
            end if

            ! Escape check
            if (y(2) > r_cam * 3.0d0) then
               image(i, j) = 0.0d0
               exit
            end if
         end do
      end do
   end do

   pixel_counts = 0
   do j = 1, ny
      do i = 1, nx
         pixel_counts(image_quantized(i, j)) = pixel_counts(image_quantized(i, j)) + 1
      end do
   end do

   write(*, *) 'Pixel distribution:'
   do i = 0, 9
      if (pixel_counts(i) > 0) then
         write(*, '("  Value ",I0,": ",I0," pixels (",F5.1,"%)")') &
         i, pixel_counts(i), 100.0 * pixel_counts(i) / real(nx * ny)
      end if
   end do

   ! --- Write Image to File ---
   call write_image_file('image.dat', image, nx, ny)
   write(*, *)
   write(*, *) 'Ray tracing complete. Image data written to image.dat'

   contains

   subroutine write_image_file(filename, img_data, w, h)
      character(len = *), intent(in) :: filename
      integer, intent(in) :: w, h
      real(kind = 8), intent(in) :: img_data(w, h)
      integer :: unit_num, i, j

      open(newunit = unit_num, file = filename, status = 'replace')
      write(unit_num, '(A,2I6)') '# nx ny', w, h
      do j = h, 1, -1  ! Flip y-axis for proper image orientation
         do i = 1, w
            write(unit_num, '(F6.4,1X)', advance = 'no') img_data(i, j)
         end do
         write(unit_num, *)  ! New line
      end do
      close(unit_num)
   end subroutine

end program