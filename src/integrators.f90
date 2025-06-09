module integrators
   ! We need mass_bh from the geometry module to calculate the step size.
   use kerr_geometry, only: mass_bh
   implicit none
   private

   public :: rk4_step, get_adaptive_step

   contains

   ! 4th-order R-K step (Unchanged)
   subroutine rk4_step(n, t, y, h, rhs, y_new)
      integer, intent(in) :: n
      real(kind = 8), intent(in) :: t, y(n), h
      real(kind = 8), intent(out) :: y_new(n)

      interface
         subroutine rhs(t_in, y_in, dydt)
            real(kind = 8), intent(in) :: t_in, y_in(:)
            real(kind = 8), intent(out) :: dydt(:)
         end subroutine
      end interface

      real(kind = 8) :: k1(n), k2(n), k3(n), k4(n), y_temp(n)

      call rhs(t, y, k1)
      k1 = h * k1

      y_temp = y + 0.5d0 * k1
      call rhs(t + 0.5d0 * h, y_temp, k2)
      k2 = h * k2

      y_temp = y + 0.5d0 * k2
      call rhs(t + 0.5d0 * h, y_temp, k3)
      k3 = h * k3

      y_temp = y + k3
      call rhs(t + h, y_temp, k4)
      k4 = h * k4

      y_new = y + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0

   end subroutine

   function get_adaptive_step(y, r_horizon) result(h)
      real(kind = 8), intent(in) :: r_horizon, y(8)

      real(kind = 8) :: h, r, theta, ur, utheta, uphi
      real(kind = 8) :: r_far, r_near, h_base, speed_factor

      ! Tuning Parameters
      real(kind = 8), parameter :: h_max = 2.0d0
      real(kind = 8), parameter :: h_min = 0.001d0
      real(kind = 8), parameter :: h_scale = 0.05d0

      r = y(2)
      theta = y(3)
      ur = y(6)
      utheta = y(7)
      uphi = y(8)
      ! Calculate region boundaries at runtime
      r_far = 50.0d0 * mass_bh
      r_near = 4.0d0 * mass_bh

      ! 1. Determine the base step size from radial position (as before)
      if (r > r_far .and. ur > 0) then
         h_base = min(h_max, 0.001d0 * r**2)
      else if (r > r_near) then
         h_base = h_min + h_scale * (r - r_near)
      else
         h_base = h_min
      end if

      ! 2. Calculate a "speed factor" to account for angular velocity.
      !    This shrinks the step size during periods of rapid angular motion.
      !    The 1/r term accounts for the fact that angular velocity is
      !    physically larger for a given 'u' at smaller radii.
      speed_factor = 1.0d0 + (abs(utheta) + abs(uphi) / max(1.d-6, sin(theta))) / r

      ! 3. Combine the two to get the final step size.
      h = h_base / speed_factor

      ! 4. Enforce bounds and direction.
      h = max(h_min, min(h_max, h)) * 0.5d0  ! Scale down for stability
      h = -h
   end function

end module