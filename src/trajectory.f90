module trajectory
    use kerr_geometry
    use integrators
    implicit none
    private
    
    public :: integrate_geodesic, write_trajectory
    
contains
    
    subroutine integrate_geodesic(y_init, t_start, t_end, n_steps, &
                                 trajectory_out, times_out, adaptive)
        real(kind=8), intent(in) :: y_init(8), t_start, t_end
        integer, intent(in) :: n_steps
        real(kind=8), intent(out) :: trajectory_out(8, n_steps+1)
        real(kind=8), intent(out) :: times_out(n_steps+1)
        logical, intent(in), optional :: adaptive
        
        real(kind=8) :: t, h, h_try, h_did, h_next
        real(kind=8) :: y_current(8), y_new(8)
        real(kind=8), parameter :: eps = 1.0d-8  ! Error tolerance
        logical :: use_adaptive
        integer :: i
        
        ! Default to adaptive
        use_adaptive = .true.
        if (present(adaptive)) use_adaptive = adaptive
        
        ! Initialize
        t = t_start
        y_current = y_init
        h_try = (t_end - t_start) / real(n_steps, 8)
        
        ! Store init conds
        trajectory_out(:, 1) = y_current
        times_out(1) = t
        
        ! Integrate using RK4 or adaptive RK4
        do i = 2, n_steps + 1
            if (use_adaptive) then
                call rk4_adaptive(8, t, y_current, h_try, eps, h_did, h_next, &
                                 kerr_geodesic_rhs, y_new)
                h_try = h_next
            else
                h = (t_end - t_start) / real(n_steps, 8)
                call rk4_step(8, t, y_current, h, kerr_geodesic_rhs, y_new)
                h_did = h
            end if
            
            ! Update
            t = t + h_did
            y_current = y_new
            
            ! Store
            trajectory_out(:, i) = y_current
            times_out(i) = t
            
            ! Check for horizon crossing (r < r_+)
            if (y_current(2) < horizon_radius()) then
                write(*,*) 'Horizon crossed at step', i
                exit
            end if
            
            ! Stop if we've reached end time
            if (t >= t_end) exit
        end do
        
    end subroutine
    
    function horizon_radius()
        real(kind=8) :: horizon_radius
        horizon_radius = mass_bh + sqrt(mass_bh**2 - spin_bh**2)
    end function
    
    subroutine write_trajectory(filename, trajectory, times, n_points)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n_points
        real(kind=8), intent(in) :: trajectory(8, n_points)
        real(kind=8), intent(in) :: times(n_points)
        
        integer :: i, unit_num
        
        open(newunit=unit_num, file=filename, status='replace')
        write(unit_num, '(A)') '# time, t, r, theta, phi, dt/dtau, dr/dtau, dtheta/dtau, dphi/dtau'
        
        do i = 1, n_points
            write(unit_num, '(9ES15.7)') times(i), trajectory(:, i)
        end do
        
        close(unit_num)
        
    end subroutine
    
end module