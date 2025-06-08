program kerr_simulator
    use kerr_geometry
    use trajectory
    implicit none
    
    real(kind=8), parameter :: pi = 4.0d0*atan(1.0d0)

    ! Sim params
    real(kind=8), parameter :: M = 1.0d0      ! BH mass
    real(kind=8), parameter :: a = 0.9d0      ! Spin
    integer, parameter :: n_steps = 10000
    logical, parameter :: adaptive = .false.  ! Use adaptive step size
    
    ! Initial conditions
    real(kind=8) :: y_init(8)
    real(kind=8) :: r_orbit = 6.0d0           ! Init radius
    real(kind=8) :: inclination = pi/3.0d0    ! 60 degrees inc
    
    
    ! Integration arrays
    real(kind=8) :: trajectory(8, n_steps+1)
    real(kind=8) :: times(n_steps+1)
    real(kind=8) :: t_start = 0.0d0, t_end = 100.0d0
    
    ! Constants of motion (for double checking)
    ! real(kind=8) :: energy, ang_mom, carter_const
    real(kind=8) :: E_init, L_init, Q_init, E_final, L_final, Q_final
    
    write(*,*) '=== Kerr Black Hole Geodesic Simulator ==='
    write(*,*) 'Black hole mass M =', M
    write(*,*) 'Spin parameter a =', a
    write(*,*) 'Initial orbit radius =', r_orbit
    write(*,*) 'Inclination =', inclination * 180.0d0/pi, 'degrees'
    write(*,*) 'Adaptive Integration = ', adaptive
    write(*,*)
    
    ! Setup BH
    call set_bh_params(M, a)
    
    call debug_circular_orbit(r_orbit)
    ! Setup circular orbit inits 
    call set_circular_orbit_ic(r_orbit, y_init)

    write(*,*) 'Starting equatorial circular orbit at r =', r_orbit
    ! Set up circular orbit initial conditions
    call set_circular_orbit_ic(r_orbit, y_init)
    
    ! DEBUG: Print initial conditions
    write(*,*) 'Initial conditions y_init:'
    write(*,*) 't0 =', y_init(1)
    write(*,*) 'r0 =', y_init(2) 
    write(*,*) 'theta0 =', y_init(3)
    write(*,*) 'phi0 =', y_init(4)
    write(*,*) 'ut =', y_init(5)
    write(*,*) 'ur =', y_init(6)
    write(*,*) 'utheta =', y_init(7)
    write(*,*) 'uphi =', y_init(8)
    write(*,*)
    
    ! Check init constants of motion
    call kerr_constants_of_motion(y_init(2), y_init(3), y_init(5:8), &
                                 E_init, L_init, Q_init)
    
    write(*,*) 'Initial constants of motion:'
    write(*,*) 'Energy E =', E_init
    write(*,*) 'Angular momentum L =', L_init  
    write(*,*) 'Carter constant Q =', Q_init
    write(*,*)
    
    ! Integrate geodesic
    write(*,*) 'Integrating geodesic...'
    call integrate_geodesic(y_init, t_start, t_end, n_steps, &
                           trajectory, times, adaptive=adaptive)
    
    ! Check final constants of motion
    call kerr_constants_of_motion(trajectory(2, n_steps+1), trajectory(3, n_steps+1), &
                                 trajectory(5:8, n_steps+1), E_final, L_final, Q_final)
    
    write(*,*) 'Final constants of motion:'
    write(*,*) 'Energy E =', E_final, '  (error:', abs(E_final - E_init), ')'
    write(*,*) 'Angular momentum L =', L_final, '  (error:', abs(L_final - L_init), ')'
    write(*,*) 'Carter constant Q =', Q_final, '  (error:', abs(Q_final - Q_init), ')'
    write(*,*)
    
    ! Write trajectory to file
    call write_trajectory('trajectory.dat', trajectory, times, n_steps+1)
    write(*,*) 'Trajectory written to trajectory.dat'
    write(*,*) 'Done!'
    
end program