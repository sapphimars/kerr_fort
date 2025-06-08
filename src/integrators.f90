module integrators
    implicit none
    private
    
    public :: rk4_step, rk4_adaptive
    
contains
    
    ! 4th-order R-K step
    subroutine rk4_step(n, t, y, h, rhs, y_new)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: t, y(n), h
        real(kind=8), intent(out) :: y_new(n)
        
        interface
            subroutine rhs(t_in, y_in, dydt)
                real(kind=8), intent(in) :: t_in, y_in(:)
                real(kind=8), intent(out) :: dydt(:)
            end subroutine
        end interface
        
        real(kind=8) :: k1(n), k2(n), k3(n), k4(n), y_temp(n)
        
        ! k1 = h * f(t, y)
        call rhs(t, y, k1)
        k1 = h * k1
        
        ! k2 = h * f(t + h/2, y + k1/2)
        y_temp = y + 0.5d0 * k1
        call rhs(t + 0.5d0*h, y_temp, k2)
        k2 = h * k2
        
        ! k3 = h * f(t + h/2, y + k2/2)
        y_temp = y + 0.5d0 * k2
        call rhs(t + 0.5d0*h, y_temp, k3)
        k3 = h * k3
        
        ! k4 = h * f(t + h, y + k3)
        y_temp = y + k3
        call rhs(t + h, y_temp, k4)
        k4 = h * k4
        
        ! Final result
        y_new = y + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4) / 6.0d0
        
    end subroutine
    
    ! Adaptive RK4 with error control
    subroutine rk4_adaptive(n, t, y, h_try, eps, h_did, h_next, rhs, y_new)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: t, y(n), h_try, eps
        real(kind=8), intent(out) :: h_did, h_next, y_new(n)
        
        interface
            subroutine rhs(t_in, y_in, dydt)
                real(kind=8), intent(in) :: t_in, y_in(:)
                real(kind=8), intent(out) :: dydt(:)
            end subroutine
        end interface
        
        real(kind=8) :: y_full(n), y_half1(n), y_half2(n)
        real(kind=8) :: error, error_max, h_temp
        real(kind=8), parameter :: safety = 0.9d0, p_grow = -0.2d0, p_shrink = -0.25d0
        real(kind=8), parameter :: err_con = 1.89d-4  ! (5/safety)^(1/p_grow)
        integer, parameter :: max_tries = 10
        integer :: i, j
        
        h_temp = h_try
        
        do i = 1, max_tries
            ! One full step
            call rk4_step(n, t, y, h_temp, rhs, y_full)
            
            ! Two half steps
            call rk4_step(n, t, y, 0.5d0*h_temp, rhs, y_half1)
            call rk4_step(n, t + 0.5d0*h_temp, y_half1, 0.5d0*h_temp, rhs, y_half2)
            
            ! Error estimate (Richardson extrapolation)
            error_max = 0.0d0
            do j = 1, n
                error = abs(y_half2(j) - y_full(j)) / 15.0d0  ! Truncation error estimate
                if (abs(y(j)) > 1.0d-10) then
                    error = error / abs(y(j))  ! Relative error
                end if
                error_max = max(error_max, error)
            end do
            
            if (error_max <= eps) then
                ! Success - use higher order result
                y_new = y_half2 + (y_half2 - y_full) / 15.0d0  ! Richardson extrapolation
                h_did = h_temp
                
                if (error_max > err_con) then
                    h_next = safety * h_temp * (error_max**p_grow)
                else
                    h_next = 5.0d0 * h_temp  ! No more than 5x increase
                end if
                return
            else
                ! Reduce stepsize, but not too much
                h_temp = safety * h_temp * (error_max**p_shrink)
                if (h_temp < 0.1d0 * h_try) then
                    h_temp = 0.1d0 * h_try 
                end if
            end if
        end do
        
        ! If we get here, we failed to converge :( sad!
        write(*,*) 'rk4_adaptive: too many step reductions'
        h_did = h_temp
        h_next = h_temp
        y_new = y_full
        
    end subroutine
    
end module