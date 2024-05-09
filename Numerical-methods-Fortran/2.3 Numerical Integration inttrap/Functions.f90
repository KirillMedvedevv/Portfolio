module functions
    use omp_lib
    implicit none
    contains
    function f(x) result(y)
        real(8), intent(in) :: x
        real(8) :: y
        y = ((((x * x) - 3.5d0) * x + 2.5d0) * x - 7d0) * x - 6.4d0
    end function
    
    
    function F1(x) result(y)
        real(8), intent(in) :: x
        real(8) :: y
        y = (((((1d0/6d0) * x * x -  3.5d0/4d0) * x + 2.5d0/3d0) * x - 7d0/2d0) * x - 6.4d0) * x
    end function
    
    
    function trap(a, b, N, h) result(int)
        integer, intent(in) :: N
        integer :: i
        real(8) :: h, int, a, b
        int = 0d0
        !$OMP PARALLEL
        !$OMP DO
        do i = 2, N
            int = int + f(a+real(i-1)*h)
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        int = (2d0*int + f(a) + f(b))*h/2d0
        
        
    end function
    
    
    subroutine inttrap(a, b, eps, err, h, N)
        integer :: N, M
        real(8) :: a, b, eps, err, h, int, prevint
        prevint = 0d0
        N = 1
        h = (b-a)/real(N,8)
        int = trap(a, b, N, h)
        do while (dabs(int - prevint)/3d0 >= eps)
            N = N*2
            prevint = int
            h = h/2d0
            int = trap(a, b, N, h)
        end do
        err = dabs((F1(b) - F1(a)) - int)
    end subroutine
end module