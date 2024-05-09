module functions
    implicit none
    contains
    function f(x) result(y)
        real(8), intent(in) :: x
        real(8)             :: y
        y = x**2d0 + 1d0 - dACOS(x)
    end function
    
    function df(x) result(y)
        real(8), intent(in) :: x
        real(8)             :: y
        y = 2d0*x + 1d0/dsqrt(1d0 - x**2d0)
    end function
    
    function uH(grid,x,n) result(y)
        integer, intent(in) :: n
        real(8), intent(in) :: x
        real(8), dimension(:,:) :: grid(n,4)
        integer :: i, j, k
        real(8) :: y, mult, insum
        y=0d0
        do j = 1,n
            mult = 1d0
            do i = 1,n
                if (i /= j) mult = mult*((x-grid(i,1))/(grid(j,1) - grid(i,1)))**2d0
            end do
            insum = 0d0
            do k = 1,n
                if (k /= j) insum = insum + ((x-grid(j,1))/(grid(j,1) - grid(k,1)))
            end do
            y = ((x - grid(j,1))*grid(j,3) + (1d0 - 2d0*insum)*grid(j,2))*mult + y
        end do 
    end function 
    
    function H(grid,ras,x,n) result(y)
        integer, intent(in) :: n
        real(8), intent(in) :: x
        real(8), dimension(:,:) :: grid(n,4), ras(2*n-1)
        integer :: i
        real(8) :: y, mult
        
        y = grid(1,2)
        mult = 1d0
        do i = 1,2*n-1
            mult = mult*(x-grid(i/2 + mod(i,2),1))
            y = y + ras(i)*mult
        end do
        write(*,*) n,x,mult
    end function
    
end module