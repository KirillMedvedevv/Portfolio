module methods
    implicit none
    contains
    function f(x, y) result(p)
         real(8), intent(in) :: x, y
         real(8) :: p
         p = (x**2 + 1d0) + 2d0*x*y/(x**2 + 1d0)
    end function
    
    function S(x, C) result(y)
        real(8), intent(in) :: x, C
        real(8) :: Y
        y = ((x + C)*x + 1d0)*x + C
    end function
    
    function euler(x, ypr, h, eps, N, p) result(y) !Вычисление уравнения в конкректной точке по предыдущей точке
    real(8), intent(in) ::  x, ypr, eps
    real(8) :: y, err, y1pr, y2pr, h
    integer :: i, N, p !Число отрезков допразбиения
    
    y1pr = ypr + h*f(x, ypr)
    if (p == 0) then
        y = ypr + h*f(x, ypr)
    else
        y1pr = ypr + h*f(x, ypr)
        N = 2
        h = h/2d0
        y2pr = ypr + h*f(x, ypr)
        y2pr = y2pr + h*f(x + h, y2pr)
        err = dabs(y2pr - y1pr)
        do while(err >= eps)
            h = h/2d0
            N = N*2
            y1pr = y2pr
            y2pr = ypr
            do i = 1, N
                y2pr = y2pr + h*f(x + h*(i-1d0), y2pr)
            end do
            err = dabs(y2pr - y1pr)
        end do
        y = y2pr
    end if
    end function
end module