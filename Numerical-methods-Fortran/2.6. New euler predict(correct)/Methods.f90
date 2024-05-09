module methods
    implicit none
    contains
    function f(x, y) result(p)
         real(8), intent(in) :: x, y
         real(8) :: p
         p = (x**2 + 1d0) + 2d0*x*y/(x**2 + 1d0)
    end function
    
    function f1(x, ypr, h) result(y) !Получаем сразу y(k)
         real(8), intent(in) :: x, ypr, h
         real(8) :: y
         y = ypr + h*f(x-h, ypr)
         y = ypr + h*f(x, y)
         !y = (ypr+h*(x**2+1d0))/(1d0 - 2d0*x*h/(x**2+1d0))
    end function
    
    function S(x, C) result(y)
        real(8), intent(in) :: x, C
        real(8) :: Y
        y = ((x + C)*x + 1d0)*x + C
    end function
    
    function nev_euler(x, ypr, h, eps, N, p) result(y) !Вычисление уравнения в конкректной точке по предыдущей точке
    real(8), intent(in) ::  x, ypr, eps
    real(8) :: y, err, y1pr, y2pr, h
    integer :: i, N, p !Число отрезков допразбиения
    
    if (p == 0) then
        y = f1(x+h, ypr, h)
    else if (p==1) then
        y1pr = f1(x+h, ypr, h)
        N = 2
        h = h/2d0
        y2pr = f1(x+h, ypr, h)
        y2pr = f1(x+2d0*h, y2pr, h)
        err = dabs(y2pr - y1pr)
        do while(err >= eps)
            h = h/2d0
            N = N*2
            y1pr = y2pr
            y2pr = ypr
            do i = 1, N
                y2pr = f1(x+(i+0d0)*h, y2pr, h)
            end do
            err = dabs(y2pr - y1pr)
        end do
        y = y2pr  
    end if
    end function
end module