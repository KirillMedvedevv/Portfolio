module functions
    implicit none
    contains
    function f(x) result(y)
        real(16) :: x, y
        y = (((((x * x) - 3.5q0) * x + 2.5q0) * x - 7q0) * x - 6.4q0)*qcos(0.4q0*x)
    end function
    
    function F1(x) result(y)
    real(16) :: p1, p2, x, y
    p1 = qsin(x*0.4q0)
    p2 = qcos(x*0.4q0)
    y = p1*(2.5q0*x**5 - 321.25q0*x**3 + 6.25q0*x**2 + 12029.375q0*x - 94.125q0) + p2*(31.25q0*x**4 - 2409.375q0*x**2 + 31.25q0*x + 30073.4375q0)
    end function
    
    function rado(a, b) result(y) !¬ычисление интеграла на отрезке разбиени€
        real(16) :: x1,x2, A0, A1, A2, y
        real(16), intent(in) :: a, b
        x1 = ((1q0 - qsqrt(6q0))/5q0)*(b-a)/2q0 + (a+b)/2q0
        x2 = ((1q0 + qsqrt(6q0))/5q0)*(b-a)/2q0 + (a+b)/2q0
        A0 = 2q0/9q0
        A1 = 1.0249716523768433q0
        A2 = 0.7528061254009345q0 
        y = A0*f(a) + A1*f(x1) + A2*f(x2)
    
    end function
    
    subroutine intrad(a, b, eps, err, h, N, prevint)
    real(16) :: a, b, err, eps, int, prevint, h
    integer :: i, N
    prevint = 0q0
    int = 0q0
    N = 1
    h = (b-a)
    do i = 1, N
        int = int + rado(a+(i-1q0)*h, a + (i + 0q0)*h)
    end do
    int = int*h/2q0
    do while(qabs(int - prevint)/31q0 >= eps)
        N = N*2
        h = h/2q0
        prevint = int
        int = 0q0
        do i = 1, N
            int = int + rado(a+(i-1q0)*h, a + (i + 0q0)*h)
        end do
        int = int*h/2q0
    end do
    err = qabs(int - F1(b) + F1(a))
    end subroutine
end module