module methods
    implicit none
    contains
    subroutine getfile()
    end subroutine
    function  TX(Fo, X, Bi, eps, T0, Te, p, mur, n) result(Theta)
        real(8) :: Fo, X, Bi, eps, T0, Te
        real(8) :: Thetan, Theta, mu
        real(8), dimension(:) :: mur
        integer :: n, p
        n = 1
        Theta = 0d0
        Thetan = 1d0
        if (p == 0) then !Узнать количество членов ряда
            do while(dabs(Thetan) >= eps) 
                mu = MPD(n, Bi, eps)
                Thetan = (2d0*dsin(mu)*dcos(mu*X)/(mu + dsin(mu)*dcos(mu)))*dexp(-Fo*(mu)**2)
                Theta = Theta + Thetan
                write(*,*) Thetan
                n = n + 1
            end do
        else if (p == 1) then !Сгенерировать mur
            do while(dabs(Thetan) >= eps) 
                mu = MPD(n, Bi, eps)
                mur(n) = mu
                Thetan = (2d0*dsin(mu)*dcos(mu*X)/(mu + dsin(mu)*dcos(mu)))*dexp(-Fo*(mu)**2)
                Theta = Theta + Thetan
                n = n + 1
            end do
        else
            do while(dabs(Thetan) >= eps) !Cчитать по имеющимся значениям
                mu = mur(n)
                Thetan = (2d0*dsin(mu)*dcos(mu*X)/(mu + dsin(mu)*dcos(mu)))*dexp(-Fo*(mu)**2)
                Theta = Theta + Thetan
                n = n + 1
            end do
        end if
        
        Theta = TtoTheta(Theta, T0+273.15d0, Te+273.15d0, -1)
    end function
    
    function MPD(n, Bi, eps) result(y)
        integer, intent(in) :: n
        real(8), parameter :: pi = 4d0*datan(1d0)
        real(8) :: a, b, c, y, Bi, eps
        
        a = pi*real((n-1), 8)
        b = pi*real(2*n - 1,8)/2d0
        do while (dabs(a-b) >= eps)
            c = (a+b)/2d0
            if(f(b, Bi)*f(c,Bi) < 0d0) then
                a = c
            else
                b = c
            end if
        end do
        y = (a+b)/2d0
    end function
    
    function f(x, Bi) result(y)
        real(8), intent(in) :: x, Bi
        real(8) :: y
        y = dcotan(x) - x/Bi
    end function
    
    function getFoMax(Bi) result(FoMax)
        real(8), intent(in) :: Bi
        real(8) :: FoMax
        if(Bi < 1.25d0) then
            FoMax = 3.11d0*((Bi)**(-0.88d0))
        else if (Bi > 20d0) then
            FoMax = 1.10d0
        else 
            FoMax = 2.76d0*((Bi)**(-0.31d0))
        end if
    end function
    
    function xtoX(x, d, m) result(X1)
        real(8), intent(in) :: x, d
        integer, intent(in) :: m
        real(8) :: X1
        if(m == 1) then
            X1 = 2d0*x/d
        else
            X1 = d*x/2d0
        end if
    end function
    
    function TtoTheta(T, T0, Te, m) result(Theta)
        real(8), intent(in) :: T, T0, Te
        integer, intent(in) :: m
        real(8) :: Theta
        if(m == 1) then
            Theta = (T - Te)/(T0 - Te)
        else
            Theta = T*(T0 - Te) + Te - 273.15d0
        end if
    end function
    
    function tautoFo(tau, l, d, c, p, m) result(Fo1)
        real(8), intent(in) :: tau, l, d, c, p
        integer, intent(in) :: m
        real(8) :: Fo1
        if(m == 1) then
            Fo1 = tau*l*4d0/(c*p*(d**2))
        else
            Fo1 = c*p*(d**2)*tau/(4d0*l)
        end if
    end function
    
    function atoBi(a,d,l,m) result(Bi)
        real(8), intent(in) :: a, d, l
        integer, intent(in) :: m
        real(8) :: Bi
        if(m == 1) then
            Bi = a*d/(2d0*l)
        else
            Bi = 2d0*l*a/d
        end if
    end function
end module