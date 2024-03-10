module methods
    implicit none
    contains
    function TX(Fo, X, Bi, eps) result(Theta)
        real(8), intent(in) :: Fo, X, Bi, eps
        real(8) :: Thetan, Theta, mu
        integer :: n
        n = 1
        Theta = 0d0
        Thetan = 1d0
        do while(dabs(Thetan) >= eps) 
            mu = MPD(n,Bi, eps)
            Thetan = (2d0*dsin(mu)*dcos(mu*X)/(mu + dsin(mu)*dcos(mu)))*dexp(-Fo*(mu)**2d0)
            Theta = Theta + Thetan
            n = n + 1
        end do
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
            Theta = T*(T0 - Te) + Te
        end if
    end function
    
    function tautoFo(tau, l, d, c, p, m) result(Fo)
        real(8), intent(in) :: tau, l, d, c, p
        integer, intent(in) :: m
        real(8) :: Fo
        if(m == 1) then
            Fo = tau*l*4d0/(c*p*(d**2d0))
        else
            Fo = c*p*(d**2d0)*tau/(4d0*l)
        end if
    end function
    
    function atoBi(a,d,l,m) result(Bi)
        real(8), intent(in) :: a, d, l
        integer, intent(in) :: m
        real(8) :: Bi
        if(m == 1) then
            Bi = a*d/2d0*l
        else
            Bi = 2d0*l*a/d
        end if
    end function
end module