module methods
    USE OMP_LIB
    implicit none
    type :: point
        real(8) :: X, Fo, Theta
    end type
    contains
    subroutine set_ics(dX, Fo, Theta, p_c, N)
        type(point), dimension(:) :: p_c
        real(8) :: dX, Fo, Theta
        integer :: i, N
        do i = 1, N+1
            p_c(i) % X = dX*(i-1d0)
            p_c(i) % Fo = 0d0
            p_c(i) % Theta = 1d0
        end do
    end subroutine
    
    subroutine add_data(d1, d2, d3, p_c, N, Fof, dFo, T0, Te, d)
        type(point), dimension(:) :: p_c, d1, d2, d3
        real(8), dimension(:) :: Fof
        real(8) dFo, T0, Te, d
        integer :: i, N
        if (dabs(p_c(1) % Fo - Fof(1)) <= dFo) then
            d1 % Fo = Fof(3)
            do i = 1, N+1
                d1(i) % Theta = TtoTheta(p_c(i) % Theta, T0, Te, -1)
                d1(i) % X = xtoX(p_c(i) % X, d, -1)
            end do
        else if (dabs(p_c(1) % Fo - Fof(2)) <= dFo) then
            d2 % Fo = Fof(2)
            do i = 1, N+1
                d2(i) % Theta = TtoTheta(p_c(i) % Theta, T0, Te, -1)
                d2(i) % X = xtoX(p_c(i) % X, d, -1)
            end do
        else if (dabs(p_c(1) % Fo - Fof(3)) <= dFo) then
            d3 % Fo = Fof(3)
            do i = 1, N+1
                d3(i) % Theta = TtoTheta(p_c(i) % Theta, T0, Te, -1)
                d3(i) % X = xtoX(p_c(i) % X, d, -1)
            end do
        end if
    end subroutine
    
    subroutine FDM_emp(dFo, dX, p_c, p_f, Bi, N) 
        type(point), dimension(:) :: p_c, p_f
        real(8) :: dFo, dX, Bi
        integer :: i, N
        p_f % Fo = p_c % Fo + dFo
        p_f % X = p_c % X
        do i = 2, N
            p_f(i) % Theta = (dFo/(dX**2))*(p_c(i+1) % Theta + p_c(i-1) % Theta) + (p_c(i) % Theta)*(1d0 - 2d0*dFo/(dX**2))
        end do
        p_f(1) % Theta = p_f(2) % Theta
        p_f(N+1) % Theta = (p_f(N) % Theta)/(1d0 + dX*Bi) 
    end subroutine
    
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