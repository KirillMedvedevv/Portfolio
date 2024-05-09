module methods
    implicit none
    type :: point
        real(8) :: X, Fo, Theta
    end type
    contains
    subroutine set_ics(dX, dFo, Fo, Theta, p_c, N, A, B, C, Bi)
        type(point), dimension(:) :: p_c
        real(8) :: dX, dFo, Fo, Theta, Bi
        real(8), dimension(:) :: A, B, C
        integer :: i, N
        do i = 1, N+1
            p_c(i) % X = dX*(i-1d0)
            p_c(i) % Fo = 0d0
            p_c(i) % Theta = 1d0
        end do
        A(1) = 2d0/(dX**2)
        A(2:N) = 1d0/(dX**2)
        A(N+1) = 0d0
        B(1:N) = -2d0/(dX**2) - 1d0/dFo
        B(N+1) = -2d0/(dX**2) - 1d0/dFo - 2d0*Bi/(dX)
        C(1) = 0d0
        C(2:N) = A(2:N)
        C(N+1) = 2d0/(dX**2)
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
    
    subroutine FDM_imp(dFo, dX, p_c, p_f, Bi, N, A, B, C) 
        real(8), allocatable :: D(:), F(:)
        real(8), dimension(:) :: A, B, C
        type(point), dimension(:) :: p_c, p_f
        real(8) :: dFo, dX, Bi
        integer :: i, N
        allocate(D(N+1), F(N+1))
        D = (p_c % theta)/dFo
        call Progonka(N+1, A, B, C, D, F)
        p_f % Fo = p_c % Fo + dFo
        p_f % X = p_c % X
        p_f % Theta = F
        deallocate(D, F)
    end subroutine
    
    subroutine Progonka(Im, A, B, C, D, F) ! Метод прогонки
        implicit none
        integer, intent(in):: Im
        real(8), dimension(1:Im), intent(in):: A, B, C, D
        real(8), dimension(1:Im), intent(out):: F
        real(8), dimension(1:Im):: alpha, beta
        real:: k0
        integer::i

        !Прямой ход
        alpha(1) = -A(1) / B(1)
        beta(1) = -D(1) / B(1)
        do i = 2, (Im-1)
            k0 = (B(i) + C(i)*alpha(i-1))
            alpha(i) = -A(i) / k0
            beta(i) = -(D(i) + C(i)*beta(i-1)) / k0
        end do

        !Обратный ход
        F = 0d0
        F(Im) = -(D(Im) + C(Im)*beta(Im-1)) / (B(Im) + C(Im)*alpha(Im-1))
        do i = (Im-1), 1, -1
            F(i) = alpha(i)*F(i+1) + beta(i)
        end do
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
        else if (m == -1) then
            Theta = T*(T0 - Te) + Te - 273.15d0
        else
            Theta = T
        end if
    end function
    
    function tautoFo(tau, l, d, c, p, m) result(Fo1)
        real(8), intent(in) :: tau, l, d, c, p
        integer, intent(in) :: m
        real(8) :: Fo1
        if(m == 1) then
            Fo1 = tau*l*4d0/(c*p*(d**2))
        else if (m==-1) then
            Fo1 = c*p*(d**2)*tau/(4d0*l)
        else
            Fo1 = tau
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