module methods
    implicit none
    contains
    
    subroutine TDMA(Im, A, B, C, D, F)
        integer, intent(in):: Im
        real(16), dimension(:) :: A, B, C, D
        real(16), dimension(:) :: F
        real(16), dimension(Im) :: alpha, beta
        real(16):: k0
        integer::i
        
        !Прямой ход
        alpha(1) = -A(1)/B(1)
        beta(1) = -D(1)/B(1)
        do i = 2, (Im-1)
            k0 = (B(i) + C(i)*alpha(i-1))
            alpha(i) = -A(i) / k0
            beta(i) = -(D(i) + C(i)*beta(i-1))/k0
        end do

        !Обратный ход
        F = 0q0
        F(Im) = -(D(Im) + C(Im)*beta(Im-1)) / (B(Im) + C(Im)*alpha(Im-1))
        do i = (Im-1), 1, -1
            F(i) = alpha(i)*F(i+1) + beta(i)
        end do
    end subroutine
    
    subroutine set_ics(a1, b1, h, ym, X, A, B, C, D, k, ord)
        integer, intent(in) :: k, ord
        real(16), intent(in) :: a1, b1, h, ym
        real(16), dimension(:) :: A, B, C, D, X
        integer :: i, j
        
        forall(i=1:k) A(i) = 1q0/h + p(X(i+1))/(2q0)
        forall(i=1:k) B(i) = (-2q0/h + q(X(i+1))*h)
        forall(i=1:k) C(i) = 1q0/h - p(X(i+1))/(2q0)
        forall(i=1:k) D(i) = -f(X(i+1))*h
        D(1) = D(1) + C(1)*a1
        if (ord == 1) then
            B(k) = h+b1
            C(k) = -b1
            D(k) = -ym*h
            
        else if (ord == 2) then
            D(k) = D(k) + A(k)*ym
        end if
    end subroutine
    
    pure function p(x) result(y)
        real(16), intent(in) :: x
        real(16) :: y
        y = 1q0 + (qsin(x))**2
    end function
    
    pure function q(x) result(y)
        real(16), intent(in) :: x
        real(16) :: y
        y = -1q0*(qcos(x))**2
    end function
    
    pure function f(x) result(y)
        real(16), intent(in) :: x
        real(16) :: y
        !y = 3d0*dexp(x)
        y = qexp(x)*(1d0 + 2q0*(qsin(x))**2)
    end function
    
end module