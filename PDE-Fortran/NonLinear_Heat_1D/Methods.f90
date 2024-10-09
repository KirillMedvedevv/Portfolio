module methods
    implicit none
    contains
    !Метод прогонки; F - вектор решения
    !A, B, C, D – векторы прогоночных коэффициентов
    !Im – число уравнений (нумерация начинается с Ib, стнандарт 1)
    subroutine Progonka(Ib, Im, A, B, C, D, F)
        integer, intent(in):: Im, Ib
        real(8), dimension(1:Im), intent(in):: A, B, C, D
        real(8), dimension(1:Im), intent(out) :: F
        real(8), dimension(1:Im) :: alpha, beta
        real(8) :: k0
        integer :: i
        F = 0d0
        !Прямой ход
        alpha(Ib) = -A(Ib) / B(Ib)
        beta(Ib) = -D(Ib) / B(Ib)
        do i = Ib+1, (Im-Ib)
            k0 = (B(i) + C(i)*alpha(i-1))
            alpha(i) = -A(i) / k0
            beta(i) = -(D(i) + C(i)*beta(i-1)) / k0
        end do

        !Обратный ход
        F(Im-Ib+1) = -(D(Im-Ib+1) + C(Im-Ib+1)*beta(Im-Ib)) / (B(Im-Ib+1) + C(Im-Ib+1)*alpha(Im-Ib))
        do i = (Im-Ib), Ib, -1
            F(i) = alpha(i)*F(i+1) + beta(i)
        end do
        
    end subroutine Progonka
    
    subroutine get_phi(Im, sig, phi) !Получаем значение потенциала на сетке
        integer :: i, Im
        real(8), dimension(Im) :: Ap, Bp, Cp, Dp, Fp, sig, phi
        
        !Вектора Слау
        Ap(2) = sig(2)*sig(3)/(sig(2)+sig(3))
        Cp(2) = 0d0
        Dp(2) = 0d0
        Bp(2) = -(sig(3)*sig(2)/(sig(3) + sig(2)) + sig(2)*sig(1)/(sig(2) + sig(1)))
        do i = 3, Im-2
            Ap(i) = sig(i)*sig(i+1)/(sig(i)+sig(i+1))
            Cp(i) = sig(i-1)*sig(i)/(sig(i-1) + sig(i))
            Bp(i) = -(Ap(i) + Cp(i))
            Dp(i) = 0d0
        end do
        Ap(Im-1) = 0d0
        Cp(Im-1) = (sig(Im-2)*sig(Im-1)/(sig(Im-2) + sig(Im-1)))
        Bp(Im-1) = -(Cp(Im-1) + sig(Im-1)*sig(Im)/(sig(Im-1) + sig(Im)))
        Dp(Im-1) = phi(Im)*sig(Im-1)*sig(Im)/(sig(Im-1) + sig(Im))
        
        !Вызов метода прогонки
        call Progonka(2,Im,Ap,Bp,Cp,Dp,Fp)
        phi(2:Im-1) = Fp(2:Im-1)
        
    end subroutine
    
    subroutine get_T(Im, At, Bt, Ct, T, phi, sig, T0, dX, d, a, kt, phim, sig0, eps)
        integer :: i, Im
        real(8) :: T0, dX, d, a, kt, phim, eps, sig0
        real(8), dimension(Im) :: At, Bt, Ct, Dt, Ft, T, phi, sig, Tprev
        T = T0
        Tprev = T0
        phi(Im) = phim
        do while(1)
            sig = sig0*(1d0 + kt*(Tprev-T0))**(-1) !Проводимость
            call get_phi(Im, sig, phi) !Потенциал
            !Вычисляем темппературу
            Dt(2) = sig(2)*((phi(3)-phi(1))/(2d0*dX))**2 + T0*Ct(2) + 4d0*a*T0/d
            do i = 3,Im-2
                Dt(i) = sig(i)*((phi(i+1)-phi(i-1))/(2d0*dX))**2 + 4d0*a*T0/d
            end do
            Dt(Im-1) = sig(Im-1)*((phi(Im)-phi(Im-2))/(2d0*dX))**2 + T0*At(Im-1) + 4d0*a*T0/d
            call Progonka(2, Im, At, Bt, Ct, Dt, Ft)
            T(2:Im-1) = Ft(2:Im-1)
            T(1) = T0
            T(Im) = T0
            if(maxval(abs(T - Tprev))/maxval(abs(T)) < eps) exit
            Tprev = T
        end do
    end subroutine
    
    function prim(Im, Q, dX) result(Q_tot)
        integer :: Im
        real(8) :: dX, Q_tot
        real(8), dimension(Im) :: Q
        
        Q_tot = sum(Q(2:Im-1))*dX
    end function
    
    subroutine get_data(input_data, L, d, T0, Tm, lambda, a, sig0, kt, logic, a1, b1, phim)
        integer :: logic
        real(8) :: L, d, T0, Tm, lambda, a, sig0, kt, a1, b1, phim
        real(8), dimension(12) :: input_data
        L = input_data(1)
        d = input_data(2)
        T0 = input_data(3)
        Tm = input_data(4)
        lambda = input_data(5)
        a = input_data(6)
        sig0 = input_data(7)
        kt = input_data(8)
        logic = int(input_data(9))
        a1 = input_data(10)
        b1 = input_data(11)
        phim = input_data(12)
    end subroutine
end module