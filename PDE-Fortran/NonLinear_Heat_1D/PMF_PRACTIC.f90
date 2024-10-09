program PMF_PRACTIC
use methods
implicit none
    integer :: i, logic
    integer, parameter :: Im = 41
    real(8), parameter :: eps = 1e-5, pi = 4d0*datan(1d0)
    real(8) :: Q_tot, potok_s, potok_t, c1, Tb, Tc, dX
    real(8) :: L, d, T0, Tm, lambda, a, sig0, kt, a1, b1, phim !Входные данные
    real(8), dimension(Im) :: At, Bt, Ct !Для рассчёта
    real(8), dimension(Im) :: T, sig, phi, Q !Выходные данные
    real(8), dimension(12) :: input_data
    
    !Считываем входные данные
    open(1, File="Input_Data.txt")
        do i = 1,12
            read(1,*) input_data(i)
        end do
    close(1)
    call get_data(input_data, L, d, T0, Tm, lambda, a, sig0, kt, logic, a1, b1, phim)
    dX = L/(Im - 1d0) !Шаг сетки
    
    !Вычисление At, Bt, Ct
    At(2:Im-1) = lambda/dX**2
    Bt(2:Im-1) = -2d0*(lambda/(dX**2) + a*2d0/d)
    Ct(2:Im-1) = At(2:Im-1)  
    
    if (logic == 1) then
        !Нахождение phim МПД
        do while(abs(b1-a1) >= eps)
            c1 = (a1+b1)/2d0
            call get_T(Im, At, Bt, Ct, T, phi, sig, T0, dX, d, a, kt, b1, sig0, eps)
            Tb = maxval(T) - Tm
            call get_T(Im, At, Bt, Ct, T,phi, sig, T0, dX, d, a, kt, c1, sig0, eps)
            Tc = maxval(T) - Tm
            if (Tb*Tc < 0d0) then
                a1 = c1
            else
                b1 = c1
            end if
        end do
        phim = (a1+b1)/2d0
    else
        !phim было введено. Ничего не делаем
    end if
    
    call get_T(Im, At, Bt, Ct, T, phi, sig, T0, dX, d, a, kt,  phim, sig0, eps) 
    
    !Формируем и интегрируем Q
    Q=0d0
    Q(1) = sig(1)*((-3d0*phi(1)+4d0*phi(2)-phi(3))/(2d0*dX))**2 
    do i = 2, Im-1
        Q(i) = sig(i)*((phi(i+1)-phi(i-1))/(2d0*dX))**2
    end do
    Q(Im) = sig(Im)*((3d0*phi(Im)-4d0*phi(Im-1)+phi(Im-2))/(2d0*dX))**2
    
    Q_tot = (prim(Im, Q, dX))*pi*(d**2)/4d0
    potok_s = pi*a*D*prim(Im, T-T0,dX)
    potok_t = 0.25d0*pi*(d**2)*lambda*(T(2)-T(1) + T(Im-1) - T(Im))/dX !First order
    
    !Вывод данных
    10  format(f16.6, a, f16.6)
    open(1, File="T(x).csv")
    open(2, File="sig(x).csv")
    open(3, File="Q(x).csv")
    write(1,*) "x", ";", "T"
    write(2,*) "x", ";", "sig"
    write(3,*) "x", ";", "Q"
        do i = 1, Im
            write(1,10) (0d0 + dX*(i-1d0)), ';', T(i)
            write(2,10) (0d0 + dX*(i-1d0)), ';', sig(i)
            write(3,10) (0d0 + dX*(i-1d0)), ';', Q(i)
        end do
        write(*,*) "phim", phim
        write(*,*) "numeric_Tm", maxval(T)
        write(*,*) "potok_t", potok_t
        write(*,*) "potok_s", potok_s
        write(*,*) "Q_tot", Q_tot
        write(*,*) "nev", (potok_s+potok_t)/Q_tot - 1d0
    close(1)
    close(2)
    close(3)
pause
end program PMF_PRACTIC

