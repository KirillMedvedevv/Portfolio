    program Lab1Furie 
    use methods
    use omp_lib
    implicit none
    integer :: i, n, m(2), k, num
    real(8) :: l, d, T0, Te, eps, p, c
    real(8) :: a(3), Bi(3), FoMax(3), dFo, X, Fof(3,3), Fo, Xof(3), dX, mur1(1)
    real(8) :: time1, time2, buffer
    real(8), allocatable :: mur(:)
    
    call CPU_TIME(time1)
    l = 55d0
    d = 0.5d0
    a(1) = 20d0
    a(2) = 300d0
    a(3) = 20000d0
    T0 = 200d0
    Te = 400d0
    dFo = 0.00001d0
    dX = 0.01d0
    p = 7900d0
    c = 500d0
    n = dint(1d0/dX)
    !Вычисление чисел Bio
    Bi(1) = atoBi(a(1), d, l, 1)
    Bi(2) = atoBi(a(2), d, l, 1)
    Bi(3) = atoBi(a(3), d, l, 1)
    FoMax(1) = getFoMax(Bi(1))
    FoMax(2) = getFoMax(Bi(2))
    FoMax(3) = getFoMax(Bi(3))
    Fof(:,1) = 0.1d0*FoMax(:)
    Fof(:,2) = 0.5d0*FoMax(:)
    Fof(:,3) = 0.9d0*FoMax(:)
    Xof(1) = xtoX(0d0, d, 1)
    Xof(2) = xtoX(d/4d0, d, 1)
    Xof(3) = xtoX(d/2d0, d, 1)
    !Вычисление шагов
    
    eps = 1e-12!Выбор точности вычисления
    write(*,*) (FoMax(i), i = 1,3)
    write(*,*) (tautoFo(FoMax(i), l, d, c, p, -1), i = 1,3)
    
    !Произведём псевдо рассчёт для определения максимального голичесвта членов ряда(Взять наименьшее био с которым считаешь)
    buffer = TX(dFo/2d0, X, Bi(1), eps, T0, Te, 0, mur1, num)
    write(*,*) num
    allocate(mur(num))
    buffer = TX(dFo/2d0, X, Bi(1), eps, T0, Te, 1, mur, num)
    write(*,*) num
    10  format(f18.9, 3(a,f18.9))
    20  format(f18.9,3(a,f18.9))
    open(1, file = "T1(x).csv")
    open(2, file = "T2(x).csv")
    open(3, file = "T3(x).csv")
        write(1,*) "X", ";", "T(a1, Fo11)", ";", "T(a1, Fo21)", ";", "T(a1, Fo31)"
        write(2,*) "X", ";", "T(a2, Fo12)", ";", "T(a2, Fo22)", ";", "T(a2, Fo32)"
        write(3,*) "X", ";", "T(a3, Fo13)", ";", "T(a3, Fo23)", ";", "T(a3, Fo33)"
        do i = 0, n
            X = dX*real(i,8)
            write(1,10) xtoX(X, d, -1), ";", TX(Fof(1,1), X, Bi(1), eps, T0, Te, 2, mur, num), ";", TX(Fof(1,2), X, Bi(1), eps, T0, Te, 2, mur, num), ";", TX(Fof(1,3), X, Bi(1), eps, T0, Te, 2, mur, num)
            write(2,10) xtoX(X, d, -1), ";", TX(Fof(2,1), X, Bi(2), eps, T0, Te, 2, mur, num), ";", TX(Fof(2,2), X, Bi(2), eps, T0, Te, 2, mur, num), ";", TX(Fof(2,3), X, Bi(2), eps, T0, Te, 2, mur, num)
            write(3,10) xtoX(X, d, -1), ";", TX(Fof(3,1), X, Bi(3), eps, T0, Te, 2, mur, num), ";", TX(Fof(3,2), X, Bi(3), eps, T0, Te, 2, mur, num), ";", TX(Fof(3,3), X, Bi(3), eps, T0, Te, 2, mur, num)
        end do
    close(1)
    close(2)
    close(3)
    
    open(4, file="Ta1(tau).csv")
    open(5, file="Ta3(tau).csv")
        write(4,*) "tau", ";", "Ta1(X1)", ";", "Ta1(X2)", ";", "Ta1(X3)"
        write(5,*) "tau", ";", "Ta3(X1)", ";", "Ta3(X2)", ";", "Ta3(X3)"
        m(1) = dint(FoMax(1)/dFo) + 1
        m(2) = dint(FoMax(3)/dFo) + 1
        write(4, 20) tautoFo(0d0, l, d, c, p, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1)
        do i = 1, m(1)
            Fo = dFo*real(i,8)
            write(4, 20) tautoFo(Fo, l, d, c, p, -1), ";", TX(Fo, Xof(1), Bi(1), eps, T0, Te, 2, mur, num), ";", TX(Fo, Xof(2), Bi(1), eps, T0, Te, 2, mur, num), ";", TX(Fo, Xof(3), Bi(1), eps, T0, Te, 2, mur, num)
        end do
        write(5, 20) tautoFo(0d0, l, d, c, p, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1), ";", TtoTheta(1d0, T0+273.15d0, Te+273.15d0, -1)
        do i = 1, m(2)
            Fo = dFo*real(i,8)
            write(5, 20) tautoFo(Fo, l, d, c, p, -1), ";", TX(Fo, Xof(1), Bi(3), eps, T0, Te, 2, mur, num), ";", TX(Fo, Xof(2), Bi(3), eps, T0, Te, 2, mur, num), ";", TX(Fo, Xof(3), Bi(3), eps, T0, Te, 2, mur, num)
        end do
    close(4)
    close(5)
    deallocate(mur)
    call CPU_TIME(time2)
    write(*,*) (time2 - time1)
    pause
    end program Lab1Furie

