    program Lab1Furie 
    use methods
    implicit none
    integer :: i, n, m(3)
    real(8) :: l, d, T0, Te, eps
    real(8) :: a(3), Bi(3), FoMax(3), dFo, X, Fof(3,3), Fo, Xof(3), dX
    l = 55d0
    d = 0.5d0
    a(1) = 20d0
    a(2) = 300d0
    a(3) = 20000d0
    T0 = 200d0
    Te = 400d0
    dFo = 0.01d0
    dX = 0.01d0
    n = dint(1d0/dX)
    !Вычисление чисел Bio
    Bi(1) = atoBi(a(1), d, l, 1)
    write(*,*) Bi(1)
    Bi(2) = atoBi(a(2), d, l, 1)
    write(*,*) Bi(2)
    Bi(3) = atoBi(a(3), d, l, 1)
    write(*,*) Bi(3)
    !Вычисление FoMax и Fof 
    FoMax(1) = getFoMax(Bi(1))
    FoMax(2) = getFoMax(Bi(1))
    FoMax(3) = getFoMax(Bi(1))
    Fof(:,1) = 0.1d0*FoMax(:)
    Fof(:,2) = 0.5d0*FoMax(:)
    Fof(:,3) = 0.9d0*FoMax(:)
    Xof(1) = xtoX(0d0, d, 1)
    Xof(2) = xtoX(d/4d0, d, 1)
    Xof(1) = xtoX(d/2d0, d, 1)
    !Вычисление шагов
    m(1) = dint(FoMax(1)/dFo) + 1
    m(2) = dint(FoMax(2)/dFo) + 1
    m(3) = dint(FoMax(3)/dFo) + 1
    
    eps = 1e-7 !Выбор точности вычисления
    
    10  format(f18.9, 3(a,f18.9))
    20  format(f18.9,3(a,f18.9))
    open(1, file = "T1(x).csv")
    open(2, file = "T2(x).csv")
    open(3, file = "T3(x).csv")
        write(1,*) "X", ";", "T(a1, Fo11)", ";", "T(a1, Fo21)", ";", "T(a1, Fo31)"
        write(2,*) "X", ";", "T(a2, Fo12)", ";", "T(a2, Fo22)", ";", "T(a2, Fo32)"
        write(3,*) "X", ";", "T(a3, Fo13)", ";", "T(a3, Fo23)", ";", "T(a3, Fo33)"
        do i = 1, n
            X = dX*real(i,8)
            write(1,10) X, ";", TX(Fof(1,1), X, Bi(1), eps), ";", TX(Fof(1,2), X, Bi(1), eps), ";", TX(Fof(1,3), X, Bi(1), eps)
            write(2,10) X, ";", TX(Fof(2,1), X, Bi(2), eps), ";", TX(Fof(2,2), X, Bi(2), eps), ";", TX(Fof(2,3), X, Bi(2), eps)
            write(3,10) X, ";", TX(Fof(3,1), X, Bi(3), eps), ";", TX(Fof(3,2), X, Bi(3), eps), ";", TX(Fof(3,3), X, Bi(3), eps)
        end do
    close(1)
    close(2)
    close(3)
    
    open(4, file="Ta1(tau).csv")
    open(5, file="Ta3(tau).csv")
        write(4,*) "tau", ";", "Ta1(X1)", ";", "Ta1(X2)", ";", "Ta1(X3)"
        write(5,*) "tau", ";", "Ta3(X1)", ";", "Ta3(X2)", ";", "Ta3(X3)"
        do i = 0, m(1)
            Fo = dFo*real(i,8)
            write(4, 20) Fo, ";", TX(Fo, Xof(1), Bi(1), eps), ";", TX(Fo, Xof(2), Bi(1), eps), ";", TX(Fo, Xof(3), Bi(1), eps)
        end do
        do i = 0, m(3)
            Fo = dFo*real(i,8)
            write(5, 20) Fo, ";", TX(Fo, Xof(1), Bi(3), eps), ";", TX(Fo, Xof(2), Bi(3), eps), ";", TX(Fo, Xof(3), Bi(3), eps)
        end do
    close(4)
    close(5)
    pause
    end program Lab1Furie

