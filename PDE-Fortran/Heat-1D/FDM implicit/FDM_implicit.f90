program FDM_implicit
    use methods
    implicit none
    integer :: i,j, N, m(2), k, Xof(3)
    real(8) :: l, d, T0, Te, eps, p, c
    real(8) :: a(3), Bi(3), FoMax(3), dFo, X, Fof(3,3), Fo, dX
    real(8) :: time1, time2, treads
    type(point), dimension(:), allocatable :: p_c, p_f, data1, data2, data3
    real(8), dimension(:), allocatable :: A1, B1, C1
    call CPU_TIME(time1)
    !Создаём объект символизирующий решение в какой-либо точке
    l = 55d0
    d = 0.5d0
    a(1) = 20d0
    a(2) = 300d0
    a(3) = 20000d0
    T0 = 200d0 + 273.15d0
    Te = 400d0 + 273.15d0
    N = 80 !Im - 1
    dX = 1d0/(N+0d0)
    dFo = (dX**2)*0.5d0
    p = 7900d0
    c = 500d0
    
    !Вычисление чисел Bio
    Bi(1) = atoBi(a(1), d, l, 1)
    Bi(2) = atoBi(a(2), d, l, 1)
    Bi(3) = atoBi(a(3), d, l, 1)
    FoMax(1) = getFoMax(Bi(1))
    FoMax(2) = getFoMax(Bi(2)) 
    FoMax(3) = getFoMax(Bi(3))
    write(*,*) FoMax(3)
    Fof(:,1) = 0.1d0*FoMax(:)
    Fof(:,2) = 0.5d0*FoMax(:)
    Fof(:,3) = 0.9d0*FoMax(:)
    Xof(1) = dint(xtoX(0d0, d, 1)/dX) + 1
    Xof(2) = dint(xtoX(d/4d0, d, 1)/dX) + 1
    Xof(3) = dint(xtoX(d/2d0, d, 1)/dX) + 1
    !Нужно будет хранить один временной срез
    10  format(f18.9, 3(a,f18.9))
    20  format(f18.9,3(a,f18.9))
    open(1, file = "T1(x).csv")
    open(2, file = "T2(x).csv")
    open(3, file = "T3(x).csv")
    open(4, file="Ta1(tau).csv")
    open(5, file="Ta3(tau).csv")
        write(1,*) "X", ";", "T(a1, Fo11)", ";", "T(a1, Fo12)", ";", "T(a1, Fo13)"
        write(2,*) "X", ";", "T(a2, Fo21)", ";", "T(a2, Fo22)", ";", "T(a2, Fo23)"
        write(3,*) "X", ";", "T(a3, Fo31)", ";", "T(a3, Fo32)", ";", "T(a3, Fo33)"
        write(4,*) "tau", ";", "Ta1(X1)", ";", "Ta1(X2)", ";", "Ta1(X3)"
        write(5,*) "tau", ";", "Ta3(X1)", ";", "Ta3(X2)", ";", "Ta3(X3)"
        allocate(p_c(N+1), p_f(N+1), data1(N+1), data2(N+1), data3(N+1), A1(N+1), B1(N+1), C1(N+1))
        
        call set_ics(dX, dFo, 0d0, 1d0, p_c, N, A1, B1, C1, Bi(1))
        do while(p_c(1) % Fo <= FoMax(1) + dFo)
            call FDM_imp(dFo, dX, p_c, p_f, Bi(1), N, A1, B1, C1) !Расчитали будущий на основе текущего
            call add_data(data1, data2, data3, p_c, N, Fof(1,:), dFo, T0, Te, d)
            write(4, 20) tautoFo(p_c(1) % Fo, l, d, c, p, -1), ';', TtoTheta(p_c(Xof(1)) % Theta, T0, Te, -1), ';',TtoTheta(p_c(Xof(2)) % Theta, T0, Te, -1), ';',TtoTheta(p_c(Xof(3)) % Theta, T0, Te, -1)
            p_c = p_f
        end do
        
        do i = 1, N+1
            write(1, 10) data1(i) % X,';', data1(i) % Theta,';', data2(i) % Theta,';', data3(i) % Theta  
        end do

        call set_ics(dX, dFo, 0d0, 1d0, p_c, N, A1, B1, C1, Bi(2))
        do while(p_c(1) % Fo <= FoMax(2) + dFo)
            CALL FDM_imp(dFo, dX, p_c, p_f, Bi(2), N, A1, B1, C1)
            call add_data(data1, data2, data3, p_c, N, Fof(2,:), dFo, T0, Te, d)
            p_c = p_f
        end do
        do i = 1, N+1
            write(2, 10) data1(i) % X,';', data1(i) % Theta,';', data2(i) % Theta,';', data3(i) % Theta  
        end do
     
        call set_ics(dX, dFo, 0d0, 1d0, p_c, N, A1, B1, C1, Bi(3))
        do while(p_c(1) % Fo <= FoMax(3) + dFo)
            CALL FDM_imp(dFo, dX, p_c, p_f, Bi(3), N, A1, B1, C1)
            call add_data(data1, data2, data3, p_c, N, Fof(3,:), dFo, T0, Te, d)
            write(5, 20) tautoFo(p_c(1) % Fo, l, d, c, p, 0), ';', TtoTheta(p_c(Xof(1)) % Theta, T0, Te, -1), ';',TtoTheta(p_c(Xof(2)) % Theta, T0, Te, -1), ';',TtoTheta(p_c(Xof(3)) % Theta, T0, Te, -1) 
            p_c = p_f
        end do
        do i = 1, N+1
            write(3, 10) data1(i) % X,';', data1(i) % Theta,';', data2(i) % Theta,';', data3(i) % Theta  
        end do
    close(1)
    close(2)
    close(3)
    deallocate(p_c, p_f, data1, data2, data3, A1, B1, C1)
    call CPU_TIME(time2)
    write(*,*) dabs(time1 - time2)
    pause
end program FDM_implicit

