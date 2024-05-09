program Lab7
use methods
    implicit none
    integer :: m, i, j, k, ord !Число узлов
    real(16) :: a0, b0, a1, b1, ym, h, ym_new
    real(16), allocatable :: A(:), B(:), C(:), D(:), X(:), Y(:), sol(:)
    m = 20
    a0 = 0q0 !Левый конец
    b0 = 1q0 !Правый конец
    a1 = 1q0 !A
    h = (b0-a0)/(m+0d0)
    ord = 1 !Выбор порядка метода
    if (ord == 1) then
        k = m
        b1 = 1q0
        ym = 2q0*qexp(1q0)
    else if (ord == 2) then
        k = m-1
        b1 = 0q0
        ym = qexp(1q0)
    end if
    allocate(A(k), B(k), C(k), D(k), Y(k), X(m+1), sol(m+1))
    do i = 1, m+1
        X(i) = a0 + (i-1)*h
    end do
    call set_ics(a1, b1, h, ym, X, A, B, C, D, k, ord)
    call TDMA(k, A, B, C, D, Y)
    if (ord == 1) then
        sol(1) = a1
        sol(2:) = Y
    else if (ord == 2) then
        sol(1) = a1
        sol(2:m) = Y
        sol(m+1) = ym
    end if
    open(1, File = "Result.csv")
        write(1, *) "X", ';', "solve", ';', "n_solve"
        do i = 1, m+1
            write(1,"(f14.10, a, f14.10, a, f14.10)") X(i), ';', qexp(X(i)), ';', sol(i)
        end do
    close(1)
    deallocate(A, B, C, D, Y, X, sol)
    
    open(2, File = "Err.csv")
        write(2, *) "h", ';', "err", ';', 'N'
        m = 2
        h = (b0-a0)/(m+0q0)
        do i = 1, 10
            m=m*2
            h = h/2q0
            if (ord == 1) then
                k = m
                b1 = 1q0
                ym = 2q0*qexp(1q0)
            else if (ord == 2) then
                k = m-1
                b1 = 0q0
                ym = qexp(1q0)
            end if
            allocate(A(k), B(k), C(k), D(k), Y(k), X(m+1), sol(m+1))
            write(*,*) i, "Hi1"
            do j = 1, m+1
                X(j) = a0 + (j-1q0)*h
            end do
            call set_ics(a1, b1, h, ym, X, A, B, C, D, k, ord)
            call TDMA(k, A, B, C, D, Y)
            sol(1) = 1q0
            if (ord == 1) then
                sol(2:) = Y
            else if (ord == 2) then
                sol(2:m) = Y
                sol(m+1) = ym
            end if
            write(2,"(f14.10, a, f14.10, a, f14.10)") qlog10(h), ';', qlog10(maxval(qabs(sol - qexp(X)))), ';', qlog(m+0q0)/qlog(2q0)
            deallocate(A, B, C, D, Y, X, sol)
        end do
    close(2)
    
    open(3, File="Err(x0).csv")
    write(3, *) "y0", ';', "err"
    m = 100
    h = (b0-a0)/(m+0d0)
    a0 = 0q0
    b0 = 1q0
    if (ord == 1) then
        k = m
        b1 = 1q0
        ym = 2q0*qexp(1q0)
    else if (ord == 2) then
        k = m-1
        b1 = 0q0
        ym = qexp(1q0)
    end if
    a1 = 1.000q0 - 0.010q0
    do i = 1, 201
        allocate(A(k), B(k), C(k), D(k), Y(k), X(m+1), sol(m+1))
        a1 = a1 + 0.0001q0
        do j = 1, m+1
            X(j) = a0 + (j-1q0)*h
        end do
        call set_ics(a1, b1, h, ym, X, A, B, C, D, k, ord)
        call TDMA(k, A, B, C, D, Y)
        sol(1) = a1
        if (ord == 1) then
            sol(2:) = Y
        else if (ord == 2) then
            sol(2:m) = Y
            sol(m+1) = ym
        end if
        write(3,"(f19.14, a, f19.14)") a1, ';', (maxval(qabs(sol - qexp(X))))
        deallocate(A, B, C, D, Y, X, sol)
    end do
    close(3)
    pause
end program Lab7

