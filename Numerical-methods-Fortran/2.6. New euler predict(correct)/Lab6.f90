program Lab6
    use methods
    implicit none
    real(8) :: x0, y0, y, x, C, a, b, h, eps, solve, maxerr, err, h1, y01
    integer :: N, steps, i, j, maxstep
    !initial conditions set
    a = 0d0
    b = 2d0
    x0 = a
    y0 = 0d0
    C = 0d0
    !Number of otrezok
    N = 20
    10 format(f16.13, a, f16.13, a, f16.13, a, f16.13)    
    30 format(f16.13, a, f16.13, a, f16.13)
    open(1, File = "Result1.csv")
    open(2, File = "Err(h).csv")
    open(3, File = "Err(y0).csv")
        write(1,*) "X", ";", "n_solve(x)", ";", "solve(x)", ";", "f_Err(x)"
        write(2, *) "h", ";", "err"
        write(3, *) "y0", ";", "err"
        !Построение грификов без заданной точности
        y = y0
        x = x0
        maxerr = 0d0
        h = (b-a)/(N+0d0)
        do i = 1, N
            x = x + h
            solve = S(x, C)
            y = nev_euler(x-h, y, h, eps, steps, 0)
            write(1,10) x, ";", y, ";", solve, ";", (dabs(solve - y))
        end do
        !Строим err(y0)
        eps = 0.1d0
        N = 2000
        h1 = (b-a)/(N+0d0) 
        C = 0d0
        do i = 1, 201
            maxerr = 0d0
            maxstep = 0
            y01 = -0.01d0 + 0.0001d0*(i-1d0)
            y = y01
            x = x0
            do j = 1, N
                h = h1
                x = x + h
                solve = S(x, C)
                y = nev_euler(x-h, y, h, eps, steps, 0)
                err = dabs(y - solve)
                if (err > maxerr) maxerr = err
            end do
            write(3,30) y01, ';', maxerr
        end do
        
        !Строим err(h)
        N = 1
        do i = 1, 0
            y = y0
            x = x0
            N = N*10
            maxerr = 0d0
            h = (b-a)/(N+0d0)
            do j = 1, N
                x = x + h
                solve = S(x, C)
                y = nev_euler(x-h, y, h, eps, steps, 0)
                err = dabs(y - solve)
                if (err > maxerr) maxerr = err
            end do
            write(*,*) "Hello there"
            write(2,30) dlog10(h), ';', dlog10(maxerr)
        end do
    close(1)
    close(2)
    close(3)
    pause
end program Lab6

