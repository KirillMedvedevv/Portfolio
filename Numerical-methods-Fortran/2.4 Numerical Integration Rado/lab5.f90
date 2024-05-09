program lab5
    use methods
    implicit none
    real(8) :: x0, y0, y, x, C, a, b, h, eps, solve, maxerr, err, h1
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
    open(3, File = "Err and N(eps).csv")
        write(1,*) "X", ";", "n_solve(x)", ";", "solve(x)", ";", "f_Err(x)"
        write(2, *) "h", ";", "err"
        write(3, *) "eps", ";", "err", ";", "N"
        !Построение грификов без заданной точности
        y = y0
        x = x0
        maxerr = 0d0
        h = (b-a)/(N+0d0)
        do i = 1, N
            x = x + h
            solve = S(x, C)
            y = euler(x-h, y, h, eps, steps, 0)
            write(1,10) x, ";", y, ";", solve, ";", dlog10(dabs(solve - y))
        end do
        !Строим err(eps)
        eps = 0.1d0
        N = 20
        h1 = (b-a)/(N+0d0) 
        C = 0d0
        do i = 1, 7
            maxerr = 0d0
            maxstep = 0
            y = y0
            x = x0
            do j = 1, N
                h = h1
                x = x + h
                solve = S(x, C)
                y = euler(x-h, y, h, eps, steps, 1)
                err = dabs(y - solve)
                if (err > maxerr) maxerr = err
                if (steps > maxstep) maxstep = steps
            end do
            write(*,*) eps
            write(3,30) dlog10(eps), ';', dlog10(maxerr), ';', dlog(steps + 0d0)/dlog(2d0)
            eps = eps*0.1d0
        end do
        
        !Строим err(h)
        N = 1
        do i = 1, 8
            y = y0
            x = x0
            N = N*10
            maxerr = 0d0
            h = (b-a)/(N+0d0)
            do j = 1, N
                x = x + h
                solve = S(x, C)
                y = euler(x-h, y, h, eps, steps, 0)
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
end program lab5

