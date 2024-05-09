    program Lab4
    use functions
    implicit none
    integer :: N, i
    real(16) :: a, b, eps, h, err, prevint
10  format(f18.10, a, f18.10, a, f18.10)
20  format(f18.10, a, f18.10)   
    a = -2.4q0
    b = -0.5q0
    eps = 0.1q0
    err = 0q0
    prevint = 0q0
    h = 0q0
    N=1
    open(1, file="Err(eps) N(eps).csv")
    open(2, file="Err(h).csv")
    write(1,*) "eps", ";", "Err", ";", "N"
    write(2,*) "h", ";", "Err"
    do i = 1,12
        write(*,*) "Hi"
        call intrad(a, b, eps, err, h, N, prevint)
        write(1,10) qlog10(eps), ";", qlog10(Err), ";", qlog(N + 0q0)/qlog(2q0)
        write(2,20) qlog10(h), ";", qlog10(Err)
        eps = eps*0.1d0
    end do
    close(1)
    close(2)
    pause
end program Lab4

