program Lab3
use functions
    implicit none
    integer :: N, i
    real(8) :: a, b, eps, h, err, time1, time2
10  format(f18.10, a, f18.10, a, f18.10)
20  format(f18.10, a, f18.10)
    call CPU_TIME(time1)
    N = 1
    a = -2.4d0
    b = -0.5d0
    eps = 0.1d0
    err = 0d0
    h = 0d0
    N = 1
    open(1, file="Err(eps) N(eps).csv")
    open(2, file="Err(h).csv")
    write(1,*) "eps", ";", "Err", ";", "N"
    write(2,*) "h", ";", "Err"
    do i = 1,14
        call inttrap(a, b, eps, err, h, N)
        write(1,10) dlog10(eps), ";", dlog10(Err), ";", dlog(N + 0d0)/dlog(2d0)
        write(2,20) dlog10(h), ";", dlog10(Err)
        eps = eps/10d0
    end do
    close(1)
    close(2)
    call CPU_TIME(time2)
    write(*,*) dabs(time1 - time2)
    pause
end program Lab3

