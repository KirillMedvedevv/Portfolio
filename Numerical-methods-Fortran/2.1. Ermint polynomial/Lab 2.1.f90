    program Lab21
    use Methods
    use functions
    implicit none
    integer :: n,m, i ,j, k ,p!n - число узлом, m - число находимых точек
    real(8), allocatable :: grid(:,:), rgrid(:,:), maxerr(:,:), griderr(:,:),ras(:)
    real(8) :: err,buff, a, b
    n = 6
    m = 10**3
    p = 17
    a = -0.5d0
    b = 0.5d0
    allocate(grid(n,3), rgrid(m+1,4),maxerr(p,2),griderr(n-1,2), ras(2*n-1))
    call gengrid(grid,griderr,ras,n,a,b)
    !Строим графики для небольшого числа узлов = m и фактичекую ошибку
    do i = 1,m+1
        rgrid(i,1) = a + (b-a)*real(i-1,8)/real(m,8)
        !rgrid(i,2) = uH(grid,rgrid(i,1),n)
        rgrid(i,2) = H(grid,ras,rgrid(i,1),n)
        rgrid(i,3) = f(rgrid(i,1))
        rgrid(i,4) = rgrid(i,2) - rgrid(i,3)
    end do 
    deallocate(griderr,grid,ras)
    !Зависимость максимальной ошибки
    do i = 2, p
        allocate(griderr(i-1,2), grid(i,3), ras(2*i-1))
        grid = 0d0
        ras = 0d0
        err = -100d0
        buff = 0d0
        griderr = 0d0
        call gengrid(grid,griderr,ras,i,a,b)
        do j = 1,i-1
            buff = dlog10(dabs(H(grid,ras,griderr(j,1),i) - griderr(j,2)))
            !buff = qlog10(qabs(uH(grid,griderr(j,1),i) - griderr(j,2)))
            if (buff>err) err = buff
        end do
        maxerr(i-1,1) = real(i,8)
        maxerr(i-1,2) = err
        deallocate(griderr,grid,ras)
    end do
    !write(*,*) rgrid(15,3)
    open(1,file="Result.csv")
    open(2,file="Err.csv")
    write(1,*) "X", ";", "Y",";","table",";", "err"
    do i =1,m+1
        write(1,fmt="(f16.12, a, f16.12, a, f16.12, a, f16.12)") rgrid(i,1),";", rgrid(i,2),";", rgrid(i,3), ";", rgrid(i,4)
    end do
    
    write(2,*) "Num", ";", "err"
    do i =2,p
        write(2,fmt="(f16.12, a, f16.12)") maxerr(i,1),";", maxerr(i,2)
    end do
    close(1)
    close(2)
    deallocate(rgrid, maxerr)
    pause
    end program Lab21

