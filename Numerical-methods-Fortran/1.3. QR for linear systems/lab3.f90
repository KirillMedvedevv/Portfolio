    program lab3
    use functions
    implicit none
    double precision, allocatable :: res(:,:), resg(:,:), resb(:,:), A(:,:), b(:), X(:), nX(:), dX(:), dA(:,:), w(:), nev(:), nX1(:), iA(:,:)
    integer :: m, n, i, j, k
    double precision :: di, a3, b3, rand
    character (len=10) format_string1, format_string2, format_string3
    open(1, File="INPUT.txt")
    read(1,*) m !Размерность СЛАУ
    read(1,*) n !Число тестов
    allocate(res(n,4), resg(n*100,3), resb(n*100,3), A(m,m), b(m), X(m), nX(m), dX(m), dA(m,m), w(m), nev(m), nX1(m), iA(m,m))
    close(1)
    
    di = 2d0
    CALL RANDOM_SEED()
    do i = 1, n
        call RANDOM_NUMBER(w)
        call RANDOM_NUMBER(X)
        X = (100d0-1d0)*X + 1d0 !Cлучайные корни от 1 до 10
        w = (10d0-1d0)*w + 1d0 !Вектор хаусхолдера
        call genA(A, w, m, di) !Генерируем А с заданым числом обусловленности
        b = MATMUL(A,X) !Получаем правую часть СЛАУ
        
        call QR(A,nX,b,m) !Получаем численное речение nX.      
        dX = nX - X !Разность решений
        nev = MATMUL(A,nX) - b !Невязка
        res(i,1) = di
        res(i,2) = vnorm(dX)
        res(i,3) = vnorm(nev)
        res(i,4) = 0d0
        !write(*,*) resb(i,1)-resb(i,2)
        !call RANDOM_NUMBER(dA)
        dA = A*0.01d0
        !w = b
        iA = A+dA
        CALL QR(iA, nX1, b, m)
        dX = nX - nX1
        a3 = norm2(dX)/norm2(X)
        b3 = ((norm2(dA)/norm2(A))*di/(1d0-di*(norm2(dA)/norm2(A))))
        if (a3 <= b3) res(i,4) = 1d0
        !write(*,*) int(res(i,4)), a3, b3, res(i,1)
        di = di + 2.07d0**(dble(i)/2d0)
    end do
    write(*,*) sum(res(:,4))
    !Проводим тесты на неравенство
    call RANDOM_SEED()
    n=n*100
    
    !Хорошая обусловлнность
    di = 10d0
    do i = 1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(w)
        call genA(A, w, m, di)
        X = (10d0-1d0)*X + 1d0
        w = (10d0-1d0)*w + 1d0  
        call genA(A, w, m, di)
        b = MATMUL(A,X)
        iA = A
        do j =1,m
            do k = 1,m
                call random_number(rand)
                rand = rand*0.000002d0*dble(i)
                iA(j,k) = iA(j,k)*(1d0 + rand)
            end do
        end do
        dA = A-iA
        CALL QR(A,nX,b,m)
        CALL QR(iA,nX1,b,m)
        dX = nX - nX1
        resg(i,1) = norm2(dA)/norm2(A)
        resg(i,2) = norm2(dX)/norm2(X)
        resg(i,3) = di*resg(i,1)/(1d0 - di*resg(i,1))
    end do
    
    !Плохая обусловленность
    di = 1e+10
    
    do i = 1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(w)
        X = (10d0-1d0)*X + 1d0
        w = (10d0-1d0)*w + 1d0  
        call genA(A, w, m, di)
        b = MATMUL(A,X) ! Получили столбец решений b
        iA = A
        do j =1,m
            do k = 1,m
                call random_number(rand)
                rand = rand*0.000002d0*dble(i)
                iA(j,k) = iA(j,k)*(1d0 + rand)
            end do
        end do
        dA = A-iA
        CALL QR(A,nX,b,m)
        CALL QR(iA,nX1,b,m)
        dX = nX - nX1
        resb(i,1) = norm2(dA)/norm2(A) 
        resb(i,2) = norm2(dX)/norm2(nX)
        resb(i,3) = di*resb(i,1)/(1d0 - di*resb(i,1))
        WRITE(*,*) i
    end do
    open(1, File="Data1.csv")
    open(2, File = "Datag.csv")
    open(3, File = "Datab.csv")
    write(1,*) "Cond(A)",";","X-X*",";", "AnX-B*"
    write(2,*) "dA/A",";","dX/X", ";", "e(dX/X)"
    write(3,*) "dA/A",";","dX/X", ";", "e(dX/X)"
    do i=1,n
        if (i<=100) write(1,*) res(i,1),";",res(i,2),";", res(i,3)
        write(2,*) resg(i,1),";",resg(i,2),";", resg(i,3)
        write(3,*) resb(i,1),";",resb(i,2),";", resb(i,3)
    end do
    close(1)
    close(2)
    close(3)
    
    deallocate(res, resg, resb, A, b, X, nX, dX, dA, w, nev,nX1,iA)
    write(*,*) "End of life"
    read(*,*)
    contains
    function vnorm(N) result (norm)
        double precision, dimension(:) :: N
        double precision :: norm
        norm = dsqrt(sum(N**2d0))
    end function
    
    function Anorm(Ai,m) result (norm)
        double precision, dimension(:,:) :: Ai
        double precision :: norm
        integer :: m
        norm = dsqrt(sum(Ai**2d0))
    end function
    end program lab3

