program lab2
use procedures
implicit none
    double precision, allocatable :: res(:,:), A(:,:), b(:), X(:), nX(:), iA(:,:), w(:), db(:), dx(:), dA(:,:)
    integer :: m, i, j, n, k
    double precision :: di, dn, a3,b3
    logical :: nan
    character (len=10) format_string1, format_string2, format_string3
    nan = .FALSE.
    open(1, File="INPUT.txt")
    read(1,*) m !Размерность СЛАУ
    read(1,*) n !Число тестов
    dn = DBLE(n)
    allocate(res(n,4), b(m), A(m,m), X(m), nX(m), iA(m,m), w(m),db(m),dx(m),dA(m,m))
    close(1)
    write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
    write(format_string2,"(a,i2,a)") "(",2*m,"f9.5)"
    write(format_string3,"(a,i2,a)") "(",4,"f10.5)"
    b = 0d0
    di = 2d0 !cond(A)
    dn = 1d0
    A = 0d0
    !Lab2 test
    call RANDOM_SEED()
    do i=1,n
        call RANDOM_NUMBER(X)
        call RANDOM_NUMBER(w)
        X = (10d0-1d0)*X + 1d0 !Пусть корни -  числа от 1 до 10
        w = DBLE(floor(w*(4d0-1d0) + 1d0))
        call genA(A, w, m, di)
        b = MATMUL(A,X) ! Получили столбец решений b
        iA = A
        CALL Jordan2(iA,m, nan)
        !На этом этапе проверим наличие NaN и Infinity
        if(nan) then
            EXIT
        end if
        
        nX = MATMUL(iA,b) !Численное решение
        w = MATMUL(A,nX) - b !Столбец невязка
        res(i,3) = vnorm(w) !Норма невязки
        res(i,1) = di !Число обусловненности
        w = nX - X !Фактическая ошибка
        res(i,2) = vnorm(w)
       !Зададим возмущение входных данных
        call RANDOM_NUMBER(w)
        db = w*b*0.01d0
        dx = nX - MATMUL(iA, b+db)
        if(vnorm(dx)/vnorm(nX) <= di*vnorm(db)/vnorm(b)) then
            res(i,4) = 1d0
        else
            res(i,4) = 0d0
        end if
        X = 0d0
        nX = 0d0
        iA = 0d0
        b=0d0
        w = 0d0
        dn = dn+1d0
        db = 0d0
        dx = 0d0
        di = (di + 2.07d0**(dn/2d0))
    end do
    !Непроизводить запись если был NaN
    if(nan) then
        write(*,*) "NAN ERROR"
    else
        open(1, File="Data.csv")
        write(1,*) "Cond(A)",";","X-X*",";", "AnX-B*",";"
        do i=1,n
            write(1,*) (res(i,1)),";",res(i,2),";", res(i,3)
            write(*,*) int(res(i,4))
        end do
        close(1)
    end if
    
    deallocate(res, A, X, b, nX, iA, w, dx, db, dA)
    read(*,*)
    contains
function vnorm(N) result (norm)
    double precision, dimension(:) :: N
    double precision :: norm
    norm = dsqrt(sum(N**2d0))
end function
end program lab2

