program lab4
    use functions
    implicit none
    double precision, allocatable :: A(:,:), b(:), X(:), nX(:), res1(:,:), res2(:,:)
    double precision :: eps, iter, alpha, det
    integer :: i,j, m
    open(1, file="Text1.txt")
    read(1,*) m
    allocate(A(m,m),b(m), X(m), nX(m), res1(100,3), res2(m-1,2))
    close(1)
    
    call RANDOM_SEED()
    det = 12d0**m
    call genA(A, m, det)
    write(*,*) A
    alpha = 2d0*((det**(1d0/m))/6d0 + 6d0*det**(1d0/m))
    !Точность, ошибка, невязка
    eps = 0.1d0
    do i = 1,100
        CALL RANDOM_NUMBER(X)
        X = X*10d0 + 5d0
        b = MATMUL(A,X)
        
        CALL Zeidel(A,b,m,nX,eps,X/X,iter,alpha)
        write(*,*) "Complete"
        res1(i,1) = eps
        res1(i,2) = norm2(X-nX)
        res1(i,3) = norm2(MATMUL(A,nX) - b)
        eps = eps**(1.027d0)
    end do
    write(*,*) "GOOD"
    iter = 0d0
    det = (12d0**m)
    eps = 1e-14
    alpha = 2d0*((det**(1d0/m))/6d0 + 6d0*det**(1d0/m))
    call genA(A,m,det)
    do i = 1,m-1
        det = det*(4d0**dble(i-1))
        alpha = 2d0*((det**(1d0/m))/6d0 + 6d0*det**(1d0/m))
        CALL RANDOM_NUMBER(X)
        X = X*4d0 + 1d0
        b = MATMUL(A,X)
        CALL Zeidel(A,b,m,nX,eps,X+10d0,iter,alpha)
        res2(i,1) = det
        res2(i,2) = iter
        A(i,i) = A(i,i)*4d0
    end do
    
    open(1, file="res1.csv")
    open(2, file="res2.csv")
    write(1,*) "Accur",";", "||X-nX||", ";","||AnX - b||"
    write(2,*) "det", ";", "iter"
    do i = 1,100
        write(1,*) res1(i,1), ";", res1(i,2), ";", res1(i,3)
    end do
    do i = 1,m-1
        write(2,*) res2(i,1), ";", res2(i,2)
    end do
    close(1)
    close(2)
    
    deallocate(A, b, X, nX,res1,res2)
    pause
end program lab4

