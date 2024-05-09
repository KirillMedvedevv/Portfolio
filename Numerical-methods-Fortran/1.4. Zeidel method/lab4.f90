    program lab4
    use functions
    implicit none
    double precision, allocatable :: A(:,:), b(:), X(:), nX(:), res1(:,:), res2(:,:), D(:)
    double precision :: eps, iter, det, opt
    integer :: i,j, m
    open(1, file="input.txt")
    read(1,*) m
    allocate(A(m,m),b(m), X(m), nX(m), res1(100,3), res2(m-1,2), D(m))
    close(1)
    
    call RANDOM_SEED()
    det = 13d0**m
    call gen2A(A,m,det,opt)
    
    !eps = 0.1d0
    !do i = 1,100
        !call gen2A(A,m,det,opt)
        !CALL RANDOM_NUMBER(X)
        !X = X*10d0 + 1d0
        !b = MATMUL(A,X)
        !CALL Zeidel(A,b,m,nX,eps,X*0d0,iter, opt)
        !res1(i,1) = eps
        !res1(i,2) = norm2(X-nX)
        !res1(i,3) = norm2(MATMUL(A,nX) - b)
        !eps = eps**(1.027d0)
    !end do
    
    iter = 0d0
    det = (10d0)
    eps = 1e-14
    call gen2A(A,m,det,opt)
    do i = 1,m-1
        CALL RANDOM_NUMBER(X)
        X = X*10d0 + 1d0
        b = MATMUL(A,X)
        forall(j=1:m) D(j)= A(j,j)
        opt = 2d0/(maxval(D) + minval(D))
        CALL Zeidel(A,b,m,nX,eps, X*0d0, iter, opt)
        res2(i,1) = det
        res2(i,2) = iter
        det = det/2d0
        A(i,i) = A(i,i)/2d0
        iter = 0d0
        D = 0d0
        write(*,*) "Complete"
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
    
    deallocate(A, b, X, nX,res1,res2,D)
    pause
    end program lab4

