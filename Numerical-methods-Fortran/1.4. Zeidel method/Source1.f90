module functions
implicit none
contains
   subroutine genA(A, m, det) !Генерация основной матрицы
        double precision, dimension(:,:) :: A
        double precision :: det
        double precision, allocatable :: E(:,:), Q(:,:), D(:), w(:)
        integer :: i, j, m
        character (len=10) format_string1
        allocate(E(m,m), Q(m,m), D(m), w(m))
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(w)
        w =  w*2d0 + 1d0
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/(norm2(w)**2d0)
        Q = E - 2d0*Q
        
        !диагональное преобладание
        forall(i=1:m) D(i) = det**(1d0/m)
        forall(i=1:4) D(i) = D(i)/2d0
        forall(i=5:8) D(i) = D(i)*2d0
        CALL RANDOM_NUMBER(A)
        do i=1,m
            do j=1,m
                if(j<i) A(i,j) = 0d0
            end do
        end do
        forall(i=1:m) A(i,i) = D(i)
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        !write(*,format_string1) (A(i,:),i=1,m)
        !write(*,*) "End"
        do i = 1,m
            do j =1,m
                if(i/=j) A(i,j) = A(i,j)/100d0
            end do
        end do
        deallocate(E,Q,D,w)
    end subroutine
    
    subroutine Zeidel(A, b, m, nX, eps, x0, iter, opt)
    double precision, dimension(:) :: nX, b, x0
    double precision, dimension(:,:) :: A
    integer :: m, i,j
    double precision :: eps, iter, nC, alp,bet,prev, opt
    character (len=10) format_string1
    double precision, allocatable :: x1(:), x2(:), dup(:), E(:,:), C(:,:), g(:)
    allocate(x1(m),x2(m),dup(m), E(m,m), C(m,m),g(m))
    E = 0d0
    forall(i=1:m) E(i,i) = 1d0
    write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
   
    !a1 = 2d0/(a1/2d0 + a1*2d0)
    
    C = E - A*opt
    g = b*opt
    !Вычислим nC
    prev = 0d0
    do i = 2,m
        alp = 0d0
        bet = 0d0
        do j = 1,i-1
            alp = alp + dabs(C(i,j))
        end do
        do j = i,m
            bet = bet + dabs(C(i,j))
        end do
        nC = bet/(1d0-alp)
        if(nC>prev) prev = nC
    end do
    
    nC = prev
    x1 = x0
    iter = 0d0
    do while(1)
        dup = x1
        do i =1,m
            x2(i) = dot_product(C(i,:),x1) + g(i)
            x1(i) = x2(i)
        end do
        if(maxval(dabs(x2 - dup)) < nC*eps) EXIT
        x1 = x2
        iter = iter + 1d0
    end do
    nX = x2
    deallocate(x1,x2,dup, E, C,g)
    end subroutine
    
    subroutine gen2A(A, m, det, opt) !Генерация основной матрицы
        double precision, dimension(:,:) :: A
        double precision :: det, opt
        double precision, allocatable :: E(:,:), Q(:,:), D(:), w(:)
        integer :: i, j, m
        character (len=10) format_string1
        allocate(E(m,m), Q(m,m), D(m), w(m))
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(w)
        w =  w*(5d0-1d0) + 1d0 !От 1 до 5
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/((norm2(w))**2d0)
        Q = E - 2d0*Q
        !Генерация диагональной матрицы
        CALL RANDOM_NUMBER(D)
        A = 0d0
        D = D*(3d0 - 1d0) + 1d0
        D(1:m-1) = D(1:m-1)*dexp((dlog(det) - sum(dlog(D)))/dble(m-1))
        opt = 2d0/(maxval(D) + minval(D))
        forall(i=1:m) A(i,i) = D(i)
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        do i = 1,m
            do j = 1,m
                if (i /= j) A(i,j) = A(i,j)/100d0
            end do
        end do
        deallocate(E,Q,D,w)
    end subroutine
end module