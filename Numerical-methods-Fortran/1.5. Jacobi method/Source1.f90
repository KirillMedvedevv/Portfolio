module functions
implicit none
    contains
    subroutine genA(A, m, det) !Генерация матрицы с заданным определителем
        double precision, allocatable :: E(:,:), Q(:,:), D(:), w(:)
        double precision, dimension(:,:) :: A
        integer :: i, j, m
        double precision :: det
        character (len=10) format_string1
        
        allocate(E(m,m), Q(m,m), D(m), w(m))
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(w)
        w =  w*2d0  + 1d0
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/(norm2(w)**2d0)
        Q = E - 2d0*Q
        
        forall(i=1:m) D(i) = det**(1d0/m)
        forall(i=1:6) D(i) = D(i)/dble(i)
        forall(i=7:12) D(i) = D(i)*dble(i-6)
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
        !диагональное преобладание
        do i = 1,m
            do j =1,m
                if(i/=j) A(i,j) = A(i,j)/10d0
            end do
        end do
        deallocate(E,Q,D,w)
    end subroutine
    
    subroutine Zeidel(A, b, m, nX, eps, x0, iter, alpha)
        character (len=10) format_string1
        integer :: i, j, m
        double precision, allocatable :: x1(:), x2(:), dup(:), E(:,:), C(:,:), g(:)
        double precision, dimension(:,:) :: A
        double precision, dimension(:) :: nX, b, x0
        double precision :: eps, iter, alpha, alp, bet, prev, k , norm
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        allocate(x1(m),x2(m), dup(m), E(m,m), C(m,m), g(m))
    
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        C = E - A*alpha
        g = b*alpha
    
        !Вычислим поправку к условию остановки
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
            k = bet/(1d0-alp)
            if(k>prev) prev = k
        end do
        k=prev
        write(*,*) k
        !Метод Зейделя
        x1 = x0
        iter = 0d0
        do while(1)
            dup = x1
            do i =1,m
                x2(i) = dot_product(C(i,:),x1) + g(i)
                write(*,*) x2(i)
                x1(i) = x2(i)
            end do
            if(maxval(dabs(x2 - dup)) < k*eps) EXIT
            x1 = x2
            iter = iter + 1d0
        end do
        nX = x2
        deallocate(x1,x2,dup, E, C,g)
    end subroutine
end module