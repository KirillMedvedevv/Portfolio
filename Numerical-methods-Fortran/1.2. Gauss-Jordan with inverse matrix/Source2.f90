module procedures
    implicit none
    contains
    subroutine genA(A, w, m, cond) !Генерация основной матрицы
        double precision, dimension(:,:) :: A
        double precision, dimension(:) :: w
        double precision :: norm, cond
        double precision, allocatable :: E(:,:), Q(:,:), D(:) 
        character (len=10) format_string1
        integer :: i, j, m
        allocate(E(m,m), Q(m,m), D(m))
        write(format_string1,"(a,i2,a)") "(",m,"f16.5)"
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        norm = 0d0
        do i = 1,m
            norm = norm + w(i)**2d0
        end do
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/norm
        Q = E - 2d0*Q
        !write(*,*) w
        !write(*,format_string1) (Q(i,:), i=1,m)
        CALL RANDOM_SEED() 
        CALL RANDOM_NUMBER(D)
        A = 1d0
        D = ((cond-1d0)*D)+1d0
        D(1) = 0d0
        D(2) = 0d0
        D(3) = cond
        do i=1, m
            do j = 1,m
                if(j<i) then 
                    A(i,j) = 0d0
                end if
            end do
        end do
        forall(i=1:m) A(i,i) = D(i)
        
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        A = A*1d0
        
        deallocate(E,Q,D)
    end subroutine
    
    subroutine Jordan2(Ai, m, nan) 
        integer :: i, j, m, k !k - координата вектора, который был выбран ведущим
        double precision, dimension(:,:) :: Ai
        double precision, allocatable :: vector(:), matrix(:,:)
        double precision :: a, b
        logical :: nan
        character (len=10) format_string1
        write(format_string1,"(a,i2,a)") "(",1*m,"f9.5)"
        allocate(vector(2*m), matrix(m,2*m))
        matrix = 0d0
        forall(j=1:m) matrix(j,:m) = Ai(j,:)
        forall(j=1:m) matrix(j,m+j) = 1d0 ! Cформировали расширенную матрицу
        
        do i = 1, m
            a = dabs(matrix(i,i))
            k = i
            do j = 1,m
                if (dabs(matrix(j,i)) > a) then
                    a = dabs(matrix(j,i))
                    k = j
                end if
            end do
            if (k>i) then 
                vector = matrix(i,:)
                matrix(i,:) = matrix(k,:)
                matrix(k,:) = vector
            end if

            matrix(i,:) = matrix(i,:)/matrix(i,i)
            if(isnan(matrix(i,i))) then
                nan = .TRUE.
                RETURN
            end if
            do j = 1, m
                if (j /= i) matrix(j,:) = matrix(j,:) - matrix(i,:)*matrix(j,i)
            end do
        end do
        forall(j=1:m) Ai(j,:) = matrix(j,m+1:)
        deallocate(vector,matrix)
    end subroutine
    
   
    end module