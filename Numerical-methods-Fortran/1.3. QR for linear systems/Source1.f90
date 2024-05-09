module functions
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
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        !norm = 0d0
        !do i = 1,m
        !    norm = norm + w(i)**2d0
        !end do
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/(norm2(w)**2d0)
        Q = E - 2d0*Q
        CALL RANDOM_SEED() 
        
        CALL RANDOM_NUMBER(D)
        CALL RANDOM_NUMBER(A)
        D = ((cond-1d0)*D)+1d0
        D(1) = 1d0
        D(2) = cond
        do i=1, m
            do j = 1,m
                if(j<i) then 
                    A(i,j) = 0d0
                end if
            end do
            !A(i,i) = D(i)
        end do
        forall(i=1:m) A(i,i) = D(i)
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        A = A*100d0
        !write(*, format_string1) (A(i,:), i=1,m)
        deallocate(E,Q,D)
    end subroutine
    subroutine QR(A,nX,b, m)
    double precision, dimension(:,:) :: A
    double precision, dimension(:) :: nX, b
    double precision :: a1, b1
    double precision, allocatable :: Q(:,:), matrix(:,:), R(:,:), vec(:), w(:)
    integer :: m, i, j,k
    character (len=10) format_string1
    write(format_string1,"(a,i2,a)") "(",m+1,"f9.5)"
    allocate(Q(m,m), matrix(m,m+1),vec(m), R(m,m),w(m))
    R = 0d0
    !Q = A
    !Q(:,1) = Q(:,1)/(dsqrt(sum(Q(:,1)**2d0)))
    !R(1,1) = dot_product(A(:,1),A(:,1))
    do j=1,m                               
        vec(:)=A(:,j)!Вектор подлежащий ортогонализации                    
        do i=1,j-1
            R(i,j) = dot_product(Q(:,i),vec(:))
            vec(:)=vec(:)-R(i,j)*Q(:,i)/dot_product(Q(:,i),Q(:,i))
        end do  
        Q(:,j)= vec(:)
        R(j,j) = (dot_product(vec,vec))
        w(j) = sum(Q(:,j)*b)
    end do
    matrix(:,m+1) = w
    forall(i=1:m) matrix(i,:m) = R(i,:)
    forall(i=1:m) matrix(i,:) = matrix(i,:)/matrix(i,i)
    do i = 1, m
        do j = i+1, m
            matrix(m - j+1,:) = matrix(m - j+1,:) - matrix(m-i+1,:)*matrix(m - j+1,m-i+1)
        end do
    end do
    !write(*, format_string1) (matrix(i,:), i=1,m)
    !write(*,*) "End of life"
    nX = matrix(:,m)
    deallocate(Q,R,matrix,vec,w)
    end subroutine
end module