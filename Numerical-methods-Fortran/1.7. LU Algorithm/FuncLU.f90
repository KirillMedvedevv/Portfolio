module funcLU
    implicit none
    contains
    
    subroutine LU(A,m,L,U)
    real(16), dimension(:,:) :: L, U, A
    integer :: i, j, k, m
    real(16) :: sum
    
    U = 0q0
    L = 0q0
    forall(i=1:m) L(i,i) = 1q0
    do i =1,m
        do j = 1,m
            if (i<=j) then
                sum = 0q0
                do k = 1,i
                    sum = sum + L(i,k)*U(k,j)
                end do
                U(i,j) = A(i,j) - sum
            else
                sum = 0q0
                do k = 1,j
                    sum = sum + L(i,k)*U(k,j)
                end do
                L(i,j) = (A(i,j) - sum)/U(j,j)
            end if
        end do
    end do
    end subroutine
    
    subroutine gen2A(A, m, D, Q) !Генерация основной матрицы
        real(16), dimension(:,:) :: A, Q
        real(16), dimension(:) :: D
        real(16), allocatable :: E(:,:), w(:)
        character(len=10) :: format_string1
        integer :: i, j, m
        allocate(E(m,m), w(m))
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        E = 0q0
        forall(i=1:m) E(i,i) = 1q0
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(w)
        A = 0q0
        Q=0q0
        w =  w*(20q0) + 20q0 !От 1 до 5 
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/((norm2(w))**2q0)
        Q = E - 2q0*Q
        
        forall(i=1:m) A(i,i) = D(i)
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        !write(*,format_string1) (A(i,:),i=1,m)
        deallocate(E, w)
    end subroutine
    
    function mnorm(A,m) result(norm)
    real(16), dimension(:,:) :: A
    real(16) :: norm
    integer :: j,m
    norm = 0q0
    do j = 1,m
        if (maxval(qabs(A(:,j))) > norm) norm = maxval(qabs(A(:,j)))
    end do
    end function
end module