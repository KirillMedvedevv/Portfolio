module methods
    use functions
    implicit none
    contains
     subroutine gen2A(A, m, D) !Генерация основной матрицы
        real(8), dimension(:,:) :: A
        real(8), dimension(:) :: D
        real(8), allocatable :: E(:,:), Q(:,:), w(:)
        character(len=10) :: format_string1
        integer :: i, j, m
        allocate(E(m,m), Q(m,m), w(m))
        write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        
        CALL RANDOM_SEED()
        CALL RANDOM_NUMBER(w)
        CALL RANDOM_NUMBER(A)
        A = 10d0*A + 3d0
        w =  w*(5d0-1d0) + 1d0 !От 1 до 5 
        forall(i=1:m) Q(i,:) = w(i)*w(:)
        Q = Q/((norm2(w))**2d0)
        Q = E - 2d0*Q
        do i = 1,m
            do j = 1,m
                A(i,i) = D(i)
                if (j < i) A(i,j) = 0d0
            end do
        end do
        A = MATMUL(MATMUL(transpose(Q),A),Q)
        !write(*,format_string1) (A(i,:),i=1,m)
        deallocate(E,Q, w)
    end subroutine
    
    subroutine Get_polyV2(A,p,m,B)
        real(8), dimension(:,:) :: A
        real(8), dimension(:) :: p
        real(8), dimension(:,:,:) :: B
        real(8), allocatable :: Ai(:,:), E(:,:)
        integer :: j, k, m,i 
        allocate(Ai(m,m), E(m,m))
        E = 0d0
        p = 0d0
        B=0d0 !Обнуляем, чтобы получать новое
        forall(i=1:m) E(i,i) = 1q0
        
        p(1) = 1d0 !p0
        p(2) = tr(A,m) ! нашли p1
        B(:,:,1) = A-p(2)*E ! нашли B1
        do k=3,m+1 ! С 2 по m
            Ai = MATMUL(A,B(:,:,k-2)) !Построили A2 = A*B1
            p(k) = tr(Ai,m)/(real(k-1,8)) !Нашли p2
            B(:,:,k-1) = Ai - p(k)*E !НАШЛИ B2
        end do
        B(:,:,m) = 0d0
        p(2:m+1) = -1d0*p(2:m+1)
        deallocate(Ai,E)
    end subroutine
    
    subroutine Vector(eval,B, EV,m)
        real(8), dimension(:,:) :: EV
        real(8), dimension(:) :: eval
        real(8), dimension(:,:,:) :: B
        real(8), allocatable :: R(:,:), E(:,:)
        integer :: i,j,m
        allocate(R(m,m),E(m,m))
        R = 0d0
        EV = 0d0 !Обнуляем, чтобы получать новое
        E = 0d0
        forall(i=1:m) E(i,i) = 1d0
        !Последовательно вычисляем собственные векторы
        do i = 1, m
            R = (eval(i)**real(m-1,8))*E
            do j = 2,m
                R = R + (eval(i)**real(m-j,8))*B(:,:,j-1)
            end do
            EV(:,i) = R(:,i)/norm2(R(:,i))
            R = 0d0
        end do
        deallocate(R, E)
    end subroutine
    end module