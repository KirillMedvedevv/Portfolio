module funcLR
    use funcLU
    implicit none
    contains
    
    subroutine vector(A,G,eval,m,eps)
    real(16), dimension(:) :: eval
    real(16), dimension(:,:) :: G, A
    real(16), allocatable :: L(:,:), U(:,:), E(:,:), y(:), x(:), xprev(:)
    real(16) :: eps, sum, err,norm
    integer :: i,j,k,m
    character(len=10) :: format_string1
    write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
    allocate(L(m,m),U(m,m),E(m,m),y(m),x(m),xprev(m))
    E = 0q0
    forall(i=1:m) E(i,i) = 1q0
    call RANDOM_SEED()
    do i = 1,m
        CALL LU(A - (eval(i) + 0.005q0)*E, m, L, U)
        call RANDOM_number(xprev)
        forall(j=1:m) xprev(j) = xprev(j)/maxval(qabs(xprev))
        err = 1000q0
        do while(err > eps)
            !Решаем СЛАУ Ly = xprev
            y = 0d0
            do j = 1,m
                sum = 0q0
                do k=1,j-1
                    sum = sum + L(j,k)*y(k)
                end do
                y(j) = (xprev(j) - sum)/L(j,j)
            end do
            !Решаем Слау Ux = y
            x = 0q0
            do j = m,1,-1
                sum = 0q0
                do k=j+1,m
                    sum = sum + U(j,k)*x(k)
                end do
                x(j) = (y(j) - sum)/U(j,j)
            end do
            
            forall(k=1:m) x(k) = x(k)/maxval(qabs(x))
            xprev = xprev*x(1)/xprev(1)
            err = maxval(qabs(xprev-x))
            xprev = x
        end do
        G(:,i) = x
    end do
    forall(i=1:m) G(:,i) = G(:,i)/norm2(G(:,i))
    deallocate(L,U,E,y,x,xprev)
    end subroutine
    
    subroutine LR(A, eg, eps, m, count)
    integer :: i, j, k, m, count, horosho
    real(16), dimension(:) :: eg
    real(16), dimension(:,:) :: A
    real(16), allocatable :: E(:,:), Ak(:,:), L(:,:), U(:,:), matrix(:), F1(:), F2(:)
    real(16) :: shift, eps, err,delta,b
    character(len=10) :: format_string1
    write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
    allocate(E(m,m), Ak(m,m), L(m,m), U(m,m),matrix(m), F1(m), F2(m))
    
    E = 0q0
    forall(i=1:m) E(i,i) = 1q0
    
    err = 1000q0
    Ak = A
    count = 0
    do while (err >= eps)
        !delta = (Ak(m-1,m-1) - Ak(m,m))/2q0
        !shift = Ak(m,m) - qsign(1q0,delta)*((Ak(m,m-1))**2q0)/(qabs(delta) + qsqrt((Ak(m,m-1))**2q0 + delta**2q0)) !Сдвиг вилкстона
        shift = Ak(m,m) ! Сдвиг 
        !shift = 0q0 ! БЕЗ СДВИГА
        
        call LU(Ak - shift*E, m, L, U)
        Ak = MATMUL(U,L) + shift*E 
        !Условие выхода
        err = 0q0
        !forall(i=1:m) F2(i) = Ak(i,i)
        !err = maxval(qabs(F2-F1))
        matrix = 0q0
        do i = 1,m
            do j = 1,i-1
                matrix(i) = matrix(i) + Ak(i,j)**2q0
            end do
        end do
        matrix = qsqrt(matrix)
        err = (sum(matrix))
        count = count + 1
    end do
    
    forall(i=1:m) eg(i) = Ak(i,i)
    deallocate(E, Ak, L, U,matrix, F1,F2)
    end subroutine
    
    subroutine SortingA(eval,B,m)
    real(16), dimension(:,:) :: B
    real(16), dimension(:) :: eval
    real(16) :: t
    integer :: i,j,m
    double precision, allocatable :: buff(:)
    allocate(buff(m))
    
    do i=m-1,1,-1
        do j=1,i
        if (eval(j).gt.eval(j+1)) then
            t=eval(j)
            buff = B(:,j)
            eval(j)=eval(j+1)
            B(:,j) = B(:,j+1)
            eval(j+1)=t
            B(:,j+1) = buff
        end if
        end do
    end do
    deallocate(buff)
    end subroutine
end module