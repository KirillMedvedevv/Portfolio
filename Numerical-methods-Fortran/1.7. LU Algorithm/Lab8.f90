    program Lab8
    use funcLR
    use funcLU
    use sorting_module 
    implicit none
    real(16), allocatable :: A(:,:),  D(:), res(:,:), G(:,:), Eg(:), Q(:,:)
    real(16) :: eps, normi
    character(len=10) :: format_string1
    integer :: i, j, m, count
    m = 13
    
    write(format_string1,"(a,i2,a)") "(",m,"f9.5)"
    allocate(A(m,m), D(m), Eg(m), G(m,m), res(14,5), Q(m,m))
    
    call RANDOM_SEED()
    
    eps = 0.1q0
    !call RANDOM_NUMBER(D)
    !D = (D*10q0 + 11q0)
    D(m) = 5q0
    do i = m-1,1,-1
        D(i) = D(i+1) + 5q0
    end do
    call sort_ascending(D)
    CALL GEN2A(A,m, D, Q)
    do i = 1, 14
        CALL LR(A, Eg,eps,m, count)
        CALL sort_ascending(Eg)
        call vector(A,G,Eg,m,eps)
        
        do j = 1, m
            if(G(1,j)*Q(1,j) < 0) G(:,j) = G(:,j)*(-1d0)
        end do
        
        res(i,1) = qlog10(eps)
        res(i,2) = qlog10(maxval(qabs(Eg(:) - D(:))))
        res(i,3) = 0q0
        do j = 1, m
            normi = maxval(qabs(matmul(A,G(:,j)) - Eg(j)*G(:,j)))
            if (normi > res(i,3)) res(i,3) = normi
        end do
        res(i,3) = qlog10(res(i,3))
        res(i,4) = 0q0
        do j = 1, m
            normi = maxval(qabs(G(:,j) - Q(:,j)))
            if (normi > res(i,4)) res(i,4) = normi
        end do
        res(i,4) = qlog10(res(i,4))
        res(i,5) = (real(count,16))
        eps = eps*0.1q0
        write(*,*) i
    end do
    
    open(1,file="Result.csv")
    write(1,*) "eps", ";", "||L-L*||",";","||Ax-eval*x||",";", "||x-x*||", ";", "N", ";", "eps1"
    do i =1,14
        write(1,fmt="(f13.5,a, f13.5,a,f13.5,a, f13.5, a, f13.7, a, f13.7)") res(i,1),";", res(i,2),";", res(i,3),";", res(i,4), ";", res(i,5), ";", res(i,1)
    end do
    close(1)
    do i = 1, 100
    end do
    deallocate(A, D, G, Eg, res, Q)
    pause
    end program Lab8

