program lab51
    use methods
    use functions
    use polyroots_module
    use sorting_module
    implicit none
        integer :: istat,i,j,m
        real(8), allocatable :: A(:,:), P(:), D(:), realsol(:), complexsol(:), res(:,:), B(:,:,:), EV(:,:)
        real(8) :: otd
        m = 10
        allocate(A(m,m),D(m),P(m+1), realsol(m),complexsol(m),res(100,3),B(m,m,m),EV(m,m))
        
        !Генерируем рандома
        CALL RANDOM_SEED()
        otd = 0.93d0
        call random_number(D)
        call sort_ascending(D)
        D = 18d0*D + 1d0
        do i = 1,100
            !Задаём отделимость
            D(3) = D(2) + (D(3) - D(2))*otd
            !D(6) = D(5) + (D(6) - D(5))*otd
            CALL gen2A(A,m,D)
            CALL get_PolyV2(A,P,m,B)
            call rpoly(P, m, realsol, complexsol, istat)
            call sort_ascending(realsol) !Сортируем полученные собственные числа
            call vector(realsol,B,EV,m)
            res(i,1) = dlog10(dabs(D(3) - D(2))) 
            res(i,2) = dlog10(dabs(realsol(3)-D(3)))
            res(i,3) = dlog10(norm2(matmul(A,EV(:,3)) - realsol(3)*EV(:,3)))
        end do
        
        open(1,file="Result.csv")
        write(1,*) "otd", ";", "||L-L*||",";","||Ax-eval*x||"
        do i =1,100
            write(1,fmt="(f9.5,a,f9.5,a,f9.5)") res(i,1),";", res(i,2),";", res(i,3)
        end do
        close(1)
        DEALLOCATE(A,D,P, realsol,complexsol,res,B, EV) 
    pause
    end program lab51

