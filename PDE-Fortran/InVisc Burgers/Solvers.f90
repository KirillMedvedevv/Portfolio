module Solvers
    implicit none
    contains
  

    subroutine exp_anti_potok(M, N, dx, dt, c, mesh)
        integer :: M, N, i, j
        real(8) :: c, mesh(N, M, 3), dx, dt, h, f1, f2, f3
        
        do i = 1, N-1
            !Вычисление значений в промежуточных узлах независит от выбора с
            do j = 2, M-1
                f1 = mesh(i, j-1, 3)**2
                f2 = mesh(i, j, 3)**2
                f3 = mesh(i, j+1, 3)**2
                c = mesh(i, j, 3)
                !Сразу учитываем знак c
                if (c > 0) then 
                    h = (dt/(2d0*dx))*(f2 - f1)
                else
                    h = (dt/(2d0*dx))*(f3 - f2)
                end if
                mesh(i+1, j, 3) = mesh(i, j, 3) - h
            end do
            
            !Обработка крайних узлов
            !Первый узел
            c = mesh(i, 1, 3)
            if (c >= 0) then
                f1 = mesh(i, M-1, 3)**2
                f2 = mesh(i, 1, 3)**2
                h = (dt/(2d0*dx))*(f2 - f1)
            else if (c<0) then
                f2 = mesh(i, 1, 3)**2
                f3 = mesh(i, 2, 3)**2
                h = (dt/(2d0*dx))*(f3 - f2)
            else
                h=0d0
            end if
            mesh(i+1, 1, 3) = mesh(i, 1, 3) - h
            
            !Последний узел
            c = mesh(i, M, 3)
            if (c >= 0) then
                f1 = mesh(i, M-1, 3)**2
                f2 = mesh(i, M, 3)**2
                h = (dt/(2d0*dx))*(f2 - f1)
            else if(c < 0) then
                f2 = mesh(i, M, 3)**2
                f3 = mesh(i, 2, 3)**2
                h = (dt/(2d0*dx))*(f3 - f2)
            else
                h=0d0
            end if
            mesh(i+1, M, 3) = mesh(i, M, 3) - h
        end do
    end subroutine
    
    subroutine imp_anti_potok(M, N, dx, dt, mesh)
        integer :: M, N, i, j
        real(8) :: VNM(M), mesh(N, M, 3), dx, dt
        !Параметры СЛАУ и dgels
        real(8) :: A(M,M), D(M)
        real(8), allocatable :: WORK(:)
        integer :: lwork, info
        !
        !Вычисление параметра LWORK
        VNM = 1d0*dt/dx
        call genA(A, VNM, M)
        D = mesh(1,:,3)
        allocate(WORK(1))
        call dgels('N', M, M, 1, A, M, D, M, WORK, -1, INFO)
        lwork = dint(work(1))
        deallocate(WORK)
        allocate(WORK(lwork))
        
        do i = 1, N-1
            VNM = mesh(i, :, 3)*dt/(2d0*dx)
            call genA(A, VNM, M)
            D = mesh(i,:,3)
            call dgels('N', M, M, 1, A, M, D, M, WORK, lwork, INFO)
            mesh(i+1, :, 3) = D 
        end do
        
        deallocate(work)
    end subroutine
    
    subroutine genA(A, VNM, M)
        integer :: M, i
        real(8) :: A(M,M), VNM(M)
        
        A = 0d0       
        forall(i=1:M) A(i,i) = 1d0+dabs(VNM(i))
        do i = 2, M-1
            if (VNM(i) >= 0d0) then
            !A(1,M) = VNM_plus
            A(i,i-1) = -VNM(i)
            else
            !A(M,1) = VNM_minus
            A(i,i+1) = VNM(i)
            end if
        end do
        
        !Учёт граничных условий
        if (VNM(1) >= 0d0) then
            A(1,M-1) = -VNM(1)
        else
            A(1, 2) = VNM(1)
        end if
        
        if (VNM(M) >= 0d0) then
            A(M, M-1) = -VNM(M)
        else
            A(M,2) = +VNM(M)
        end if
    end subroutine
end module