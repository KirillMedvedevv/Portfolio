module Solvers
    implicit none
    contains
    
    subroutine exp_anti_potok(M, N, dx, dt, c, mesh)
        integer :: M, N, i, j
        real(8) :: c, mesh(N, M, 3), dx, dt, h
        
        do i = 1, N-1
            !���������� �������� � ������������� ����� ��������� �� ������ �
            do j = 2, M-1
                h = (dt/(2d0*dx))*(-(c + dabs(c))*mesh(i, j-1, 3) + 2d0*dabs(c)*mesh(i, j, 3) + (c - dabs(c))*mesh(i, j+1, 3))
                mesh(i+1, j, 3) = mesh(i, j, 3) - h
            end do
            
            if (c > 0) then
                !��������� ��������� ����
                h = (dt/(2d0*dx))*(-(c + dabs(c))*mesh(i, M-1, 3) + 2d0*dabs(c)*mesh(i, M, 3))
                mesh(i+1, M, 3) = mesh(i, M, 3) - h
                !�� ������������� ��� ������� ����
                mesh(i+1, 1, 3) = mesh(i, 1, 3) - c*(dt/dx)*(mesh(i, 1, 3) - mesh(i, M-1, 3)) !�� �������������
            else
                !��������� ������ ����
                h = (dt/(2d0*dx))*(2d0*dabs(c)*mesh(i, 1, 3) + (c - dabs(c))*mesh(i, 2, 3))
                mesh(i+1, 1, 3) = mesh(i, 1, 3) - h
                !�� ������������� ��� ���������� ����
                h = (dt/(2d0*dx))*(2d0*dabs(c)*mesh(i, M, 3) + (c - dabs(c))*mesh(i, 2, 3))
                mesh(i+1, M, 3) = mesh(i, M, 3) - h
            end if
        end do
    end subroutine
    
    subroutine imp_anti_potok(M, N, dx, dt, c, mesh)
        integer :: M, N, i, j
        real(8) :: c, mesh(N, M, 3), dx, dt, VNM
        !��������� ���� � dgels
        real(8) :: A(M,M), D(M)
        real(8), allocatable :: WORK(:)
        integer :: lwork, info
        VNM = c*dt/dx
        write(*,*) VNM
        
        !��������� ������� A
        call genA(A, VNM, M)
        
        !���������� ��������� LWORK
        D = mesh(1,:,3)
        allocate(WORK(1))
        call dgels('N', M, M, 1, A, M, D, M, WORK, -1, INFO)
        lwork = dint(work(1))
        deallocate(WORK)
        allocate(WORK(lwork))
        
        do i = 1, N-1
            call genA(A, VNM, M)
            D = mesh(i,:,3)
            call dgels('N', M, M, 1, A, M, D, M, WORK, lwork, INFO)
            mesh(i+1, :, 3) = D 
        end do
        
        deallocate(work)
    end subroutine
    
    subroutine genA(A, VNM, M)
        integer :: M, i
        real(8) :: A(M,M), VNM, VNM_plus, VNM_minus
        
        VNM_plus = -(VNM+dabs(VNM))/2d0
        VNM_minus = (VNM - dabs(VNM))/2d0
        A = 0d0       
        forall(i=1:M) A(i,i) = 1d0+dabs(VNM)
        if (VNM > 0) then
            A(1,M-1) = VNM_plus
            forall(i=1:M-1) A(i+1,i) = VNM_plus
        else
            A(M,1) = VNM_minus
            forall(i=2:M) A(i-1,i) = VNM_minus
        end if
    end subroutine
end module