program NGGD1
    use solvers
    implicit none
    integer :: M, N, info, k, p
    real(8), allocatable :: mesh(:,:,:)
    real(8) :: dx, dt, VNM, xm, tm, eps !параметры схемы
    real(8) :: a, w, U0, time1, time ! Параметры задачи
    real(8) :: M_3, M_6, N_5 
    
    a = 1d0
    w = 1d0
    U0 = 1d0
    
    VNM = 1d0
    tm = 60d0
    eps=1d-2
    xm =10d0 !оценка высоты пограничного слоя
    
    M=100
    dx = xm/(M-1)
    !dt = 0.5d0*VNM*(xm/100)**2
    !dt=xm/100
    dt=dx
    N = dint(tm/dt) + 1
    
    M_3 = dint(3d0/dx) + 1
    M_6 = dint(6d0/dx) + 1
    N_5 = dint(10d0/dt) + 1
    allocate(mesh(M,N,3))
    
    call gen_mesh(mesh, M, N, U0, w, a, dt, dx)
    call implicit3(mesh, M, N, a, dt, dx)
    !call explicit(mesh, M, N, U0, w, a, dt, dx, p)
    !write(*,*) 'M3 N5'
    !write(*,*) mesh(M_3, N_5, :), 100d0*dabs(U(U0, w, a, mesh(M_3, N_5, 1), mesh(M_3, N_5, 2)) - mesh(M_3, N_5, 3))/U(U0, w, a, mesh(M_3, N_5, 1), mesh(M_3, N_5, 2))
    !write(*,*) 'M6 N5'
    !write(*,*) mesh(M_6, N_5, :), 100d0*dabs(U(U0, w, a, mesh(M_6, N_5, 1), mesh(M_6, N_5, 2)) - mesh(M_6, N_5, 3))/U(U0, w, a, mesh(M_6, N_5, 1), mesh(M_6, N_5, 2))
    
    call writefile(mesh, M, N)
    
    deallocate(mesh)
    pause
    contains 
    
    subroutine writefile(mesh, M, N)
        integer :: i,j, N, M
        real(8) :: mesh(M,N,3)
        
        
        open(1, File='An.csv')
        open(2, File='Num.csv')
        !open(3, File='res3.csv')
        write(1, '(5a)') 'x', ';', 't', ';', 'u'
        write(2, '(5a)') 'x', ';', 't', ';', 'u'
        do i = 1, M
            do j = 1, N
                !write(1,*) mesh(i,j, 1), ';', mesh(i,j,2), ';', U(mesh(i,j,1), mesh(i,j,2))
                write(1,*) mesh(i,j, 1), ';', mesh(i,j,2), ';', U(U0, w, a, mesh(i,j,1), mesh(i,j,2))
                write(2,*) mesh(i,j, 1), ';', mesh(i,j,2), ';', mesh(i,j,3)
                !write(3,*) mesh(i,j,3), ';', U(U0, w, a, mesh(i,j,1), mesh(i,j,2))
            end do
        end do
        close(1)
        close(2)
        !close(3)
    end subroutine
end program NGGD1
