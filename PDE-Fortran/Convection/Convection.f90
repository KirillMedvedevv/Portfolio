program Convection
use Solvers
implicit none
    real(8), parameter :: pi = 4d0*datan(1d0)
    real(8) :: m1, L, c, tm
    integer :: M, N
    real(8) :: dx, dt, VNM
    real(8), allocatable :: mesh(:,:,:)
    c=1d0
    L = 1d0
    
    M = 201
    dx = L/(M - 1d0)
    
    VNM = 0.5d0
    
    tm = 15d0
    dt = dx*VNM/dabs(c)
    N = dint(tm/dt) + 1
    
    !Ğåøàåì Óğàâíåíèÿ
    allocate(mesh(N, M, 3)) ! (t, x, u)
    
    !sin m=5, 14 imp, exp
    
    L=2d0/5d0
    dx = L/(M-1d0)
    dt = dx*VNM/dabs(c)
    call genmesh(0, 5)
    call exp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='sin_5_EXP.csv')
    call writefile(1)
    close(1)
    
    call genmesh(0, 5)
    call imp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='sin_5_IMP.csv')
    call writefile(1)
    close(1)
    
    L=1d0
    dx = L/(M-1d0)
    dt = dx*VNM/dabs(c)
    call genmesh(0, 14)
    call exp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='sin_14_EXP.csv')
    call writefile(1)
    close(1)
    
    call genmesh(0, 14)
    call imp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='sin_14_IMP.csv')
    call writefile(1)
    close(1)
    
    !Ïğÿìîóãîëüíûé ñèãíàë

    call genmesh(1, 14)
    call imp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='prim_IMP.csv')
    call writefile(1)
    
    call genmesh(1, 14)
    call exp_anti_potok(M, N, dx, dt, c, mesh)
    open(1, file='prim_EXP.csv')
    call writefile(1)
    close(1)
    
    close(1)
    deallocate(mesh)
    pause
    contains
    
    subroutine genmesh(logic, m1)
        integer :: i, j, logic, m1
        !Ñåòêà
        do i = 1, N
            do j = 1, M
                mesh(i, j, 1) = dt*(i-1d0)
                mesh(i, j, 2) = dx*(j-1d0)
                mesh(i, j, 3) = 0d0
            end do
        end do
        
        !ÍÓ
        do j = 1, M
            if (logic==0) then
                mesh(1, j, 3) = sig1(dble(m1), L, mesh(1, j, 2))
            else
                mesh(1, j, 3) = sig2(L, mesh(1, j, 2))
            end if
        end do
    end subroutine
    
    function sig1(m, L, x) result(y)
        real(8) :: m, L, x, y
        y = dsin(pi*m*x)
    end function 
    
    function sig2(L, x) result(y)
        real(8) :: L, x, y
        if (x >= 0d0 .and. x<=0.2d0*L) then
            y = 1d0
        else if (x > 0.2d0*L .and. x<=0.4d0*L) then
            y = 3d0
        else if (x > 0.4d0*L .and. x<=0.6d0*L) then
            y = 1d0
        else if (x > 0.6d0*L .and. x<=0.8d0*L) then
            y = 2d0
        else if (x > 0.8d0*L .and. x<=1d0*L) then
            y = 1d0
        end if
    end function
    
    subroutine writefile(ifo)
        integer :: ifo, i, j
        
        write(ifo, '(5a)') 't', ';', 'x', ';', 'u'
        do i = 1, N
            do j = 1, M
                write(ifo, *) mesh(i,j, 1), ';', mesh(i,j, 2), ';', mesh(i,j, 3)
            end do
        end do
    end subroutine
    
end program Convection

