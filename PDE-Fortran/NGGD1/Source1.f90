module solvers
    use omp_lib
    implicit none
    contains
    subroutine gen_mesh(mesh, M, N, U0, w, a, dt, dx)
        integer, intent(in) :: M, N
        real(8), intent(in) :: U0, w, a, dt, dx
        integer :: i, j
        real(8) :: mesh(M, N, 3)
        
        do i = 1, M
            do j=1, N
                mesh(i, j, 1) = dx*(i-1d0)
                mesh(i,j,2) = dt*(j-1d0)
                mesh(i,j,3) = 0d0
            end do
        end do
        
        !ГУ
        do j = 1, N
            mesh(1,j,3) = func(U0, w, mesh(1,j,2))
        end do
        
        !НУ
        !do i = 1, M
            !mesh(i,1,3) = U(U0, w, a, mesh(i,1,1), mesh(i,1,2))
        !end do
    end subroutine
    
    function U(U0, w, a, x, t) result(p) !Аналитическое решение
        real(8), intent(in) :: U0, w, x, a, t
        real(8) :: p
        p = U0*dexp(-dsqrt(w/(2d0*a))*x)*dsin(w*t - dsqrt(w/(2d0*a))*x)
    end function
    
    function func(U0, w, t) result(U) ! НУ
        real(8), intent(in) :: U0, w, t
        real(8) :: U
        U = U0*dsin(w*t)
    end function
    
    subroutine explicit(mesh, M, N, U0, w, a, dt, dx, p) !Явная схема
        integer, intent(in) :: M, N, p
        real(8), intent(in) :: dt, dx, a, U0, w
        integer :: i, j, k
        real(8) :: mesh(M, N, 3)
        do j = 2, N
            !$omp parallel do num_threads(p) shared(mesh) private(i)
            do i = 2, M - 1 !координата
                mesh(i,j,3) = (a*dt/(dx**2))*(mesh(i+1,j-1,3) + mesh(i-1, j-1,3)) + mesh(i,j-1, 3)*(1d0 - 2d0*a*dt/(dx**2))
            end do
            !$omp end parallel do
            mesh(M,j,3) = (4d0*mesh(M-1, j, 3) - mesh(M-1, j, 3))/3d0
        end do
    end subroutine
    
    subroutine implicit2(mesh, M, N, a, dt, dx) !Двухслойная неявная
        integer :: i, j, M, N, info
        real(8) :: mesh(M, N, 3), dt, dx, a
        real(8) :: A1(M-2), B1(M-1), C1(M-2), D(M-1)
        real(8) :: alpha, beta, time

        alpha = a * dt / dx**2
        beta = 1d0 / dt
        do j = 1, N-1
            ! Заполнение матрицы
            A1 = -alpha
            B1 = 1d0 + 2d0 * alpha
            C1 = -alpha
            B1(M-1) = B1(M-1) - alpha  

            ! Правая часть
            do i = 1, M-1
                D(i) = mesh(i+1, j, 3)
            end do
            D(1) = D(1) + alpha * mesh(1, j+1, 3)

            ! Решение СЛАУ
            call dgtsv(M-1, 1, A1, B1, C1, D, M-1, info)
            mesh(2:M, j+1, 3) = D
    end do
    end subroutine
    
    subroutine implicit3(mesh, M, N, a, dt, dx) !Трехслойная неявная
        integer :: i, j, k, M, N, nrhs, ldb, info
        real(8) :: mesh(M, N, 3), dt, dx, a
        real(8) :: A1(M-2), B1(M-1), C1(M-2), D(M-1)
        real(8) :: alpha, beta, time
        
    
        !Один шаг делаем двухслойной неявной схемой
        A1 = -a * dt / dx**2
        B1 = 1d0 + 2d0 * a * dt / dx**2
        C1 = -a * dt / dx**2
        B1(M-1) = 1d0 + a * dt / dx**2
        
        do i = 1, M-1
            D(i) = mesh(i+1, 1, 3)
        end do
        D(1) = D(1) + A1(1)*mesh(1, 2, 3)
        
        call dgtsv(M-1, 1, A1, B1, C1, D, M-1, info)
        mesh(2:M,2, 3) = D
        
        !Трехслойная схема
        alpha = 2d0*a * dt / dx**2
        beta = 3d0 + 4d0*a*dt/dx**2
        
        time = omp_get_wtime()
        do j = 2, N-1 !Время
            !Создание матрицы прогоночных коэффициентов
            A1 = -alpha
            B1 = beta
            C1 = -alpha
            B1(M-1) = B1(M-1) + 4d0*C1(M-2)/3d0
            A1(M-2) = A1(M-2) - C1(M-2)/3d0
    
            !Cоставляем матрицу свободных коэффициентов
            do i = 1, M-1 !Координата
                D(i) = 4d0*mesh(i+1, j, 3) - mesh(i+1,j-1, 3)
            end do
            D(1) = D(1) - A1(1)*mesh(1, j+1, 3)
     
            !Вызываем решатель MKL
            call dgtsv(M-1, 1, A1, B1, C1, D, M-1, info)
            mesh(2:M,j+1, 3) = D
        end do
        !write(*,*) omp_get_wtime() - time
    end subroutine
    
    subroutine p_dgtsv(N, nrhs, A, B, C, D, ldb, info, p)
        integer, intent(in) :: N, nrhs, ldb, p
        integer :: id, start, endd, i, j, k, info, box
        real(8) :: A(N), B(N), C(N), D(N), coef, cord
        real(8) :: Ap(p), Bp(p), Cp(p), Dp(p)
        
        !Проверка
        if (mod(N, p) /= 0) then
            !info = -1
            !return
        else
            box = N/p
        end if
        
        call omp_set_num_threads(p)
        !$omp parallel shared(A, B, C, D, box) private(id, start, endd, i, j, coef, cord)
        !Назначение итераций
        id = omp_get_thread_num()
        start = id*box + 1
        endd = (id+1d0)*box
        !write(*,*) id, start, endd, box
        
        !Прямой ход
        do i = start+1, endd
            coef = A(i)/B(i-1)
            A(i) = -coef*A(i-1)
            B(i) =  B(i) - coef**C(i-1)
            D(i) = D(i) - coef*D(i-1)
        end do
        
        !Обратный ход
        do i = endd-1, start, -1
            coef = C(i)/B(i)
            C(i) = -coef*C(i+1)
            A(i) = A(i) -coef*A(i+1)
            D(i) = D(i) - coef*D(i+1)
        end do
        !$omp barrier !Барьер синхранизации
        !$omp single
        !Обмен данными
        do i = 1,p-1
            cord = i*box+1
            coef = C(cord-1)/B(cord)
            C(cord-1) = -coef*C(cord)
            B(cord-1) = B(cord-1) - coef*A(cord)
            D(cord-1) = D(cord-1) - coef*D(cord)          
            Ap(i) = A(cord-1)
            Bp(i) = B(cord-1) 
            Cp(i) = C(cord-1) 
            Dp(i) = D(cord-1) 
        end do
        Ap(p) = A(N)
        Bp(p) = B(N) 
        Cp(p) = 0d0
        Dp(p) = D(N)
        !Параллельное решение слау меньшего размера
        call dgtsv(box, 1, Ap(2:box), Bp, Cp(1:box-1), Dp, box, info)
        
        !$omp end single
        
        !Обратный ход
        do i = start, endd
            cord = (id+1)*box !Номер известного корня Номер в массиве = id+1
            if (id == 0 .or. id==p-1)  then
                
            else
                
            end if
        end do
        !$omp end parallel
    end subroutine
end module