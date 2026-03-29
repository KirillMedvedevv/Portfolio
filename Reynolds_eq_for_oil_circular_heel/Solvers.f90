module solvers
    implicit none
    real(8), parameter :: pi = 4d0*datan(1d0)
    contains
    
    subroutine genmesh(p, Nr, Nr1, Nphi, Rs)
        integer, intent(in) :: Nr, Nr1, Nphi
        real(8), intent(in) :: Rs
        integer :: i, j
        real(8), dimension(Nr, Nphi, 3) :: p
        real(8) :: dr, dphi
        
        dr = Rs/real(Nr1 - 1,8)
        dphi = pi/(3d0*(Nphi-1d0))
        
        do i = 1, Nr
            do j = 1, Nphi
                p(i,j,1) = dr*real(i-1,8)
                p(i,j,2) = dphi*real(j-1,8)
                p(i,j,3) = 0d0
            end do
        end do
    end subroutine
    
    subroutine get_solution(p, Nr, Nr1, Nphi, Rs, m, kappa, h, order, eps, maxiter, name, omega, console_log)
        integer, intent(in) :: Nr, Nr1, Nphi, order, maxiter, console_log
        real(8), dimension(Nr, Nphi, 3) :: p
        integer :: i, j, iter, logic
        character(len=*) :: name
        real(8) :: Rs, dr, dphi, m, kappa, h, nev, nev_local, nev_1, force, eps, omega
        
        
        dr = Rs/(Nr1 - 1d0)
        dphi = pi/(3d0*(Nphi-1d0)) 
        iter = 0
        
        logic = 1
        open(2, file=(Trim(name) // '_log.csv'))
        open(3, file='force.csv')
        write(2, '(21a)') 'Iteration', ',', 'residual', ',', 'force'
        write(3, '(15a)') 'Iteration', ',', 'force', ','
        nev = 0d0
        nev_local = 0d0
        do while(1)
            !Обработка границ
            if (order == 1) then
                call first_order_bc(p, Nr, Nr1, Nphi, dr, dphi, m, kappa, h, omega)
            else if (order == 2) then
                call second_order_bc(p, Nr, Nr1, Nphi, dr, dphi, m, kappa, h, omega)
            else
                write(*,*) 'WRONG ORDER'
                exit
            end if
            
            !Внутренний домен 
            
            nev = 0d0
            nev_local = 0d0
            do j = 2, Nphi-1
                do i=2, Nr1-1
                    call new_p(p, Nr, Nphi, dr, dphi, i, j, nev_local, omega)
                    if (nev_local > nev) nev = nev_local
                end do
                
                do i=Nr1+1, Nr-1
                    call new_p(p, Nr, Nphi, dr, dphi, i, j, nev_local, omega)
                    if (nev_local > nev) nev = nev_local
                end do
            end do
            
            if (logic == 1) then
                nev_1 = nev
                logic = 0
            end if
            
            nev = nev/nev_1
            force = 0d0
            force = calculate_force(p, Nr, Nr1, Nphi, dr, dphi)
            if (console_log == 1) write(*,*) iter, nev, force
            
            write(2, '(I8, a, E15.7, a, E15.7)') iter, ',', nev, ',', force
    
            iter = iter + 1
            
            if (iter == maxiter .or. nev < eps) exit
        end do
        close(2)
        close(3)
    end subroutine
    
    subroutine first_order_bc(p, Nr, Nr1, Nphi, dr, dphi, m, kappa, h, omega)
        integer, intent(in) :: Nr, Nr1, Nphi
        real(8), intent(in) :: dr, dphi, m, kappa, h, omega
        real(8) :: const, r_c, p_old, p_new
        real(8), dimension(Nr, Nphi, 3) :: p
        integer :: i, j
        !центральная точка
        
        p_old = p(1, 1, 3)
        p_new = m*dr*theta(p(1,1,3))/(3d0*kappa*h**3) + p(2,1,3) 
        p(1, :, 3) = omega*p_new + (1d0 - omega)*p_old
        
        !Радиальная канавка
        do i = 2, Nr1-1
            r_c = p(i, 1, 1)
            const = 2d0*(dr**2)/(kappa*dphi*r_c)
            p_old = p(i,1,3)
            p_new = (const*p(i,2,3) + p(i+1,1,3) + p(i-1, 1,3))/(2d0 + const)
            p(i,1,3) = omega*p_new + (1d0 - omega)*p_old
        end do
        
        r_c = p(Nr1, 1, 1)
        !Точка поворота
        const = 2d0*dr/(r_c*dphi)
        p_old = p(Nr1, 1, 3)
        p_new = (const*p(Nr1, 2, 3) + p(Nr1-1,1,3))/(1d0+const)
        p(Nr1, 1, 3) = omega*p_new + (1d0 - omega)*p_old
        
        !Кольцевая канавка
        const = -kappa*dr/((r_c*dphi)**2)
        do j = 2, Nphi-1
            p_old = p(Nr1,j,3)
            p_new = (const*(p(Nr1,j+1, 3) + p(Nr1,j-1,3)) - (p(Nr1+1,j,3) - p(Nr1-1,j,3)))/(2d0*const)
            p(Nr1,j,3) = omega*p_new + (1d0 - omega)*p_old
        end do
        
        !Симметрия 1
        p(Nr1+1:Nr-1, 1, 3) = (p(Nr1+1:Nr-1, 2, 3))*omega + (1d0-omega)*p(Nr1+1:Nr-1, 1, 3)
        
        !Симметрия 2
        p(2:Nr-1, Nphi, 3) = (p(2:Nr-1, Nphi-1, 3))*omega + (1d0-omega)*p(2:Nr-1, Nphi, 3)
        
        !Выход
        p(Nr, :, 3) = 0d0
        
    end subroutine
    
    subroutine second_order_bc(p, Nr, Nr1, Nphi, dr, dphi, m, kappa, h, omega)
        integer, intent(in) :: Nr, Nr1, Nphi
        real(8), intent(in) :: dr, dphi, m, kappa, h, omega
        real(8), dimension(Nr, Nphi, 3) :: p
        
        real(8) :: r_c, const, p_old, p_new
        
        integer :: i, j
        
        !центральная точка
        p_old = p(1, 1, 3)
        p_new = ((2d0*dr*m/(3d0*kappa*(h**3)))*theta(p(1, 1, 3)) + 4d0*p(2,1,3) - p(3,1,3))/3d0 
        p(1, :, 3) = p_new*omega + (1d0-omega)*p_old
        
        !Радиальная канавка
        do i = 2, Nr1-1
            r_c = p(i, 1, 1)
            const = (dr**2)/(kappa*dphi*r_c)
            p_old = p(i,1,3)
            p_new = (const*(4d0*p(i,2,3) - p(i, 3, 3)) + p(i+1,1,3) + p(i-1, 1,3))/(2d0 + 3d0*const)
            p(i,1,3) = omega*p_new + (1d0-omega)*p_old
        end do
        
        r_c = p(Nr1, 1, 1)
        
        !Точка поворота
        const = 2d0*dr/(r_c*dphi)
        p_old = p(Nr1, 1, 3)
        p_new = (const*(4d0*p(Nr1,2,3) - p(Nr1,3,3)) + 4d0*p(Nr1-1,1,3) - p(Nr1-2, 1, 3))/(3d0*(1d0+const))
        p(Nr1, 1, 3) = omega*p_new + (1d0-omega)*p_old
        
        !Кольцевая канавка HARD с разрывами
        const = -2d0*kappa*dr/((r_c*dphi)**2)
        do j = 2, Nphi-1
            p_old = p(Nr1,j,3)
            p_new = (const*(p(Nr1,j+1, 3) + p(Nr1,j-1,3)) - 4d0*(p(Nr1+1,j,3) - p(Nr1-1,j,3)) + (p(Nr1+2,j,3) - p(Nr1-2,j,3)))/(2d0*const)
            p(Nr1,j,3) = omega*p_new + (1d0-omega)*p_old
        end do
        
        !Симметрия 1
        p(Nr1+1:Nr-1, 1, 3) = ((4d0*p(Nr1+1:Nr-1, 2, 3) - p(Nr1+1:Nr-1,3,3))/3d0)*omega + (1d0-omega)*p(Nr1+1:Nr-1, 1, 3)
        
        !Симметрия 2
        p(2:Nr-1, Nphi, 3) = ((4d0*p(2:Nr-1, Nphi-1, 3) - p(2:Nr-1,Nphi-2,3))/3d0)*omega + (1d0-omega)*p(2:Nr-1, Nphi, 3)
        
        !Выход
        p(Nr,:, 3) = 0d0
    end subroutine
    
    function theta(p) result (p1)
        real(8) p, p1
        p1 = 0d0
        if (p <= 1d0) p1 = dsqrt(1d0 - p)
    end function
    
    subroutine new_p(p, Nr, Nphi, dr, dphi, i, j, nev, omega) !Тут все правильно
        integer, intent(in) :: Nr, Nphi, i, j
        real(8), intent(in) :: dr, dphi, omega
        real(8), dimension(Nr, Nphi, 3) :: p
        
        real(8) :: r_p, r_c, r_n, p_p_r, p_n_r, p_p_p, p_n_p, nev, p_old, p_new        
        r_p = p(i-1,j,1) ! prev
        r_c = p(i,j,1) ! current
        r_n = p(i+1,j,1) ! next
                    
        p_p_r = p(i-1,j,3) ! prev по r
        p_n_r = p(i+1,j,3) ! next по r
        p_p_p = p(i,j-1,3)
        p_n_p = p(i,j+1,3)
        
        nev = p(i,j,3)
        
        p_old = p(i,j,3)
        p_new = (p_p_r*(r_c + r_p)/(2d0*(dr**2)) + p_n_r*(r_c + r_n)/(2d0*(dr**2)) + (1/(r_c*(dphi**2)))*(p_n_p + p_p_p))/((r_p + 2d0*r_c + r_n)/(2d0*(dr**2)) + 2d0/(r_c*(dphi**2)))
        p(i,j,3) = omega*p_new + (1d0 - omega)*p_old
        nev = dabs(p(i,j,3) - nev)
        !write(*,*) p(i,j,3)
    end subroutine
    
    
    function calculate_force(p, Nr, Nr1, Nphi, dr, dphi) result(force_sum)
        integer, intent(in) :: Nr, Nphi, Nr1
        real(8), intent(in) :: dr, dphi
        real(8), dimension(Nr, Nphi, 3), intent(in) :: p 
        real(8) :: force_sum
        integer :: i, j
        real(8) :: r_val, p_val, w_phi, sum_phi(Nr)
        real(8) :: sum_total, sector_sum

        sum_total = 0d0
        sector_sum = 0d0

        !Интегрируем по углу phi (вдоль каждой кольцевой линии i)
        do i = 1, Nr
            sector_sum = 0.0d0
            do j = 1, Nphi
                if (j == 1 .or. j == Nphi) then
                    w_phi = 1.0d0
                else if (mod(j, 2) == 0) then
                    w_phi = 4.0d0
                else
                    w_phi = 2.0d0
                end if
                sector_sum = sector_sum + p(i, j, 3) * w_phi
            end do
            sum_phi(i) = sector_sum * (dphi / 3.0d0)
        end do

        !Интегрируем по радиусу r с учетом излома в Nr1
        sum_total = sum_total + simpson(sum_phi, p(:,1,1), 1, Nr1)
        sum_total = sum_total + simpson(sum_phi, p(:,1,1), Nr1, Nr)
    

        force_sum = 6d0 * sum_total
    end function

    function simpson(f, r, start_idx, end_idx) result(res)
        integer, intent(in) :: start_idx, end_idx
        real(8), dimension(:), intent(in) :: f, r
        real(8) :: res, dr_local
        integer :: n_points, i
    
        res = 0.0d0
        n_points = end_idx - start_idx + 1
        dr_local = r(start_idx + 1) - r(start_idx)
    
        if (n_points < 2) return

        if (mod(n_points, 2) /= 0) then
            ! Стандартный Симпсон для нечетного числа узлов
            do i = start_idx, end_idx
                if (i == start_idx .or. i == end_idx) then
                    res = res + f(i) * r(i) * 1.0d0
                else if (mod(i - start_idx + 1, 2) == 0) then
                    res = res + f(i) * r(i) * 4.0d0
                else
                    res = res + f(i) * r(i) * 2.0d0
                end if
            end do
            res = res * (dr_local / 3.0d0)
        else
            ! Если число узлов четное, используем стандартного Симпсона на (n-3) узлах и правило 3/8
            do i = start_idx, end_idx - 3
                if (i == start_idx .or. i == end_idx - 3) then
                    res = res + f(i) * r(i) * 1.0d0
                else if (mod(i - start_idx + 1, 2) == 0) then
                    res = res + f(i) * r(i) * 4.0d0
                else
                    res = res + f(i) * r(i) * 2.0d0
                end if
            end do
            res = res * (dr_local / 3.0d0)
            ! Добавляем последний сегмент 3/8
            res = res + (3d0*dr_local/8d0)*(f(end_idx-3)*r(end_idx-3) + 3d0*f(end_idx-2)*r(end_idx-2) + 3d0*f(end_idx-1)*r(end_idx-1) + f(end_idx)*r(end_idx))
        end if
    end function
    
    function analytic_p(r, phi, m, Rs) result(p)
        real(8) :: p, r, m, Rs, phi, A
        
        A = 0.5d0*(-m**2/(4d0*pi**2)*dlog(1d0/Rs) + dsqrt((m**2/(4d0*pi**2)*dlog(1d0/Rs))**2 + (m/pi)**2))
        
        if (r < Rs) then
            p = A*dlog(1d0/Rs)
        else
            p = A*dlog(1d0/r)
        end if
    end function
end module