module methods
    implicit none
    contains
    
    function f(x) result(y) !Функция исходная(с погрешностями) 
        real(8), intent(in) ::  x
        real(8) :: y
        y = x**2d0 + 1d0 - (4d0*datan(1d0)/2d0 - x - (x**3d0)/6d0)
        !y = x**2d0 + 1d0 - dacos(x)
    end function
    
    function e(x,n) result(y) !Базисная функция 
        integer, intent(in) :: n
        real(8), intent(in) ::  x
        real(8) :: y
        y = x**real(n-1,8)
        
    end function
    
    function pf(res,x,n) result(y) !Вычисление приближения в точке 
        integer, intent(in) :: n
        real(8), intent(in) ::  x
        integer :: i
        real(8), dimension(:) :: res(n)
        real(8) :: y
        y=0d0
        do i = 1,n
            y = y + res(i)*e(x,i)
        end do
    end function
    
    subroutine gengrid(grid,n,a,b) !Генерация равномерной сетки 
        integer, intent(in) :: n
        real(8),intent(in) :: a, b
        real(8), dimension(:,:) :: grid(n,4)
        real(8) :: AA
        integer :: i
        !Генерация сетки
        grid = 0d0
        do i = 1,n
            grid(i,1) = (b-a)*real(i,8)/real(n-1,8) + a - (b-a)/real(n-1,8)
            grid(i,2) = f(grid(i,1))
        end do
        grid(:,3) = grid(:,2)
        
        !Выбросы(промахи)
        AA = dabs(maxval(grid(:,2)) - minval(grid(:,2))) !Амплитуда
        grid(2,3) = grid(2,3)*2d0*AA
        grid(10,3) = grid(10,3)*2d0*AA
        !grid(6,3) = dabs(grid(6,3)*2d0*AA)
    end subroutine
    
    subroutine MNK(grid,w,n,m,res) !метод наименьших квадратов с конкректным весами. 
        integer, intent(in) :: n, m !n - число узлов m-1 - степень строящегося полинома
        integer :: i,j,k,info
        real(8), dimension(:,:) :: grid(n,4)
        real(8), dimension(:) :: w(n), res(m) !ВЕСА и решение
        real(8), allocatable :: A(:,:), b(:), y1(:), L(:,:)
        allocate(A(m,m),b(m),y1(m), L(m,m))
        !Составление матрицы СЛАУ с весами
        b = 0d0
        A = 0d0
        y1 = 0d0
        do k = 1, m
            do i = 1,n
                b(k) = b(k) + w(i)*grid(i,3)*e(grid(i,1),k)
            end do
            
            do j = 1,m
                do i = 1,n
                    A(k,j) = A(k,j) + w(i)*e(grid(i,1),j)*e(grid(i,1),k)
                end do
            end do
        end do
        L = A
        write(*,*) "Hi"
        write(*,*) b
        write(*,"(4f12.8)") (A(i,:), i = 1,m)
        call DPOTRF("L", m, A, m,info)
        
        do i = 1,m-1
            A(i,i+1:) = 0d0
        end do
        !Решаем Ly = b
        do i = 1,m
            do j = 1, i-1
                b(i) = b(i) - A(i,j)*y1(j)
            end do
            y1(i) = b(i)/A(i,i)
        end do
        !Решаем (L**T)x = y
        A = transpose(A)
        res = 0d0
        do i = m,1,-1
            do j = m, i+1,-1
                y1(i) = y1(i) - A(i,j)*res(j)
            end do
            res(i) = y1(i)/A(i,i)
        end do
        write(*,*) "res"
        write(*,*) res
        deallocate(A, b, y1, L)
        
    end subroutine
    
    subroutine MNK2(grid,n,m,res) !Взвешенный МНК
        integer, intent(in) :: n, m !n - число узлов m - степень строящегося полинома
        integer :: i,j,k
        real(8), dimension(:,:) :: grid(n,4)
        real(8), dimension(:) :: w(n), res(m)
        real(8), allocatable :: epsi(:), wprev(:)
        real(8) :: eps
        allocate(epsi(n),wprev(n))
        eps = 1e-10
        epsi = 1d0
        res = 0d0
        wprev = 0d0
        w = 1d0
        w = w/sum(w)
        k = 0
        do while(sum(abs(wprev-w))>=eps)
            wprev = w
            write(*,*) "w"
            write(*,*) w
            call MNK(grid,w,n,m,res)
            do i = 1,n
                epsi(i) = (grid(i,3) - pf(res,grid(i,1),m))
                w(i) = 1d0/(epsi(i)**2d0 + eps)
            end do
            w = w/sum(w)
            k = k + 1
        end do
        do i = 1,n
            grid(i,4) = pf(res,grid(i,1),m)
        end do
        deallocate(epsi,wprev)
        write(*,*) k
    end subroutine
end module