    program Lab22
    use methods
        implicit none
        integer :: i, j, k, n, m, count
        real(8) :: a, b, sum
        real(8), allocatable :: grid(:,:), mapp(:,:), w(:), res(:)
        integer :: count1, count_rate, count_max
        call SYSTEM_CLOCK(count1, count_rate, count_max)
        n = 12 !Эксперименты
        m = 4 !Степень полинома - 1
        k = 0 !Тесты
        count = 1000
        allocate(grid(n,4),mapp(count,3), res(m))
        !grid 1 - сетка 2 - значения функции 3 - значения с выбросами 4 - МНК значения
        !mapp 1 - аргумент 2 - значения функции
        a = -0.5d0
        b = 0.5d0
        call gengrid(grid,n,a,b)
        !Запись экспериментальных данных
        open(3, file = "ferr.csv") 
        write(3,*) "X", ";", "Y"
        do i = 1, n
            write(3,"(f18.14,a,f18.14)") grid(i,1),";", grid(i,3)
        end do
        close(3)
        call MNK2(grid,n,m,res)
        write(*,"(4f14.8)") (grid(i,:), i = 1,n)
        !Построение полинома и вычисление фактической ошибки
        do i = 1,count
            mapp(i,1) = a + (b-a)*real(i-1,8)/real(count,8)
            mapp(i,2) = pf(res,mapp(i,1),m)
            mapp(i,3) = (mapp(i,2) - f(mapp(i,1)))
        end do
        open(1, file = "Result.csv")
        write(1,*) "x", ";", "y", ";", "py", ";", "error"
        do i = 1, count
            write(1, "(f16.12, a, f16.12, a, f16.12, a, f16.12)") mapp(i,1),";", f(mapp(i,1)),";",mapp(i,2), ";", mapp(i,3)
        end do
        close(1)
        deallocate(grid, mapp)
        open(2, file= "Err.csv")
        write(2,*) "n", ";", "maxerr"
        do i = 12,k
            allocate(grid(i,4))
            w = 1d0
            call gengrid(grid,i,a,b)
            call MNK2(grid,i,m,res)
            sum = 0d0
            do j= 1,i-1
                sum = sum + (grid(i,4) - grid(i,3))**2d0
            end do
            sum = dsqrt(sum/real(i-1,8))
            write(2,"(f16.12,a, f16.12)") (real(i,8)), ";", dlog10(sum)
            deallocate(grid)
        end do
        close(2)
        deallocate(res)
        write(*,*) count1, count_rate, count_max
        pause
    end program Lab22

