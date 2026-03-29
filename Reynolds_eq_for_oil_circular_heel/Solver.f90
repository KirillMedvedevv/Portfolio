    program Solver
    use omp_lib
    use solvers
    implicit none
    integer :: Nr, Nphi, Nr1, maxiter, order, console_log
    real(8) :: kappa, m, hs, h, Rs, force, kappa_true, eps, time, omega, dr, dphi
    character(len=25) :: name
    real(8), allocatable :: p(:,:,:)
    
    
    call read_parameters('Input.txt')
    
    
    kappa = kappa_true/(h**3)
    Nr1 = dnint((Nr-1)*Rs) + 1
    dr = Rs/(Nr1 - 1d0)
    dphi = 4d0*datan(1d0)/(3d0*(Nphi-1d0))
    
    write(*,*) 'Nr1 = ', Nr1
    pause
    
    
    write(*,*) 'Running'
    
    time = omp_get_wtime()
    
    allocate(p(Nr, Nphi, 3))
    
    call genmesh(p, Nr, Nr1, Nphi, Rs)
    call get_solution(p, Nr, Nr1, Nphi, Rs, m, kappa, h, order, eps, maxiter, name, omega, console_log)
    call write_tecplot_file(p, Nr, Nphi, m, Rs, name)
    
    force = calculate_force(p, Nr, Nr1, Nphi, dr, dphi)
    time = omp_get_wtime() - time
    
    open(1, file=(trim(name) // '_report.txt'))
    write(1,*) 'force = ', force
    write(1,*) 'perepad = ', (p(1,1,3)-p(Nr1, 1,3))*100d0/p(1,1,3)
    write(1,*) 'time = ', time
    close(1)
    
    write(*,*) 'Main Results:'
    write(*,*) '  force = ', force
    write(*,*) '  perepad = ', (p(1,1,3)-p(Nr1, 1,3))*100d0/p(1,1,3)
    write(*,*) '  time = ', time
    
    deallocate(p)
    pause
    contains
    
    subroutine read_parameters(filename)
        character(*), intent(in) :: filename
        integer :: unit, io_stat
        character(len=256) :: line
        
        unit = 10
        
        open(unit=unit, file=filename, status='old', action='read', iostat=io_stat)
        
        if (io_stat /= 0) then
            print *, "ERROR:", filename
            stop
        end if

        read(unit, *, IOSTAT=io_stat) Rs
        read(unit, *, IOSTAT=io_stat) h
        read(unit, *, IOSTAT=io_stat) m
        read(unit, *, IOSTAT=io_stat) kappa_true
        read(unit, *, IOSTAT=io_stat) Nr
        read(unit, *, IOSTAT=io_stat) Nphi
        read(unit, *, IOSTAT=io_stat) order
        read(unit, *, IOSTAT=io_stat) omega
        read(unit, *, IOSTAT=io_stat) eps
        read(unit, *, IOSTAT=io_stat) maxiter   
        read(unit, *, IOSTAT=io_stat) name
        read(unit, *, IOSTAT=io_stat) console_log
        close(unit)
        
        ! ¬ывод считанных параметров дл€ проверки
        write(*,*) "Task parameters"
        write(*,*), "  Rs:", Rs
        write(*,*), "  h:", h
        write(*,*), "  m:", m
        write(*,*), "  kappa:", kappa_true
        
        write(*,*) "Mesh parameters"
        write(*,*), "  Nr:", Nr
        write(*,*), "  Nphi:", Nphi
        
        write(*,*) "Solver parameters"
        write(*,*), "  order:", order
        write(*,*), "  omega:", omega
        write(*,*) "  eps:", eps
        write(*,*), "  maxiter:", maxiter
        
        write(*,*) "Other parameters"
        write(*,*), "  Path to save:", name
        write(*,*), "  Console_log:", console_log
    end subroutine
        
    subroutine write_csv_file(p, Nr, Nphi, m, Rs)
        integer, intent(in) :: Nr, Nphi
        integer :: i, j
        real(8), dimension(Nr, Nphi, 3) :: p
        real(8) :: m, Rs
        
        open(1, file='Res.csv')
        write(1,'(11a)') 'r', ',', 'phi', ',', 'p',',','A_p'
        do i = 1, Nr
            do j = 1, Nphi
                write(1,'(f15.3, a, f15.3, a, f15.3, a, f15.3)') p(i,j,1), ',', p(i,j,2), ',', p(i,j,3), ',', analytic_p(p(i,j,1), p(i,j,2), m, Rs)
            end do
        end do
        close(1)
    end subroutine
    
    subroutine write_tecplot_file(p, Nr, Nphi, m, Rs, filename)
        integer, intent(in) :: Nr, Nphi
        real(8), dimension(Nr, Nphi, 3), intent(in) :: p
        character(len=*), intent(in) :: filename
        integer :: i, j
        real(8) :: x, y, m, Rs
    
        open(1, file=trim(filename)//'.plt', status='replace')
    
        write(1, '(A)') 'TITLE = "Results"'
        write(1, '(A)') 'VARIABLES = "X", "Y", "P", "A_p"'
        write(1, '(A,I5,A,I5,A)') 'ZONE I=', Nr, ', J=', Nphi, ', F=POINT'
    
        do j = 1, Nphi
            do i = 1, Nr
                x = p(i,j,1) * dcos(p(i,j,2))
                y = p(i,j,1) * dsin(p(i,j,2))
                write(1, '(4ES15.6)') x, y, p(i,j,3), analytic_p(p(i,j,1), p(i,j,2), m, Rs)
            end do
        end do
    
        close(1)
    end subroutine write_tecplot_file
    end program Solver

