module functions
    implicit none
    contains
    
    function getDet(A, m) result(det)
        real(8), dimension(:,:) :: A
        real(8), allocatable :: Ai(:,:)
        real(8) :: det
        integer :: i, j, m
        allocate(Ai(m,m))
        Ai = A
        do i = 1, m-1
            do j = i+1, m
                Ai(j,:) = Ai(j,:) - Ai(i,:)*Ai(j,i)/Ai(i,i)
            end do
        end do
        
        det = 1d0
        do i =1,m
            det = det*Ai(i,i)
        end do
        deallocate(Ai)
    end function
    
    function tr(A,m) result (p)
        real(8), dimension(:,:) ::  A
        integer i, m
        real(8) :: p
        p = 0d0
        do i = 1,m
            p = p + A(i,i)
        end do
    end function
    
    subroutine mpow(A,k,m, Ai) 
    real(8), dimension(:,:)  :: A, Ai
    integer :: i, k, m
    Ai = A
    do i =1,k-1
        Ai = MATMUL(Ai,A)
    end do
    end subroutine
end module