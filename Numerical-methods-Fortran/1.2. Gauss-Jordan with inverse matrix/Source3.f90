module old
    implicit none
contains
    subroutine Jordan(Ai, m) 
        integer :: i, j, m, k !k - координата вектора, который был выбран ведущим
        double precision, dimension(:,:) :: Ai
        double precision, allocatable :: matrix(:,:)
        double precision :: a
        character (len=10) format_string1
        write(format_string1,"(a,i2,a)") "(",2*m,"f9.5)"
        allocate(matrix(m,2*m))
        matrix = 0d0
        forall(j=1:m) matrix(j,:m) = Ai(j,:)
        forall(j=1:m) matrix(j,m+j) = 1d0 ! Cформировали расширенную матрицу
        
        do i = 1, m
            a = dabs(matrix(i,i))
            k = i
            do j = 1,m
                if (dabs(matrix(j,i)) > a) then
                    a = dabs(matrix(j,i))
                    k = j
                end if
            end do
        
            matrix(i,:) = matrix(i,:)/matrix(i,i)
            do j = i+1, m
                matrix(j,:) = matrix(j,:) - matrix(i,:)*matrix(j,i)
                matrix(j,:) = matrix(j,:)/matrix(j,j)
            end do
            !do j = i+1, m
             !   matrix(m - j+1,:) = matrix(m - j+1,:) - matrix(m-i+1,:)*matrix(m - j+1,m-i+1)
            !end do
            !write(*,format_string1) (matrix(j,:), j=1,m)
            !write(*,*) "Delitel"
            do j = 1, i-1
                matrix(j,:) = matrix(j,:) - matrix(i,:)*matrix(j,i)/matrix(i,i)
            end do
            !write(*,format_string1) (matrix(j,:), j=1,m)
            !write(*,*) "Delitel"
        end do
        forall(j=1:m) Ai(j,:) = matrix(j,m+1:)
        deallocate (matrix)
    end subroutine
end module