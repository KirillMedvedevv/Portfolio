    SUBROUTINE MC(f, eps, a, b, buffer)
        DOUBLE PRECISION, EXTERNAL :: f
        DOUBLE PRECISION :: eps, x1, x2, x3, a, b, i
        DOUBLE PRECISION, DIMENSION(2) :: buffer
        buffer(1) = 0d0
        buffer(2) = 0d0
        x1 = a
        x2 = b
        i = 0d0
        write(*,*) x1
        write(*,*) x2
        DO
            i = i + 1d0
            x3 = x2 - f(x2)*((x2 - x1)/(f(x2) - f(x1)))
            IF(dABS((x3-x2)*(x2-x1))<eps) THEN
                EXIT
            END IF
            WRITE(*,*) x3
            x1 = x2
            x2 = x3
            
        END DO
        buffer(1) = i
        buffer(2) = x3
    END SUBROUTINE MC
    
    SUBROUTINE MPD(f, eps, a1, b1, buffer)
        DOUBLE PRECISION, EXTERNAL :: f
        DOUBLE PRECISION :: eps, a, b, c, a1, b1, i
        DOUBLE PRECISION, DIMENSION(2) :: buffer
        buffer(1) = 0d0
        buffer(2) = 0d0
        a = a1
        b = b1
        i = 0d0
        DO WHILE(dABS(b - a)/2d0 >= eps)
            i = i+1d0
            c = (a+b)/2d0
            IF(f(a)*f(c) < 0d0) THEN
                b = c
            ELSE
                a = c
            END IF
        END DO
        
        buffer(1) = i
        buffer(2) = (a+b)/2d0
    END SUBROUTINE MPD
    
    PROGRAM lab1
    IMPLICIT NONE
        DOUBLE PRECISION :: eps, ideal1, ideal2, j, c1,c2, macheps
        DOUBLE PRECISION, DIMENSION(2,2) :: mat
        DOUBLE PRECISION, DIMENSION(2) :: buffer
        DOUBLE PRECISION, DIMENSION(3, 200) :: tableMC1, tableMC2, tableMPD1, tableMPD2
        integer :: i
        ideal1 = -12.491851485282d0 !Точный корен полинома
        ideal2 = 1.10900581714465d0 !Точный корен трансцендентного уравнения
        c1 = 1080d0/(2d0*2547d0)
        c2 = 12.5114158d0
        macheps = 1.110223024625157E-018 !машинный эпсилон
        
    
        mat(1,1) = -3.8394d0
        mat(1,2) = -13.0d0
        mat(2,1) = 1.0d0
        mat(2,2) = 2.0d0
        
        eps = 1e-14
        WRITE(*,*) "Secant Method:"
        CALL MC(eq1, eps/c1, mat(1,1), mat(1,2), buffer)
        WRITE(*,*) buffer(1)
        WRITE(*,*) buffer(2)
        CALL MC(eq2, eps/c2, mat(2,1), mat(2,2), buffer)
        WRITE(*,*) buffer(2)
        
        
        j=0d0
        eps = 1d0
        DO i=1,200
            j = j+1d0 !номер итерации
            eps = 1/(2**(j/5)) !eps на i итерации
            
            CALL MC(eq1, eps/c1, mat(1,1), mat(1,2), buffer)
            tableMC1(1,i) = eps
            tableMC1(2,i) = buffer(1)
            tableMC1(3,i) = dABS(ideal1 - buffer(2))
            
            CALL MC(eq2, eps/c2, mat(2,1), mat(2,2), buffer)
            tableMC2(1,i) = eps
            tableMC2(2,i) = buffer(1)
            tableMC2(3,i) = dABS(ideal2 - buffer(2))
            
        END DO
        
        
        OPEN(1, file = "TableMC1.csv")
        OPEN(2, file = "TableMC2.csv")
        
        WRITE(1,*) "EPS",";", "N", ";", "Delta"
        WRITE(2,*) "EPS",";", "N", ";", "Delta"
        
        
        DO i = 1,200
            WRITE(1,*) tableMC1(1,i),";", tableMC1(2,i), ";", tableMC1(3,i)
            WRITE(2,*) tableMC2(1,i),";", tableMC2(2,i), ";", tableMC2(3,i)
        END DO
        
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        CLOSE(4)
        
        
    READ(*,*)   
    CONTAINS 
     FUNCTION eq1(x) RESULT(y)
        DOUBLE PRECISION :: x, y
        y = (1d0)*x**(4d0) + (12d0)*x**3d0 - (6d0)*x**2d0 + (1d0)*x - 10d0
    END FUNCTION
    
    FUNCTION eq2(x) RESULT(y)
        DOUBLE PRECISION :: x, y
        y = (7d0)**x - (6d0)*x - 2d0 
    END FUNCTION
    END PROGRAM lab1

