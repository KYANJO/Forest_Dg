program LagrangePolynomial
implicit none
    DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE    :: xsample,xroots
    DOUBLE PRECISION                :: x,p,xl,xi,xj, Pi = 3.14
    INTEGER                         :: k,i,j,N

    OPEN(unit=10,FILE='dataCheb.txt',ACTION='READWRITE')
    OPEN(unit=1,FILE='outputc1.txt',ACTION='READWRITE')
    OPEN(unit=2,FILE='outputc2.txt',ACTION='READWRITE')
    OPEN(unit=3,FILE='outputc3.txt',ACTION='READWRITE')

    READ(10,*) N
    ALLOCATE(xroots(N))
    ALLOCATE(xsample(N))
    DO i=1,N
        READ(10,*) xroots(i), xsample(i)
	
    END DO
    
    

!   CALL LagrangeBasis(N,xsample,xroots,xl,p)

    DO i=1,100
        xl=cos(((2.0*i+1.0)/(2.0*100+2.0))*Pi)
        CALL LagrangeBasis(N,xsample,xroots,xl,p)
        WRITE(1,*)xl,p
	
    END DO!i
    

    

    DEALLOCATE(xroots)
    DEALLOCATE(xsample)

    STOP

    CONTAINS
    

    SUBROUTINE LagrangeBasis(N,xsample,xroots,xl,output)
    
    IMPLICIT NONE
    INTEGER, INTENT(IN)                         :: N
    DOUBLE PRECISION, DIMENSION(N),INTENT(IN)   :: xsample,xroots
    DOUBLE PRECISION, INTENT(IN)                :: xl
    DOUBLE PRECISION, INTENT(OUT)               :: output
    DOUBLE PRECISION                            :: L, S
    INTEGER                                     :: i,j

    output=0.0D0
    
	     do i = 1, N 
		xi = xroots(i)
		L = 1.0D0
		
		S = 0.0D0
		do j = 1, N 
		     xj = xroots(j)
		     if (j .ne. i) then
			L = L*(xl-xj)/(xi-xj)
			S = S+abs(L)
		     endif
		
		enddo
		write(3,*)S
		write(2,*)i,L
	        output=output+L*xsample(i)
	      enddo
     END SUBROUTINE LagrangeBasis
end program LagrangePolynomial
