SUBROUTINE fdjac(x,fvec,df)
	USE nrtype; USE nrutil, ONLY :assert_eq
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: fvec
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: df
	INTERFACE
		FUNCTION funcv(x)
			USE nrtype
			IMPLICIT NONE
			REAL(DP), DIMENSION(:), INTENT(IN) :: x
			REAL(DP), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	REAL(DP), PARAMETER :: EPS=1.0e-6_dp
	! Computes forward-difference approximation to Jacobian. On input, x is the point at which
	! the Jacobian is to be evaluated, and fvec is the vector of function values at the point,
	! both arrays of length N. df is the N × N output Jacobian. FUNCTION funcv(x) is a
	! fixed-name, user-supplied routine that returns the vector of functions at x.
	! Parameter: EPS is the approximate square root of the machine precision.
	INTEGER(I4B) :: i,j,n
	REAL(DP), DIMENSION(size(x)) :: xsav,xph,h
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
	xsav=x
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h ! Trick to reduce finite precision error.
	h=xph-xsav
	do j=1,n
		x(j)=xph(j)
		df(:,j)=(funcv(x)-fvec(:))/h(j) ! Forward difference formula.
		x(j)=xsav(j)
    end do
    
    !write(*,*) ''
    !do i=1,n
    !    write(*,'(6ES30.16E3)') (df(i,j),j=1,n)
    !enddo
    !write(*,*) ''
    !continue
    
    return
    
END SUBROUTINE fdjac	

