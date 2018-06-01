SUBROUTINE polint(xa,ya,x,y,dy)
	USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
	REAL(DP), INTENT(IN) :: x
	REAL(DP), INTENT(OUT) :: y,dy
	! Given arrays xa and ya of length N, and given a value x, this routine returns a value y,
	! and an error estimate dy. If P(x) is the polynomial of degree N - 1 such that P(xa_i) =
	! ya_i, i = 1,...,N, then the returned value y = P(x).
	INTEGER(I4B) :: m,n,ns
	REAL(DP), DIMENSION(size(xa)) :: c,d,den,ho
	n=assert_eq(size(xa),size(ya),'polint')
	c=ya ! Initialize the tableau of c's and d's.
	d=ya
	ho=xa-x
	ns=iminloc(abs(x-xa))   ! Find index ns of closest table entry.
	y=ya(ns)                ! This is the initial approximation to y.
	ns=ns-1
	do m=1,n-1                          ! For each column of the tableau,
		den(1:n-m)=ho(1:n-m)-ho(1+m:n)  ! we loop over the current c's and d's and update them.
		if (any(den(1:n-m) == 0.0)) &
		    call nrerror('polint: calculation failure')
		    ! This error can occur only if two input xa's are (to within roundoff) identical.
		den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
		d(1:n-m)=ho(1+m:n)*den(1:n-m) ! Here the c's and d's are updated.
		c(1:n-m)=ho(1:n-m)*den(1:n-m)
		if (2*ns < n-m) then    ! After each column in the tableau is completed, we decide
		    dy=c(ns+1) 		    ! which correction, c or d, we want to add to our accumulating
		else 			        ! value of y, i.e., which path to take through
			dy=d(ns) 	        ! the tableau-forking up or down. We do this in such a
			ns=ns-1 	        ! way as to take the most "straight line" route through the
		end if 			        ! tableau to its apex, updating ns accordingly to keep track
		y=y+dy 			        ! of where we are. This route keeps the partial approximations
	end do 				        ! centered (insofar as possible) on the target x. The
						        ! last dy added is thus the error indication.
END SUBROUTINE polint