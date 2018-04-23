SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
	USE nrtype; USE nrutil, ONLY :assert_eq, nrerror,vabs
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(DP), INTENT(IN) :: fold,stpmax
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x
	REAL(DP), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	INTERFACE
		FUNCTION func(x)
			USE nrtype
			IMPLICIT NONE
			REAL(DP) :: func
			REAL(DP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	REAL(DP), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
	! Given an N-dimensional point xold, the value of the function and gradient there, fold
	! and g, and a direction p, finds a new point x along the direction p from xold where the
	! function func has decreased “sufficiently.” xold, g, p, and x are all arrays of length N.
	! The new function value is returned in f. stpmax is an input quantity that limits the length
	! of the steps so that you do not try to evaluate the function in regions where it is undefined
	! or subject to overflow. p is usually the Newton direction. The output quantity check is
	! false on a normal exit. It is true when x is too close to xold. In a minimization algorithm,
	! this usually signals convergence and can be ignored. However, in a zero-finding algorithm
	! the calling program should check whether the convergence is spurious.
	! Parameters: ALF ensures sufficient decrease in function value; TOLX is the convergence
	! criterion on Dx.
	INTEGER(I4B) :: ndum
	REAL(DP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs ! Scale if attempted step is too big.
	slope=dot_product(g,p)
	if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_sp)) ! Compute lamda-min.
	alam=1.0                                               ! Always try full Newton step first.
	do                                                     ! Start of iteration loop.
		x(:)=xold(:)+alam*p(:)
		f=func(x)
		if (alam < alamin) then                            ! Convergence on Dx. For zero finding,
			x(:)=xold(:)                                   ! the calling program should
			check=.true.                                   ! verify the convergence.
		RETURN
		else if (f <= fold+ALF*alam*slope) then            ! Sufficient function decrease.
		RETURN
		else                                               ! Backtrack.
			if (alam == 1.0) then                          ! First time.
			tmplam=-slope/(2.0_sp*(f-fold-slope))
			else                                           ! Subsequent backtracks.
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/ &
				(alam-alam2)
				if (a == 0.0) then
				tmplam=-slope/(2.0_sp*b)
				else
					disc=b*b-3.0_sp*a*slope
					if (disc < 0.0) then
					tmplam=0.5_sp*alam
					else if (b <= 0.0) then
					tmplam=(-b+sqrt(disc))/(3.0_sp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam ! lamda ≤ 0.5*lamda1.
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_sp*alam) ! lamda ≥ 0.1*lamda1.
	end do ! Try again.
END SUBROUTINE lnsrch







