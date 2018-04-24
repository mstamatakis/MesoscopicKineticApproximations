SUBROUTINE newt(x,check)
    
    use global_constants
	
	USE nrtype; USE nrutil, ONLY :nrerror,vabs
	USE nr, ONLY :fdjac,lnsrch,lubksb,ludcmp
	USE fminln ! Communicates with fmin.
	
	IMPLICIT NONE
	
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(DP), PARAMETER :: TOLF=1.0e-12_dp,TOLMIN=1.0e-12_dp,TOLX=epsilon(x), &
	STPMX=100.0

    INTERFACE
        SUBROUTINE jacobian(x,jaca)
            USE nrtype
            IMPLICIT NONE
            REAL(DP), DIMENSION(:), INTENT(IN) :: x
            REAL(DP), DIMENSION(size(x),size(x)) :: jaca
		END SUBROUTINE jacobian
	END INTERFACE
	
	! Given an initial guess x for a root in N dimensions, find the root by a globally convergent
	! Newton's method. The length N vector of functions to be zeroed, called fvec in the routine
	! below, is returned by a user-supplied routine that must be called funcv and have the
	! declaration FUNCTION funcv(x). The output quantity check is false on a normal return
	! and true if the routine has converged to a local minimum of the function fmin defined
	! below. In this case try restarting from a different initial guess.
	! Parameters: MAXITS is the maximum number of iterations; TOLF sets the convergence
	! criterion on function values; TOLMIN sets the criterion for deciding whether spurious convergence
	! to a minimum of fmin has occurred; TOLX is the convergence criterion on dx;
	! STPMX is the scaled maximum step length allowed in line searches.
	
	INTEGER(I4B) :: its
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(DP) :: d,f,fold,stpmax
	REAL(DP), DIMENSION(size(x)) :: g,p,xold
	REAL(DP), DIMENSION(size(x)), TARGET :: fvec
	REAL(DP), DIMENSION(size(x),size(x)) :: fjac
	
	fmin_fvecp=>fvec ! fvec is also computed by this call.
	f=fmin(x)
	
	if (maxval(abs(fvec(:))) < 0.01_dp*TOLF) then ! Test for initial guess being a root.
		check=.false.                             ! Use more stringent test than
		RETURN                                    ! simply TOLF.
	end if
	
	stpmax=STPMX*max(vabs(x(:)),real(size(x),sp)) ! Calculate stpmax for line searches.
    
    write(*,'(8x,"Iteration   Func-count       |f(x)|")')
    
	do its=1,MAXITS ! Start of iteration loop.
        
        fcncounts = 0 ! reset function counts
        
		! If analytic Jacobian is available, you can replace the routine fdjac below with your own
		! routine.
        !call fdjac(x,fvec,fjac)
		call jacobian(x,fjac)
		g(:)=matmul(fvec(:),fjac(:,:)) ! Compute Del-f for the line search.
		xold(:)=x(:)                   ! Store x,
		fold=f                         ! and f.
		p(:)=-fvec(:)                  ! Right-hand side for linear equations.
		call ludcmp(fjac,indx,d)       ! Solve linear equations by LU decomposition.
		call lubksb(fjac,indx,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin)
		! lnsrch returns new x and f. It also calculates fvec at the new x when it calls fmin.

        write(*,'(8x,I9,1x,I12,2x,E11.5E2)') its, fcncounts, vabs(fvec)
        
        if (maxval(abs(fvec(:))) < TOLF) then ! Test for convergence on function values.
			check=.false. 
			RETURN
		end if
		
		if (check) then ! Check for gradient of f zero, i.e., spurious
			check=(maxval(abs(g(:))*max(abs(x(:)),1.0_dp) / & ! convergence.
			max(f,0.5_dp*size(x))) < TOLMIN)
			RETURN ! Test for convergence on dx.
		end if
		if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_dp)) < TOLX) &
		RETURN
	end do
	
	call nrerror('MAXITS exceeded in newt')

END SUBROUTINE newt