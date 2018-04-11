	real*8 function fmin(x,enerv)
	! See chapter 9 of Numerical Recipes in Fortran by Press et al.
        use approx_inst
	integer n, np
	real*8 x(*), fvec
	parameter(np=40)
	common /newtv/ fvec(np), n
	save /newtv/
	integer i
        real*8, dimension(2**obj_approx%nsites) :: enerv
        
	real*8 sum

	call funcv(fvec,enerv)
	sum=0.d0
	do i=1,n
	 sum=sum+fvec(i)**2
	end do
	fmin=0.5*sum
	return
	end 
