	real*8 function fmin(x)
	! See chapter 9 of Numerical Recipes in Fortran by Press et al.
	integer n, np
	real*8 x(*), fvec
	parameter(np=40)
	common /newtv/ fvec(np), n
	save /newtv/
	integer i
	real*8 sum

	call funcv(n,x,fvec)
	sum=0.d0
	do i=1,n
	 sum=sum+fvec(i)**2
	end do
	fmin=0.5*sum
	return
	end 
