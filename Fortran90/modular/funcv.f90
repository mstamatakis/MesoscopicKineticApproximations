	subroutine funcv(n,x,fvec)
	implicit none
	integer n, v1(1), v2(2)
	real*8 x(n), fvec(n), corr	
	
	v1=1
	fvec(1)=log(corr(v1,1,n,x))
	v1=2
	fvec(1)=fvec(1)-log(corr(v1,1,n,x))
	v1=1
	fvec(2)=log(corr(v1,1,n,x))
	v1=8
	fvec(2)=fvec(2)-log(corr(v1,1,n,x))
	v2(1)=1
	v2(2)=2
	fvec(3)=log(corr(v2,2,n,x))
	v2(1)=2
	v2(2)=3
	fvec(3)=fvec(3)-log(corr(v2,2,n,x))
	v2(1)=1
	v2(2)=2
	fvec(4)=log(corr(v2,2,n,x))
	v2(1)=2
	v2(2)=8
	fvec(4)=fvec(4)-log(corr(v2,2,n,x))
	end subroutine funcv
