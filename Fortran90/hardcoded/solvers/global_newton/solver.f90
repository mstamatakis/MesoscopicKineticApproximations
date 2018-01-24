	subroutine solver(x,n,check)
	integer n, nn, np, maxits
	logical check
	real*8 x(n), fvec, tolf, tolmin, tolx, stpmx
	parameter(np=40,maxits=200,tolf=1.0d-4,tolmin=1.0d-6,tolx=1.0d-7,stpmx=100.d0)
	common /newtv/fvec(np),nn
	save /newtv/
	integer i, its, j, indx(np)
	real*8 d, den, f, fold, stpmax, sum, temp, test, fjac(np,np)
	real*8 g(np), p(np), xold(np), fmin
	external fmin

	nn=n
	f=fmin(x)
	test=0.d0
	do i=1,n
	 if(abs(fvec(i)).gt.test)test=abs(fvec(i))
	end do
	if(test.lt.0.01d0*tolf)then
	 check=.false.
	 return
	end if
	sum=0.d0
	do i=1,n
	 sum=sum+x(i)**2
	end do	
	stpmax=stpmx*max(sqrt(sum),float(n))
	do its=1,maxits
	 call fdjac(n,x,fvec,np,fjac)
	 do i=1,n
	  sum=0.d0
	  do j=1,n
	   sum=sum+fjac(j,i)*fvec(j)
	  end do
 	  g(i)=sum
	 end do
	 do i=1,n
	  xold(i)=x(i)
	 end do
	 fold=f
	 do i=1,n
	  p(i)=-fvec(i)
	 end do
	 call ludcmp(fjac,n,np,indx,d)
	 call lubksb(fjac,n,np,indx,p)
	 call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
	 test=0.d0
	 do i=1,n
	  if(abs(fvec(i)).gt.test)test=abs(fvec(i))
	 end do
	 if(test.lt.tolf)then
	  test=0.d0
	  den=max(f,0.5d0*n)
	  do i=1,n
	   temp=abs(g(i))*max(abs(x(i)),1.d0)/den
	   if(temp.gt.test)test=temp 
	  end do 
	  if(test.lt.tolmin)then
	   check=.true.
	  else
	   check=.false.
	  end if
	  return
	 end if
	 test=0.d0
	 do i=1,n
	  temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
	  if(temp.gt.test)test=temp
	 end do
	 if(test.lt.tolx)return
	end do

	contains
	 subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
	 ! See chapter 9 of Numerical Recipes in Fortran by Press et al.
	  integer n
	  logical check
	  real*8 f, fold, stpmax, g(n), p(n), x(n), xold(n), func, alf, tolx
	  parameter (alf=1.d-4,tolx=1.d-7)
	  external func
	  integer i
	  real*8 a, alam, alam2, alamin, b, disc, f2, rhs1, rhs2
	  real*8 slope, sum, temp, test, tmplam

	  check=.false.
	  sum=0.d0
	  do i=1,n
	   sum=sum+p(i)*p(i)
	  end do
	  sum=sqrt(sum)
	  if(sum.gt.stpmax)then
	   do i=1,n
	    p(i)=p(i)*stpmax/sum
	   end do
	  end if
	  slope=0.d0
	  do i=1,n
	   slope=slope+g(i)*p(i)
	  end do
	  test=0.d0
	  do i=1,n
	   temp=abs(p(i))/max(abs(xold(i)),1.d0)
	   if(temp.gt.test)test=temp
	  end do
	  alamin=tolx/test
	  alam=1.d0
1	  continue
	   do i=1,n
	    x(i)=xold(i)+alam*p(i)
	   end do
	   f=func(x)
	   if(alam.lt.alamin)then
	    do i=1,n
	     x(i)=xold(i)
	    end do
	    check=.true.
	    return
	   else if(f.le.fold+alf*alam*slope)then
	    return
	   else
	    if(alam.eq.1.d0)then
	     tmplam=-slope/(2.d0*(f-fold-slope))
	    else
	     rhs1=f-fold-alam*slope
	     rhs2=f2-fold-alam2*slope
	     a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
	     b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
	     if(a.eq.0.d0)then
	      tmplam=-slope/(2.d0*b)
	     else
	      disc=b*b-3.0*a*slope
	      if(disc.lt.0.d0)then
	       tmplam=0.5d0*alam
	      else if(b.le.0.d0)then
	       tmplam=(-b+sqrt(disc))/(3.d0*a)
	      else
	       tmplam=-slope/(b+sqrt(disc))
	      end if 
	     end if
	      if(tmplam.gt.0.5d0*alam)tmplam=0.5d0*alam
	    end if
	   end if
	   alam2=alam
	   f2=f
	   alam=max(tmplam,0.1d0*alam)
	  goto 1
	 end subroutine	lnsrch 
	
	 subroutine ludcmp(a,n,np,indx,d)
	 integer n, np, indx(n), nmax
	 real*8 d, a(np,np), tiny
	 parameter(nmax=500,tiny=1.0d-20)
	 integer i, imax, j, k
	 real*8 aamax, dum, sum, vv(nmax)

	 d=1.d0
	 do i=1,n
	  aamax=0.d0
	  do j=1,n
	   if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
	  end do
	  vv(i)=1.d0/aamax
	 end do
	 do j=1,n
	  do i=1, j-1
	   sum=a(i,j) 
	   do k=1,i-1
	    sum=sum-a(i,k)*a(k,j)
	   end do
	   a(i,j)=sum
	  end do
	  aamax=0.d0
	  do i=j,n
	   sum=a(i,j)
	   do k=1,j-1
	    sum=sum-a(i,k)*a(k,j)
	   end do
	   a(i,j)=sum
	   dum=vv(i)*abs(sum)
	   if(dum.ge.aamax)then
	    imax=i
	    aamax=dum
	   end if
	  end do
	  if(j.ne.imax)then
	   do k=1,n
	    dum=a(imax,k)
	    a(imax,k)=a(j,k)
	    a(j,k)=dum
	   end do
	   d=-d
	   vv(imax)=vv(j)
	  end if
	  indx(j)=imax
	  if(a(j,j).eq.0.d0)a(j,j)=tiny
	  if(j.ne.n)then
	   dum=1.d0/a(j,j)
	   do i=j+1,n
	    a(i,j)=a(i,j)*dum
	   end do
	  end if
	 end do
	 return
	 end subroutine ludcmp
	 
	 subroutine lubksb(a,n,np,indx,b)
	 integer n, np, indx(n)
	 real*8 a(np,np), b(n)	
	 integer i, ii, j, ll
	 real*8 sum

	 ii=0
	 do i=1,n
	  ll=indx(i)
	  sum=b(ll)
	  b(ll)=b(i)
	  if(ii.ne.0)then
	   do j=ii,i-1
	    sum=sum-a(i,j)*b(j)
	   end do
	  else if(sum.ne.0.d0)then
	    ii=i
	  end if
	  b(i)=sum
	 end do
	 do i=n,1,-1
	  sum=b(i)
	  do j=i+1,n
	   sum=sum-a(i,j)*b(j)
	  end do
	  b(i)=sum/a(i,i)
	 end do
	 return
	 end subroutine lubksb 	

	 subroutine fdjac(n,x,fvec,np,df)
	 integer n, np, nmax
	 real*8 df(np,np),fvec(n),x(n),eps
	 parameter (nmax=40,eps=1.0d-4)
	 integer i, j
	 real*8 h, temp, f(nmax)

	 do j=1,n
	  temp=x(j)
	  h=eps*abs(temp)
	  if(h.eq.0)h=eps
	  x(j)=temp+h
	  h=x(j)-temp
	  call funcv(n,x,f)
	  x(j)=temp
	  do i=1,n
	   df(i,j)=(f(i)-fvec(i))/h
	  end do
	 end do
	 end subroutine fdjac	

	 subroutine funcv(n,x,fvec)
	 implicit none
	 integer n
	 real*8 x(n), fvec(n)	

	 end subroutine funcv
	end 
