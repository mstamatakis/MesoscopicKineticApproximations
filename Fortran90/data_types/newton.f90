	subroutine solver(x,n,check)
	! See chapter 9 of Numerical Recipes in Fortran by Press et al.
        use lu 
        use approx_inst
        implicit none
	integer n, nn, np, maxits, d, c
	logical check
	real*8 x(n), fvec, tolf, tolmin, tolx, stpmx
	parameter(np=40,maxits=200,tolf=1.0d-4,tolmin=1.0d-6,tolx=1.0d-7,stpmx=100.d0)
	common /newtv/fvec(np),nn
	save /newtv/
	integer i, its, j, indx(np)
	real*8 den, f, fold, stpmax, sum, tmp, test, fjac(np,np)
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

 	 call fdjac(fjac)
         write(*,*) 'ANALYTIC'
         do i=1,n
           write(*,*) (fjac(i,j), j=1,n)
         end do 
         call fdjac2(n,x,fvec,np,fjac)
         write(*,*) 'NUMERICAL'
         do i=1,n
           write(*,*) (fjac(i,j), j=1,n)
         end do

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
	 call ludcmp(fjac,n,np,indx,d,c)
	 call lubksb(fjac,n,np,indx,p)
         !call gelim(fjac,-fvec,n,np,p)
         write(*,*) its,sqrt(fvec(1)**2+fvec(2)**2+fvec(3)**2+fvec(4)**2), check
         write(18,*) its
	 call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
	 test=0.d0
	 do i=1,n
	  if(abs(fvec(i)).gt.test)test=abs(fvec(i))
	 end do
	 if(test.lt.tolf)then
          check=.false.
          return
         end if
         if(check)then
	  test=0.d0
	  den=max(f,0.5d0*n)
	  do i=1,n
	   tmp=abs(g(i))*max(abs(x(i)),1.d0)/den
	   if(temp.gt.test)test=tmp 
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
	  tmp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
	  if(tmp.gt.test)test=tmp
	 end do
	 if(test.lt.tolx)return
	end do
        contains
	subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
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
1	 continue
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
	
	subroutine fdjac(df)
        use commons
        use meso_approx
        use approx_inst
	integer i, j, k, l, m
        integer p, q, s, t
        integer, dimension(obj_approx%nsites) :: state
        real*8 corf, totl, totr
        real*8, dimension(obj_approx%eqn%neqns,obj_approx%eqn%neqns) :: df

	do i=1,obj_approx%eqn%neqns
         do j=1,obj_approx%eqn%neqns
          totl=0.d0
          do k=1,2**obj_approx%nsites
           call confs(state,k)
           p=1
           corf=0.d0
           do l=1,obj_approx%eqn%ncol(i)
            p=p*state(obj_approx%eqn%lhs(i,l))
           end do
           s=0
           do l=1,obj_approx%hamilt%corr%nrows
            q=1
            if(obj_approx%hamilt%corr%intmap(l).eq.j) then
             do m=1,obj_approx%hamilt%corr%ncol(l)
              q=q*state(obj_approx%hamilt%corr%term(l,m))
             end do
            end if
            s=s+q
           end do
           t=0
           do l=1,obj_approx%nsites
            t=t+state(l)
           end do
           totl=totl+p*s*exp((chemp*t-obj_approx%ener(obj_approx,state)-h0)/(kb*temp))
           corf=corf+p*exp((chemp*t-obj_approx%ener(obj_approx,state)-h0)/(kb*temp))
          end do
          totl=-totl/(kb*temp)

          totr=0.d0
          do k=1,2**obj_approx%nsites
           call confs(state,k)
           p=1
           corf=0.d0
           do l=1,obj_approx%eqn%ncol(i)
            p=p*state(obj_approx%eqn%rhs(i,l))
           end do
           s=0
           do l=1,obj_approx%hamilt%corr%nrows
            q=1
            if(obj_approx%hamilt%corr%intmap(l).eq.j) then
             do m=1,obj_approx%hamilt%corr%ncol(l)
              q=q*state(obj_approx%hamilt%corr%term(l,m))
             end do
            end if
            s=s+q
           end do
           t=0
           do l=1,obj_approx%nsites
            t=t+state(l)
           end do
           totr=totr+p*s*exp((chemp*t-obj_approx%ener(obj_approx,state)-h0)/(kb*temp))
           corf=corf+p*exp((chemp*t-obj_approx%ener(obj_approx,state)-h0)/(kb*temp))
          end do
          totr=-totr/(kb*temp)
          df(j,i)=totl-totr
          !write(*,*) corf, df(i,j)
         end do
        end do
	end subroutine fdjac

	subroutine fdjac2(n,x,fvec,np,df)
	integer n, np, nmax
	real*8 df(np,np),fvec(n),x(n),eps
	parameter (nmax=40,eps=1.0d-6)
	integer i, j
	real*8 h, temp, f(nmax)

	do j=1,n
	 temp=x(j)
	 h=0.1d0*eps*abs(temp)
	 if(h.eq.0)h=eps
	 x(j)=temp+h
	 call funcv(f)
	 x(j)=temp
	 do i=1,n
	  df(i,j)=(f(i)-fvec(i))/h
	 end do
	end do
	end subroutine fdjac2	
        end
