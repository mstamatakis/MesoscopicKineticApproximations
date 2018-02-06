	program bpec
	use constants
	implicit none
	logical check
	integer i, j
	integer, dimension(nsi) :: state
	integer, allocatable, dimension(:) :: vec
	real*8 hamiltonian, part, corr
	real*8 cov, k
	real*8, dimension(2) :: x, fvec, tmp
	real*8 jac(2,2)

	x(1)=0.72657562404266485d0
	x(2)=-9.4768471953067229d-2
	jac=0.d0
	chemp=0.07d0
	call jacobian(2,x,jac)
	do i=1,2
	 write(*,*) (jac(i,j),j=1,2)
	end do
	jac=0.d0
	call jbn(2,x,jac)
	write(*,*)
	do i=1,2
	 write(*,*) (jac(i,j),j=1,2)
	end do
stop
	check=.false.
	state=0
	allocate(vec(1))
	chemp=-1.40
	do i=1,240
	 chemp=chemp+0.01d0
	 call newton(x,2) 
	 cov=0.d0
	 do j=1,nsi
	  vec(1)=j
	  cov=cov+corr(vec,1,1,x)/part(2,x)
	 end do
	 k=nsi*1.d0
	 cov=cov/k
	 write(16,*) chemp, cov
	end do
	end program

	subroutine newton(x,p)
	implicit none
	integer i, k, p, d, c
	integer, dimension(p) :: indx
	real*8, dimension(p) :: x, f, v
	real*8, dimension(p,p) :: jac
	real*8 corr, hamiltonian, part, tol
	
	write(*,*) 'NEWTON ALGORITHM STARTS NOW'
	jac=0.d0
	tol=0.d0
	write(*,*) 'Fs'
	do i=1,p
	 call funcv(p,x,f)
	 tol=tol+f(i)**2
	 write(*,*) f(i)
	end do
	tol=sqrt(tol)
	write(*,*) 'TOLERANCE:', tol
	do while(tol.gt.0.001d0)
	 call jbn(p,x,jac)
	 write(*,*) 'JACOBIAN'
	 do i=1,p	
	  write(*,*) (jac(i,k),k=1,p)
	 end do
	 v=-f
	 call ludcmp(jac,p,indx,d,c) 
	 call lubksb(jac,p,indx,v)
	 do i=1,p
	  x(i)=x(i)+v(i)
	  write(*,*) x(i)
	 end do
	 tol=0.d0
	 call funcv(p,x,f)
         do i=1,p
          tol=tol+f(i)**2
         end do
	 tol=sqrt(tol)
	 write(*,*) 'TOLERANCE:', tol
	end do
	end subroutine

	subroutine jacobian(n,x,jac)
	use constants
	implicit none
	integer i, k, n 
	real*8 d
	real*8, dimension(n) :: x, xfwd, xbwd, fvec
	real*8, dimension(n,n) :: jac
	real*8 func, corr, hamiltonian, part
	
	fvec=0.d0
	jac=0.d0
	d=sqrt(2.2d-18)
	do i=1,n
	 do k=1,n
	  xfwd=x
	  xfwd(k)=xfwd(k)*(1.d0+d)
	  xbwd=x
	  xbwd(k)=xbwd(k)*(1.d0-d)
	  call funcv(n,xfwd,fvec)
	  jac(i,k)=fvec(i)/(2.d0*x(k)*d)
	  call funcv(n,xbwd,fvec)
	  jac(i,k)=jac(i,k)-fvec(i)/(2.d0*x(k)*d)
	 end do
	end do
	end subroutine

	subroutine funcv(n,x,fvec)
	implicit none
	integer n, v1(1), v2(2)
	real*8 x(n), fvec(n), corr	

    
	v1=1
	fvec(1)=log(corr(v1,1,n,x))
	v1=2
	fvec(1)=fvec(1)-log(corr(v1,1,n,x))
	v2(1)=1
	v2(2)=2
	fvec(2)=log(corr(v2,2,n,x))
	v2(1)=2
	v2(2)=3
	fvec(2)=fvec(2)-log(corr(v2,2,n,x))
    
    
	!v1=1
	!fvec(1)=corr(v1,1,n,x)
	!v1=2
	!fvec(1)=fvec(1)-corr(v1,1,n,x)
	!v2(1)=1
	!v2(2)=2
	!fvec(2)=corr(v2,2,n,x)
	!v2(1)=2
	!v2(2)=3
	!fvec(2)=fvec(2)-corr(v2,2,n,x)
	end subroutine funcv

	subroutine confs(state,l)
	use constants
	implicit none
	integer cls
	integer, dimension(nsi) :: state
	integer i, k, l, mod, f

	do i=1,nsi
	 f=2**(nsi-i)
	 k=l
	 if(k.gt.f) then
	  do while(k.gt.f)
	   k=k-2*f
	  end do
	 end if
	 if(k.gt.0) then
	  state(i)=1
	 else
	  state(i)=0
	 end if
	end do
	end subroutine

	real*8 function part(n,x)
	! This calculates the partition function the brute-force way
	use constants
	implicit none
	integer i, n
	real*8 h, j
	integer, dimension(nsi) :: state
	real*8 hamiltonian, b 
	real*8, dimension(n) :: x

	state=0
	part=0.d0
        do i=1,2**nsi
         call confs(state,i)
         part=part+exp((chemp*sum(state)-hamiltonian(n,x,state))/(kb*temp))
        end do
	end function

        real*8 function corr(v,m,n,x)
	use constants
        implicit none
        integer i, k, m, n, s
        integer, dimension(m) :: v
        integer, dimension(nsi) :: state
        real*8 hamiltonian
	real*8, dimension(n) :: x

        state=0
        corr=0.d0
        do i=1,2**nsi
         call confs(state,i)
	 s=1
	 do k=1,m
	  s=s*state(v(k))
	 end do
         corr=corr+s*exp((chemp*sum(state)-hamiltonian(n,x,state))/(kb*temp))
        end do
        end function

        subroutine ludcmp(a,n,indx,d,code)
        ! By J-P Moreau, Paris
        parameter(nmax=100,tiny=1.5D-16)
        real*8  amax, dum, sum, a(n,n), vv(nmax)
        integer code, d, indx(n)
       
        d=1; code=0
       
        do i=1,n
          amax=0.d0
          do j=1,n
            if (dabs(a(i,j)).gt.amax) amax=dabs(a(i,j))
          end do ! j loop
          if(amax.lt.tiny) then
            code = 1
            return
          end if
          vv(i) = 1.d0/amax
        end do ! i loop
       
        do j=1,n
          do i=1,j-1
            sum=a(i,j)
            do k=1,i-1
              sum=sum-a(i,k)*a(k,j) 
            end do ! k loop
            a(i,j)=sum
          end do ! i loop
          amax=0.d0
          do i=j,n
            sum=a(i,j)
            do k=1,j-1
              sum=sum-a(i,k)*a(k,j) 
            end do ! k loop
            a(i,j)=sum
            dum=vv(i)*dabs(sum)
            if(dum.ge.amax) then
              imax=i
              amax=dum
            end if
          end do ! i loop  
         
          if(j.ne.imax) then
            do k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
            end do ! k loop
            d=-d
            vv(imax)=vv(j)
          end if
        
          indx(j)=imax
          if(dabs(a(j,j))<tiny) a(j,j)=tiny
        
          if(j.ne.n) then
            dum=1.d0/a(j,j)
            do i=j+1,n
              a(i,j) = a(i,j)*dum
            end do ! i loop
          end if
        end do ! j loop
       
        return
        end subroutine ludcmp
       
        subroutine lubksb(a,n,indx,b)
        ! By J-P Moreau, Paris
        real*8  sum, a(n,n), b(n)
        integer indx(n)
       
        ii=0
       
        do i=1,n
          ll=indx(i)
          sum=b(ll)
          b(ll)=b(i)
          if(ii.ne.0) then
            do j=ii,i-1
              sum=sum-a(i,j)*b(j)
            end do ! j loop
          else if(sum.ne.0.d0) then
            ii=i
          end if
          b(i)=sum
        end do ! i loop
       
        do i=n,1,-1
          sum=b(i)
          if(i<n) then
            do j=i+1,n
              sum=sum-a(i,j)*b(j)
            end do ! j loop
          end if
          b(i)=sum/a(i,i)
        end do ! i loop
       
        return
        end subroutine lubksb
	
	real*8 function hamiltonian(n,x,state)
	use constants
	implicit none
	integer i, k, n
	integer, dimension(nsi) :: state
	real*8, dimension(n) :: x
	
	hamiltonian=hads*state(1)
	do i=2,7
	 hamiltonian=hamiltonian+(hads+x(1))*state(i)+jint*state(i)*state(1)
	end do
	do i=2,6
	 hamiltonian=hamiltonian+(jint+x(2))*state(i)*state(i+1)
	end do
	hamiltonian=hamiltonian+(jint+x(2))*state(7)*state(2)
	end function
	
	subroutine jbn(n,x,jac)
	use constants
	implicit none
	real*8 corr, hamiltonian
	integer i, v1(1), v2(2), v3(3), v4(4), n
	real*8 x(n), jac(n,n), t1, t2

	t1=0.0d0
	do i=2,7
	 v2(1)=1
	 v2(2)=i
	 t1=t1+corr(v2,2,n,x)
	end do
	v1=1
	t1=-t1/((kb*temp)*corr(v1,1,n,x))
	t2=0.0d0
	do i=2,7
	 v2(1)=2
	 v2(2)=i
	 t2=t2+corr(v2,2,n,x)
	end do
	v1=2
	t2=-t2/((kb*temp)*corr(v1,1,n,x))
	jac(1,1)=t1-t2
	
        t1=0.0d0
        do i=2,6
         v3(1)=1
         v3(2)=i
         v3(3)=i+1
         t1=t1+corr(v3,3,n,x)
        end do
        v3(1)=1
        v3(2)=2
        v3(3)=7
        t1=t1+corr(v3,3,n,x)
        v1=1
        t1=-t1/((kb*temp)*corr(v1,1,n,x))
        t2=0.0d0
        do i=2,6
         v3(1)=2
         v3(2)=i
         v3(3)=i+1
         t2=t2+corr(v3,3,n,x)
        end do
        v3(1)=2
        v3(2)=2
        v3(3)=7
        t2=t2+corr(v3,3,n,x)
        v1=2
        t2=-t2/((kb*temp)*corr(v1,1,n,x))
        jac(1,2)=t1-t2

	t1=0.0d0
	do i=2,7
	 v3(1)=1
	 v3(2)=2
	 v3(3)=i
	 t1=t1+corr(v3,3,n,x)
	end do
	v2(1)=1
	v2(2)=2
	t1=-t1/((kb*temp)*corr(v2,2,n,x))
	t2=0.0d0
	do i=2,7
	 v3(1)=2
	 v3(2)=3
	 v3(3)=i
	 t2=t2+corr(v3,3,n,x)
	end do
	v2(1)=2
	v2(2)=3
	t2=-t2/((kb*temp)*corr(v2,2,n,x))
	jac(2,1)=t1-t2

	t1=0.d0
	do i=2,6
	 v4(1)=1
	 v4(2)=2
	 v4(3)=i
	 v4(4)=i+1
	 t1=t1+corr(v4,4,n,x)
	end do
	v4(1)=1
	v4(2)=2
	v4(3)=2
	v4(4)=7
	t1=t1+corr(v4,4,n,x)
        v2(1)=1
        v2(2)=2
        t1=-t1/((kb*temp)*corr(v2,2,n,x))
        t2=0.d0
        do i=2,6
         v4(1)=2
         v4(2)=3
         v4(3)=i
         v4(4)=i+1
         t2=t2+corr(v4,4,n,x)
        end do
        v4(1)=2
        v4(2)=3
        v4(3)=2
        v4(4)=7
        t2=t2+corr(v4,4,n,x)
        v2(1)=2
        v2(2)=3
        t2=-t2/((kb*temp)*corr(v2,2,n,x))
	jac(2,2)=t1-t2
	end subroutine
