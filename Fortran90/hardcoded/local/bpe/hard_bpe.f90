	program bp
	implicit none
	integer i
	integer, dimension(13) :: state
	integer, allocatable, dimension(:) :: vec
	real*8 h, j, hamiltonian, part, corr, func
	real*8 temp, chemp, kb, cov
	real*8, dimension(1) :: x

	! Initialisation
	write(*,*) '---------------------------------'
	write(*,*) '           ENERGETICS            '
	write(*,*) '---------------------------------'
	write(*,*) 'Please, insert the adsorption energy below'
	read(*,*) h
	write(*,*) 'And now, the Ising two-body interaction'
	read(*,*) j
	write(*,*) '---------------------------------'
	write(*,*) '          THERMODYNAMICS         '
	write(*,*) '---------------------------------'
	write(*,*) 'What is the chemical potential of the system?'
	read(*,*) chemp
	write(*,*) 'Last parameter to be inserted is the temperature of the system'
	read(*,*) temp
	write(*,*) '---------------------------------'
	write(*,*)
	write(*,*) 'Thanks!'
	write(*,*)

	! Calculating the MF correction parameters
	call newton(h,j,x,1,chemp,temp)
	write(*,*) x(1)	
	! Calculating the coverage
	cov=0.d0
	allocate(vec(1))
	do i=1,7
	 vec=i
	 cov=cov+corr(vec,1,h,j,x,1,chemp,temp)/part(h,j,x,1,chemp,temp)
	end do
	write(*,*) cov/7
	end program

	subroutine newton(h,j,x,p,chemp,temp)
	implicit none
	integer i, k, p, d, c
	integer, dimension(p) :: indx
	real*8, dimension(p) :: x, f, v
	real*8, dimension(p,p) :: jac
	real*8 tol, h, j, chemp, temp, func, corr, hamiltonian, part
	
	x=0.50d0
	jac=0.d0
	tol=0.d0
	write(*,*) 'Fs'
	do i=1,p
	 f(i)=func(i,h,j,x,p,chemp,temp)
	 tol=tol+f(i)**2
	 write(*,*) f(i)
	end do
	tol=sqrt(tol)
	write(*,*) 'TOLERANCE:', tol
	do while(tol.gt.0.001d0)
	 call jacobian(h,j,x,p,temp,chemp,jac)
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
         do i=1,p
          f(i)=func(i,h,j,x,p,chemp,temp)
          tol=tol+f(i)**2
         end do
	 tol=sqrt(tol)
	 write(*,*) 'TOLERANCE:', tol
	end do
	end subroutine

	subroutine jacobian(h,j,x,p,temp,chemp,jac)
	implicit none
	integer i, k, p
	real*8 d, h, j, temp, chemp
	real*8, dimension(p) :: x, xfwd, xbwd
	real*8, dimension(p,p) :: jac
	real*8 func, corr, hamiltonian, part
	real*8 test
	
	jac=0.d0
	d=sqrt(2.2d-18)
	do i=1,p
	 do k=1,p
	  xfwd=x
	  xfwd(k)=xfwd(k)*(1.d0+d)
	  xbwd=x
	  xbwd(k)=xbwd(k)*(1.d0-d)
	  jac(i,k)=(func(i,h,j,xfwd,p,chemp,temp)-func(i,h,j,xbwd,p,chemp,temp))/(2.d0*x(k)*d)
	 end do
	end do
	end subroutine

	subroutine confs(state,l)
	implicit none
	integer cls
	integer, dimension(13) :: state
	integer i, k, l, mod, f

	do i=1,13
	 f=2**(13-i)
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

	real*8 function hamiltonian(h,j,x,p,state)
	implicit none
	integer i, k, p
	real*8 h, j, a1, a2, b1, b2
	integer, dimension(13) :: state
	real*8, dimension(p) :: x
	
	hamiltonian=h*state(1)
	do i=2,7
	 hamiltonian=hamiltonian+(h+x(1))*state(i)+j*state(i)*state(1)
	end do
	do i=2,6
	 hamiltonian=hamiltonian+j*state(i)*state(i+1)
	end do
	hamiltonian=hamiltonian+j*state(7)*state(2)
	end function

	real*8 function part(h,j,x,p,chemp,temp)
	! This calculates the partition function the brute-force way
	implicit none
	integer i, p
	real*8 h, j
	integer, dimension(13) :: state
	real*8 hamiltonian, chemp, kb, temp
	real*8, dimension(p) :: x

	state=0
	kb=8.6173303d-5
	kb=kb*temp
	part=0.d0
        do i=1,2**13
         call confs(state,i)
         part=part+exp((chemp*sum(state)-hamiltonian(h,j,x,p,state))/kb)
        end do
	end function

        real*8 function corr(vec,n,h,j,x,p,chemp,temp)
        implicit none
        integer i, k, s, n, p
        real*8 h, j
        integer, dimension(n) :: vec
        integer, dimension(13) :: state
        real*8 hamiltonian, chemp, kb, temp
	real*8, dimension(p) :: x

        state=0
        kb=8.6173303d-5
        kb=kb*temp
        corr=0.d0
        do i=1,2**13
         call confs(state,i)
	 s=1
	 do k=1,n
	  s=s*state(vec(k))
	 end do
         corr=corr+s*exp((chemp*sum(state)-hamiltonian(h,j,x,p,state))/kb)
        end do
        end function
	
	real*8 function func(i,h,j,x,p,chemp,temp)
	implicit none
        integer s, p, i
        real*8 h, j
        integer, dimension(1) :: v1
	integer, dimension(1) :: v2
        integer, dimension(2) :: v3
	integer, dimension(2) :: v4
        integer, dimension(13) :: state
        real*8 hamiltonian, chemp, kb, temp, corr
	real*8, dimension(p) :: x

	if(i.eq.1) then
	 v1(1)=1
	 v2(1)=2
	 func=log(corr(v1,1,h,j,x,p,chemp,temp))-log(corr(v2,1,h,j,x,p,chemp,temp))
	end if
	if(i.eq.2) then
	 v1(1)=1
	 v2(1)=8
	 func=log(corr(v1,1,h,j,x,p,chemp,temp))-log(corr(v2,1,h,j,x,p,chemp,temp))
	end if
	if(i.eq.3) then
	 v3(1)=1
	 v3(2)=2
	 v4(1)=2
	 v4(2)=3
	 func=log(corr(v3,2,h,j,x,p,chemp,temp))-log(corr(v4,2,h,j,x,p,chemp,temp))
	end if
	if(i.eq.4) then
	 v3(1)=1
	 v3(2)=2
	 v4(1)=2
	 v4(2)=8
	 func=log(corr(v3,2,h,j,x,p,chemp,temp))-log(corr(v4,2,h,j,x,p,chemp,temp))
	end if
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
