	program cluster_mean_field
	implicit none
	integer i, lev, cls, alpha, beta
	real*8 r, rmax, latcon, pi, t, xt, yt
	real*8, allocatable, dimension(:) :: x, y 
	integer, allocatable, dimension(:) :: state, v1, v2
	real*8 h, j, g, hamiltonian, part, corr
	real*8 energy, temp, chemp, kb
	real*8 s1, s2, diff
	integer c, imax, sgn
	real*8 ga, fa, gb, fb, gc, fc

	! Initialisation
	write(*,*) '---------------------------------'
	write(*,*) '       LATTICE GENERATION        '
	write(*,*) '---------------------------------'
	write(*,*) 'Please, insert the lattice constant below'
	read(*,*) latcon
	write(*,*) 'What is the degree of approximation of the simulation? (BP=1, KNN=2, ...)'
	read(*,*) lev
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

	! Writing down the lattice
	cls=((2*lev+1)**2)-(lev**2)-lev
	allocate(x(cls),y(cls))
	call lat_gen(lev,x,y,cls,latcon,c)

	! Calculating the MF correction parameter
	allocate(v1(1),v2(1))
	v1=c
	v2=cls
	!! Bisector method to solve for g
	ga=-1.d0
	gb=1.d0
	g=0.d0
        fa=log(corr(v1,1,lev,cls,h,j,ga,x,y,chemp,temp))-log(corr(v2,1,lev,cls,h,j,ga,x,y,chemp,temp))
	fb=log(corr(v1,1,lev,cls,h,j,gb,x,y,chemp,temp))-log(corr(v2,1,lev,cls,h,j,gb,x,y,chemp,temp))	
	if((fa.lt.0.d0).and.(fb.gt.0.d0)) then
	 i=0
	 do while(i.ne.1000)
	  i=i+1
	  gc=(gb+ga)/2.d0
	  fc=log(corr(v1,1,lev,cls,h,j,gc,x,y,chemp,temp))-log(corr(v2,1,lev,cls,h,j,gc,x,y,chemp,temp))
	  if((abs(fc).le.0.00001d0).or.(((gb-ga)/2.d0).le.0.00001d0)) then
	   g=gc
	   exit
	  end if
	  if(sgn(fc).eq.sgn(fa)) then
	   ga=gc
	  else
	   gb=gc
	  end if
	 end do
	end if
	write(*,*) g
	end program

        subroutine lat_gen(lev,x,y,cls,latcon,c)
	implicit none
	real*8 pi, t, rmax, latcon, xt, yt, r
	integer i, c, lev, cls, alpha, beta
	real*8, dimension(cls) :: x, y

        pi=3.14159265359d0
        t=pi/3.d0
        rmax=lev*latcon
        x=0.d0
        y=0.d0
        write(13,*) cls
        write(13,*)
        i=0
        do alpha=-lev,lev
         do beta=-lev,lev
          xt=alpha*latcon+beta*latcon*cos(t)
          yt=beta*latcon*sin(t)
          r=sqrt(xt*xt+yt*yt)
          if(r.le.(rmax+0.001d0)) then
           i=i+1
           x(i)=xt
           y(i)=yt
           write(13,*) 'H ', x(i), y(i), 0.d0
           if((alpha.eq.0).and.(beta.eq.0)) then
            c=i
           end if
          end if
         end do
        end do
        call system('mv fort.13 lattice.xyz')
	end subroutine

	subroutine confs(state,cls,l)
	implicit none
	integer cls
	integer, dimension(cls) :: state
	integer i, k, l, mod, f

	do i=1,cls
	 f=2**(cls-i)
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

	real*8 function hamiltonian(lev,cls,x,y,latcon,h,j,g,state)
	implicit none
	integer i, k, lev, cls
	real*8 h, j, g, latcon, ri, rk
	real*8, dimension(cls) :: x, y
	integer, dimension(cls) :: state
	
	hamiltonian=0.d0
	do i=1,cls
	 ri=sqrt(x(i)*x(i)+y(i)*y(i))
	 if(abs(ri-latcon*(lev-1)).le.0.001d0) then
	  do k=1,cls
	   rk=sqrt(x(k)*x(k)+y(k)*y(k))
	   if(abs(rk-latcon*(lev-1)).le.0.001d0) then
	    if((abs(ri-rk).le.(latcon+0.001d0)).and.(i.ne.k)) then
	     hamiltonian=hamiltonian+state(i)*state(k)*j/2.d0
	    end if
	   else
	     hamiltonian=hamiltonian+state(i)*state(k)*j
	   end if
	  end do
	  hamiltonian=hamiltonian+state(i)*h
	 else
	  hamiltonian=hamiltonian+state(i)*(h+g)
	 end if
	end do
	end function

	real*8 function part(lev,cls,h,j,g,x,y,chemp,temp)
	! This calculates the partition function the brute-force way
	implicit none
	integer i, cls, lev
	real*8 latcon, h, j, g
	real*8, dimension(cls) :: x, y
	integer, dimension(cls) :: state
	real*8 hamiltonian, chemp, kb, temp

	state=0
	kb=8.6173303*10E-5
	kb=kb*temp
	part=0.d0
        do i=1,2**cls
         call confs(state,cls,i)
         part=part+exp(chemp*sum(state)-hamiltonian(lev,cls,x,y,latcon,h,j,g,state)/kb)
        end do
	end function

        real*8 function corr(vec,n,lev,cls,h,j,g,x,y,chemp,temp)
        ! This calculates the partition function the brute-force way
        implicit none
        integer i, k, cls, lev, s, n
        real*8 latcon, h, j, g, tot
        real*8, dimension(cls) :: x, y
        integer, dimension(n) :: vec
        integer, dimension(cls) :: state
        real*8 hamiltonian, chemp, kb, temp

        state=0
        kb=8.6173303*10E-5
        kb=kb*temp
        corr=0.d0
        do i=1,2**cls
         call confs(state,cls,i)
	 s=1
	 do k=1,n
	  s=s*state(vec(k))
	 end do
         corr=corr+s*exp(chemp*sum(state)-hamiltonian(lev,cls,x,y,latcon,h,j,g,state)/kb)
        end do
        end function
	
	integer function sgn(arg)
	implicit none
	real*8 arg

	if(arg.lt.0.0) then
	 sgn=-1
	else
	 sgn=1
	end if
	end function
