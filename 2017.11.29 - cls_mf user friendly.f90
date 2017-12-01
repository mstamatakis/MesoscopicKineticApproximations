	program cluster_mean_field
	implicit none
	integer i, lev, cls, alpha, beta
	real*8 r, rmax, latcon, pi, t, xt, yt
	real*8, allocatable, dimension(:) :: x, y 
	integer, allocatable, dimension(:) :: state
	real*8 h, j, g, hamiltonian
	real*8 energy, temp, part, chemp, kb
	real*8 s1, s2, diff
	integer c

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
	pi=3.14d0
	t=pi/3.d0
	cls=((2*lev+1)**2)-(lev**2)-lev
	rmax=lev*latcon
	allocate(x(cls),y(cls))
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

	! Calculating the MF correction parameter
	allocate(state(cls))
	kb=8.6173303*10E-5
	kb=kb*temp
	diff=10E5
	do while(diff.ge.10E-5)
		s1=0.d0
		s2=0.d0
		part=0.d0
		do i=1,2**cls
			call confs(state,cls,i)
			part=part+exp(chemp*sum(state)-hamiltonian(lev,cls,x,y,latcon,h,j,g,state)/kb)
			s1=s1+state(c)*exp(chemp*sum(state)-hamiltonian(lev,cls,x,y,latcon,h,j,g,state)/kb)
			s2=s2+state(cls)*exp(chemp*sum(state)-hamiltonian(lev,cls,x,y,latcon,h,j,g,state)/kb)
		end do
		s1=s1/part
		s2=s2/part
		diff=abs(s1-s2)
		if(diff.gt.10E-05) then
			if((s1-s2).lt.0) then
				g=g+0.0001d0
			else
				g=g-0.0001d0
			end if
		else
			write(*,*) 'field correction parameter=', g, 'eV'
		end if
	end do
	end program

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
