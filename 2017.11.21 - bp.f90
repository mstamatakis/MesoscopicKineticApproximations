	program bethe_peierls
	implicit none
	real*8 chemp,temp,h0,h,j,g, hamiltonian
	integer, dimension(7) :: state
	
	h0=-0.027d0
	h=-1.2d0
	j=0.3d0
	g=0.d0
	chemp=1.d0
	temp=500.d0
	call partfun(h0,h,j,g,chemp,temp)
	end program

	real*8 function hamiltonian(state,h0,h,j,g)
	! The legend here is: 'g' stands for 'correction'
	implicit none
	integer i, k
	integer, dimension(7) :: state
	real*8 h0, h, j, g

	hamiltonian=state(1)*h+h0
	do i=1,6
	 hamiltonian=hamiltonian+(g+h)*state(i+1)
	end do
	do i=2,7
	  hamiltonian=hamiltonian+j*state(1)*state(i)
	end do
	end function

	subroutine partfun(h0,h,j,g,chemp,temp)
	implicit none
	integer, dimension(2**7,7) :: sts
	integer, dimension(7) :: state
	integer i, k
	real*8 h0, h, j, g, kb, temp, chemp, part, s1, s2, diff
	real*8 hamiltonian

	kb=8.6173303*10E-5
	! Inverse of beta below.
	kb=kb*temp
	call conf(sts)
	diff=10E5
	do while(diff.ge.10E-5)
	 s1=0.d0
	 s2=0.d0
	 part=0.d0
	 do i=1,2**7
	  do k=1,7
	   state(k)=sts(i,k)
 	  end do
	  part=part+exp(chemp*sum(state)-hamiltonian(state,h0,h,j,g)/kb)
	  s1=s1+state(1)*exp(chemp*sum(state)-hamiltonian(state,h0,h,j,g)/kb)
	  s2=s2+state(2)*exp(chemp*sum(state)-hamiltonian(state,h0,h,j,g)/kb)
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
	end subroutine

	subroutine conf(v)
	implicit none
	integer, dimension(2**7,7) :: v
	integer i, j, k, l, f, mod

	do i= 1, 7
	 f= 2**(7-i)
	 do j= 1, 2**i
	  do k= 1, f
	   l= k+f*(j-1)
	   v(l,i)= mod(j-1,2)
          end do
	 end do
	end do
	end subroutine

	integer function mod(j,m)
	implicit none
	integer j, m

	if(j.lt.m) then
	 j= j
	else
	 j= j-m
	end if
	end function
