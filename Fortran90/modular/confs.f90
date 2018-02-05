	subroutine confs(state,l)
        ! Confs is a subroutine which generates all the possible cluster
        ! configurations and stores them into a matrix called "states".
        ! On top of that, the state vector "state" is updated with
        ! values as taken from the "l"-th row of the matrix "states".
        ! When called with argument l=0, confs calculates and stores
        ! "states"; else, i.e. l>0, "state" is updated.
	use constants
	implicit none
	integer cls
	integer, dimension(nsi) :: state
	integer, dimension(2**nsi,nsi) :: states
	integer i, j, k, l, mod, f
	save states

	if(l.eq.0) then
	 do i=1,2**nsi
	  do j=1,nsi
        ! "f" here stands for a frequency. The idea is to build "states"
        ! by columns. Each column is made of 1s and 0s which alternate
        ! with a given frequency, e.g. the first column is the less
        ! frequent and only alternates 1s and 0s once (11111...00000),
        ! while the last has the highest frequency (trades 1s and 0s 
        ! at each row number, 101010101010...). 
	   f=2**(nsi-j)
	   k=i
	   if(k.gt.f) then
	    do while(k.gt.f)
	     k=k-2*f
	    end do
	   end if
	   if(k.gt.0) then
	    states(i,j)=1
	   else
	    states(i,j)=0
	   end if
	  end do
 	 end do
	else
	 do i=1,nsi
	  state(i)=states(l,i)
	 end do
	end if
	end subroutine
