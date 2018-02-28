	subroutine confs(state,l)
        ! Confs is a subroutine which generates all the possible cluster
        ! configurations and stores them into a matrix called "states".
        ! On top of that, the state vector "state" is updated with
        ! values as taken from the "l"-th row of the matrix "states".
        ! When called with argument l=0, confs calculates and stores
        ! "states"; else, i.e. l>0, "state" is updated.
        use meso_approx
	implicit none
	integer, dimension(nsites) :: state
	integer, allocatable, dimension(:,:) :: states
	integer i, j, k, l, f
	save states

	if(l.eq.0) then
         allocate(states(2**nsites,nsites))
        !$OMP PARALLEL
        !$OMP DO
	 do i=1,2**nsites
	  do j=1,nsites
        ! "f" here stands for a frequency. The idea is to build "states"
        ! by columns. Each column is made of 1s and 0s which alternate
        ! with a given frequency, e.g. the first column is the less
        ! frequent and only alternates 1s and 0s once (11111...00000),
        ! while the last has the highest frequency (trades 1s and 0s 
        ! at each row number, 101010101010...). 
	   f=2**(nsites-j)
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
        !$OMP END DO
        !$OMP END PARALLEL
	else
        !$OMP PARALLEL
        !$OMP DO
	 do i=1,nsites
	  state(i)=states(l,i)
	 end do
        !$OMP END DO
        !$OMP END PARALLEL
	end if
	end subroutine
