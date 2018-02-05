	subroutine confs(state,l)
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
