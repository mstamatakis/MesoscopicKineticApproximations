	real*8 function energy(appr,state)
        ! Input: the number of parameters "n"; the vector of parameters
        ! "x"; and the state vector "state".
        use commons
	use meso_approx
	implicit none
	integer i, j, s
	integer, dimension(nsites) :: state
        type(approximation_type) :: appr 
	
        energy=0.d0
        ! Hamiltonian original block
        do i=1,2**nsites
         s=1
         do j=1,nsites
          if(appr%hamilt%orig%term(i,j).ne.0)then
           s=s*state(appr%hamilt%orig%term(i,j))
          end if
         end do
         energy=energy+s*appr%hamilt%orig%value(i)
        end do
        ! Hamiltonian correction block
        do i=1,2**nsites
         s=1
         do j=1,nsites
          if(appr%hamilt%corr%term(i,j).ne.0)then
           s=s*state(appr%hamilt%corr%term(i,j))
          end if
         end do
         if(appr%hamilt%corr%intmap(i).ne.0) then
          energy=energy+s*appr%hamilt%corr%value(appr%hamilt%corr%intmap(i))
         end if
        end do        
	end function
