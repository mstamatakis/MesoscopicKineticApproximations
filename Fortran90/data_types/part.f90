	real*8 function part(appr)
	! This calculates the partition function the brute-force way.
        ! Input: number of parameters, vector of parameters
	use commons
        use meso_approx
	implicit none
	integer i
	integer, dimension(nsites) :: state
	real*8 energy
        type (approximation_type) :: appr

	state=0
	part=0.d0
        do i=1,2**nsites
         call confs(state,i)
         part=part+exp((chemp*sum(state)-energy(appr,state))/(kb*temp))
        end do
	end function
