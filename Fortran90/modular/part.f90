	real*8 function part(n,x)
	! This calculates the partition function the brute-force way
	use constants
	implicit none
	integer i, n
	real*8 h, j
	integer, dimension(nsi) :: state
	real*8 hamiltonian 
	real*8, dimension(n) :: x

	state=0
	part=0.d0
        do i=1,2**nsi
         call confs(state,i)
         part=part+exp((chemp*sum(state)-hamiltonian(n,x,state))/(kb*temp))
        end do
	end function
