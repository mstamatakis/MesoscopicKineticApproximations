	real*8 function hamiltonian(n,x,state)
        ! Kikuchi 2nn C1 Hamiltonian
        ! Input: the number of parameters "n"; the vector of parameters
        ! "x"; and the state vector "state".
	use constants
	implicit none
	integer i, k, n
	integer, dimension(nsi) :: state
	real*8, dimension(4) :: x
	
	hamiltonian=hads*state(1)
	do i=2,7
	 hamiltonian=hamiltonian+(hads+x(1))*state(i)+jint*state(i)*state(1)
	end do
        do i=8,13
         hamiltonian=hamiltonian+(hads+x(2))*state(i)
        end do
	do i=2,6
	 hamiltonian=hamiltonian+(jint+x(3))*state(i)*state(i+1)
	end do
	hamiltonian=hamiltonian+(jint+x(3))*state(2)*state(7)
        do i=8,12
         hamiltonian=hamiltonian+(jint+x(4))*state(i)*(state(i-5)+state(i-6))
        end do	
	hamiltonian=hamiltonian+(jint+x(4))*state(13)*state(7)
	hamiltonian=hamiltonian+(jint+x(4))*state(13)*state(2)
	end function
