        real*8 function corr(v,m,n,x)
        ! This calculates correlation functions of whatever order.
        ! Input: a vector of variables (even 1-dimensional); its
        ! dimension "m"; the number of parameters in the Hamiltonian
        ! "n"; and the vector of parameters "x".
	use constants
        implicit none
        integer i, k, m, n, s
        integer, dimension(m) :: v
        integer, dimension(nsi) :: state
        real*8 hamiltonian
	real*8, dimension(n) :: x

        state=0
        corr=0.d0
        do i=1,2**nsi
         call confs(state,i)
	 s=1
	 do k=1,m
	  s=s*state(v(k))
	 end do
         corr=corr+s*exp((chemp*sum(state)-hamiltonian(n,x,state))/(kb*temp))
        end do
        end function
