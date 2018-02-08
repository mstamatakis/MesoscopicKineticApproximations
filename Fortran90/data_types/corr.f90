        real*8 function corr(v,m,appr)
        ! This calculates correlation functions of whatever order.
        ! Input: a vector of variables (even 1-dimensional); its
        ! dimension "m"; the number of parameters in the Hamiltonian
        ! "n"; and the vector of parameters "x".
	use commons
        use meso_approx
        implicit none
        integer i, k, m, s
        integer, dimension(m) :: v
        integer, dimension(nsites) :: state
        real*8 energy
        type (approximation_type) :: appr

        state=0
        corr=0.d0
        do i=1,2**nsites
         call confs(state,i)
	 s=1
	 do k=1,m
	  s=s*state(v(k))
	 end do
         corr=corr+s*exp((chemp*sum(state)-energy(appr,state))/(kb*temp))
        end do
        end function
