subroutine approx_initialise_bpec(obj_approx)

    use meso_approx
    use global_constants
    
    implicit none
    type (approximation) :: obj_approx
    integer nsites, nterms, norig, ncorc, nbodymax, i, ntm
        
    obj_approx%approxname = 'BPEC'
	nsites = 7
	obj_approx%nsites = nsites
	call obj_approx%populate_allstates()
	allocate(obj_approx%allenergs(2**nsites),source=0.d0)
	
	! Populate Hamiltonian data-structures
	nterms = 7 + 6 + 6 ! 7 single body, 6 pairwise (center-1NN), and 6 more pairwise (between 1NNs)
	norig = 2 ! adsorption energy and 1NN
	ncorc = 2 ! corrections to adsorption energies of edge sites and 1NN interactions of edge sites
	nbodymax = 2
	obj_approx%hamilt%nterms = nterms
	obj_approx%hamilt%norig = norig
	obj_approx%hamilt%ncorc = ncorc
	obj_approx%hamilt%nbodymax = nbodymax
	
	allocate(obj_approx%hamilt%interaction(nterms,nbodymax),source=0) ! encoding interaction terms
	allocate(obj_approx%hamilt%internbody(nterms),source=0) ! number of bodies for each interaction term
	obj_approx%hamilt%H0 = 0.d0 ! constant term in the Hamiltonian
	allocate(obj_approx%hamilt%origpars(0:norig),source=(/0.d0,hads,Jint/)) ! value of an original interaction term in the Hamiltonian
	allocate(obj_approx%hamilt%corcpars(0:ncorc),source=0.d0) ! value of a correction interaction term in the Hamiltonian
	allocate(obj_approx%hamilt%origterms(nterms),source=0) ! will be used as a pointer to an original interaction term in the Hamiltonian
	allocate(obj_approx%hamilt%corcterms(nterms),source=0) ! will be used as a pointer to a correction interaction term in the Hamiltonian
	
	ntm = 0 ! counter of terms
	do i = 1,nsites
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1) = i ! single body terms
		obj_approx%hamilt%internbody(ntm) = 1
	enddo
	do i = 2,7
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/1,i/) ! central-1NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	do i = 2,6
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+1/) ! 1NN-1NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/7,2/) ! final 1NN-1NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2

	! Integer arrays to be used as pointers to terms and corrections:
	obj_approx%hamilt%origterms(1:nsites) = 1
	obj_approx%hamilt%origterms(nsites+1:nterms) = 2		 
	obj_approx%hamilt%corcterms(2:7) = 1
	obj_approx%hamilt%corcterms(14:19) = 2
        
    ! Populate equations data-structures
    obj_approx%eqns%neqns = 2
    obj_approx%eqns%nterms = 4
    obj_approx%eqns%nbodymax = 2
        
    allocate(obj_approx%eqns%residual(obj_approx%eqns%neqns),source=1.d+10) ! initialise residuals to something large
    allocate(obj_approx%eqns%lhs(obj_approx%eqns%neqns),source=0)
    allocate(obj_approx%eqns%rhs(obj_approx%eqns%neqns),source=0) ! null info for left and right hand side

    allocate(obj_approx%eqns%correlation(obj_approx%eqns%nterms,obj_approx%eqns%nbodymax),source=0)
    allocate(obj_approx%eqns%corrlvalue(obj_approx%eqns%nterms),source=0.d0) ! all correlation values initialised to zero
    allocate(obj_approx%eqns%corrlnbody(obj_approx%eqns%nterms),source=0)
        
	obj_approx%eqns%corrlnbody(1:2) = 1
	obj_approx%eqns%correlation(1,1) = 1
	obj_approx%eqns%correlation(2,1) = 2 ! the terms of the 1st equation sigma(1) = sigma(2)
	obj_approx%eqns%lhs(1) = 1
	obj_approx%eqns%rhs(1) = 2 ! and finally we encode the equation
        
	obj_approx%eqns%corrlnbody(3:4) = 2
	obj_approx%eqns%correlation(3,1:2) = (/1,2/)
	obj_approx%eqns%correlation(4,1:2) = (/2,3/) ! the terms of the 2nd equation sigma(1)*sigma(2) = sigma(2)*sigma(3)
	obj_approx%eqns%lhs(2) = 3
	obj_approx%eqns%rhs(2) = 4 ! and finally we encode the equation
        
    return

end subroutine approx_initialise_bpec