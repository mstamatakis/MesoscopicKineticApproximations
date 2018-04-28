subroutine approx_initialise_k3nnc2(obj_approx)

    use meso_approx
    use global_constants
    
    implicit none
    type (approximation) :: obj_approx
    integer nsites, nterms, norig, ncorc, nbodymax, ntm, i
        
    obj_approx%approxname = 'K3NNC2'
	nsites = 19
	obj_approx%nsites = nsites
	call obj_approx%populate_allstates()
	allocate(obj_approx%allenergs(2**nsites),source=0.d0)
	
	! Populate Hamiltonian data-structures
	nterms = 19 + & ! 19 single body
        6 + & ! 6 pairwise 1NN (center-1NN)
        6 + & ! 6  pairwise 1NN (1NN-1NN)
       12 + & ! 12  pairwise 1NN (1NN-2NN)
        6 + & ! 6  pairwise 1NN (1NN-3NN)
       12 + & ! 12  pairwise 1NN (2NN-3NN)
        6 + & ! 6  pairwise 2NN (1NN-1NN) which are ghost interactions (no original terms or correction for center-2NN)
        6 + & ! 6  pairwise 2NN (2NN-2NN) which are also ghost interactions 
       12     ! 12  pairwise 2NN (2NN-3NN) which are also ghost interactions 
	norig = 2 ! adsorption energy and 1NN
	ncorc = 3 + & ! corrections to adsorption energies of edge 1NN, 2NN, 3NN sites 
        3 + & ! 1NN interactions of edge sites (1NN-1NN, 1NN-2NN, 2NN-3NN)
        4     ! 2NN interactions of edge sites (1NN-1NN, 2NN-2NN, 2NN-3NN, 3NN-3NN)
	nbodymax = 2
	obj_approx%hamilt%nterms = nterms
	obj_approx%hamilt%norig = norig
	obj_approx%hamilt%ncorc = ncorc
	obj_approx%hamilt%nbodymax = nbodymax
	
	allocate(obj_approx%hamilt%interaction(nterms,nbodymax),source=0) ! encoding interaction terms
	allocate(obj_approx%hamilt%internbody(nterms),source=0) ! number of bodies for each interaction term
	obj_approx%hamilt%H0 = 0.d0 ! constant term in the Hamiltonian
	allocate(obj_approx%hamilt%origpars(0:norig),source=0.d0)
	obj_approx%hamilt%origpars(1) = hads
	obj_approx%hamilt%origpars(2) = Jint ! value of an original interaction term in the Hamiltonian
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
    ! Encode all the 1NN interaction terms
	do i = 2,6
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+1/) ! 1NN-1NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/7,2/) ! final 1NN-1NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
	do i = 8,12
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i-6/) ! 1NN-2NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i-6+1/)
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/13,7/) ! pre-final 1NN-2NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/13,2/) ! final 1NN-2NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
	do i = 2,7
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+12/) ! 1NN-3NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
    do i = 8,12
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+6/) ! 2NN-3NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+6+1/)
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/13,19/) ! pre-final 2NN-3NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/13,14/) ! final 2NN-3NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
    ! Encode all the 2NN interaction terms
	do i = 2,5
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+2/) ! 1NN-1NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
        if (i+4 > 7) cycle ! focusing on 1NN ring of sites
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+4/) ! 1NN-1NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	do i = 8,12
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i+1/) ! 2NN-2NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/13,8/) ! final 2NN-2NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2
    
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/14,7/) ! first 2NN-3NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2    
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/14,3/) ! second 2NN-3NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2    
	do i = 15,18
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i-13/) ! 2NN-3NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
		ntm = ntm + 1
		obj_approx%hamilt%interaction(ntm,1:2) = (/i,i-13+2/) ! 2NN-3NN pairwise terms
		obj_approx%hamilt%internbody(ntm) = 2
	enddo
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/19,6/) ! pre-final 2NN-2NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2    
	ntm = ntm + 1
	obj_approx%hamilt%interaction(ntm,1:2) = (/19,2/) ! final 2NN-2NN pairwise term
	obj_approx%hamilt%internbody(ntm) = 2    
    
	! Integer arrays to be used as pointers to terms and corrections:
	obj_approx%hamilt%origterms(1:nsites) = 1
	obj_approx%hamilt%origterms(nsites+1:61) = 2		 
	obj_approx%hamilt%corcterms(2:7) = 1 ! adsorption energy correction for sites in 1NN ring
	obj_approx%hamilt%corcterms(8:13) = 2 ! adsorption energy correction for sites in 2NN ring
	obj_approx%hamilt%corcterms(14:19) = 3 ! adsorption energy correction for sites in 3NN ring
	obj_approx%hamilt%corcterms(26:31) = 4 ! 1NN pairwise correction for 1NN-1NN interactions
    obj_approx%hamilt%corcterms(32:43) = 5 ! 1NN pairwise correction for 1NN-2NN interactions
    obj_approx%hamilt%corcterms(44:49) = 6 ! 1NN pairwise correction for 1NN-3NN interactions
    obj_approx%hamilt%corcterms(50:61) = 7 ! 1NN pairwise correction for 2NN-3NN interactions
	obj_approx%hamilt%corcterms(62:67) = 8 ! 2NN pairwise correction for 1NN-1NN interactions
	obj_approx%hamilt%corcterms(68:73) = 9 ! 2NN pairwise correction for 2NN-2NN interactions
	obj_approx%hamilt%corcterms(74:85) = 10 ! 2NN pairwise correction for 2NN-3NN interactions
        
    ! Populate equations data-structures
    obj_approx%eqns%neqns = 10
    obj_approx%eqns%nterms = 13
    obj_approx%eqns%nbodymax = 2
        
    allocate(obj_approx%eqns%residual(obj_approx%eqns%neqns),source=1.d+10) ! initialise residuals to something large
    allocate(obj_approx%eqns%jacobian(obj_approx%eqns%neqns,obj_approx%eqns%neqns),source=0.d+0) ! initialise jacobian to zero
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

	obj_approx%eqns%corrlnbody(3) = 1
	obj_approx%eqns%correlation(3,1) = 8 ! need one more term for the 2nd equation sigma(1) = sigma(8)
	obj_approx%eqns%lhs(2) = 1
	obj_approx%eqns%rhs(2) = 3 ! and finally we encode the equation
    
	obj_approx%eqns%corrlnbody(4) = 1
	obj_approx%eqns%correlation(4,1) = 14 ! need one more term for the 3rd equation sigma(1) = sigma(14)
	obj_approx%eqns%lhs(3) = 1
	obj_approx%eqns%rhs(3) = 4 ! and finally we encode the equation

    obj_approx%eqns%corrlnbody(5:6) = 2
	obj_approx%eqns%correlation(5,1:2) = (/1,2/)
	obj_approx%eqns%correlation(6,1:2) = (/2,3/) ! the terms of the 4th equation sigma(1)*sigma(2) = sigma(2)*sigma(3)
	obj_approx%eqns%lhs(4) = 5
	obj_approx%eqns%rhs(4) = 6 ! and finally we encode the equation

	obj_approx%eqns%corrlnbody(7) = 2
	obj_approx%eqns%correlation(7,1:2) = (/2,8/) ! need one more term for the 5th equation sigma(1)*sigma(2) = sigma(2)*sigma(8)
	obj_approx%eqns%lhs(5) = 5
	obj_approx%eqns%rhs(5) = 7 ! and finally we encode the equation

	obj_approx%eqns%corrlnbody(8) = 2
	obj_approx%eqns%correlation(8,1:2) = (/2,14/) ! need one more term for the 6th equation sigma(1)*sigma(2) = sigma(2)*sigma(14)
	obj_approx%eqns%lhs(6) = 5
	obj_approx%eqns%rhs(6) = 8 ! and finally we encode the equation

	obj_approx%eqns%corrlnbody(9) = 2
	obj_approx%eqns%correlation(9,1:2) = (/8,14/) ! need one more term for the 7th equation sigma(1)*sigma(2) = sigma(8)*sigma(14)
	obj_approx%eqns%lhs(7) = 5
	obj_approx%eqns%rhs(7) = 9 ! and finally we encode the equation

    obj_approx%eqns%corrlnbody(10:11) = 2
	obj_approx%eqns%correlation(10,1:2) = (/1,8/)
	obj_approx%eqns%correlation(11,1:2) = (/2,4/) ! the terms of the 8th equation sigma(1)*sigma(8) = sigma(2)*sigma(4)
	obj_approx%eqns%lhs(8) = 10
	obj_approx%eqns%rhs(8) = 11 ! and finally we encode the equation
    
	obj_approx%eqns%corrlnbody(12) = 2
	obj_approx%eqns%correlation(12,1:2) = (/8,9/) ! need one more term for the 9th equation sigma(1)*sigma(8) = sigma(8)*sigma(9)
	obj_approx%eqns%lhs(9) = 10
	obj_approx%eqns%rhs(9) = 12 ! and finally we encode the equation    
    
	obj_approx%eqns%corrlnbody(13) = 2
	obj_approx%eqns%correlation(13,1:2) = (/3,14/) ! need one more term for the 9th equation sigma(1)*sigma(8) = sigma(3)*sigma(14)
	obj_approx%eqns%lhs(10) = 10
	obj_approx%eqns%rhs(10) = 13 ! and finally we encode the equation    

    return

end subroutine approx_initialise_k3nnc2