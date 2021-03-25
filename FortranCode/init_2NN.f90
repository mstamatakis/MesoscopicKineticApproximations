subroutine approx_initialise_2NN(obj_approx)

        use global_constants
        use constants_module, only: nnam0, inttermsmax
        use calculation_setup_module
	use meso_approx
		
	implicit none
	class (approximation) :: obj_approx
        integer nsites, nterms, norig, ncorc, nbodymax, ntm, i
        character(nnam0), allocatable :: cor_names(:), orig_names(:)
        real(8), allocatable :: int_values(:)
        integer, allocatable :: int_term(:,:)

        interface
           subroutine  sigma1_sigma2_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface

        interface
           subroutine  sigma1_sigma8_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface

        interface
           subroutine  sigma1sigma2_sigma2sigma3_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface

        interface
           subroutine  sigma1sigma2_sigma2sigma8_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface

        interface
           subroutine  sigma1sigma8_sigma2sigma4_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface

        interface
           subroutine  sigma1sigma8_sigma8sigma9_2NN(obj_approx)
               import approximation
               class(approximation) :: obj_approx
           end subroutine
        end interface
        
	nsites =  cal_parser%get_nsites()
	obj_approx%nsites = nsites
	call obj_approx%populate_allstates()
	allocate(obj_approx%allenergs(2**nsites),source=0.d0)
		
	! Populate Hamiltonian data-structures
        nterms = cal_parser%get_nterms()  ! number of terms in the Hamiltonian
	norig = cal_parser%get_norigterms() ! adsorption energy and 1NN 
        ncorc = cal_parser%get_n_corr() ! corrections 
	nbodymax = 2
	obj_approx%hamilt%nterms = nterms
	obj_approx%hamilt%norig = norig
	obj_approx%hamilt%ncorc = ncorc
        obj_approx%hamilt%nbodymax = nbodymax
 
	allocate(obj_approx%hamilt%interaction(nterms,nbodymax),source = 0) ! encoding interaction terms
	allocate(obj_approx%hamilt%internbody(nterms),source = 0) ! number of bodies for each interaction term
        obj_approx%hamilt%H0 = 0.d0 ! constant term in the Hamiltonian
 
	allocate(obj_approx%hamilt%origpars(0:norig),source = 0.d0)
        allocate(int_values(0:norig), source = cal_parser%get_int_values())
        do i = 1, norig
             obj_approx%hamilt%origpars(i) = int_values(i) ! value of an original interaction term in the Hamiltonian
        enddo!
 
	allocate(obj_approx%hamilt%corcpars(0:ncorc),source = 0.d0) ! value of a correction interaction term in the Hamiltonian
	allocate(obj_approx%hamilt%origterms(nterms),source = 0) ! will be used as a pointer to an original interaction term in the Hamiltonian
	allocate(obj_approx%hamilt%corcterms(nterms),source = 0) ! will be used as a pointer to a correction interaction term in the Hamiltonian
		
	ntm = 0 ! counter of terms
	do i = 1,nsites
            ntm = ntm + 1
            obj_approx%hamilt%interaction(ntm,1) = i ! single body terms
            obj_approx%hamilt%internbody(ntm) = 1
        enddo

        allocate(int_term(1:inttermsmax,1:2), source = cal_parser%get_int_term())  
        do i = 1,cal_parser%get_n_int()
            ntm = ntm + 1
            obj_approx%hamilt%interaction(ntm,1:2) = int_term(i,1:2) ! two body pairwise terms
            obj_approx%hamilt%internbody(ntm) = 2
        enddo
        
        if (cal_parser%get_cluster_depth().eq.'2NN') then ! for the momment it is not really needed
           
	! Integer arrays to be used as pointers to terms and corrections:
            allocate(orig_names(0:norig), source = cal_parser%get_orig_terms_names())
            do i = 1,norig
                if (orig_names(i).eq.'SB') then 
                    obj_approx%hamilt%origterms(1:nsites) = 1
                elseif (orig_names(i).eq.'1NN') then
                    if ((ncorc.eq.2).or.(ncorc.eq.4)) then 
                        obj_approx%hamilt%origterms(nsites+1:nterms) = 2
                    elseif (ncorc.gt.4) then
                        obj_approx%hamilt%origterms(nsites+1:37) = 2
                    endif   
                endif
            enddo    
        
            allocate(cor_names(0:cal_parser%get_ncorrterms()), source = cal_parser%get_corr_terms_names())
            do i = 1,cal_parser%get_ncorrterms()        
                if ((cor_names(i).eq.'SB').or.(cor_names(i).eq.'SBWI')) then 
	            obj_approx%hamilt%corcterms(2:7) = 1 ! adsorption energy correction for sites in 1NN ring
                    obj_approx%hamilt%corcterms(8:13) = 2 ! adsorption energy correction for sites in 2NN ring
                elseif (cor_names(i).eq.'1NN') then
	            obj_approx%hamilt%corcterms(20:25) = 3 ! 1NN pairwise correction for 1NN-1NN interactions
	            obj_approx%hamilt%corcterms(26:37) = 4 ! 1NN pairwise correction for 1NN-2NN interactions
                elseif (cor_names(i).eq.'2NN') then
                    obj_approx%hamilt%corcterms(38:43) = 5 ! 2NN pairwise correction for 1NN-1NN interactions
	            obj_approx%hamilt%corcterms(44:49) = 6 ! 2NN pairwise correction for 1NN-1NN interactions
                endif
            enddo

          ! Populate equations data-structures
            obj_approx%eqns%neqns = cal_parser%get_n_eqns()
	    obj_approx%eqns%nterms = cal_parser%get_n_terms()
	    obj_approx%eqns%nbodymax = nbodymax 

	    allocate(obj_approx%eqns%residual(obj_approx%eqns%neqns),source = 1.d+10) ! initialise residuals to something large
	    allocate(obj_approx%eqns%jacobian(obj_approx%eqns%neqns,obj_approx%eqns%neqns),source = 0.d+0) ! initialise jacobian to zero
	    allocate(obj_approx%eqns%lhs(obj_approx%eqns%neqns),source = 0)
	    allocate(obj_approx%eqns%rhs(obj_approx%eqns%neqns),source = 0) ! null info for left and right hand side

	    allocate(obj_approx%eqns%correlation(obj_approx%eqns%nterms,obj_approx%eqns%nbodymax),source = 0)
	    allocate(obj_approx%eqns%corrlvalue(obj_approx%eqns%nterms),source = 0.d0) ! all correlation values initialised to zero
	    allocate(obj_approx%eqns%corrlnbody(obj_approx%eqns%nterms),source = 0)

            if (ncorc.eq.2) then !SB and SBWI
                call sigma1_sigma2_2NN(obj_approx)
                call sigma1_sigma8_2NN(obj_approx)  
            elseif (ncorc.eq.4) then  !SB plus 1NN 		
                call sigma1_sigma2_2NN(obj_approx)
                call sigma1_sigma8_2NN(obj_approx)
                call sigma1sigma2_sigma2sigma3_2NN(obj_approx)
                call sigma1sigma2_sigma2sigma8_2NN(obj_approx)
            elseif (ncorc.eq.6) then  !SB plus 1NN and ghost 2NN
                call sigma1_sigma2_2NN(obj_approx)
                call sigma1_sigma8_2NN(obj_approx)
                call sigma1sigma2_sigma2sigma3_2NN(obj_approx)
                call sigma1sigma2_sigma2sigma8_2NN(obj_approx)
                call sigma1sigma8_sigma2sigma4_2NN(obj_approx)
                call sigma1sigma8_sigma8sigma9_2NN(obj_approx)      
            endif
    
        endif
        
        return

end subroutine approx_initialise_2NN

subroutine sigma1_sigma2_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx
    
        obj_approx%eqns%corrlnbody(1:2) = 1
        obj_approx%eqns%correlation(1,1) = 1
        obj_approx%eqns%correlation(2,1) = 2 ! the terms of the 1st equation sigma(1) = sigma(2)
        obj_approx%eqns%lhs(1) = 1
        obj_approx%eqns%rhs(1) = 2 ! and finally we encode the equation

        return          
 end subroutine sigma1_sigma2_2NN

subroutine sigma1_sigma8_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx

        obj_approx%eqns%corrlnbody(3) = 1
        obj_approx%eqns%correlation(3,1) = 8 ! need one more term for the 2nd equation sigma(1) = sigma(8)
        obj_approx%eqns%lhs(2) = 1
        obj_approx%eqns%rhs(2) = 3 ! and finally we encode the equation   

        return          
end subroutine sigma1_sigma8_2NN

subroutine sigma1sigma2_sigma2sigma3_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx

        obj_approx%eqns%corrlnbody(4:5) = 2
        obj_approx%eqns%correlation(4,1:2) = (/1,2/)
        obj_approx%eqns%correlation(5,1:2) = (/2,3/) ! the terms of the 3rd equation sigma(1)*sigma(2) = sigma(2)*sigma(3)
        obj_approx%eqns%lhs(3) = 4
        obj_approx%eqns%rhs(3) = 5 ! and finally we encode the equation
    
        return          
end subroutine sigma1sigma2_sigma2sigma3_2NN

subroutine sigma1sigma2_sigma2sigma8_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx

        obj_approx%eqns%corrlnbody(6) = 2
        obj_approx%eqns%correlation(6,1:2) = (/2,8/) ! the terms of the 4th equation sigma(1)*sigma(2) = sigma(2)*sigma(8)
        obj_approx%eqns%lhs(4) = 4
        obj_approx%eqns%rhs(4) = 6 ! and finally we encode the equation
    
        return          
end subroutine sigma1sigma2_sigma2sigma8_2NN

subroutine sigma1sigma8_sigma2sigma4_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx

        obj_approx%eqns%corrlnbody(7:8) = 2
        obj_approx%eqns%correlation(7,1:2) = (/1,8/)
        obj_approx%eqns%correlation(8,1:2) = (/2,4/) ! the terms of the 5th equation sigma(1)*sigma(8) = sigma(2)*sigma(4)
        obj_approx%eqns%lhs(5) = 7
        obj_approx%eqns%rhs(5) = 8 ! and finally we encode the equation
    
        return          
end subroutine sigma1sigma8_sigma2sigma4_2NN

subroutine sigma1sigma8_sigma8sigma9_2NN(obj_approx)
        use meso_approx
        implicit none
        class (approximation) :: obj_approx

        obj_approx%eqns%corrlnbody(9) = 2
        obj_approx%eqns%correlation(9,1:2) = (/8,9/) ! the terms of the 5th equation sigma(1)*sigma(8) = sigma(8)*sigma(9)
        obj_approx%eqns%lhs(6) = 7
        obj_approx%eqns%rhs(6) = 9 ! and finally we encode the equation
    
        return          
end subroutine sigma1sigma8_sigma8sigma9_2NN 
  
  
  
