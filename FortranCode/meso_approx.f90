module meso_approx

	implicit none
	private

	type :: hamiltonian
		real H0 ! constant term in the Hamiltonian
		integer nterms     ! number of terms (original Hamiltonian terms or corrections)
		integer nbodymax   ! consider up to nbodymax-body terms (1-body, 2-body pairwise etc)
		integer, allocatable, dimension(:) :: internbody ! number of bodies for each interaction term (in array right below),    e.g. 1, 1, 2, ...
		integer, allocatable, dimension(:,:) :: interaction ! intended as a (nterms x nbodymax) array encoding interaction terms e.g. 1, 2, 1-2, ...
		integer, allocatable, dimension(:) :: origterms ! integers pointing to origpars (see below) to return original Hamiltonian parameters for each of the interaction terms in vector "interaction" above
		integer, allocatable, dimension(:) :: corcterms ! integers pointing to corcpars (see below) to return correction parameters for each of the interaction terms in vector "interaction" above
        
		integer norig      ! number of original parameters in the Hamiltonian, e.g. = 2 if only adsorption energy and 1NN interaction energy are used
		integer ncorc      ! number of correction terms in the Hamiltonian, e.g. = 1 if only one correction for adsorption energy is imposed (e.g. in BP, BPE)
		real(8), allocatable, dimension(:) :: origpars ! value of an original interaction term in the Hamiltonian
		real(8), allocatable, dimension(:) :: corcpars ! value of a correction interaction term in the Hamiltonian		
		! The sizes of the above are intended to be:
		!    internbody(1:nterms)
		!    interaction(1:nterms,1:nbodymax)
		!    origterms(1:nterms)
		!    corcterms(1:nterms)
		!    origpars(0:norig)
		!    corcpars(0:ncorc) 
        ! In origpars and corcpars the zero-th element is there on purpose and must always be equal to zero. This is because for example, not all terms of corcpars are 
        ! referenced in corcterms(1:nterms), as some interactions are not corrected. Thus corcterms may contain zeros.
		! The intention is to do a fast addition: origpars(origterms) + corcpars(corcterms) and then evaluate the Hamiltonian
	end type
	 
	type :: equation
		integer neqns     ! number of equations
		integer nterms    ! number of distinct terms
		integer nbodymax  ! consider up to nbodymax-body correlation terms (1-body, 2-body pairwise etc)
		integer, allocatable, dimension(:,:) :: correlation ! intended as a (nterms x nbodymax) array encoding correlation terms e.g. 1, 2, 1-2, ...
		integer, allocatable, dimension(:) :: corrlnbody ! number of bodies for each correlation term (in array right above),    e.g. 1, 1, 2, ...
		integer, allocatable, dimension(:) :: corrlvalue ! number of bodies for each correlation term (in array right above),    e.g. 1, 1, 2, ...
		integer, allocatable, dimension(:) :: lhs  ! points to array "correlation", encoding the left hand side of each equation
		integer, allocatable, dimension(:) :: rhs  ! points to array "correlation", encoding the right hand side of each equation
		real(8), allocatable, dimension(:) :: residual  ! residual of each equation
		! The sizes of the above are intended to be:
		!    correlation(1:nterms,1:nbodymax)
		!    corrlnbody(1:nterms)
		!    corrlvalue(1:nterms)
		!    lhs(1:neqns)
		!    rhs(1:neqns)
		!    residual(1:neqns)
		! The number of equations must be the same as the number of correction terms in the Hamiltonian (neqns = ncorc)
		! to be able to solve the approximation		
    end type 

    type, public :: approximation
        character*10 :: approxname
        integer nsites
        real(8) temp ! temperature
        real(8) mu ! chemical potential
        integer, allocatable, dimension(:,:) :: allstates ! intended size (2^nsites,nsites) - we focus on 1-adsorbate approximations
        real(8), allocatable, dimension(:) :: allenergs ! intended size (2^nsites) 
        type (hamiltonian) :: hamilt
        type (equation) :: eqn
    contains
        !procedure, nopass :: energy => calc_energy
        !procedure :: partfcn => calc_partitionfcn
        !procedure, nopass :: corfun => calc_correlation
        procedure :: init => approx_initialise
        procedure :: prnt => approx_print
        procedure :: populate_allstates
        procedure :: calc_energ => calculate_energies
        !procedure :: residuals => residuals
    end type
    
    contains

    subroutine populate_allstates(this)
    
        implicit none
        class (approximation) :: this
        integer i, j, dec, count
        
        allocate(this%allstates(2**this%nsites,this%nsites),source=0)
        
        do i = 0,2**this%nsites-1
            dec = i
            count = 0
            do j = 1,this%nsites
                if (mod(dec,2)==0) then
                    this%allstates(i+1,j) = 0
                else
                    this%allstates(i+1,j) = 1
                end if
                dec = dec/2
                count = count + 1
                if (dec == 0) then
                    exit  
                endif
            enddo
    
        enddo
    
        return
    
    end subroutine populate_allstates

    
    subroutine calculate_energies(this)
    
        implicit none
        class (approximation) :: this
        integer i, j, dec, count
        integer, allocatable, save :: temptermvec(:,:) ! #MSTAM: this is a temporary array storing terms like sigma(1)*sigma(2) etc, for each state.
        ! These are saved for computational efficiency, but for large approximations it can be memory intensive!
        
        if (.not. allocated(temptermvec)) then
            allocate(temptermvec(2**this%nsites,this%hamilt%nterms),source=1)
            do i = 1,this%hamilt%nterms
                do j = 1,this%hamilt%internbody(i)
                    temptermvec(:,i) = temptermvec(:,i)*this%allstates(:,this%hamilt%interaction(i,j))
                enddo
            enddo
        endif

        this%allenergs = this%hamilt%H0
        do i = 1,this%hamilt%nterms                        
            this%allenergs = this%allenergs + & 
                (this%hamilt%origpars(this%hamilt%origterms(i)) + &
                 this%hamilt%corcpars(this%hamilt%corcterms(i)))*temptermvec(:,i)
        enddo
        
        return
        
    end subroutine calculate_energies
    
	subroutine approx_initialise(this)

	    use global_constants
        implicit none
        class (approximation) :: this
        integer nsites, nterms, norig, ncorc, nbodymax, ntm, i, j, k, s1, s2
        real*8 tot1, tot2
	
	    this%mu = mu0
	    this%temp = temp

        
	    this%approxname = 'BPEC'
	    nsites = 7
	    this%nsites = nsites
	    call this%populate_allstates()
	    allocate(this%allenergs(2**nsites),source=0.d0)
	
	    ! Populate Hamiltonian data-structures
	    nterms = 7 + 6 + 6 ! 7 single body, 6 pairwise (center-1NN), and 6 more pairwise (between 1NNs)
	    norig = 2 ! adsorption energy and 1NN
	    ncorc = 2 ! corrections to adsorption energies of edge sites and 1NN interactions of edge sites
	    nbodymax = 2
	    this%hamilt%nterms = nterms
	    this%hamilt%norig = norig
	    this%hamilt%ncorc = ncorc
	    this%hamilt%nbodymax = 2
	
	    allocate(this%hamilt%interaction(nterms,nbodymax),source=0) ! encoding interaction terms
	    allocate(this%hamilt%internbody(nterms),source=0) ! number of bodies for each interaction term
	    this%hamilt%H0 = 0.d0 ! constant term in the Hamiltonian
	    allocate(this%hamilt%origpars(0:norig),source=(/0.d0,hads,Jint/)) ! value of an original interaction term in the Hamiltonian
	    allocate(this%hamilt%corcpars(0:ncorc),source=0.d0) ! value of a correction interaction term in the Hamiltonian
	    allocate(this%hamilt%origterms(nterms),source=0) ! will be used as a pointer to an original interaction term in the Hamiltonian
	    allocate(this%hamilt%corcterms(nterms),source=0) ! will be used as a pointer to a correction interaction term in the Hamiltonian
	
	    ntm = 0 ! counter of terms
	    do i = 1,nsites
		    ntm = ntm + 1
		    this%hamilt%interaction(ntm,1) = i ! single body terms
		    this%hamilt%internbody(ntm) = 1
	    enddo
	    do i = 2,7
		    ntm = ntm + 1
		    this%hamilt%interaction(ntm,1:2) = (/1,i/) ! central-1NN pairwise terms
		    this%hamilt%internbody(ntm) = 2
	    enddo
	    do i = 2,6
		    ntm = ntm + 1
		    this%hamilt%interaction(ntm,1:2) = (/i,i+1/) ! 1NN-1NN pairwise terms
		    this%hamilt%internbody(ntm) = 2
	    enddo
	    ntm = ntm + 1
	    this%hamilt%interaction(ntm,1:2) = (/7,2/) ! final 1NN-1NN pairwise term
	    this%hamilt%internbody(ntm) = 2

	    ! Integer arrays to be used as pointers to terms and corrections:
	    this%hamilt%origterms(1:nsites) = 1
	    this%hamilt%origterms(nsites+1:nterms) = 2		 
	    this%hamilt%corcterms(2:7) = 1
	    this%hamilt%corcterms(14:19) = 2
	 
        return
	
	end subroutine approx_initialise
	
	
	subroutine approx_print(this)

	    use parser_module
	    class (approximation) :: this
	    integer i, j

	    write(*,*) 'Approximation name:',this%approxname
	    write(*,*) 'Number of sites:',this%nsites
	    write(*,*) 'HAMILTONIAN'
	    write(*,*) 'original parameters:   ',(this%hamilt%origpars(j),j=1,this%hamilt%norig)
	    write(*,*) 'correction parameters: ',(this%hamilt%corcpars(j),j=1,this%hamilt%ncorc)
	    write(*,*) 'nterms:',this%hamilt%nterms
	    write(*,*) 'Original par pointers, values; correction par pointers, values; table of terms:'
	    do i = 1,this%hamilt%nterms
		    write(*,'(I3,a,2(I3,F16.8),' // trim(int2str(this%hamilt%internbody(i))) // 'I6)') &
                i,')', &
                this%hamilt%origterms(i), &
			    this%hamilt%origpars(this%hamilt%origterms(i)), &
                this%hamilt%corcterms(i), &
                this%hamilt%corcpars(this%hamilt%corcterms(i)), &
			    (this%hamilt%interaction(i,j),j=1,this%hamilt%internbody(i))
        enddo
        
        write(*,*) ''
        write(*,*) 'All states:'
	    do i = 1,2**this%nsites
		    write(*,'(' // int2str(this%nsites) // 'I3)') (this%allstates(i,j), j=this%nsites,1,-1)
        enddo
	
	    return
	
	end subroutine approx_print

end module