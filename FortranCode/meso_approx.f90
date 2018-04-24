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
		integer, allocatable, dimension(:,:) :: stateprods ! array storing terms that appear in the Hamiltonian, like sigma(1)*sigma(2) etc, for each state.
		integer, allocatable, dimension(:,:) :: sumstateprods ! array storing sums of the above, which make it faster to calculate dHamiltonian/dcorcpar, for each state.
        ! Both the above are stored as members of the class for computational efficiency, but for large approximations it can be memory intensive!
        
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
        !    stateprods(1:2**nsites,1:nterms)
        !    sumstateprods(1:2**nsites,ncorc)
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
		real(8), allocatable, dimension(:) :: corrlvalue ! value for each correlation term
		integer, allocatable, dimension(:) :: lhs  ! points to array "correlation", encoding the left hand side of each equation
		integer, allocatable, dimension(:) :: rhs  ! points to array "correlation", encoding the right hand side of each equation
		real(8), allocatable, dimension(:) :: residual  ! residual of each equation
		real(8), allocatable, dimension(:,:) :: jacobian  ! jacobian of each equation
        integer, allocatable, dimension(:,:) ::  stateprods(:,:) ! array storing terms that appear in the correlations, like sigma(1)*sigma(2) etc, for each state.
        ! These are saved for computational efficiency, but for large approximations it can be memory intensive!
        
		! The sizes of the above are intended to be:
		!    correlation(1:nterms,1:nbodymax)
		!    corrlnbody(1:nterms)
		!    corrlvalue(1:nterms)
		!    lhs(1:neqns)
		!    rhs(1:neqns)
		!    residual(1:neqns)
        !    stateprods(1:2**nsites,1:nterms)
		! The number of equations must be the same as the number of correction terms in the Hamiltonian (neqns = ncorc)
		! to be able to solve the approximation		
    end type 

    type, public :: approximation
        character*10 :: approxname
        integer nsites
        real(8) temp ! temperature
        real(8) mu ! chemical potential
        real(8) partfcn ! partition function
        integer, allocatable, dimension(:,:) :: allstates ! intended size (2^nsites,nsites) - we focus on 1-adsorbate approximations
        integer, allocatable, dimension(:) :: nparticles ! intended size (2^nsites,nsites) - we focus on 1-adsorbate approximations
        real(8), allocatable, dimension(:) :: allenergs ! intended size (2^nsites) 
        type (hamiltonian) :: hamilt
        type (equation) :: eqns
    contains
        procedure :: init => approx_initialise
        procedure :: prnt => approx_print
        procedure :: populate_allstates
        procedure :: calc_energ => calculate_energies
        procedure :: calc_resid => calculate_residuals
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
    
        allocate(this%nparticles(2**this%nsites),source=sum(this%allstates,2))
        
        return
    
    end subroutine populate_allstates

    
    subroutine calculate_energies(this)
    
        implicit none
        class (approximation) :: this
        integer i, j
        
        ! Preparatory steps: allocate and precompute stateprods and sumstateprods
        if (.not. allocated(this%hamilt%stateprods)) then
            allocate(this%hamilt%stateprods(2**this%nsites,this%hamilt%nterms),source=1)
            do i = 1,this%hamilt%nterms
                do j = 1,this%hamilt%internbody(i)
                    this%hamilt%stateprods(:,i) = &
                        this%hamilt%stateprods(:,i)*this%allstates(:,this%hamilt%interaction(i,j))
                enddo
            enddo
        endif
        if (.not. allocated(this%hamilt%sumstateprods)) then
            allocate(this%hamilt%sumstateprods(2**this%nsites,this%hamilt%ncorc),source=0)
            do i = 1,2**this%nsites
                do j = 1,this%hamilt%ncorc
                    this%hamilt%sumstateprods(i,j) = sum(this%hamilt%stateprods(i,:),MASK=this%hamilt%corcterms(1:)==j)
                enddo
            enddo
        endif

        this%allenergs = this%hamilt%H0
        do i = 1,this%hamilt%nterms                        
            this%allenergs = this%allenergs + & 
                (this%hamilt%origpars(this%hamilt%origterms(i)) + &
                 this%hamilt%corcpars(this%hamilt%corcterms(i)))*this%hamilt%stateprods(:,i)
        enddo
        
        return
        
    end subroutine calculate_energies
    
    
    subroutine calculate_residuals(this)
    
        use global_constants

        implicit none
        
        class (approximation) :: this
        integer i, j, k
        real(8) lhsderivativeterm, rhsderivativeterm
        
        ! Preparatory steps: allocate and precompute stateprods
        if (.not. allocated(this%eqns%stateprods)) then
            allocate(this%eqns%stateprods(2**this%nsites,this%eqns%nterms),source=1)
            do i = 1,this%eqns%nterms
                do j = 1,this%eqns%corrlnbody(i)
                    this%eqns%stateprods(:,i) = &
                        this%eqns%stateprods(:,i)*this%allstates(:,this%eqns%correlation(i,j))
                enddo
            enddo
        endif

        call this%calc_energ() ! Note that we calculate the energies here, so if a program unit is calling the correlations subroutine,
        ! it would be unnecessary (and a waste of time) to compute the energies in the calling program unit
        
        this%eqns%corrlvalue = 0.d0

        this%partfcn = sum(exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp)))
        do i = 1,this%eqns%nterms
            ! The following two expressions should give the same results (numerical accuracy issues excluded)
            ! In the Matlab code the first expression is used, i.e. not the actual correlation function, but the 
            ! non-normalised partial sum that corresponds to that correlation
            this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp)))
            ! this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp)))/this%partfcn
        enddo        
    
        this%eqns%residual = 0.d0
        do i = 1,this%eqns%neqns
            ! Again, two options. In Matlab we have used the version of the equations with the logarithms
            this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i))) - log(this%eqns%corrlvalue(this%eqns%rhs(i)))
            ! this%eqns%residual(i) = this%eqns%corrlvalue(this%eqns%lhs(i)) - this%eqns%corrlvalue(this%eqns%rhs(i))
        enddo
        
        this%eqns%jacobian = 0.d0
        do i = 1,this%eqns%neqns
            do j = 1,this%eqns%neqns
                lhsderivativeterm = sum(this%eqns%stateprods(:,this%eqns%lhs(i)) &
                    *exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp)) &
                    *this%hamilt%sumstateprods(:,j))    ! some part of this expression has been computed before - could be saved for efficiency
                rhsderivativeterm = sum(this%eqns%stateprods(:,this%eqns%rhs(i)) &
                    *exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp)) &
                    *this%hamilt%sumstateprods(:,j))
                this%eqns%jacobian(i,j) = 1.d0/this%eqns%corrlvalue(this%eqns%lhs(i))*lhsderivativeterm &
                    - 1.d0/this%eqns%corrlvalue(this%eqns%rhs(i))*rhsderivativeterm
            enddo
        enddo
        this%eqns%jacobian = -1/(kboltz*this%temp)*this%eqns%jacobian

        return
    
    end subroutine calculate_residuals
    
    
	subroutine approx_initialise(this,approx_name)

	    use global_constants
        implicit none
        class (approximation) :: this
        character(*) approx_name
	
	    this%mu = mu0
	    this%temp = temp

        select case (approx_name)
            
            case ('BPEC')
                call approx_initialise_bpec(this)

            case ('K2NNC2')
                call approx_initialise_k2nnc2(this)
                
            case default
                write(*,*) 'Unknown approximation',approx_name
                
        end select
            
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
        
     !   write(*,*) ''
     !   write(*,*) 'All states:'
	    !do i = 1,2**this%nsites
		   ! write(*,'(' // int2str(this%nsites) // 'I3)') (this%allstates(i,j), j=this%nsites,1,-1)
     !   enddo
        
	    write(*,*) 'EQUATIONS'
	    write(*,*) 'number of distinct terms:   ',this%eqns%nterms
	    write(*,*) 'number of equations:        ',this%eqns%neqns
	    write(*,*) 'list of equations:   ',this%eqns%nterms
        do i = 1,this%eqns%neqns
            write(*,'(I3,a,5x,' // trim(int2str(this%eqns%corrlnbody(this%eqns%lhs(i)))) // '("s(",I3,")"),' &
                             // ', " = ", ' // trim(int2str(this%eqns%corrlnbody(this%eqns%rhs(i)))) // '("s(",I3,")"))' ) &
                i,')', &
                (this%eqns%correlation(this%eqns%lhs(i),j),j=1,this%eqns%corrlnbody(this%eqns%lhs(i))), &
                (this%eqns%correlation(this%eqns%rhs(i),j),j=1,this%eqns%corrlnbody(this%eqns%rhs(i)))
            continue
        enddo
        	
        return
	
	end subroutine approx_print

end module