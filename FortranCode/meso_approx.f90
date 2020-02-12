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
		integer, allocatable, dimension(:,:) :: sumstateprodsorig ! array storing sums of the above, which make it faster to calculate dHamiltonian/dcorcpar, for each state.
		integer, allocatable, dimension(:) :: sumstateprodscorc ! array storing sums of the above, which make it faster to calculate dHamiltonian/dcorcpar, for each state.
        ! All three of the above are stored as members of the class for computational efficiency, but for large approximations it can be memory intensive!
        integer, allocatable, dimension(:) :: stateprods_new ! array storing terms that appear in the Hamiltonian, like sigma(1)*sigma(2) etc, for each state.

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
        !    sumstateprodsorig(1:2**nsites,norig)
        !    sumstateprodscorc(1:2**nsites,ncorc)
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
        integer, allocatable, dimension(:) ::  stateprods(:) ! array storing terms that appear in the correlations, like sigma(1)*sigma(2) etc, for each state.
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
        integer, allocatable, dimension(:) :: nparticles ! intended size (2^nsites) - we focus on 1-adsorbate approximations
        real(8), allocatable, dimension(:) :: allenergs ! intended size (2^nsites) 
        real(8), allocatable, dimension(:) :: expenergies ! intended size (2^nsites) - stores expressions exp(-(E-mu*N)/(kB*T)) for efficiency
        type (hamiltonian) :: hamilt
        type (equation) :: eqns
    contains
        procedure :: init => approx_initialise
        procedure :: prnt => approx_print
        procedure :: populate_allstates
        procedure :: calc_energ => calculate_energies
        procedure :: calc_resid => calculate_residuals
		procedure :: approx_initialise_bpec
		procedure :: approx_initialise_k2nnc2
		procedure :: approx_initialise_k3nnc2
    end type

	! Interfaces for approximation-specific initialization subroutines
	
    interface
        subroutine approx_initialise_bpec(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface
    
    interface
        subroutine approx_initialise_k2nnc2(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface
    
    interface
        subroutine approx_initialise_k3nnc2(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface

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
    
        allocate(this%nparticles(2**this%nsites))
        do i = 1,2**this%nsites
            this%nparticles(i) = sum(this%allstates(i,:))
        enddo
        
        return
    
    end subroutine populate_allstates

    
    subroutine calculate_energies(this)
    
        implicit none
        class (approximation) :: this
        integer i, j, k, upper_range
        real(4) t1, t2, t3, t4
        real(8) tmp_val
        logical is_even_ncorc

        upper_range = 2**this%nsites

        ! Preparatory steps: allocate and precompute stateprods, sumstateprodsorig, sumstateprodscorc
        if (.not. (allocated(this%hamilt%sumstateprodsorig) .and. allocated(this%hamilt%sumstateprodscorc))) then 
            call cpu_time(t1) ! function for calculating elapsed CPU time
            if (.not. allocated(this%hamilt%stateprods)) then
                allocate(this%hamilt%stateprods(upper_range,this%hamilt%nterms),source=1)
                
                !!!!!!
                ! allocate(this%hamilt%stateprods_new(upper_range * this%hamilt%nterms), source = 1)
                !!!!!!

                do i = 1,this%hamilt%nterms
                    do j = 1,this%hamilt%internbody(i)
                        this%hamilt%stateprods(:,i) = &
                            this%hamilt%stateprods(:,i)*this%allstates(:,this%hamilt%interaction(i,j))
                        ! do k = 1, upper_range
                        !     this%hamilt%stateprods_new(i + k * this%hamilt%nterms) = this%hamilt%stateprods(k,i)
                        ! enddo
                    enddo
                enddo
            endif
            if (.not. allocated(this%hamilt%sumstateprodsorig)) then
                allocate(this%hamilt%sumstateprodsorig(upper_range,this%hamilt%norig),source=0)
                do i = 1,upper_range
                    do j = 1,this%hamilt%norig
                        this%hamilt%sumstateprodsorig(i,j) = sum(this%hamilt%stateprods(i,:),MASK=this%hamilt%origterms(1:)==j)
                    enddo
                enddo
            endif
            if (.not. allocated(this%hamilt%sumstateprodscorc)) then
                ! allocate(this%hamilt%sumstateprodscorc(upper_range,this%hamilt%ncorc),source=0)
                ! do i = 1,upper_range
                !     do j = 1,this%hamilt%ncorc
                !         this%hamilt%sumstateprodscorc(i,j) = sum(this%hamilt%stateprods(i,:),MASK=this%hamilt%corcterms(1:)==j)
                !     enddo
                ! enddo
                allocate(this%hamilt%sumstateprodscorc(upper_range * this%hamilt%ncorc), source=0)
                do i = 1,upper_range
                    do j = 1,this%hamilt%ncorc
                        this%hamilt%sumstateprodscorc(i + (j - 1) * upper_range) = &
                            sum(this%hamilt%stateprods(i,:),MASK=this%hamilt%corcterms(1:)==j)
                    enddo
                enddo
            endif
            deallocate(this%hamilt%stateprods) ! free up some memory; this array not needed anymore, unless we use the first method to calculate energies
            call cpu_time(t2) ! function for calculating elapsed CPU time
            write(*,*) 'time (energy 1)',t2-t1
        endif
        
        !this%allenergs = this%hamilt%H0
        !do i = 1,this%hamilt%nterms                        
        !    do j = 1,2**this%nsites
        !        this%allenergs(j) = this%allenergs(j) + & 
        !            (this%hamilt%origpars(this%hamilt%origterms(i)) + &
        !             this%hamilt%corcpars(this%hamilt%corcterms(i)))*this%hamilt%stateprods(j,i)
        !    enddo
        !enddo
        
        ! ! call cpu_time(t1) ! function for calculating elapsed CPU time
        ! this%allenergs = this%hamilt%H0
        ! do i = 1,this%hamilt%norig
        !     tmp_val = this%hamilt%origpars(i)
        !     do j = 1,2**this%nsites
        !         this%allenergs(j) = this%allenergs(j) + & 
        !             tmp_val*this%hamilt%sumstateprodsorig(j,i)
        !     enddo
        ! enddo
        ! ! call cpu_time(t2) ! function for calculating elapsed CPU time
        ! ! write(*,*) 'Elapsed CPU-time (residual 2)',t2-t1, this%hamilt%norig, 2**this%nsites

        ! ! call cpu_time(t3) ! function for calculating elapsed CPU time
        ! do i = 1,this%hamilt%ncorc
        !     tmp_val = this%hamilt%corcpars(i)
        !     do j = 1,2**this%nsites
        !         this%allenergs(j) = this%allenergs(j) + & 
        !             tmp_val * this%hamilt%sumstateprodscorc(j,i)
        !     enddo
        ! enddo
        ! ! call cpu_time(t4) ! function for calculating elapsed CPU time
        ! ! write(*,*) 'Elapsed CPU-time (residual 3)',t4-t3, this%hamilt%ncorc, 2**this%nsites

        ! ! write(*,*) 'Elapsed CPU-time (residual full)',t4-t3 + t2-t1


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! call cpu_time(t1) ! function for calculating elapsed CPU time
        this%allenergs = this%hamilt%H0
        ! upper_range = 2**this%nsites
        do i = 1,this%hamilt%norig
            tmp_val = this%hamilt%origpars(i)
            do j = 1,upper_range
                this%allenergs(j) = this%allenergs(j) + & 
                    tmp_val*this%hamilt%sumstateprodsorig(j,i)
            enddo
        enddo
        ! call cpu_time(t2) ! function for calculating elapsed CPU time
        ! write(*,*) 'Elapsed CPU-time (residual 2)',t2-t1, this%hamilt%norig, 2**this%nsites

        ! call cpu_time(t3) ! function for calculating elapsed CPU time
        do i = 1,this%hamilt%ncorc
            tmp_val = this%hamilt%corcpars(i)
            do j = 1,upper_range
                this%allenergs(j) = this%allenergs(j) + & 
                    tmp_val * this%hamilt%sumstateprodscorc(j + (i - 1) * upper_range)
                    ! tmp_val * this%hamilt%sumstateprodscorc(j,i)
            enddo
        enddo
        ! call cpu_time(t4) ! function for calculating elapsed CPU time
        ! write(*,*) 'Elapsed CPU-time (residual 3)',t4-t3, this%hamilt%ncorc, 2**this%nsites

        ! write(*,*) 'Elapsed CPU-time (residual full)',t4-t3 + t2-t1
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! call cpu_time(t3) ! function for calculating elapsed CPU time
        ! this%allenergs = this%hamilt%H0
        ! upper_range = 2**this%nsites
        ! if (mod(this%hamilt%ncorc,2) .eq. 0) then
        !     is_even_ncorc = .TRUE.
        ! else
        !     is_even_ncorc = .FALSE.
        ! endif

        ! do j = 1,upper_range
        !     tmp_val = 0.
        !     ! first loop
        !     ! if (this%hamilt%norig == 2) then
        !     !     tmp_val = tmp_val + & 
        !     !     this%hamilt%origpars(1)*this%hamilt%sumstateprodsorig(j,1) + &
        !     !     this%hamilt%origpars(2)*this%hamilt%sumstateprodsorig(j,2)
        !     ! else
        !         do i = 1,this%hamilt%norig
        !             tmp_val = tmp_val + & 
        !                     this%hamilt%origpars(i)*this%hamilt%sumstateprodsorig(j,i)
        !         enddo
        !     ! endif

        !     ! second loop
        !     ! if (is_even_ncorc) then
        !     !     do i = 1,this%hamilt%ncorc,2
        !     !         tmp_val = tmp_val + & 
        !     !                 this%hamilt%corcpars(i) * this%hamilt%sumstateprodscorc(j,i) + &
        !     !                 this%hamilt%corcpars(i + 1) * this%hamilt%sumstateprodscorc(j,i + 1)
        !     !     enddo
        !     ! else
        !         do i = 1,this%hamilt%ncorc
        !             tmp_val = tmp_val + & 
        !                     this%hamilt%corcpars(i) * this%hamilt%sumstateprodscorc(j,i)
        !         enddo
        !     ! endif
        !     this%allenergs(j) = this%allenergs(j) + tmp_val
        ! enddo
        ! call cpu_time(t4) ! function for calculating elapsed CPU time
        ! write(*,*) 'Elapsed CPU-time (residual full)',t4-t3, this%hamilt%ncorc, 2**this%nsites

        return
        
    end subroutine calculate_energies
    
    
    subroutine calculate_residuals(this,calcjac)
    
        use global_constants
        use omp_lib

        implicit none
        
        class (approximation) :: this
        integer i, j, k, upper_range
        real(8) lhsderivativeterm, rhsderivativeterm
        logical, intent(in), optional :: calcjac
        logical :: calculatejacobian
        real(4) t1, t2, t3, t4, t5, t6, t7, t8
        integer lhs_i, rhs_i, id, tmp_id
        real(8) tmp_var1, tmp_var2

        upper_range = 2**this%nsites

        if (present(calcjac)) then
            calculatejacobian = calcjac
        else
            calculatejacobian = .true.
        endif
        
        ! Preparatory steps: allocate and precompute stateprods
        if (.not. allocated(this%eqns%stateprods)) then
            allocate(this%eqns%stateprods(upper_range * this%eqns%nterms), source = 1)

            do i = 1,this%eqns%nterms
                do j = 1,this%eqns%corrlnbody(i)
                    tmp_id = this%eqns%correlation(i, j)
                    do k = 1, upper_range
                        id = k + (i - 1) * upper_range
                        this%eqns%stateprods(id) = &
                                    this%eqns%stateprods(id) * this%allstates(k, tmp_id)
                    enddo
                enddo
            enddo
        endif
        if (.not. allocated(this%expenergies)) then
            allocate(this%expenergies(upper_range),source=0.d0)
        endif

        call cpu_time(t1) !!!
        call this%calc_energ() ! Note that we calculate the energies here, so if a program unit is calling the correlations subroutine,
        ! it would be unnecessary (and a waste of time) to compute the energies in the calling program unit

        ! this%expenergies = exp(-(this%allenergs-this%mu*this%nparticles)/(kboltz*this%temp))
        ! tmp_var1 = 1.d0 / (kboltz*this%temp)
        do i = 1, upper_range
            this%expenergies(i) = exp(-(this%allenergs(i) - this%mu * this%nparticles(i)) / (kboltz*this%temp))
        enddo
        call cpu_time(t2) !!!
        write(*,*) 'time (residual 1)',t2-t1

        this%eqns%corrlvalue = 0.d0

        call cpu_time(t3) !!!
        ! this%partfcn = sum(this%expenergies)
        this%partfcn = 0.d0
        do i = 1,upper_range
            this%partfcn = this%partfcn + this%expenergies(i)
        enddo

        do i = 1,this%eqns%nterms
            ! The following two expressions should give the same results (numerical accuracy issues excluded)
            ! In the Matlab code the first expression is used, i.e. not the actual correlation function, but the 
            ! non-normalised partial sum that corresponds to that correlation
            ! this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*this%expenergies)
            tmp_var1 = 0.d0
            do k = 1,upper_range,1
                ! tmp_var1 = tmp_var1 &
                !                     + this%eqns%stateprods(k,i) * this%expenergies(k)

                tmp_var1 = tmp_var1 &
                                    + this%eqns%stateprods(k + (i - 1) * upper_range) &
                                    * this%expenergies(k)
            enddo
            this%eqns%corrlvalue(i) = tmp_var1
            ! this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*this%expenergies)/this%partfcn
        enddo  
        call cpu_time(t4) !!!
        write(*,*) 'time (residual 2)',t4-t3, this%eqns%nterms
    
        call cpu_time(t5) !!!
        ! this%eqns%residual = 0.d0     ! no point to flush it, since it will be rewritten in the loop
        do i = 1,this%eqns%neqns
            ! Again, two options. In Matlab we have used the version of the equations with the logarithms
            
            ! mathematically, log(a/b) = log(a) - log(b). however, if implemented in this way the roundoff errors change the
            ! results. so, keep it as it was
            ! this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i)) / this%eqns%corrlvalue(this%eqns%rhs(i)))    ! One log is better than two
            this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i))) - log(this%eqns%corrlvalue(this%eqns%rhs(i)))    ! One log is better than two
            
            ! this%eqns%residual(i) = this%eqns%corrlvalue(this%eqns%lhs(i)) - this%eqns%corrlvalue(this%eqns%rhs(i))
        enddo
        call cpu_time(t6) !!!
        write(*,*) 'time (residual 3)',t6-t5, this%eqns%neqns

        this%eqns%jacobian = 0.d0
        
        if (.not.calculatejacobian) return
        
        call cpu_time(t7) !!!
        ! !$OMP PARALLEL
        do i = 1,this%eqns%neqns
            lhs_i = this%eqns%lhs(i)
            rhs_i = this%eqns%rhs(i)

            do j = 1,this%eqns%neqns
                lhsderivativeterm = 0.d0
                rhsderivativeterm = 0.d0
                tmp_var1 = 0.d0
                tmp_var2 = 0.d0
                ! !$OMP SIMD
                do k = 1,upper_range,2
                    tmp_var1 = this%expenergies(k    ) &
                            * this%hamilt%sumstateprodscorc(k + (j - 1) * upper_range)
                            ! * this%hamilt%sumstateprodscorc(k,     j)
                    tmp_var2 = this%expenergies(k + 1) &
                            * this%hamilt%sumstateprodscorc(k + (j - 1) * upper_range)
                            ! * this%hamilt%sumstateprodscorc(k + 1, j)
                    ! tmp_var1 = this%expenergies(k    ) * this%hamilt%sumstateprodscorc(k,     j)
                    ! tmp_var2 = this%expenergies(k + 1) * this%hamilt%sumstateprodscorc(k + 1, j)
                    lhsderivativeterm = lhsderivativeterm &
                                      + this%eqns%stateprods(k + (lhs_i - 1) * upper_range) * tmp_var1 &
                                      + this%eqns%stateprods(k + 1 + (lhs_i - 1) * upper_range) * tmp_var2
                    rhsderivativeterm = rhsderivativeterm &
                                      + this%eqns%stateprods(k + (rhs_i - 1) * upper_range) * tmp_var1 &
                                      + this%eqns%stateprods(k + 1 + (rhs_i - 1) * upper_range) * tmp_var2
                enddo
                ! !$OMP END SIMD
                this%eqns%jacobian(i,j) = &
                        1.d0/this%eqns%corrlvalue(lhs_i)*lhsderivativeterm &
                      - 1.d0/this%eqns%corrlvalue(rhs_i)*rhsderivativeterm
            enddo
        enddo
        ! !$OMP END PARALLEL
        this%eqns%jacobian = -1.d0/(kboltz*this%temp)*this%eqns%jacobian
        call cpu_time(t8) !!!
        write(*,*) 'time (residual 4)',t8-t7, this%eqns%neqns

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
                call this%approx_initialise_bpec()

            case ('K2NNC2')
                call this%approx_initialise_k2nnc2()
                
            case ('K3NNC2')
                call this%approx_initialise_k3nnc2()
                
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
            write(*,'(I3,a,5x,' // trim(int2str(this%eqns%corrlnbody(this%eqns%lhs(i)))) // '("s(",I3,")")' &
                             // ', " = ", ' // trim(int2str(this%eqns%corrlnbody(this%eqns%rhs(i)))) // '("s(",I3,")"))' ) &
                i,')', &
                (this%eqns%correlation(this%eqns%lhs(i),j),j=1,this%eqns%corrlnbody(this%eqns%lhs(i))), &
                (this%eqns%correlation(this%eqns%rhs(i),j),j=1,this%eqns%corrlnbody(this%eqns%rhs(i)))
            continue
        enddo
        	
        return
	
	end subroutine approx_print
	
end module