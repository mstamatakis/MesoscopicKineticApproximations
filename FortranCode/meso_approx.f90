module meso_approx
  use calculation_setup_module
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
		integer, allocatable, dimension(:) :: stateprods ! array storing terms that appear in the Hamiltonian, like sigma(1)*sigma(2) etc, for each state.
		integer, allocatable, dimension(:) :: sumstateprodsorig ! array storing sums of the above, which make it faster to calculate dHamiltonian/dcorcpar, for each state.
        integer, allocatable, dimension(:) :: sumstateprodscorc ! array storing sums of the above, which make it faster to calculate dHamiltonian/dcorcpar, for each state.
        ! All three of the above are stored as members of the class for computational efficiency, but for large approximations it can be memory intensive!
        
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
		procedure :: approx_initialise_1NN
		procedure :: approx_initialise_2NN
		procedure :: approx_initialise_3NN
        procedure :: allocate_memory
    end type

	! Interfaces for approximation-specific initialization subroutines
	
    interface
        subroutine approx_initialise_1NN(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface
    
    interface
        subroutine approx_initialise_2NN(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface
    
    interface
        subroutine approx_initialise_3NN(obj_approx)
            import approximation
            class(approximation) :: obj_approx
        end subroutine
    end interface

    contains

    subroutine populate_allstates(this)
    
        implicit none
        class (approximation) :: this
        integer i, j, dec, allstatescount
        
        allstatescount = 2**this%nsites

        allocate(this%allstates(allstatescount,this%nsites),source=0)
        
        !$OMP PARALLEL DO PRIVATE(i, j, dec) 
        do i = 0,allstatescount-1
            dec = i
            do j = 1,this%nsites
                if (mod(dec,2)==0) then
                    this%allstates(i+1,j) = 0
                else
                    this%allstates(i+1,j) = 1
                end if
                dec = dec * 0.5
                if (dec == 0) then
                    exit  
                endif
            enddo
        enddo
        !$OMP END PARALLEL DO
    
        allocate(this%nparticles(allstatescount))
        !$OMP PARALLEL DO PRIVATE(i, dec) 
        do i = 1,allstatescount
            this%nparticles(i) = sum(this%allstates(i,:))
        enddo
        !$OMP END PARALLEL DO
        
        return
    
    end subroutine populate_allstates

    
    subroutine allocate_memory(this)
        implicit none
        class (approximation) :: this
        integer i, j, k, allstatescount, tmp_id1, tmp_id2, tmp_id3
        real(8) tmp_val
        allstatescount = 2**this%nsites

        ! Preparatory steps: allocate and precompute stateprods
        if (.not. allocated(this%eqns%stateprods)) then
            allocate(this%eqns%stateprods(allstatescount * this%eqns%nterms), source = 1)
            do i = 1,this%eqns%nterms
                tmp_id2 = (i - 1) * allstatescount
                do j = 1,this%eqns%corrlnbody(i)
                    tmp_id1 = this%eqns%correlation(i, j)
                    !$OMP PARALLEL PRIVATE(k, tmp_id3)
                    !$OMP DO SIMD
                    do k = 1, allstatescount
                        tmp_id3 = k + tmp_id2
                        this%eqns%stateprods(tmp_id3) = &
                            this%eqns%stateprods(tmp_id3) * this%allstates(k, tmp_id1)
                    enddo
                    !$OMP END DO SIMD
                    !$OMP END PARALLEL
                enddo
            enddo
        endif
        if (.not. allocated(this%expenergies)) then
            allocate(this%expenergies(allstatescount),source=0.d0)
        endif

        ! Preparatory steps: allocate and precompute stateprods, sumstateprodsorig, sumstateprodscorc
        if (.not. (allocated(this%hamilt%sumstateprodsorig) .and. allocated(this%hamilt%sumstateprodscorc))) then 
            if (.not. allocated(this%hamilt%stateprods)) then
                allocate(this%hamilt%stateprods(allstatescount * this%hamilt%nterms),source=1)
                do i = 1,this%hamilt%nterms
                    tmp_id3 = (i - 1) * allstatescount
                    do j = 1,this%hamilt%internbody(i)
                        tmp_id2 = this%hamilt%interaction(i, j)    
                        !$OMP PARALLEL PRIVATE(k, tmp_id1) 
                        !$OMP DO SIMD
                        do k = 1, allstatescount
                            tmp_id1 = k + tmp_id3
                            this%hamilt%stateprods(tmp_id1) = &
                                        this%hamilt%stateprods(tmp_id1) * this%allstates(k, tmp_id2)
                        enddo
                        !$OMP END DO SIMD
                        !$OMP END PARALLEL
                    enddo
                enddo
            endif
            if (.not. allocated(this%hamilt%sumstateprodsorig)) then
                allocate(this%hamilt%sumstateprodsorig(allstatescount * this%hamilt%norig),source=0)
                !$OMP PARALLEL PRIVATE(i, j, k, tmp_id1) 
                !$OMP DO SIMD
                do i = 1,allstatescount
                    do j = 1,this%hamilt%norig
                        tmp_id1 = 0
                        do k = 1, this%hamilt%nterms
                            if (this%hamilt%origterms(k) .eq. j) &
                                tmp_id1 = tmp_id1 + this%hamilt%stateprods(i + (k - 1) * allstatescount)
                        enddo
                        this%hamilt%sumstateprodsorig(i + (j - 1) * allstatescount) = tmp_id1
                    enddo
                enddo
                !$OMP END DO SIMD
                !$OMP END PARALLEL
            endif
            if (.not. allocated(this%hamilt%sumstateprodscorc)) then
                allocate(this%hamilt%sumstateprodscorc(allstatescount * this%hamilt%ncorc), source=0)
                !$OMP PARALLEL PRIVATE(i, j, k, tmp_id1)
                !$OMP DO SIMD
                do i = 1,allstatescount
                    do j = 1,this%hamilt%ncorc
                        tmp_id1 = 0
                        do k = 1, this%hamilt%nterms
                            if (this%hamilt%corcterms(k) .eq. j) &
                                tmp_id1 = tmp_id1 + this%hamilt%stateprods(i + (k - 1) * allstatescount)
                        enddo
                        this%hamilt%sumstateprodscorc(i + (j - 1) * allstatescount) = tmp_id1
                    enddo
                enddo
                !$OMP END DO SIMD
                !$OMP END PARALLEL
            endif
            deallocate(this%hamilt%stateprods) ! free up some memory; this array not needed anymore, unless we use the first method to calculate energies
        endif
    end subroutine allocate_memory


    subroutine calculate_energies(this)
    
        implicit none
        class (approximation) :: this
        integer i, j, k, allstatescount, tmp_int, tmp_int2, tmp_int3
        real(8) tmp_val
        allstatescount = 2**this%nsites
        
        if ((this%hamilt%norig == 2) .and. (mod(this%hamilt%ncorc,2) == 0)) then
            !$OMP PARALLEL PRIVATE(i, j, tmp_val)
            !$OMP DO SIMD
            do j = 1,allstatescount
                tmp_val = 0.d0
                tmp_int = j - allstatescount
                tmp_val = tmp_val &
                        + this%hamilt%sumstateprodsorig(tmp_int + allstatescount) * this%hamilt%origpars(1) &
                        + this%hamilt%sumstateprodsorig(tmp_int + 2 * allstatescount) * this%hamilt%origpars(2)
                do i = 1,this%hamilt%ncorc
                    tmp_val = tmp_val &
                            + this%hamilt%sumstateprodscorc(tmp_int + i * allstatescount) * this%hamilt%corcpars(i) !&
                            ! + this%hamilt%sumstateprodscorc(tmp_int + (i + 1) * allstatescount) * this%hamilt%corcpars(i + 1)
                enddo
                this%allenergs(j) = this%hamilt%H0 + tmp_val
            enddo
            !$OMP END DO SIMD
            !$OMP END PARALLEL
        else
            !$OMP PARALLEL DO PRIVATE(i, j, tmp_val) 
            do j = 1,allstatescount
                tmp_val = 0.d0
                tmp_int = j - allstatescount
                do i = 1,this%hamilt%norig
                    tmp_val = tmp_val + this%hamilt%sumstateprodsorig(tmp_int + i * allstatescount) * this%hamilt%origpars(i)
                enddo
                do i = 1,this%hamilt%ncorc
                    tmp_val = tmp_val + this%hamilt%sumstateprodscorc(tmp_int + i * allstatescount) * this%hamilt%corcpars(i)
                enddo
                this%allenergs(j) = this%hamilt%H0 + tmp_val
            enddo
            !$OMP END PARALLEL DO
        end if
        
        return
        
    end subroutine calculate_energies
    
    
    subroutine calculate_residuals(this,calcjac)
    
      use constants_module

        implicit none
        
        class (approximation) :: this
        integer i, j, k, kk, allstatescount
        ! real(8) lhsderivativeterm, rhsderivativeterm
        real(8) derivativeterm_lhs, derivativeterm_rhs
        logical, intent(in), optional :: calcjac
        logical :: calculatejacobian
        integer id, tmp_id1, tmp_id2, tmp_id3
        real(8) tmp_var1, tmp_var2, tmp_var3
        real(8) kboltz_times_temp

        allstatescount = 2**this%nsites
        kboltz_times_temp = kboltz*this%temp
        
        if (present(calcjac)) then
            calculatejacobian = calcjac
        else
            calculatejacobian = .true.
        endif

        call this%calc_energ() ! Note that we calculate the energies here, so if a program unit is calling the correlations subroutine,
        ! it would be unnecessary (and a waste of time) to compute the energies in the calling program unit
        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO SIMD
        do i = 1, allstatescount
            this%expenergies(i) = exp(-(this%allenergs(i) - this%mu * this%nparticles(i)) / kboltz_times_temp)
        enddo
        !$OMP END DO SIMD
        !$OMP END PARALLEL
        
        this%eqns%corrlvalue = 0.d0

        tmp_var1 = 0.d0
        !$OMP PARALLEL
        !$OMP DO SIMD REDUCTION(+:tmp_var1) 
        do i = 1,allstatescount
            tmp_var1 = tmp_var1 + this%expenergies(i)
        enddo
        !$OMP END DO SIMD
        !$OMP END PARALLEL
        this%partfcn = tmp_var1

        do i = 1,this%eqns%nterms
            ! The following two expressions should give the same results (numerical accuracy issues excluded)
            ! In the Matlab code the first expression is used, i.e. not the actual correlation function, but the 
            ! non-normalised partial sum that corresponds to that correlation
!            this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*this%expenergies)
            tmp_var1 = 0.d0
            !$OMP PARALLEL PRIVATE(k)
            !$OMP DO SIMD REDUCTION(+:tmp_var1)
            do k = 1,allstatescount,1
                tmp_var1 = tmp_var1 &
                                    + this%eqns%stateprods(k + (i - 1) * allstatescount) &
                                    * this%expenergies(k)
            enddo
            !$OMP END DO SIMD
            !$OMP END PARALLEL

            this%eqns%corrlvalue(i) = tmp_var1
            ! this%eqns%corrlvalue(i) = sum(this%eqns%stateprods(:,i)*this%expenergies)/this%partfcn
        enddo
    
        this%eqns%residual = 0.d0
        !$OMP PARALLEL DO PRIVATE(i)
        do i = 1,this%eqns%neqns
            ! Again, two options. In Matlab we have used the version of the equations with the logarithms
!            this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i))) - log(this%eqns%corrlvalue(this%eqns%rhs(i)))
            ! this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i))) - log(this%eqns%corrlvalue(this%eqns%rhs(i)))
            this%eqns%residual(i) = log(this%eqns%corrlvalue(this%eqns%lhs(i)) / this%eqns%corrlvalue(this%eqns%rhs(i)))    ! One log is better than two
            ! this%eqns%residual(i) = this%eqns%corrlvalue(this%eqns%lhs(i)) - this%eqns%corrlvalue(this%eqns%rhs(i))
        enddo
        !$OMP END PARALLEL DO
        
        this%eqns%jacobian = 0.d0
        
        if (.not.calculatejacobian) return

        do i = 1,this%eqns%neqns
            tmp_var1 = 1.d0 / this%eqns%corrlvalue(this%eqns%lhs(i))
            tmp_var2 = 1.d0 / this%eqns%corrlvalue(this%eqns%rhs(i))
            tmp_id2 = (this%eqns%lhs(i) - 1) * allstatescount
            tmp_id3 = (this%eqns%rhs(i) - 1) * allstatescount
            do j = 1,this%eqns%neqns
                derivativeterm_lhs = 0.d0
                derivativeterm_rhs = 0.d0
                tmp_id1 = (j - 1) * allstatescount

                !$OMP PARALLEL PRIVATE(k)
                !$OMP DO SIMD REDUCTION(+:derivativeterm_lhs, derivativeterm_rhs)
                do k = 1, allstatescount
                    tmp_var3 = this%hamilt%sumstateprodscorc(k + tmp_id1) &
                                * this%expenergies(k)
                    derivativeterm_lhs = derivativeterm_lhs &
                                        + this%eqns%stateprods(k + tmp_id2) &
                                        * tmp_var3
                    derivativeterm_rhs = derivativeterm_rhs &
                                        + this%eqns%stateprods(k + tmp_id3) &
                                        * tmp_var3
                enddo
                !$OMP END DO SIMD
                !$OMP END PARALLEL

                this%eqns%jacobian(i, j) = derivativeterm_lhs * tmp_var1 &
                                            - derivativeterm_rhs * tmp_var2

                this%eqns%jacobian(i, j) = -this%eqns%jacobian(i, j) / kboltz_times_temp
            enddo
        enddo
        
        return
    
    end subroutine calculate_residuals
    
    
	subroutine approx_initialise(this,cluster_name)

	    use global_constants
        implicit none
        class (approximation) :: this
        character(*) cluster_name
	
        this%mu = cal_parser% get_muIni()
       
	this%temp = cal_parser%get_temp()

        select case (cluster_name)
            
            case ('1NN')
               !call this%approx_initialise_bpec()
                call this%approx_initialise_1NN()
            case ('2NN')
               !call this%approx_initialise_k2nnc2()
                call this%approx_initialise_2NN()
                
            case ('3NN')
                !call this%approx_initialise_k3nnc2()
                call this%approx_initialise_3NN()
            case default
                write(*,*) 'Unknown approximation',cluster_name
                
        end select
        
        call this%allocate_memory()

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
