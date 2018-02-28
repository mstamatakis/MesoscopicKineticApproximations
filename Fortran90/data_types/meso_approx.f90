    module meso_approx

    implicit none
    private

    integer, public :: nsites
    integer, public :: npar
    integer, public :: nrows
            
    type :: original
     integer, allocatable, dimension(:,:) :: term
     real*8, allocatable, dimension(:) :: value 
    end type 
    type :: correction
     integer, allocatable, dimension(:,:) :: term
     integer, allocatable, dimension(:) :: intmap
     real*8, allocatable, dimension(:) :: value 
    end type 
    type :: equation
     integer, allocatable, dimension(:,:) :: lhs
     integer, allocatable, dimension(:,:) :: rhs 
    end type 
    type :: hamiltonian
     type (original) :: orig
     type (correction) :: corr
    end type
    
    type, public :: approximation_type
     character*10 :: approx
     type (hamiltonian) :: hamilt
     type (equation) :: eqn
    contains
     procedure, nopass :: ener => energy
     procedure :: part => partition
     procedure, nopass :: corfun => correlation
     procedure :: init => approx_initialisation
    end type

    contains
    real*8 function energy(appr,state)
    use commons
    implicit none
    integer i, j, s
    integer, dimension(nsites), intent(in) :: state
    class(approximation_type), intent(in) :: appr 
    
    energy=0.d0
    ! Hamiltonian original block
    !$OMP PARALLEL
    !$OMP DO
    do i=1,nrows
     s=1
     do j=1,nsites
      if(appr%hamilt%orig%term(i,j).ne.0) then
       s=s*state(appr%hamilt%orig%term(i,j))
      end if
     end do
     energy=energy+s*appr%hamilt%orig%value(i)
    end do
    !$OMP END DO
    ! Hamiltonian correction block
    !$OMP DO
    do i=1,nrows
     s=1
     do j=1,nsites
      if(appr%hamilt%corr%term(i,j).ne.0) then
       s=s*state(appr%hamilt%corr%term(i,j))
      end if
     end do
     if(appr%hamilt%corr%intmap(i).ne.0) then
      energy=energy+s*appr%hamilt%corr%value(appr%hamilt%corr%intmap(i))
     end if
    end do        
    !$OMP END DO
    !$OMP END PARALLEL
    end function
 
    real*8 function partition(appr)
    use commons
    implicit none
    integer i
    integer, dimension(nsites) :: state
    real*8 energy
    class (approximation_type), intent(in) :: appr
    
    state=0
    partition=0.d0
    !$OMP PARALLEL
    !$OMP DO
    do i=1,2**nsites
     call confs(state,i)
     partition=partition+exp((chemp*sum(state)-appr%ener(appr,state)-h0)/(kb*temp))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    end function
 
    real*8 function correlation(v,m,appr)
    use commons
    implicit none
    integer i, k, m, s
    integer, dimension(m), intent(in) :: v
    integer, dimension(nsites) :: state
    real*8 energy
    class (approximation_type), intent(in) :: appr
    
    state=0
    correlation=0.d0
    !$OMP PARALLEL
    !$OMP DO
    do i=1,2**nsites
     call confs(state,i)
     s=1
     do k=1,m
      s=s*state(v(k))
     end do
     correlation=correlation+s*exp((chemp*sum(state)-appr%ener(appr,state)-h0)/(kb*temp))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    end function
 
    subroutine approx_initialisation(appr)
    use commons
    class(approximation_type), intent(inout) :: appr
    integer i, j
 
    if(appr%approx.eq.'BP') then
     nsites=7
     npar=1
     nrows=13
     allocate(appr%hamilt%orig%term(nrows,nsites))
     allocate(appr%hamilt%orig%value(nrows))
     allocate(appr%hamilt%corr%term(nrows,nsites))
     allocate(appr%hamilt%corr%intmap(nrows))
     allocate(appr%hamilt%corr%value(npar))
     allocate(appr%eqn%lhs(npar,nsites))
     allocate(appr%eqn%rhs(npar,nsites))
     appr%hamilt%orig%term=0
     appr%hamilt%orig%value=0.d0
     appr%hamilt%corr%term=0
     appr%hamilt%corr%intmap=0
     appr%hamilt%corr%value=0.1d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0
     do i=1,7
      appr%hamilt%orig%term(i,1)=i
     end do
     appr%hamilt%orig%value(1)=hads
     j=7
     do i=1,6
      appr%hamilt%orig%term(j+i,1)=1
      appr%hamilt%orig%term(j+i,2)=i+1
      appr%hamilt%orig%value(i+1)=hads
     end do
     do i=8,13
      appr%hamilt%orig%value(i)=jint
     end do
     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=1
     end do
     appr%eqn%lhs(1,1)=1
     appr%eqn%rhs(1,1)=2
    end if

    if(appr%approx.eq.'BPE') then
     nsites=7
     npar=1
     nrows=19
     allocate(appr%hamilt%orig%term(nrows,nsites))
     allocate(appr%hamilt%orig%value(nrows))
     allocate(appr%hamilt%corr%term(nrows,nsites))
     allocate(appr%hamilt%corr%intmap(nrows))
     allocate(appr%hamilt%corr%value(npar))
     allocate(appr%eqn%lhs(npar,nsites))
     allocate(appr%eqn%rhs(npar,nsites))
     appr%hamilt%orig%term=0
     appr%hamilt%orig%value=0.d0
     appr%hamilt%corr%term=0
     appr%hamilt%corr%intmap=0
     appr%hamilt%corr%value=0.1d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0
     do i=1,7
      appr%hamilt%orig%term(i,1)=i
     end do
     appr%hamilt%orig%value(1)=hads
     j=7
     do i=1,6
      appr%hamilt%orig%term(j+i,1)=1
      appr%hamilt%orig%term(j+i,2)=i+1
      appr%hamilt%orig%value(i+1)=hads
     end do
     j=13
     do i=1,5
      appr%hamilt%orig%term(j+i,1)=i+1
      appr%hamilt%orig%term(j+i,2)=i+2
     end do
     appr%hamilt%orig%term(19,1)=2
     appr%hamilt%orig%term(19,2)=7
     do i=8,19
      appr%hamilt%orig%value(i)=jint
     end do
     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=1
     end do
     appr%eqn%lhs(1,1)=1
     appr%eqn%rhs(1,1)=2
    end if

    if(appr%approx.eq.'BPEC') then
     nsites=7
     npar=2
     nrows=19
     allocate(appr%hamilt%orig%term(nrows,nsites))
     allocate(appr%hamilt%orig%value(nrows))
     allocate(appr%hamilt%corr%term(nrows,nsites))
     allocate(appr%hamilt%corr%intmap(nrows))
     allocate(appr%hamilt%corr%value(npar))
     allocate(appr%eqn%lhs(npar,nsites))
     allocate(appr%eqn%rhs(npar,nsites))
     appr%hamilt%orig%term=0
     appr%hamilt%orig%value=0.d0
     appr%hamilt%corr%term=0
     appr%hamilt%corr%intmap=0
     appr%hamilt%corr%value=0.1d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0
     do i=1,7
      appr%hamilt%orig%term(i,1)=i
     end do
     j=7
     appr%hamilt%orig%value(1)=hads
     do i=1,6
      appr%hamilt%orig%term(j+i,1)=1
      appr%hamilt%orig%term(j+i,2)=i+1
      appr%hamilt%orig%value(i+1)=hads
     end do
     j=13
     do i=1,5
      appr%hamilt%orig%term(j+i,1)=i+1
      appr%hamilt%orig%term(j+i,2)=i+2
     end do 
     appr%hamilt%orig%term(19,1)=2
     appr%hamilt%orig%term(19,2)=7
     do i=8,19
      appr%hamilt%orig%value(i)=jint
     end do
     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
     end do
     j=6
     do i=1,5
      appr%hamilt%corr%term(j+i,1)=i+1
      appr%hamilt%corr%term(j+i,2)=i+2
     end do
     appr%hamilt%corr%term(12,1)=2
     appr%hamilt%corr%term(12,2)=7
 
     do i=1,6
      appr%hamilt%corr%intmap(i)=1
     end do
     do i=7,12
      appr%hamilt%corr%intmap(i)=2
     end do
     appr%hamilt%corr%value(1)=0.1d0
     appr%hamilt%corr%value(2)=0.1d0
     appr%eqn%lhs=0        
     appr%eqn%rhs=0        
     appr%eqn%lhs(1,1)=1
     appr%eqn%lhs(2,1)=1
     appr%eqn%lhs(2,2)=2
     appr%eqn%rhs(1,1)=2
     appr%eqn%rhs(2,1)=2
     appr%eqn%rhs(2,2)=3
    end if

    if(appr%approx.eq.'K2NNC1') then
     nsites=13
     npar=4
     nrows=37
     allocate(appr%hamilt%orig%term(nrows,nsites))
     allocate(appr%hamilt%orig%value(nrows))
     allocate(appr%hamilt%corr%term(nrows,nsites))
     allocate(appr%hamilt%corr%intmap(nrows))
     allocate(appr%hamilt%corr%value(npar))
     allocate(appr%eqn%lhs(npar,nsites))
     allocate(appr%eqn%rhs(npar,nsites))
     appr%hamilt%orig%term=0
     appr%hamilt%orig%value=0.d0
     appr%hamilt%corr%term=0
     appr%hamilt%corr%intmap=0
     appr%hamilt%corr%value=0.1d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0
     ! Original terms
     do i=1,13
      appr%hamilt%orig%term(i,1)=i
      appr%hamilt%orig%value(i)=hads
     end do
     do i=2,7
      appr%hamilt%orig%term(12+i,1)=1
      appr%hamilt%orig%term(12+i,2)=i
      appr%hamilt%orig%value(12+i)=jint
     end do 
     do i=2,6
      appr%hamilt%orig%term(18+i,1)=i
      appr%hamilt%orig%term(18+i,2)=i+1
      appr%hamilt%orig%value(18+i)=jint
     end do
      appr%hamilt%orig%term(25,1)=2
      appr%hamilt%orig%term(25,2)=7
      appr%hamilt%orig%value(25)=jint
     do i=8,12
      appr%hamilt%orig%term(18+i,1)=i
      appr%hamilt%orig%term(18+i,2)=i-6
      appr%hamilt%orig%value(18+i)=jint
     end do
     do i=8,12
      appr%hamilt%orig%term(23+i,1)=i
      appr%hamilt%orig%term(23+i,2)=i-6+1
      appr%hamilt%orig%value(23+i)=jint
     end do
     appr%hamilt%orig%term(36,1)=13
     appr%hamilt%orig%term(36,2)=2
     appr%hamilt%orig%term(37,1)=13
     appr%hamilt%orig%term(37,2)=7
     appr%hamilt%orig%value(36)=jint
     appr%hamilt%orig%value(37)=jint
    ! Correction terms
     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=1
     end do
     do i=8,13
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=2
     end do
     do i=2,6
      appr%hamilt%corr%term(11+i,1)=i   
      appr%hamilt%corr%term(11+i,2)=i+1   
      appr%hamilt%corr%intmap(11+i)=3  
     end do
     appr%hamilt%corr%term(18,1)=2
     appr%hamilt%corr%term(18,2)=7
     appr%hamilt%corr%intmap(18)=3
     do i=8,12
      appr%hamilt%corr%term(11+i,1)=i
      appr%hamilt%corr%term(11+i,2)=i-6
      appr%hamilt%corr%intmap(11+i)=4
     end do
     do i=8,12
      appr%hamilt%corr%term(16+i,1)=i
      appr%hamilt%corr%term(16+i,2)=i-6+1
      appr%hamilt%corr%intmap(16+i)=4
     end do
     appr%hamilt%corr%term(29,1)=13
     appr%hamilt%corr%term(29,2)=2
     appr%hamilt%corr%term(30,1)=13
     appr%hamilt%corr%term(30,2)=7
     appr%hamilt%corr%intmap(29)=4
     appr%hamilt%corr%intmap(30)=4
    ! Equations
     appr%eqn%lhs(1,1)=1
     appr%eqn%lhs(2,1)=1
     appr%eqn%lhs(3,1)=1
     appr%eqn%lhs(3,2)=2
     appr%eqn%lhs(4,1)=1
     appr%eqn%lhs(4,2)=2
     appr%eqn%rhs(1,1)=2
     appr%eqn%rhs(2,1)=8
     appr%eqn%rhs(3,1)=2     
     appr%eqn%rhs(3,2)=3     
     appr%eqn%rhs(4,1)=2     
     appr%eqn%rhs(4,2)=8     
    end if

    if(appr%approx.eq.'K2NNC2') then
     nsites=13
     npar=6
     nrows=42
     allocate(appr%hamilt%orig%term(nrows,nsites))
     allocate(appr%hamilt%orig%value(nrows))
     allocate(appr%hamilt%corr%term(nrows,nsites))
     allocate(appr%hamilt%corr%intmap(nrows))
     allocate(appr%hamilt%corr%value(npar))
     allocate(appr%eqn%lhs(npar,nsites))
     allocate(appr%eqn%rhs(npar,nsites))
     appr%hamilt%orig%term=0
     appr%hamilt%orig%value=0.d0
     appr%hamilt%corr%term=0
     appr%hamilt%corr%intmap=0
     appr%hamilt%corr%value=0.1d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0
     ! Original terms
     do i=1,13
      appr%hamilt%orig%term(i,1)=i
      appr%hamilt%orig%value(i)=hads
     end do
     do i=2,7
      appr%hamilt%orig%term(12+i,1)=1
      appr%hamilt%orig%term(12+i,2)=i
      appr%hamilt%orig%value(12+i)=jint
     end do 
     do i=2,6
      appr%hamilt%orig%term(18+i,1)=i
      appr%hamilt%orig%term(18+i,2)=i+1
      appr%hamilt%orig%value(18+i)=jint
     end do
      appr%hamilt%orig%term(25,1)=2
      appr%hamilt%orig%term(25,2)=7
      appr%hamilt%orig%value(25)=jint
     do i=8,12
      appr%hamilt%orig%term(18+i,1)=i
      appr%hamilt%orig%term(18+i,2)=i-6
      appr%hamilt%orig%value(18+i)=jint
     end do
     do i=8,12
      appr%hamilt%orig%term(23+i,1)=i
      appr%hamilt%orig%term(23+i,2)=i-6+1
      appr%hamilt%orig%value(23+i)=jint
     end do
     appr%hamilt%orig%term(36,1)=13
     appr%hamilt%orig%term(36,2)=2
     appr%hamilt%orig%term(37,1)=13
     appr%hamilt%orig%term(37,2)=7
     appr%hamilt%orig%value(36)=jint
     appr%hamilt%orig%value(37)=jint
    ! Correction terms
     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=1
     end do
     do i=8,13
      appr%hamilt%corr%term(i-1,1)=i
      appr%hamilt%corr%intmap(i-1)=2
     end do
     do i=2,6
      appr%hamilt%corr%term(11+i,1)=i   
      appr%hamilt%corr%term(11+i,2)=i+1   
      appr%hamilt%corr%intmap(11+i)=3  
     end do
     appr%hamilt%corr%term(18,1)=2
     appr%hamilt%corr%term(18,2)=7
     appr%hamilt%corr%intmap(18)=3
     do i=8,12
      appr%hamilt%corr%term(11+i,1)=i
      appr%hamilt%corr%term(11+i,2)=i-6
      appr%hamilt%corr%intmap(11+i)=4
     end do
     do i=8,12
      appr%hamilt%corr%term(16+i,1)=i
      appr%hamilt%corr%term(16+i,2)=i-6+1
      appr%hamilt%corr%intmap(16+i)=4
     end do
     appr%hamilt%corr%term(29,1)=13
     appr%hamilt%corr%term(29,2)=2
     appr%hamilt%corr%term(30,1)=13
     appr%hamilt%corr%term(30,2)=7
     appr%hamilt%corr%intmap(29)=4
     appr%hamilt%corr%intmap(30)=4
     appr%hamilt%corr%term(31,1)=2
     appr%hamilt%corr%term(31,2)=4
     appr%hamilt%corr%term(32,1)=2
     appr%hamilt%corr%term(32,2)=6
     appr%hamilt%corr%term(33,1)=3
     appr%hamilt%corr%term(33,2)=5
     appr%hamilt%corr%term(34,1)=3
     appr%hamilt%corr%term(34,2)=7
     appr%hamilt%corr%term(35,1)=4
     appr%hamilt%corr%term(35,2)=6
     appr%hamilt%corr%term(36,1)=5
     appr%hamilt%corr%term(36,2)=7
     do i=1,6
      appr%hamilt%corr%intmap(30+i)=5
     end do
     do i=8,12
      appr%hamilt%corr%term(29+i,1)=i
      appr%hamilt%corr%term(29+i,2)=i+1
      appr%hamilt%corr%intmap(29+i)=6
     end do
     appr%hamilt%corr%term(42,1)=8
     appr%hamilt%corr%term(42,2)=13
     appr%hamilt%corr%intmap(42)=6
    ! Equations
     appr%eqn%lhs(1,1)=1
     appr%eqn%lhs(2,1)=1
     appr%eqn%lhs(3,1)=1
     appr%eqn%lhs(3,2)=2
     appr%eqn%lhs(4,1)=1
     appr%eqn%lhs(4,2)=2
     appr%eqn%lhs(5,1)=1
     appr%eqn%lhs(5,2)=8
     appr%eqn%lhs(6,1)=1
     appr%eqn%lhs(6,2)=8
     appr%eqn%rhs(1,1)=2
     appr%eqn%rhs(2,1)=8
     appr%eqn%rhs(3,1)=2     
     appr%eqn%rhs(3,2)=3     
     appr%eqn%rhs(4,1)=2     
     appr%eqn%rhs(4,2)=8     
     appr%eqn%rhs(5,1)=2     
     appr%eqn%rhs(5,2)=4     
     appr%eqn%rhs(6,1)=8     
     appr%eqn%rhs(6,2)=9     
    end if
    return
    end subroutine approx_initialisation
    end module 
