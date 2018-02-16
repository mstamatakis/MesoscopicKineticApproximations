    module meso_approx

    implicit none
    private

    integer, public :: nsites
    integer, public :: npar
            
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
    do i=1,2**nsites
     s=1
     do j=1,nsites
      if(appr%hamilt%orig%term(i,j).ne.0) then
       s=s*state(appr%hamilt%orig%term(i,j))
      end if
     end do
     energy=energy+s*appr%hamilt%orig%value(i)
    end do
    ! Hamiltonian correction block
    do i=1,2**nsites
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
    do i=1,2**nsites
     call confs(state,i)
     partition=partition+exp((chemp*sum(state)-appr%ener(appr,state))/(kb*temp))
    end do
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
    do i=1,2**nsites
     call confs(state,i)
     s=1
     do k=1,m
      s=s*state(v(k))
     end do
     correlation=correlation+s*exp((chemp*sum(state)-appr%ener(appr,state))/(kb*temp))
    end do
    end function
 
    subroutine approx_initialisation(appr)
    use commons
    class(approximation_type), intent(inout) :: appr
    integer i, j
 
    if(appr%approx.eq.'bp') then
     nsites=7
     npar=1
     allocate(appr%hamilt%orig%term(2**nsites,nsites))
     allocate(appr%hamilt%orig%value(2**nsites))
     allocate(appr%hamilt%corr%term(2**nsites,nsites))
     allocate(appr%hamilt%corr%intmap(2**nsites))
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

    if(appr%approx.eq.'bpe') then
     nsites=7
     npar=1
     allocate(appr%hamilt%orig%term(2**nsites,nsites))
     allocate(appr%hamilt%orig%value(2**nsites))
     allocate(appr%hamilt%corr%term(2**nsites,nsites))
     allocate(appr%hamilt%corr%intmap(2**nsites))
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

    if(appr%approx.eq.'bpec') then
     nsites=7
     npar=2
     allocate(appr%hamilt%orig%term(2**nsites,nsites))
     allocate(appr%hamilt%orig%value(2**nsites))
     allocate(appr%hamilt%corr%term(2**nsites,nsites))
     allocate(appr%hamilt%corr%intmap(2**nsites))
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
    return
    end subroutine approx_initialisation
    end module 
