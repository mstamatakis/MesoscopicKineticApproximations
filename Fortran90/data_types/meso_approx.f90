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
         procedure :: init => approx_initialisation
        end type

    contains
    
    subroutine approx_initialisation(appr)
    use commons
    class(approximation_type), intent(inout) :: appr
    integer i

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
     appr%hamilt%corr%value=0.d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0

     do i=1,7
      appr%hamilt%orig%term(i,1)=1
     end do
     appr%hamilt%orig%term(8,1)=1
     appr%hamilt%orig%term(8,2)=2
     appr%hamilt%orig%term(9,1)=1
     appr%hamilt%orig%term(9,2)=3
     appr%hamilt%orig%term(10,1)=1
     appr%hamilt%orig%term(10,2)=4
     appr%hamilt%orig%term(11,1)=1
     appr%hamilt%orig%term(11,2)=5
     appr%hamilt%orig%term(12,1)=1
     appr%hamilt%orig%term(12,2)=6
     appr%hamilt%orig%term(13,1)=1
     appr%hamilt%orig%term(13,2)=7
     appr%hamilt%orig%term(14,1)=2
     appr%hamilt%orig%term(14,2)=3
     appr%hamilt%orig%term(15,1)=3
     appr%hamilt%orig%term(15,2)=4
     appr%hamilt%orig%term(16,1)=4
     appr%hamilt%orig%term(16,2)=5
     appr%hamilt%orig%term(17,1)=5
     appr%hamilt%orig%term(17,2)=6
     appr%hamilt%orig%term(18,1)=6
     appr%hamilt%orig%term(18,2)=7
     appr%hamilt%orig%term(19,1)=2
     appr%hamilt%orig%term(19,2)=7

     appr%hamilt%orig%value(1)=hads
     appr%hamilt%orig%value(2)=hads
     appr%hamilt%orig%value(3)=hads
     appr%hamilt%orig%value(4)=hads
     appr%hamilt%orig%value(5)=hads
     appr%hamilt%orig%value(6)=hads
     appr%hamilt%orig%value(7)=hads
     appr%hamilt%orig%value(8)=jint
     appr%hamilt%orig%value(9)=jint
     appr%hamilt%orig%value(11)=jint
     appr%hamilt%orig%value(12)=jint
     appr%hamilt%orig%value(13)=jint
     appr%hamilt%orig%value(14)=jint
     appr%hamilt%orig%value(15)=jint
     appr%hamilt%orig%value(16)=jint
     appr%hamilt%orig%value(17)=jint
     appr%hamilt%orig%value(18)=jint
     appr%hamilt%orig%value(19)=jint

     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
     end do
     appr%hamilt%corr%term(7,1)=2
     appr%hamilt%corr%term(7,2)=3
     appr%hamilt%corr%term(8,1)=3
     appr%hamilt%corr%term(8,2)=4
     appr%hamilt%corr%term(9,1)=4
     appr%hamilt%corr%term(9,2)=5
     appr%hamilt%corr%term(10,1)=5
     appr%hamilt%corr%term(10,2)=6
     appr%hamilt%corr%term(11,1)=6
     appr%hamilt%corr%term(11,2)=7
     appr%hamilt%corr%term(12,1)=2
     appr%hamilt%corr%term(12,2)=7

     appr%hamilt%corr%intmap(1)=1
     appr%hamilt%corr%intmap(2)=1
     appr%hamilt%corr%intmap(3)=1
     appr%hamilt%corr%intmap(4)=1
     appr%hamilt%corr%intmap(5)=1
     appr%hamilt%corr%intmap(6)=1
     appr%hamilt%corr%intmap(7)=2
     appr%hamilt%corr%intmap(8)=2
     appr%hamilt%corr%intmap(9)=2
     appr%hamilt%corr%intmap(10)=2
     appr%hamilt%corr%intmap(11)=2
     appr%hamilt%corr%intmap(12)=2
     appr%hamilt%corr%intmap(13)=2

     appr%hamilt%corr%value(1)=0.d0
     appr%hamilt%corr%value(2)=0.d0

     appr%eqn%lhs=0        
     appr%eqn%rhs=0        
     appr%eqn%lhs(1,1)=1
     appr%eqn%lhs(2,1)=1
     appr%eqn%lhs(2,2)=2
     appr%eqn%rhs(1,1)=2
     appr%eqn%rhs(2,1)=2
     appr%eqn%rhs(2,2)=3
    end if

    if(appr%approx.eq.'k2nnc1') then
     nsites=13
     npar=4
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
     appr%hamilt%corr%value=0.d0
     appr%eqn%lhs=0
     appr%eqn%rhs=0

     do i=1,7
      appr%hamilt%orig%term(i,1)=1
     end do
     appr%hamilt%orig%term(8,1)=1
     appr%hamilt%orig%term(8,2)=2
     appr%hamilt%orig%term(9,1)=1
     appr%hamilt%orig%term(9,2)=3
     appr%hamilt%orig%term(10,1)=1
     appr%hamilt%orig%term(10,2)=4
     appr%hamilt%orig%term(11,1)=1
     appr%hamilt%orig%term(11,2)=5
     appr%hamilt%orig%term(12,1)=1
     appr%hamilt%orig%term(12,2)=6
     appr%hamilt%orig%term(13,1)=1
     appr%hamilt%orig%term(13,2)=7
     appr%hamilt%orig%term(14,1)=2
     appr%hamilt%orig%term(14,2)=3
     appr%hamilt%orig%term(15,1)=3
     appr%hamilt%orig%term(15,2)=4
     appr%hamilt%orig%term(16,1)=4
     appr%hamilt%orig%term(16,2)=5
     appr%hamilt%orig%term(17,1)=5
     appr%hamilt%orig%term(17,2)=6
     appr%hamilt%orig%term(18,1)=6
     appr%hamilt%orig%term(18,2)=7
     appr%hamilt%orig%term(19,1)=2
     appr%hamilt%orig%term(19,2)=7

     appr%hamilt%orig%value(1)=hads
     appr%hamilt%orig%value(2)=hads
     appr%hamilt%orig%value(3)=hads
     appr%hamilt%orig%value(4)=hads
     appr%hamilt%orig%value(5)=hads
     appr%hamilt%orig%value(6)=hads
     appr%hamilt%orig%value(7)=hads
     appr%hamilt%orig%value(8)=jint
     appr%hamilt%orig%value(9)=jint
     appr%hamilt%orig%value(11)=jint
     appr%hamilt%orig%value(12)=jint
     appr%hamilt%orig%value(13)=jint
     appr%hamilt%orig%value(14)=jint
     appr%hamilt%orig%value(15)=jint
     appr%hamilt%orig%value(16)=jint
     appr%hamilt%orig%value(17)=jint
     appr%hamilt%orig%value(18)=jint
     appr%hamilt%orig%value(19)=jint

     do i=2,7
      appr%hamilt%corr%term(i-1,1)=i
     end do
     appr%hamilt%corr%term(7,1)=2
     appr%hamilt%corr%term(7,2)=3
     appr%hamilt%corr%term(8,1)=3
     appr%hamilt%corr%term(8,2)=4
     appr%hamilt%corr%term(9,1)=4
     appr%hamilt%corr%term(9,2)=5
     appr%hamilt%corr%term(10,1)=5
     appr%hamilt%corr%term(10,2)=6
     appr%hamilt%corr%term(11,1)=6
     appr%hamilt%corr%term(11,2)=7
     appr%hamilt%corr%term(12,1)=2
     appr%hamilt%corr%term(12,2)=7

     appr%hamilt%corr%intmap(1)=1
     appr%hamilt%corr%intmap(2)=1
     appr%hamilt%corr%intmap(3)=1
     appr%hamilt%corr%intmap(4)=1
     appr%hamilt%corr%intmap(5)=1
     appr%hamilt%corr%intmap(6)=1
     appr%hamilt%corr%intmap(7)=2
     appr%hamilt%corr%intmap(8)=2
     appr%hamilt%corr%intmap(9)=2
     appr%hamilt%corr%intmap(10)=2
     appr%hamilt%corr%intmap(11)=2
     appr%hamilt%corr%intmap(12)=2
     appr%hamilt%corr%intmap(13)=2

     appr%hamilt%corr%value(1)=0.d0
     appr%hamilt%corr%value(2)=0.d0

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
