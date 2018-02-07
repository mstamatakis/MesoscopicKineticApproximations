        module meso_approx
        implicit none
        !character approx
        integer, parameter :: nsites=7
        integer, parameter :: npar=2 
        type :: original
         integer term(2**nsites,nsites)
         real*8 value(2**nsites) 
        end type 
        type :: correction
         integer term(2**nsites,nsites)
         integer intmap(2**nsites)
         real*8 value(npar) 
        end type 
        type :: eqn
         integer lhs(npar,nsites)
         integer rhs(npar,nsites) 
        end type 
        type :: hamiltonian
         type (original) :: list1 
         type (correction) :: list2
        end type
        type (hamiltonian) :: model
        type (eqn) :: syst_eqn
        real*8, parameter :: hads=-1.171955613174250d0 
        real*8, parameter :: jint=0.3d0
        real*8, parameter :: kb=8.6173303d-5
        real*8 chemp, temp
        integer i

        !odel%list1%term=0
        !odel%list1%value=0.d0
        !odel%list2%term=0
        !odel%list2%intmap=0
        !odel%list2%value=0.d0
        !o i=1,7
        !model%list1%term(i,1)=1
        !nd do
        !odel%list1%term(8,1)=1
        !odel%list1%term(8,2)=2
        !odel%list1%term(9,1)=1
        !odel%list1%term(9,2)=3
        !odel%list1%term(10,1)=1
        !odel%list1%term(10,2)=4
        !odel%list1%term(11,1)=1
        !odel%list1%term(11,2)=5
        !odel%list1%term(12,1)=1
        !odel%list1%term(12,2)=6
        !odel%list1%term(13,1)=1
        !odel%list1%term(13,2)=7
        !odel%list1%term(14,1)=2
        !odel%list1%term(14,2)=3
        !odel%list1%term(15,1)=3
        !odel%list1%term(15,2)=4
        !odel%list1%term(16,1)=4
        !odel%list1%term(16,2)=5
        !odel%list1%term(17,1)=5
        !odel%list1%term(17,2)=6
        !odel%list1%term(18,1)=6
        !odel%list1%term(18,2)=7
        !odel%list1%term(19,1)=2
        !odel%list1%term(19,2)=7

        !odel%list1%value(1)=hads
        !odel%list1%value(2)=hads
        !odel%list1%value(3)=hads
        !odel%list1%value(4)=hads
        !odel%list1%value(5)=hads
        !odel%list1%value(6)=hads
        !odel%list1%value(7)=hads
        !odel%list1%value(8)=jint
        !odel%list1%value(9)=jint
        !odel%list1%value(11)=jint
        !odel%list1%value(12)=jint
        !odel%list1%value(13)=jint
        !odel%list1%value(14)=jint
        !odel%list1%value(15)=jint
        !odel%list1%value(16)=jint
        !odel%list1%value(17)=jint
        !odel%list1%value(18)=jint
        !odel%list1%value(19)=jint

        !o i=2,7
        !model%list2%term(i-1,1)=i
        !nd do
        !odel%list2%term(7,1)=2
        !odel%list2%term(7,2)=3
        !odel%list2%term(8,1)=3
        !odel%list2%term(8,2)=4
        !odel%list2%term(9,1)=4
        !odel%list2%term(9,2)=5
        !odel%list2%term(10,1)=5
        !odel%list2%term(10,2)=6
        !odel%list2%term(11,1)=6
        !odel%list2%term(11,2)=7
        !odel%list2%term(12,1)=2
        !odel%list2%term(12,2)=7

        !odel%list2%intmap(1)=1
        !odel%list2%intmap(2)=1
        !odel%list2%intmap(3)=1
        !odel%list2%intmap(4)=1
        !odel%list2%intmap(5)=1
        !odel%list2%intmap(6)=1
        !odel%list2%intmap(7)=2
        !odel%list2%intmap(8)=2
        !odel%list2%intmap(9)=2
        !odel%list2%intmap(10)=2
        !odel%list2%intmap(11)=2
        !odel%list2%intmap(12)=2
        !odel%list2%intmap(13)=2

        !odel%list2%value(1)=0.d0
        !odel%list2%value(2)=0.d0

        !odel_eqn%lhs=0        
        !odel_eqn%rhs=0        
        !odel_eqn%lhs(1,1)=1
        !odel_eqn%lhs(2,1)=1
        !odel_eqn%lhs(2,2)=2
        !odel_eqn%rhs(1,1)=2
        !odel_eqn%rhs(2,1)=2
        !odel_eqn%rhs(2,2)=3

        !contains
        !real*8 function energy(state,approx,nsites)
        !implicit none
        !integer nsites, s, i, j
        !integer, dimension(nsites) :: state
        !type (hamiltonian) :: approx
        ! 
        !energy=0.d0
        !! Hamiltonian "original" block
        !do i=1,2**nsites
        ! s=1
        ! do j=1,nsites
        !  if(approx%list1%term(i,j).gt.0) then
        !   s=s*state(approx%list1%term(i,j))
        !  end if
        ! end do
        !! The variable "s" is used as a correlation variable. If zero, "value" is not included in the sum
        ! energy=energy+s*approx%list1%value(i)
        !end do          
        !! Hamiltonian "correction" block
        !do i=1,2**nsites
        ! s=1
        ! do j=1,nsites
        !  if(approx%list2%term(i,j).gt.0) then
        !   s=s*state(approx%list2%term(i,j))
        !  end if
        ! end do
        ! energy=energy+approx%list2%value(approx%list2%intmap(i))
        !end do
        !end function 
        end module
