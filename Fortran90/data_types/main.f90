        program main 
        use commons
        use meso_approx
        implicit none
        integer, allocatable :: state(:)
        type (approximation_type) :: obj_approx
        real*8 energy

        ! Initialising the data structure of the model
        write(*,*) 'Approximation?'
        read(*,*) obj_approx%approx
        call obj_approx%init
        allocate(state(nsites)) 
        state=1
        write(*,*) energy(obj_approx,state)

        end program 
