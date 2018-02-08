        program main 
        use commons
        use meso_approx
        implicit none
        integer, allocatable :: state(:)
        type (approximation_type) :: obj_approx
        real*8 energy, part

        ! Initialising the data structure of the model
        write(*,*) 'Approximation?'
        read(*,*) obj_approx%approx
        write(*,*) 'chemp?'
        read(*,*) chemp
        write(*,*) 'temp?'
        read(*,*) temp
        call obj_approx%init
        allocate(state(nsites)) 
        call confs(state,0)
        write(*,*) part(obj_approx)

        end program 
