        program main 
        use commons
        use meso_approx
        implicit none
        logical check
        !real*8 energy, part, corr
        integer, allocatable :: state(:)
        type (approximation_type) :: obj_approx

        real*8, allocatable, dimension(:) :: fvec

        ! Initialising the data structure of the model
        write(*,*) '--------------------------------------------'
        write(*,*) '|              APPROXIMATION               |'
        write(*,*) '--------------------------------------------'
        read(*,*) obj_approx%approx
        write(*,*) '--------------------------------------------'
        write(*,*) '|              THERMODYNAMICS              |'
        write(*,*) '--------------------------------------------'
        write(*,*) 'chemp?'
        read(*,*) chemp
        write(*,*) 'temp?'
        read(*,*) temp
        call obj_approx%init

        ! Initialising the state vectors
        allocate(state(nsites),fvec(npar))
        call confs(state,0)
        call funcv(fvec)

        ! Calculating the self-consistent correction fields
        check=.false.
        !call solver(obj_approx%hamilt%corr%value,npar,check)
        end program 
