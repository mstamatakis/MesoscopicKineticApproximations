        program main 
        use commons
        use meso_approx
        use approx_inst
        implicit none
        logical check
        integer i, j
        real*8 cov
        integer, dimension(1) :: v
        integer, allocatable :: state(:)


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
        allocate(state(nsites))
        call confs(state,0)

        ! Calculating the self-consistent correction fields
        check=.false.
        call solver(obj_approx%hamilt%corr%value,npar,check)
    
        ! Coverage vs Chemical Potential Plot
        chemp=-1.40d0
        do i=1,240
         chemp=chemp+0.01d0
         call solver(obj_approx%hamilt%corr%value,npar,check)
         cov=0.d0
         do j=1,nsites
          v(1)=j
          cov=cov+obj_approx%corfun(v,1,obj_approx)/obj_approx%part()
         end do
         cov=cov/nsites
         write(16,*) chemp, cov
        end do
        end program 
