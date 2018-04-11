        program main 
        use commons
        use meso_approx
        use approx_inst
        implicit none
        logical check
        integer i, j
        real*8 cov, start, finish
        integer, dimension(1) :: v
        integer, allocatable :: state(:)
        real*8, allocatable :: enerv(:)

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

        call cpu_time(start)
        call obj_approx%init

        ! Initialising the state vectors
        allocate(state(obj_approx%nsites),enerv(2**obj_approx%nsites))
        enerv=0.d0
        call confs(state,0)

        ! Coverage vs Chemical Potential Plot
        chemp=-1.40d0
        h0=0.d0
        do i=1,240
         chemp=chemp+0.01d0
         call solver(obj_approx%hamilt%corr%value,obj_approx%eqn%neqns,check,enerv)
         cov=0.d0
         do j=1,obj_approx%nsites
          v(1)=j
          cov=cov+obj_approx%corfun(v,1,obj_approx,enerv)/obj_approx%part(obj_approx,enerv)
         end do
         cov=cov/obj_approx%nsites
         write(16,*) chemp, cov, obj_approx%part(obj_approx,enerv)
         h0=h0+kb*temp*log(obj_approx%part(obj_approx,enerv))
        end do
        call cpu_time(finish)
        print '("Time = ", f10.3, " seconds")',finish-start
        end program 
