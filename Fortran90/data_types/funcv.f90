	subroutine funcv(fvec)
        ! This subroutine applies self-consistent constraints to the
        ! parameters. Its input are the number of parameters "n" and the
        ! vector of parameters "x", while "fvec" is the output. "fvec"
        ! is a function to be zeroed by the solver, which will find the
        ! self-consistent parameters this way.
        use commons
        use approx_inst
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8 fvec(npar), corr
	
        do i=1,npar
         total=0
         do j=1,nsites
          if(obj_approx%eqn%lhs(i,j).ne.0) then
            total=total+1
          end if
         end do
         allocate(v1(total),v2(total)) 
         do j=1,total
          v1(j)=obj_approx%eqn%lhs(i,j)
          v2(j)=obj_approx%eqn%rhs(i,j)
         end do
         fvec(i)=log(corr(v1,total,obj_approx))-log(corr(v2,total,obj_approx))
         deallocate(v1,v2)
        end do
        do i=1,npar
         write(*,*) fvec(i)
        end do
	end subroutine funcv
