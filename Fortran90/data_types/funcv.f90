	subroutine funcv(appr,fvec)
        ! This subroutine applies self-consistent constraints to the
        ! parameters. Its input are the number of parameters "n" and the
        ! vector of parameters "x", while "fvec" is the output. "fvec"
        ! is a function to be zeroed by the solver, which will find the
        ! self-consistent parameters this way.
        use commons
        use meso_approx
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8 fvec(npar), corr
        type (approximation_type) :: appr
	
        do i=1,npar
         total=0
         do j=1,nsites
          if(appr%eqn%lhs(i,j).ne.0) then
            total=total+1
          end if
         end do
         allocate(v1(total),v2(total)) 
         do j=1,total
          v1(j)=appr%eqn%lhs(i,j)
          v2(j)=appr%eqn%rhs(i,j)
         end do
         fvec(i)=log(corr(v1,total,appr))-log(corr(v2,total,appr))
         deallocate(v1,v2)
        end do
        do i=1,npar
         write(*,*) fvec(i)
        end do
	end subroutine funcv
