	subroutine funcv(fvec)
        use commons
        use approx_inst
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8, intent(out) :: fvec(npar)
	
        do i=1,npar
        ! "total" accumulates the number of non-zero elements into the
        ! row vector of "obj_approx%eqn%lhs". This is equal to
        ! obj_approx%eqn%rhs, so we just count it once.
         total=0
         do j=1,nsites
          if(obj_approx%eqn%lhs(i,j).ne.0) then
            total=total+1
          end if
         end do
        ! Knowing that, we allocate the vectors to calculate
        ! correlations.
         allocate(v1(total),v2(total)) 
         do j=1,total
          v1(j)=obj_approx%eqn%lhs(i,j)
          v2(j)=obj_approx%eqn%rhs(i,j)
         end do
         fvec(i)=obj_approx%corfun(v1,total,obj_approx)-obj_approx%corfun(v2,total,obj_approx)
         deallocate(v1,v2)
        end do
	end subroutine funcv
