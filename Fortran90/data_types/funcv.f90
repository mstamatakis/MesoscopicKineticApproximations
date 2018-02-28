	subroutine funcv(fvec)
        use commons
        use meso_approx
        use approx_inst
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8, intent(out) :: fvec(npar)
	
        !$OMP PARALLEL
        !$OMP DO
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
         fvec(i)=log(obj_approx%corfun(obj_approx%eqn%lhs(i,1:total),total,obj_approx))-log(obj_approx%corfun(obj_approx%eqn%rhs(i,1:total),total,obj_approx))
        end do
        !$OMP END DO
        !$OMP END PARALLEL
	end subroutine funcv
