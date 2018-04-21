	subroutine funcv(fvec)
        use commons
        use meso_approx
        use approx_inst
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8, intent(out) :: fvec(obj_approx%eqn%neqns)
	
        
        call obj_approx%residuals()
        do i=1,obj_approx%eqn%neqns
         fvec(i)=obj_approx%res(i)
        end do
	end subroutine funcv
