	subroutine funcv(fvec)
        use commons
        use meso_approx
        use approx_inst
	implicit none
        integer i, j, total
	integer, allocatable, dimension(:) ::  v1, v2
	real*8, intent(out) :: fvec(npar)
	
        call obj_approx%sigma()
        do i=1,npar
         fvec(i)=obj_approx%res(i)
        end do
	end subroutine funcv
