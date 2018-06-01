function funcv(x)

use nrtype
use meso_approx_inst
use global_constants

implicit none

real(dp), dimension(:), intent(in) :: x
real(dp), dimension(size(x)) :: funcv

obj_approx%hamilt%corcpars(1:obj_approx%hamilt%ncorc) = x(1:obj_approx%hamilt%ncorc)

call obj_approx%calc_resid(calcjac=.false.)
funcv = obj_approx%eqns%residual

fcncounts = fcncounts + 1

return

end function funcv
    
subroutine jacobian(x,jaca)
    
    use meso_approx_inst

    USE nrtype
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: x
    REAL(DP), DIMENSION(size(x),size(x)) :: jaca
    
    ! #MSTAM: This can be wasteful: if funcv has been called the jacobian could have already been computed
    ! We'll improve speed at a later stage after making sure the code is correct
    ! For now we have switched off Jacobian calculations in funcv. Worse case scenario we duplicate
    ! the calculatios of the right hand side...
    call obj_approx%calc_resid()
    jaca = obj_approx%eqns%jacobian
    
    return
    
end subroutine jacobian    
