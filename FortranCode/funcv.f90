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