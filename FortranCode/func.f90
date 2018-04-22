real*8 function func(p)

use meso_approx_inst

real*8 p(obj_approx%hamilt%ncorc)

obj_approx%hamilt%corcpars(1:obj_approx%hamilt%ncorc) = P(1:obj_approx%hamilt%ncorc)

call obj_approx%calc_resid()
func = sum(obj_approx%eqns%residual**2)

return

end function func