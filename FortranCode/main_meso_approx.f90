program main
    
    use global_constants
    use meso_approx
    
    implicit none

    type (approximation) :: obj_approx
    
    
    call obj_approx%init()
    call obj_approx%prnt()

    !obj_approx%hamilt%origpars(1) = -1.5d0
    !obj_approx%hamilt%corcpars(1) =  0.25d0
    !obj_approx%hamilt%corcpars(2) =  0.10d0
    call obj_approx%prnt()
    
    call obj_approx%calc_energ()
    call obj_approx%calc_energ()
    
    pause
    continue    
  
stop    
    
end program