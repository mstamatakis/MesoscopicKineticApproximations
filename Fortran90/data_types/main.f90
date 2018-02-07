        program la_forge        
        
        use commons
        use meso_approx
        
        implicit none
        
        integer, allocatable :: state(:)
        type (approximation_type) :: obj_approx1

        state=1
        
        call obj_approx1%init
        
        continue
        
        stop
        
        end program 
