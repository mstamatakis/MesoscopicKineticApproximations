        program la_forge
        use data
        implicit none
        integer, dimension(nsites) :: state
        
        state=1
        write(*,*) energy(state,bpec,nsites)
        end program 
