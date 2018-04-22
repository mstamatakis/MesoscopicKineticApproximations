program main
    
    use global_constants
    use meso_approx_inst
    use parser_module
    
    implicit none
    integer i
    real(8) mu
    
    integer N, ITER
    real(8), allocatable, dimension(:) :: P
    real(8), allocatable, dimension(:,:) :: XI
    real(8) FTOL, FRET    
    
    !call obj_approx%init('BPEC')
    call obj_approx%init('K2NNC2')
    call obj_approx%prnt()
    
    !obj_approx%hamilt%origpars(1) = -1.5d0
    obj_approx%hamilt%corcpars(1) =  0.025d0
    obj_approx%hamilt%corcpars(2) =  0.010d0
    obj_approx%hamilt%corcpars(3) =  0.035d0
    obj_approx%hamilt%corcpars(4) =  0.040d0
    obj_approx%hamilt%corcpars(5) =  0.035d0
    obj_approx%hamilt%corcpars(6) =  0.040d0
    call obj_approx%prnt()
    
    call obj_approx%calc_energ()
    call obj_approx%calc_energ()
    
    !do i = 1,2**7
    !    write(*,'(a,f22.16)') trim(int2str(i)), obj_approx%allenergs(i)
    !enddo

    call obj_approx%calc_resid()
    
    N=obj_approx%hamilt%ncorc         !number of variables

    allocate(P(N),source=0.001d0)
    !allocate(P(N),source=(/0.300202951383935d0, -0.106778775188085d0/))
    allocate(XI(N,N),source=0.0d0)
    do i = 1,N
        XI(i,i) = 1.d0
    enddo

    FTOL=1.D-12

    CALL POWELL(P,XI,N,N,FTOL,ITER,FRET)

    print *,' '
    print *,' Number of iterations:', ITER
    print *,' '
    print *,' Minimum value:', FRET
    print *,' '
    print *,' at point:',P(1),' ',P(2),' ',P(3),' ',P(4),' ',P(5),' ',P(6)
    print *,' '

    open(unit=101,file='K2NNC2_Fortran_Theta_vs_Mu.txt')
    do mu = mu0,mu1,Dmu
        
        obj_approx%mu = mu
        
        call powell(p,xi,n,n,ftol,iter,fret)
        print *,' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'        
        print *,' Chemical potential:', mu
        print *,' Coverage:          ', obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn
        print *,' Number of iterations:', ITER
        print *,' Minimum value:', FRET
        print *,' at point:',P(1),' ',P(2),' ',P(3),' ',P(4),' ',P(5),' ',P(6)
       
        write(101,'(' // trim(int2str(N+2)) // 'F32.13)') mu, obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn, (P(i), i = 1,N)
        
        ! Reset direction for the next solution
        XI = 0.d0
        do i = 1,N
            XI(i,i) = 1.d0
        enddo

        
        ! Re-reference energy levels to avoid large values of the partition function sums
        obj_approx%hamilt%h0 = obj_approx%hamilt%h0 + kboltz*obj_approx%temp*log(obj_approx%partfcn)
        
        continue
        
    enddo
    
    pause
    continue    
  
stop    
    
end program