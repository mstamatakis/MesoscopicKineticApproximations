program main
    
    use global_constants
    use meso_approx_inst
    use parser_module
    use nrtype
    use nr
    
    implicit none
    integer i
    real(8) mu
    
    integer N, ITER
    real(DP), allocatable, dimension(:) :: P
    real(8), allocatable, dimension(:,:) :: XI
    real(8) FTOL, FRET    
    logical check
    character(10) approx

    approx = 'BPEC'
    call obj_approx%init(trim(approx))
    call obj_approx%prnt()
    
    !obj_approx%hamilt%corcpars(1) =  0.025d0
    !obj_approx%hamilt%corcpars(2) =  0.010d0
    !obj_approx%hamilt%corcpars(3) =  0.035d0
    !obj_approx%hamilt%corcpars(4) =  0.040d0
    !obj_approx%hamilt%corcpars(5) =  0.035d0
    !obj_approx%hamilt%corcpars(6) =  0.040d0
    !call obj_approx%prnt()
    
    !call obj_approx%calc_energ()    
    !do i = 1,2**7
    !    write(*,'(a,f22.16)') trim(int2str(i)), obj_approx%allenergs(i)
    !enddo

    call obj_approx%calc_resid()
    
    N=obj_approx%hamilt%ncorc         !number of variables

    allocate(P(N),source=0.d0)

    !P(1) = 0.0000148405758d0
    !P(2) = 0.0000593569692d0
    !P(3) = 0.0000000000000d0
    !P(4) = -0.0000148300566d0
    !P(5) = -0.0000000000000d0
    !P(6) = -0.0000148300585d0
    !
    !P(1) = 0.0857037076391d0                 
    !P(2) = 0.6409559144920d0                 
    !P(3) = 0.1061338496444d0                
    !P(4) = -0.0424423610796d0                
    !P(5) = -0.0279416468428d0                
    !P(6) = -0.1775037458119d0
    
    P(1) = 0.300202951383935d0
    P(2) = -0.106778775188085d0
    allocate(XI(N,N),source=0.0d0)
    do i = 1,N
        XI(i,i) = 1.d0
    enddo

    FTOL=1.D-12

    call powell(p,xi,n,n,ftol,iter,fret)
    !call newt(P,check)

    print *,' '
    print *,' -------------------------------'        
    print *,' Chemical potential:', mu
    print *,' Coverage:          ', obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn
    print *,' Number of iterations:', ITER
    print *,' Minimum value:', FRET
    print *,' at point:',P(1),' ',P(2) !,' ',P(3),' ',P(4),' ',P(5),' ',P(6)
    print *,' -------------------------------'        
    print *,' '
    print *,' '
    print *,' '

    open(unit=101,file=trim(approx) // '_Fortran_Theta_vs_Mu.txt')
    do mu = mu0,mu1,Dmu
    !do mu = -0.5d0,mu1,Dmu
        
        obj_approx%mu = mu
        
        call powell(p,xi,n,n,ftol,iter,fret)
        !call newt(P,check)
        print *,' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'        
        print *,' Chemical potential:', mu
        print *,' Coverage:          ', obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn
        print *,' Number of iterations:', ITER
        print *,' Minimum value:', FRET
        print *,' at point:',P(1),' ',P(2) !,' ',P(3),' ',P(4),' ',P(5),' ',P(6)
       
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