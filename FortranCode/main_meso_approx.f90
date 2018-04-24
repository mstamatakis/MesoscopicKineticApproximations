program main
    
    use global_constants
    use meso_approx_inst
    use parser_module
    use nrtype
    use nr
    
    implicit none
    integer i,j,ndeg
    real(8) mu, munext
    
    integer N, ITER
    real(DP), allocatable, dimension(:) :: P
    real(DP), allocatable, dimension(:) :: Pnewguess
    real(DP) H0, H0newguess
    real(8), allocatable, dimension(:,:) :: XI
	REAL(DP), allocatable, DIMENSION(:) :: fvec
    real(8) FTOL, FRET    
    logical check
    character(10) approx

    !real(8) x,y,dy
    !real(8) xa(2),ya(2)
    !xa = (/1.d0,2.d0/)
    !ya = (/2.d0,4.d0/)
    !x = 3.d0    
    !call polint(xa,ya,x,y,dy)
    !continue
    !stop
    
    
    !approx = 'BPEC'
    approx = 'K2NNC2'
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
    
    P(1) = 0.0857037076391d0                 
    P(2) = 0.6409559144920d0                 
    P(3) = 0.1061338496444d0                
    P(4) = -0.0424423610796d0                
    P(5) = -0.0279416468428d0                
    P(6) = -0.1775037458119d0
    
    !P(1) = 0.135d0
    !P(2) = -0.108d0

    !P(1) = 0.300202951383935d0
    !P(2) = -0.106778775188085d0
    do i = 1,N
        obj_approx%hamilt%corcpars(i) = P(i)
    enddo
    
    call obj_approx%calc_resid()

    allocate(XI(N,N),source=0.0d0)
    do i = 1,N
        XI(i,i) = 1.d0
    enddo
    allocate(fvec(N),source=0.0d0)
    fvec = obj_approx%eqns%residual
    
    call fdjac(P,fvec,XI)  
    write(*,*) '-----------------------'
    write(*,*) 'Numerical Jacobian'
    do i = 1,N
        write(*,'(' // trim(int2str(N)) // 'ES15.5E3)') (XI(i,j),j=1,N)
    enddo
    write(*,*) 'Analytical Jacobian'
    do i = 1,N
        write(*,'(' // trim(int2str(N)) // 'ES15.5E3)') (obj_approx%eqns%jacobian(i,j),j=1,N)
    enddo
    write(*,*) '-----------------------'

    !pause
    !continue
    !stop
    
    FTOL=1.D-12

    !call powell(p,xi,n,n,ftol,iter,fret)
    call newt(P,check)

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

    ! Degree of polynomial to be used in continuation of initial guesses
    ndeg = 1
    allocate(Pnewguess(N),source=0.d0)
    
    open(unit=101,file=trim(approx) // '_Fortran_Theta_vs_Mu.txt')
    do mu = mu0,mu1,Dmu
    !do mu = -0.5d0,mu1,Dmu
        
        obj_approx%mu = mu
        
        write(*,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'        
        write(*,*) ' Chemical potential:', mu

        !call powell(p,xi,n,n,ftol,iter,fret)
        call newt(P,check)
        write(*,*) ' Coverage:           ', obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn
        write(*,*) ' Partition function: ', obj_approx%partfcn
        write(*,*) ' Residual norm:      ', sqrt(sum(obj_approx%eqns%residual**2))
        write(*,'(1x," at point: "' // trim(int2str(N)) // 'F32.13)') (P(i), i = 1,N)
       
        write(101,'(' // trim(int2str(N+2)) // 'F32.13)') mu, obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn, (P(i), i = 1,N)
        
        ! Reset direction for the next solution
        XI = 0.d0
        do i = 1,N
            XI(i,i) = 1.d0
        enddo

        H0 = obj_approx%hamilt%h0 + kboltz*obj_approx%temp*log(obj_approx%partfcn) ! "absolute" energies
        munext = mu + Dmu
        call InitialGuessesContinuation(mu,munext,N,ndeg,H0,H0newguess,P,Pnewguess)
        
        ! Set parameter values to new initial guesses
        P = Pnewguess
        
        ! Re-reference energy levels to avoid large values of the partition function sums
        obj_approx%hamilt%h0 = H0newguess
        
        continue
        
    enddo
    
    pause
    continue    
  
    stop    
    
    contains
    
        subroutine InitialGuessesContinuation(x,xnext,N,ndeg,H0,H0newguess,P,Pnewguess)
            
            use nrtype
        
            implicit none
            
            integer ndeg, N
            integer, save :: ipopul = 0
            
            real(DP) H0, H0newguess, x, xnext, errestim
            real(DP), dimension(N) :: P, Pnewguess
            real(DP), dimension(:), allocatable, save :: xrange
            real(DP), dimension(:), allocatable, save :: H0prevguess
            real(DP), dimension(:,:), allocatable, save :: Pprevguess

            if (.not.allocated(H0prevguess)) then
                allocate(xrange(ndeg+1),source=0.d0)
                allocate(H0prevguess(ndeg+1),source=0.d0)
                allocate(Pprevguess(N,ndeg+1),source=0.d0)
            endif
            
            ! Create good initial guesses via the continuation/interpolation scheme
            if (ipopul < ndeg) then ! First few iterations where we need to populate the guess vectors
                ipopul = ipopul + 1
                ! The guess vector will be set to the latest solution (0th
                ! order continuation scheme) and the H0 constant to minus the free energy
                xrange(ipopul) = x
                H0prevguess(ipopul) = H0
                Pprevguess(:,ipopul) = P
                H0newguess = H0
                Pnewguess = P
                
            elseif (ipopul == ndeg) then ! we have enough information to generate the first initial guess with continuation
                
                ipopul = ipopul + 1
                xrange(ipopul) = x
                H0prevguess(ipopul) = H0
                Pprevguess(:,ipopul) = P
                call polint(xrange,H0prevguess,xnext,H0newguess,errestim)
                do i = 1,N
                    call polint(xrange,Pprevguess(i,:),xnext,Pnewguess(i),errestim)
                enddo
                continue
    
            elseif (ipopul > ndeg) then ! by now ipopul = ndeg+1 and will remain so for the rest of the run

                ! shifting back and deleting old entries - one could use linked lists for this,
                ! but we are not dealing with large arrays so probably it's not going to make a difference
                continue
                do i = 1,ndeg
                    xrange(i) = xrange(i+1)
                    H0prevguess(i) = H0prevguess(i+1)
                    Pprevguess(:,i) = Pprevguess(:,i+1)
                enddo
                xrange(ipopul) = x
                H0prevguess(ipopul) = H0
                Pprevguess(:,ipopul) = P
                continue
                call polint(xrange,H0prevguess,xnext,H0newguess,errestim)
                do i = 1,N
                    call polint(xrange,Pprevguess(i,:),xnext,Pnewguess(i),errestim)
                enddo
                continue
                
            endif
        
        end subroutine InitialGuessesContinuation

end program