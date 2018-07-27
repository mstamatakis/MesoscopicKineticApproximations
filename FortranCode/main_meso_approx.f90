program main
    
    use global_constants
    use meso_approx_inst
    use parser_module
    use nrtype
    use nr
    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
	
    implicit none
    integer i,j,imu,Nmu,ndeg
    real(8) mu, munext
    
    integer N, ITER
    real(DP), allocatable, dimension(:) :: P
    real(DP), allocatable, dimension(:) :: Pnewguess
    real(DP) H0, H0newguess
    real(8), allocatable, dimension(:,:) :: XI
	REAL(DP), allocatable, DIMENSION(:) :: fvec
    real(8) FTOL, FRET
    real(4) t1, t2
    logical check
    character(10) approx    
    
	open(unit=501,file="calculation_input.dat",form="formatted",status="old",action="read")
	read(501,"(a)") approx
	close(501)
    !approx = 'BPEC'
    !approx = 'K2NNC2'
    !approx = 'K3NNC2'
    call obj_approx%init(trim(approx))
    call obj_approx%prnt()
    
    call obj_approx%calc_resid()
    
    N=obj_approx%hamilt%ncorc         !number of variables

    allocate(P(N),source=0.d0)
    
    !allocate(XI(N,N),source=0.0d0)
    !do i = 1,N
    !    XI(i,i) = 1.d0
    !enddo
    
    ! Comment test...
    
    FTOL=1.D-12
    ITER = 0
    FRET = 0.d0
    
    !call powell(p,xi,n,n,ftol,iter,fret)
    call newt(P,check)

    write(*,*) ' '
    write(*,*) ' -------------------------------'        
    write(*,*) ' Chemical potential:', obj_approx%mu
    write(*,*) ' Coverage:           ', obj_approx%eqns%corrlvalue(1)/obj_approx%partfcn
    write(*,*) ' Partition function: ', obj_approx%partfcn
    write(*,*) ' Residual norm:      ', sqrt(sum(obj_approx%eqns%residual**2))
    write(*,'(1x," at point: "' // trim(int2str(N)) // 'F32.13)') (P(i), i = 1,N)

    ! Degree of polynomial to be used in continuation of initial guesses
    ndeg = 2
    allocate(Pnewguess(N),source=0.d0)
    
    open(unit=101,file=trim(approx) // '_Fortran_Theta_vs_Mu.txt')
    
    call cpu_time(t1) ! function for calculating elapsed CPU time
	Nmu = nint((mu1-mu0)/Dmu)
    do imu = 1,Nmu+1
        
		mu = mu0 + (imu-1)*Dmu
		
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
        !XI = 0.d0
        !do i = 1,N
        !    XI(i,i) = 1.d0
        !enddo

        H0 = obj_approx%hamilt%h0 + kboltz*obj_approx%temp*log(obj_approx%partfcn) ! "absolute" energies
        munext = mu + Dmu
        call InitialGuessesContinuation(mu,munext,N,ndeg,H0,H0newguess,P,Pnewguess)
        
        ! Set parameter values to new initial guesses
        P = Pnewguess
        
        ! Re-reference energy levels to avoid large values of the partition function sums
        obj_approx%hamilt%h0 = H0newguess
        
        continue
        
    enddo
    call cpu_time(t2) ! function for calculating elapsed CPU time
    
    write(*,*) 'Elapsed CPU-time',t2-t1
    
	!print *, 'Program paused, press Enter to continue.'
	!read(stdin,*) 
    !continue    
  
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