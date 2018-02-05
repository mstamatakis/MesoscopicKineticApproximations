!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine gelim(a,b,n,n0,x)	

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
								
! Subroutine for solving system the linear system A*x = b
! with Gauss Elimination

! Input variables:
! a = matrix A
! b = vector b
! n = order of the system

! Output variables:
! x = vector, solution of the system 

! Other variables:
! sa, sb = helping variables that are used in interchanging
!          the rows (for pivoting)
! fa = variable that contains the value of the 1st non-zero 
!	   element of every row, used in elimination procedure
! amax = value of the max element of the current row
! imax = s/n (index) of the row that contains the max element
!        of the current column
! prs = variable used in back substitution

! Parameters:
! n0 = maximum allowed order of system

! Subroutines called by gelim:
! augmprnt (for printing the augmented matrix)

! BEWARE: n0 should be the same in the main program and this 
!         subroutine otherwise errors occur.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

implicit none
Integer n0
Double Precision a(n0,n0), b(n0), sa(n0), x(n0), fa, amax, sb, prs
Integer i, j, k, n, imax

!write(*,*) '[gelim] Matrix as given from main.'
!call augmprnt(a,b,n,n0)
!pause
 
do i = 1,n					   

! Find max element for correct pivoting
 amax = 0.D0					   
 imax = 0.D0
 do j = i,n
  if ((abs(a(j,i))).gt.(amax)) then
   amax = abs(a(j,i))
   imax = j
  endif
 enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Check for singularity: If all elements of a column = 0 then STOP!
! (Practically " = 0" means " < threshold" = 1e-40 here)
 if (amax.lt.1D-40) then	   
   write(*,*) "[gelim] Singular or ill-conditioned matrix. EXITING!"
   pause
   stop
 endif

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Pivoting: Interchange rows so as max element is used in elimination
! The row that has the desired pivotal element is temporarily stored in sa
! so that the interchange can be performed.
 if (imax.ne.i) then
  do j = 1,n					  
   sa(j) = a(imax,j)
   a(imax,j) = a(i,j)
   a(i,j) = sa(j)
  enddo
   sb = b(imax)
   b(imax) = b(i)
   b(i) = sb
 endif

! write(*,*) '[gelim] Matrix after pivoting.'
! call augmprnt(a,b,n,n0)
! pause

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Elimination Proccess for the next rows. fa is the first non zero element of each
! row. 
 do j = i+1,n			! j = Row counter
  fa = a(j,i)												
  do k = i,n			! k = column counter
   a(j,k) = a(j,k)-(fa/(a(i,i)))*a(i,k)
  enddo
   b(j) = b(j)-(fa/(a(i,i)))*b(i)
 enddo

! write(*,*) '[gelim] matrix after elimination]'
! call augmprnt(a,b,n,n0)
! pause

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

enddo
 
! Perform back substitution
x(n) = b(n)/a(n,n)
do i = 1,n-1
 k = n-i
 prs = 0.D0
 do j = k+1,n
  prs = prs+a(k,j)*x(j)
 enddo
 x(k) = (b(k)-prs)/a(k,k)
enddo

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

return
end ! End of subroutine gelim


!###################################################################################


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!subroutine augmprnt(a,b,n,n0)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Subroutine for printing A, b as an augmented matrix

! Input variables:
! A = matrix
! b = vector
! n = order of the system

! Parameters:
! m = maximum allowed order of system

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!Implicit none
!Integer n0
!Double Precision a(n0,n0), b(n0)
!Integer n, k, j

!do k = 1,n
!   write(*,8) (a(k,j), j = 1,n), b(k)
!8  format(<n0+1>(F10.5,1x))
!   write(*,*) "  "
! enddo
! write(*,*) " "
! write(*,*) " "
!return
!end ! End of subroutine prnt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

