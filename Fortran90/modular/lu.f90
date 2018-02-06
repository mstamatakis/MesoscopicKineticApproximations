!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     * 
!*******************************************************
module lu

    contains

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
!  MS 2018-02-06: Defined matrix and vectors to have max dimension
!  of n0. However, only the first n elements are manipulated by the
!  algorithm.
 subroutine ludcmp(a,n,n0,indx,d,code)
 parameter(tiny=1.5d-16)
 real*8  amax,dum, sum, a(n0,n0),vv(n0)
 integer code, d, indx(n0)

 d=1; code=0

 do i=1,n
   amax=0.d0
   do j=1,n
     if (dabs(a(i,j)).gt.amax) amax=dabs(a(i,j))
   end do ! j loop
   if(amax.lt.tiny) then
     code = 1
     return
   end if
   vv(i) = 1.d0 / amax
 end do ! i loop

 do j=1,n
   do i=1,j-1
     sum = a(i,j)
     do k=1,i-1
       sum = sum - a(i,k)*a(k,j) 
     end do ! k loop
     a(i,j) = sum
   end do ! i loop
   amax = 0.d0
   do i=j,n
     sum = a(i,j)
     do k=1,j-1
       sum = sum - a(i,k)*a(k,j) 
     end do ! k loop
     a(i,j) = sum
     dum = vv(i)*dabs(sum)
     if(dum.ge.amax) then
       imax = i
       amax = dum
     end if
   end do ! i loop  
   
   if(j.ne.imax) then
     do k=1,n
       dum = a(imax,k)
       a(imax,k) = a(j,k)
       a(j,k) = dum
     end do ! k loop
     d = -d
     vv(imax) = vv(j)
   end if

   indx(j) = imax
   if(dabs(a(j,j)) < tiny) a(j,j) = tiny

   if(j.ne.n) then
     dum = 1.d0 / a(j,j)
     do i=j+1,n
       a(i,j) = a(i,j)*dum
     end do ! i loop
   end if 
 end do ! j loop

 return
 end subroutine ludcmp


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 subroutine lubksb(a,n,n0,indx,b)
 real*8  sum, a(n0,n0),b(n0)
 integer indx(n0)

 ii = 0

 do i=1,n
   ll = indx(i)
   sum = b(ll)
   b(ll) = b(i)
   if(ii.ne.0) then
     do j=ii,i-1
       sum = sum - a(i,j)*b(j)
     end do ! j loop
   else if(sum.ne.0.d0) then
     ii = i
   end if
   b(i) = sum
 end do ! i loop

 do i=n,1,-1
   sum = b(i)
   if(i < n) then
     do j=i+1,n
       sum = sum - a(i,j)*b(j)
     end do ! j loop
   end if
   b(i) = sum / a(i,i)
 end do ! i loop

 return
 end subroutine lubksb

end module lu

! end of file lu.f90