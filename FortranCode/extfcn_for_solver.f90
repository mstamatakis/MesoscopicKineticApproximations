      SUBROUTINE FCN(N,X,FVEC,IFLAG)
      INTEGER N,IFLAG
      DOUBLE PRECISION X(N),FVEC(N)
      INTEGER K
      DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
      DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/

      IF (IFLAG .NE. 0) GO TO 5

!     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.

      RETURN
    5 CONTINUE
      DO 10 K = 1, N
         TEMP = (THREE - TWO*X(K))*X(K)
         TEMP1 = ZERO
         IF (K .NE. 1) TEMP1 = X(K-1)
         TEMP2 = ZERO
         IF (K .NE. N) TEMP2 = X(K+1)
         FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
   10    CONTINUE
      RETURN
    END  
    
    
         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
         INTEGER N,LDFJAC,IFLAG
         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
!         ----------
!         Calculate the Jacobian at X and return this
!         matrix in FJAC.  FVEC contains the function
!         values at X and should not be altered.
!         ----------
         RETURN
         END
    