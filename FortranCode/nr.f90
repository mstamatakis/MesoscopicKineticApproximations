MODULE nr
	
	INTERFACE
		SUBROUTINE fdjac(x,fvec,df)
			USE nrtype
			REAL(DP), DIMENSION(:), INTENT(IN) :: fvec
			REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
			REAL(DP), DIMENSION(:,:), INTENT(OUT) :: df
		END SUBROUTINE fdjac
	END INTERFACE

	INTERFACE
		SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func)
			USE nrtype
			REAL(DP), DIMENSION(:), INTENT(IN) :: xold,g
			REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
			REAL(DP), INTENT(IN) :: fold,stpmax
			REAL(DP), DIMENSION(:), INTENT(OUT) :: x
			REAL(DP), INTENT(OUT) :: f
			LOGICAL(LGT), INTENT(OUT) :: check
			INTERFACE
				FUNCTION func(x)
					USE nrtype
					REAL(DP) :: func
					REAL(DP), DIMENSION(:), INTENT(IN) :: x
				END FUNCTION func
			END INTERFACE
		END SUBROUTINE lnsrch
	END INTERFACE

	INTERFACE
		SUBROUTINE lubksb(a,indx,b)
			USE nrtype
			REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
			INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
			REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
		END SUBROUTINE lubksb
	END INTERFACE
	
	INTERFACE
		SUBROUTINE ludcmp(a,indx,d)
			USE nrtype
			REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
			INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
			REAL(DP), INTENT(OUT) :: d
		END SUBROUTINE ludcmp
	END INTERFACE

	INTERFACE
		SUBROUTINE newt(x,check)
			USE nrtype
			REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
			LOGICAL(LGT), INTENT(OUT) :: check
		END SUBROUTINE newt
    END INTERFACE
    
    INTERFACE
        SUBROUTINE polint(xa,ya,x,y,dy)
            USE nrtype
            REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya
            REAL(DP), INTENT(IN) :: x
            REAL(DP), INTENT(OUT) :: y,dy
        END SUBROUTINE polint
    END INTERFACE

END MODULE nr
	
				