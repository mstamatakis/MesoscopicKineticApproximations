MODULE fminln
	USE nrtype; USE nrutil, ONLY :nrerror
	REAL(DP), DIMENSION(:), POINTER :: fmin_fvecp
	CONTAINS
	FUNCTION fmin(x)
		IMPLICIT NONE
		REAL(DP), DIMENSION(:), INTENT(IN) :: x
		REAL(DP) :: fmin
		! Returns f = 1/2F*F at x. FUNCTION funcv(x) is a fixed-name, user-supplied routine that
		! returns the vector of functions at x. The pointer fmin vecp communicates the function
		! values back to newt.
		INTERFACE
			FUNCTION funcv(x)
				USE nrtype
				IMPLICIT NONE
				REAL(DP), DIMENSION(:), INTENT(IN) :: x
				REAL(DP), DIMENSION(size(x)) :: funcv
			END FUNCTION funcv
		END INTERFACE
		if (.not. associated(fmin_fvecp)) call &
		nrerror('fmin:problem with pointer for returned values')
		fmin_fvecp=funcv(x)
		fmin=0.5_sp*dot_product(fmin_fvecp,fmin_fvecp)
	END FUNCTION fmin
END MODULE fminln	

