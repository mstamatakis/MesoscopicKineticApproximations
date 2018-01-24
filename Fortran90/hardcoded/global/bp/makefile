prog.x: constants.o hard_bp.o fmin.o
	gfortran -o prog.x constants.o hard_bp.o fmin.o
constants.mod: constants.o constants.f90
	gfortran -c constants.f90
constants.o: constants.f90
	gfortran -c constants.f90
hard_bp.o: constants.mod hard_bp.f90
	gfortran -c hard_bp.f90
fmin.o: hard_bp.o fmin.f90
	gfortran -c fmin.f90
clean:
	rm constants.mod constants.o hard_bp.o fmin.o prog.x
