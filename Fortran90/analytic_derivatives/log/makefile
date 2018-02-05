prog.x: constants.o main.o 
	gfortran -o prog.x constants.o main.o
constants.mod: constants.o constants.f90
	gfortran -c constants.f90
constants.o: constants.f90
	gfortran -c constants.f90
main.o: constants.mod main.f90
	gfortran -c main.f90
clean:
	rm constants.mod constants.o main.o prog.x
