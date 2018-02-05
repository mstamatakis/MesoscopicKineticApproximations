COMPILER = gfortran

OUTPUT := ./prog.x

FLAGS := -fcheck=all -ffree-line-length-1024
LIBS := 

SOURCE = constants confs part corr funcv newton k2nnc1_hamiltonian fmin gelim

EXECUTABLE = main

OBJS = $(SOURCE:%=%.o)
ALL_OBJS = $(SOURCE:%=%.o) $(EXECUTABLE:%=%.o)

.PHONY : all
all: $(ALL_OBJS) $(EXECUTABLE)

.PHONY : clean
clean:
	@echo "Cleaning directory"
	@\rm -rf *.o *.mod *.exe $(OUTPUT)

%.o: ../%.f90
	@echo "Compiling: " $*
	$(COMPILER) -c  ${FLAGS} $< -o $@

%:  %.o $(ALL_OBJS)
	@echo "    Linking program:   " $@
	$(COMPILER) $(FLAGS) $(ALL_OBJS) $(LIBS) -o $(OUTPUT)

