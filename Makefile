DEBUG?=0
LOG?=0
D?=3
FP?=32
FLAGS=-std=c++20 -fopenmp -DDIM_SIZE=$(D) -DFP_SIZE=$(FP)
INCS=-I./misc -I./ptraits -I./mpienv -I./

EXE?=D$(D).FP$(FP)

ifeq ($(shell uname -s),Linux)
COMPILER=CC
MPI_COMPILER=CC
else
COMPILER=clang++
MPI_COMPILER=mpic++
endif

ifeq ($(DEBUG),1)
FLAGS+=-O0 -g -fsanitize=address -fno-omit-frame-pointer -DDEBUG
else
FLAGS+=-O2
endif

ifeq ($(LOG),1)
FLAGS+=-DLOG
endif

all: omp_test_driver ptgen

.PHONY: version.h

version.h:
	@echo "#define GIT_COMMIT \"$(shell git describe --always --dirty --match 'NOT A TAG')\"" > version.h

ptgen: src/ptgen.cpp misc ptraits version.h
	$(COMPILER) -o ptgen $(FLAGS) $(INCS) $<

omp_test_driver: testing/omp_test_driver.cpp misc ptraits ctree version.h
	$(MPI_COMPILER) -o omp_test_driver $(FLAGS) $(INCS) -I./ctree $<

clean:
	rm -rf ptgen omp_test_driver version.h *.dSYM
