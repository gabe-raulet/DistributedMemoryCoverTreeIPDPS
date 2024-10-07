DEBUG?=0
LOG?=0
D?=16
FP?=32
FLAGS=-std=c++20 -fopenmp -DDIM_SIZE=$(D) -DFP_SIZE=$(FP)
INCS=-I./misc -I./ptraits -I./mpienv -I./

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

all: main_mpi main ptgen

.PHONY: version.h

version.h:
	@echo "#define GIT_COMMIT \"$(shell git describe --always --dirty --match 'NOT A TAG')\"" > version.h

main_mpi: src/main_mpi.cpp misc ptraits ctree version.h
	$(MPI_COMPILER) -o main_mpi $(FLAGS) $(INCS) -I./ctree $<

main: src/main.cpp misc ptraits ctree version.h
	$(MPI_COMPILER) -o main $(FLAGS) $(INCS) -I./ctree $<

ptgen: src/ptgen.cpp misc ptraits version.h
	$(COMPILER) -o ptgen $(FLAGS) $(INCS) $<

clean:
	rm -rf main_mpi main ptgen version.h *.dSYM
