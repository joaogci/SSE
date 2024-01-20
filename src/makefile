SOURCES = main.c sse/sse.c sampling/sampling.c io/io.c sampling/observables.c hamiltonian/hamiltonian.c hamiltonian/lattice.c rng/pcg_basic.c 
INCLUDES = sse/sse.h sampling/sampling.h io/io.h sampling/observables.h hamiltonian/hamiltonian.h hamiltonian/lattice.h rng/pcg_basic.h
FLAGS = -g -fopenmp -O2 
DEBUG_FLAGS = -Wall -Wshadow
LIBS = -lm -lgfortran
PROGRAM = main

ifeq ($(OS),Windows_NT)
	$(Windows is not supported. Use a UNIX OS.)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CC = gcc
		FORT = gfortran
	endif
	ifeq ($(UNAME_S),Darwin)
		CC = gcc-13
		FORT = gfortran
	endif
endif

$(PROGRAM): $(SOURCES) $(INCLUDES) io/read_parameters.f90
	$(FORT) -o io/read_parameters.o -c io/read_parameters.f90
	$(CC) -o $(PROGRAM) $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o $(SOURCES) $(LIBS)

clean: 
	rm io/read_parameters.o 
	rm $(PROGRAM)
