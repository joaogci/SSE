SOURCES = sse/sse.c sampling/sampling.c io/io.c sampling/observables.c hamiltonian/hamiltonian.c hamiltonian/lattice.c rng/pcg_basic.c analysis/analysis.c
INCLUDES = sse/sse.h sampling/sampling.h io/io.h sampling/observables.h hamiltonian/hamiltonian.h hamiltonian/lattice.h rng/pcg_basic.h analysis/analysis.h
FLAGS = -g -fopenmp -O2 
DEBUG_FLAGS = -Wall -Wshadow
LIBS = -lm -lgfortran


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

all: main ana

main: read_parameters $(SOURCES) $(INCLUDES) main.c
	$(CC) -o main $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o main.c $(SOURCES) $(LIBS)

ana: read_parameters $(SOURCES) $(INCLUDES) ana.c
	$(CC) -o ana $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o ana.c $(SOURCES) $(LIBS)

read_parameters: io/read_parameters.f90
	$(FORT) -o io/read_parameters.o -c io/read_parameters.f90

clean: 
	rm io/read_parameters.o 
	rm main
	rm ana
