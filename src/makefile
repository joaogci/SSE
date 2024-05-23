INCLUDES_SSE = sse/sse.h sampling/observables.h sampling/sampling.h \
	sampling/transport.h hamiltonian/hamiltonian.h hamiltonian/lattice.h \
	io/io.h analysis/analysis.h rng/pcg_basic.h
SOURCES_SSE = main.c sse/sse.c sampling/observables.c sampling/sampling.c \
	sampling/transport.c hamiltonian/hamiltonian.c hamiltonian/lattice.c io/io.c \
	analysis/analysis.c rng/pcg_basic.c
INCLUDES_ANA = hamiltonian/lattice.h io/io.h analysis/analysis.h sampling/observables.h
SOURCES_ANA = ana.c hamiltonian/lattice.c io/io.c analysis/analysis.c sampling/observables.c

FLAGS = -fopenmp -O2 
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

main: read_parameters $(SOURCES_SSE) $(INCLUDES_SSE)
	$(CC) -o main $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o io/read_parameters_analysis.o $(SOURCES_SSE) $(LIBS)

ana: read_parameters $(SOURCES_ANA) $(INCLUDES_ANA) 
	$(CC) -o ana $(FLAGS) $(DEBUG_FLAGS) io/read_parameters.o io/read_parameters_analysis.o $(SOURCES_ANA) $(LIBS)

read_parameters: io/read_parameters.f90 io/read_parameters_analysis.f90
	$(FORT) -o io/read_parameters.o -c io/read_parameters.f90
	$(FORT) -o io/read_parameters_analysis.o -c io/read_parameters_analysis.f90

clean: 
	rm io/read_parameters.o 
	rm io/read_parameters_analysis.o 
	rm main
	rm ana
