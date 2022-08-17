SOURCES = main.c sse/sse.c sampling/sampling.c hamiltonian/hamiltonian.c vtx/vtx_type.c io/io.c
INCLUDES = sse/sse.h sampling/sampling.h hamiltonian/hamiltonian.h vtx/vtx_type.h rng/xorshiro256++.h io/io.h
FLAGS = -fopenmp -O3
LIBS = -lm
PROGRAM = main

$(PROGRAM): $(SOURCES) $(INCLUDES)
	gcc-11 -o $(PROGRAM) $(SOURCES) $(FLAGS) $(LIBS)

clean: 
	rm $(PROGRAM)
