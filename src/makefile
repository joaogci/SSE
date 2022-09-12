SOURCES = main.c sse/sse.c sampling/sampling.c io/io.c
INCLUDES = sse/sse.h sampling/sampling.h vtx/vtx_type.h rng/xorshiro256++.h io/io.h
FLAGS = -fopenmp -O3 -Wall
LIBS = -lm
PROGRAM = main

$(PROGRAM): $(SOURCES) $(INCLUDES)
	gcc-12 -o $(PROGRAM) $(SOURCES) $(FLAGS) $(LIBS)

test: $(SOURCES) $(INCLUDES)
	gcc-12 -o $(PROGRAM) $(SOURCES) $(FLAGS) $(LIBS) -DTESTING

series: $(SOURCES) $(INCLUDEs) sampling/autocorrelation.h sampling/autocorrelation.c
	gcc-12 -o $(PROGRAM) $(SOURCES) sampling/autocorrelation.c $(FLAGS) $(LIBS) -DM_AUTOCORRELATION

clean: 
	rm $(PROGRAM)
