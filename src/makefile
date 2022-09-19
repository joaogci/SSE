SOURCES = main.c sse/sse.c sampling/sampling.c io/io.c
INCLUDES = sse/sse.h sampling/sampling.h vtx/vtx_type.h rng/xorshiro256++.h io/io.h
FLAGS = -fopenmp -O3 -Wall
LIBS = -lm
PROGRAM = main

ifeq ($(OS),Windows_NT)
	$(Windows is not supported. Use a UNIX OS.)
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CC = gcc
	endif
	ifeq ($(UNAME_S),Darwin)
		CC = gcc-12
	endif
endif

$(PROGRAM): $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(SOURCES) $(FLAGS) $(LIBS)

test: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(SOURCES) $(FLAGS) $(LIBS) -DTESTING

series: $(SOURCES) $(INCLUDEs) sampling/autocorrelation.h sampling/autocorrelation.c
	$(CC) -o $(PROGRAM) $(SOURCES) sampling/autocorrelation.c $(FLAGS) $(LIBS) -DM_AUTOCORRELATION

clean: 
	rm $(PROGRAM)
