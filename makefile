SOURCES = src/main.c src/sse/sse.c src/sampling/sampling.c src/io/io.c src/sampling/sampling_memory.c src/rng/pcg_basic.c
INCLUDES = src/sse/sse.h src/sampling/sampling.h src/vtx/vtx_type.h src/io/io.h src/rng/pcg_basic.h
FLAGS = -fopenmp -O2
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
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS)

test: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) -Wall -Wshadow $(SOURCES) $(LIBS) -DTESTING


both: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DSPIN_CONDUCTANCE -DHEAT_CONDUCTANCE -larb

spin: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DSPIN_CONDUCTANCE -larb

heat: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DHEAT_CONDUCTANCE -larb


clean: 
	rm $(PROGRAM)
