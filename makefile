SOURCES = src/main.c src/sse/sse.c src/sampling/sampling.c src/io/io.c src/sampling/sampling_memory.c 
INCLUDES = src/sse/sse.h src/sampling/sampling.h src/vtx/vtx_type.h src/rng/xorshiro256++.h src/io/io.h
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

cond: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DSPIN_COND -larb

clean: 
	rm $(PROGRAM)
