SOURCES = main.c sse/sse.c sampling/sampling.c io/io.c sampling/sampling_memory.c 
INCLUDES = sse/sse.h sampling/sampling.h vtx/vtx_type.h rng/xorshiro256++.h io/io.h
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
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DCOND

clean: 
	rm $(PROGRAM)
