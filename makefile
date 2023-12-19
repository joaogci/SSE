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
		CC = gcc-13
	endif
endif

$(PROGRAM): $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS)

test: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) -Wall -Wshadow $(SOURCES) $(LIBS) -DTESTING

# "full", "ss", "hh", "diag", "offd", "sssh", "hhsh"
full: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_SS -DL_HH -DL_SH -DL_HS

ss: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_SS

hh: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_HH

diag: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_SS -DL_HH

offd: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_SH -DL_HS

sssh: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_SS -DL_SH -DL_HS 

hhsh: $(SOURCES) $(INCLUDES)
	$(CC) -o $(PROGRAM) $(FLAGS) $(SOURCES) $(LIBS) -DKINETIC -DL_HH -DL_SH -DL_HS

clean: 
	rm $(PROGRAM)
