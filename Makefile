EXECUTABLE=eigen
SOURCES=alloc.c output.c boundary.c init.c algo.c derivatives.c main.c matrix.c profiles.c selfgravity.c 
HEADER=eigen.h defines.h

LDFLAGS=-llapack -lblas -lm -lgomp -lgsl

CFLAGS=-c -fopenmp -Wall -O3  -g

INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib

BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/


UNAME := $(shell uname)

ifeq ($(UNAME),Linux)
CC=gcc
endif

ifeq ($(UNAME),Darwin)
CC=gcc-4.9
endif

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(COBJECTS) 
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER) 
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@

	
clean:
	rm $(COBJECTS) $(EXECUTABLE) 
