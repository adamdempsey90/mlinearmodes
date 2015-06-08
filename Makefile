EXECUTABLE=eigen
SOURCES=readinputs.c alloc.c output.c boundary.c init.c algo.c derivatives.c main.c matrix.c profiles.c selfgravity.c planets.c
HEADER=eigen.h defines.h

LAPACKLIB=-llapack -lblas
OPENMPLIB=-lgomp
MATHLIB=-lm
GSLLIB=-lgsl -lgslcblas

LDFLAGS=$(GSLLIB) $(LAPACKLIB) $(OPENMPLIB) $(MATHLIB)

CFLAGS=-c -fopenmp -Wall -O3  -g


INCLIB=
LDLIB=


BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/


UNAME := $(shell echo $(USER))

ifeq ($(UNAME),apollo)
CC=gcc
endif

ifeq ($(UNAME),jupiter)
CC=gcc-4.9
endif

ifeq ($(UNAME),helios)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif

ifeq ($(UNAME),amd616)
CC=gcc
LDLIB=-L/software/lapack/3.4.0/lib
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
