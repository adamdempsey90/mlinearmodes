EXECUTABLE=eigen
SOURCES=readinputs.c alloc.c output.c boundary.c init.c algo.c derivatives.c main.c matrix.c profiles.c selfgravity.c planets.c hdf5.c
HEADER=eigen.h defines.h

LAPACKLIB=-llapack -lblas
OPENMPLIB=-lgomp
MATHLIB=-lm
GSLLIB=-lgsl -lgslcblas
HDF5LIB=-lhdf5

LDFLAGS=$(GSLLIB) $(LAPACKLIB) $(OPENMPLIB) $(MATHLIB) $(HDF5LIB)

CFLAGS=-c -fopenmp -Wall -O3  -DH5_USE_16_API -g


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
ifeq ($(UNAME),zeus)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif


ifeq ($(UNAME),helios)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif

ifeq ($(UNAME),amd616)
CC=gcc
LDLIB=-L/software/lapack/3.4.0/lib -L/software/gsl/1.16-gcc4.8.3/lib/ -L/software/hdf5/1.8.12-serial/lib/
INCLIB=-I/software/hdf5/1.8.12-serial/include/
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
