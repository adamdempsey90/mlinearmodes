EXECUTABLE=eigen
SOURCES=eigenvalues.c
HEADER=eigen.h defines.h

LDFLAGS=-llapack -lblas -lm -lgomp -lgsl

CFLAGS=-c -fopenmp -Wall -O3 -g 


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




all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC)  $(OBJECTS) $(LDFLAGS) -o $@

%.o: %.c $(HEADER) 
	$(CC) $(CFLAGS) $< -o $@

	
clean:
	rm $(OBJECTS) $(EXECUTABLE) 
