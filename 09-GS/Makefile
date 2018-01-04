# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -Ofast
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-mpi.o askparams.o 

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-mpi

partdiff-mpi: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *~ partdiff-mpi

partdiff-mpi.o: partdiff-mpi.c Makefile

askparams.o: askparams.c Makefile


