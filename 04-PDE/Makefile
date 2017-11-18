# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-seq partdiff-openmp
OBJS = partdiff-seq.o askparams.o
OMPOBJS = partdiff-openmp.o askparams.o

# Targets ...
all: partdiff-seq partdiff-openmp

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-seq.o: partdiff-seq.c Makefile

partdiff-openmp: $(OMPOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OMPOBJS) $(LIBS)

partdiff-openmp.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(OMPOBJS)
	$(RM) $(TGTS)
