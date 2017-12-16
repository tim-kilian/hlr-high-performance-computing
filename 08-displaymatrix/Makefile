CC = mpicc

CFLAGS = -Wall -Wextra -fopenmp
LIBS   = -lm

TARGETS = partdiff-seq partdiff-mpi
SEQOBJS = partdiff-seq.o askparams.o
MPIOBJS = partdiff-mpi.o askparams.o

all: partdiff-seq partdiff-mpi

partdiff-seq: $(SEQOBJS) Makefile
	$(CC) $(CFLAGS) -o $@ $(SEQOBJS) $(LIBS)

partdiff-seq.o: partdiff-seq.c Makefile

partdiff-mpi: $(MPIOBJS) Makefile
	$(CC) $(CFLAGS) -o $@ $(MPIOBJS) $(LIBS)

partdiff-mpi.o: partdiff-mpi.c Makefile

askparams.o: askparams.c Makefile

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(TARGETS)
	$(RM) $(SEQOBJS)
	$(RM) $(MPIOBJS)
