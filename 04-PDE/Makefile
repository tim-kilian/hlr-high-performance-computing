# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-seq partdiff-openmp partdiff-element partdiff-spalten partdiff-zeilen
OBJS = partdiff-seq.o askparams.o
OMPOBJS = partdiff-openmp.o askparams.o
ELEOBJS = partdiff-element.o askparams.o
SPAOBJS = partdiff-spalten.o askparams.o
ZEIOBJS = partdiff-zeilen.o askparams.o

# Targets ...
all: partdiff-seq partdiff-openmp partdiff-element partdiff-zeilen partdiff-spalten

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-seq.o: partdiff-seq.c Makefile

partdiff-openmp: $(OMPOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OMPOBJS) $(LIBS)

partdiff-openmp.o: partdiff-seq.c Makefile

partdiff-element: $(ELEOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(ELEOBJS) $(LIBS)

partdiff-element.o: partdiff-seq.c Makefile

partdiff-spalten: $(SPAOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(SPAOBJS) $(LIBS)

partdiff-spalten.o: partdiff-seq.c Makefile

partdiff-zeilen: $(ZEIOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(ZEIOBJS) $(LIBS)

partdiff-zeilen.o: partdiff-seq.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(OMPOBJS)
	$(RM) $(ELEOBJS)
	$(RM) $(SPAOBJS)
	$(RM) $(ZEIOBJS)
	$(RM) $(TGTS)
