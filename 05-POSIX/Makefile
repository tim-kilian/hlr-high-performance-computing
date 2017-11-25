# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -pthread
LFLAGS = $(CFLAGS)
LIBS   = -lm

TGTS = partdiff-seq partdiff-posix
OBJS = partdiff-seq.o askparams.o
POSIXOBJS = partdiff-posix.o askparams.o

# Targets ...
all: partdiff-seq partdiff-posix

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-seq.o: partdiff-seq.c Makefile

partdiff-posix: $(POSIXOBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(POSIXOBJS) $(LIBS)

partdiff-posix.o: partdiff-posix.c Makefile

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(POSIXOBJS)
	$(RM) $(TGTS)
