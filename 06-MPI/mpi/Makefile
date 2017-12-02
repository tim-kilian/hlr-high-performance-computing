CC = mpicc

CFLAGS = -Wall -g

TARGETS = timempi timempi2
OBJS = timempi.o timempi2.o

all: timempi timempi2


timempi: $(OBJS)
	$(CC) $(CFLAGS) -o $@ timempi.o

timempi.o: timempi.c Makefile

timempi2: $(OBJS)
	$(CC) $(CFLAGS) -o $@ timempi2.o

timempi2.o: timempi2.c Makefile

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(TARGETS)
	$(RM) $(OBJS)

