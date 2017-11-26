HEADERS =
default: timempi

timempi.o: timempi.c $(HEADERS)
	mpicc -c timempi.c -o timempi.o

timempi: timempi.o
	mpicc timempi.o -o timempi

clean:
	-rm -f timempi.o
	-rm -f timempi
