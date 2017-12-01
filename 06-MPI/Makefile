HEADERS =
default: timempi timempi2

timempi.o: timempi.c $(HEADERS)
	mpicc -c timempi.c -o timempi.o

timempi2.o: timempi2.c $(HEADERS)
	mpicc -c timempi2.c -o timempi2.o

timempi: timempi.o
	mpicc timempi.o -o timempi

timempi2: timempi2.o
	mpicc timempi2.o -o timempi2

clean:
	-rm -f timempi.o
	-rm -f timempi
	-rm -f timempi2.o
	-rm -f timempi2

