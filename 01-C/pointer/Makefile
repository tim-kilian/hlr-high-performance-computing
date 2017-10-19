all: pointer
	
pointer: pointer.o
	gcc -o pointer pointer.o
	
pointer.o: pointer.c
	gcc -c pointer.c

clean:
	rm -f pointer
	rm -f *.o
	rm -f *~
