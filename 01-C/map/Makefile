all: map
	
map: map.o
	gcc -o map map.o
	
map.o: map.c
	gcc -c map.c

clean:
	rm -f map
	rm -f *.o
	rm -f *~
