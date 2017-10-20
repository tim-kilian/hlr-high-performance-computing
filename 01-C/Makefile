objects = map.o pointer.o
executable = map pointer

all: $(executable)

%: %.o
	gcc -o $@ $<

%.o: %.c
	gcc -c $<



clean:
		rm -f *.o
		rm -f $(executable)