CC = gcc
EXEC = map pointer
OBJ = map.o pointer.o

all: $(EXEC)

%: %.o
	$(CC) -o $@ $<

%.o: %c
	$(CC) -c $<

clean:
	rm -f $(EXEC)
	rm -f *.o
	rm -f *~
