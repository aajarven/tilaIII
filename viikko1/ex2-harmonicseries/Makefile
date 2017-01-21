CFLAGS = -g -Wall -Wextra -pedantic

all: main

clean: 
	rm main *.o

main: main.o harmonic.o
	gcc $(CFLAGS) -o main main.o harmonic.o

harmonic.o: harmonic.c harmonic.h
	gcc $(CFLAGS) -c harmonic.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
