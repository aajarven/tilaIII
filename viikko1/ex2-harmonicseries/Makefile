CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o harmonic.o
	gcc $(CFLAGS) -o bin/main main.o harmonic.o

harmonic.o: harmonic.c harmonic.h
	gcc $(CFLAGS) -c harmonic.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
