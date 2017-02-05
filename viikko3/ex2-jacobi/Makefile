CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o jacobi.o
	gcc $(CFLAGS) -o bin/main main.o jacobi.o

jacobi.o: harmonic.c harmonic.h
	gcc $(CFLAGS) -c jacobi.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
