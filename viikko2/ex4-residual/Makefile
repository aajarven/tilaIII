CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o harmonic_kahan.o
	gcc $(CFLAGS) -o bin/main main.o harmonic_kahan.o

harmonic.o: harmonic_kahan.c harmonic_kahan.h
	gcc $(CFLAGS) -c harmonic_kahan.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
