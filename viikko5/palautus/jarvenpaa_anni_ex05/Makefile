CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o broyden.o
	gcc $(CFLAGS) -o bin/main main.o broyden.o -lm -llapack -lblas

harmonic.o: broyden.c lapack.h
	gcc $(CFLAGS) -lm -llapack -lblas -c broyden.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
