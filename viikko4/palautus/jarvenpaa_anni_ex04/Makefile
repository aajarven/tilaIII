CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o roots.o
	gcc $(CFLAGS) -o bin/main main.o roots.o -lm -llapack -lblas

harmonic.o: roots.c roots.h lapack.h
	gcc $(CFLAGS) -lm -llapack -lblas -c roots.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
