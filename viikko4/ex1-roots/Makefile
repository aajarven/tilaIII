CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o roots.o
	gcc $(CFLAGS) -o bin/main main.o roots.o -lm

harmonic.o: roots.c roots.h
	gcc $(CFLAGS) -lm -c roots.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
