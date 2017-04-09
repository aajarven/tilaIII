CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o boundatyproblem.o
	gcc $(CFLAGS) -o bin/main main.o boundatyproblem.o -lm

boundatyproblem.o: boundatyproblem.c boundatyproblem.h
	gcc $(CFLAGS) -lm -c boundatyproblem.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
