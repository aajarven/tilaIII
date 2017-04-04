CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o odes.o
	gcc $(CFLAGS) -o bin/main main.o odes.o -lm

odes.o: odes.c odes.h
	gcc $(CFLAGS) -lm -c odes.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
