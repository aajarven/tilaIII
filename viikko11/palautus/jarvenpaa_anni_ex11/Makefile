CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o boundaryproblem.o
	gcc $(CFLAGS) -o bin/main main.o boundaryproblem.o -lm

boundaryproblem.o: boundaryproblem.c boundaryproblem.h
	gcc $(CFLAGS) -lm -c boundaryproblem.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
