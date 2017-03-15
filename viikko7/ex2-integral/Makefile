CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o integral.o
	gcc $(CFLAGS) -o bin/main main.o integral.o -lm

integral.o: integral.c integral.h
	gcc $(CFLAGS) -lm -c integral.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
