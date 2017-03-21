CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o forking.o
	gcc $(CFLAGS) -o bin/main main.o forking.o -lm

forking.o: forking.c forking.h
	gcc $(CFLAGS) -lm -c forking.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
