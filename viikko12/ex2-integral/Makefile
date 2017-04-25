CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o integration.o
	gcc $(CFLAGS) -o bin/main main.o integration.o -lm

integration.o: integration.c integration.h
	gcc $(CFLAGS) -lm -c integration.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
