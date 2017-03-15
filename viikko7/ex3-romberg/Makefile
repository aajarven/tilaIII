CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o romberg.o
	gcc $(CFLAGS) -o bin/main main.o romberg.o -lm

romberg.o: romberg.c romberg.h
	gcc $(CFLAGS) -lm -c romberg.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
