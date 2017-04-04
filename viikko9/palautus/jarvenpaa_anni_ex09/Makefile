CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o correlation.o
	gcc $(CFLAGS) -o bin/main main.o correlation.o -lm

correlation.o: correlation.c correlation.h
	gcc $(CFLAGS) -lm -c correlation.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
