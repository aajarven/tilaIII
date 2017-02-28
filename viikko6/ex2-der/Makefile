CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o re_der.o
	gcc $(CFLAGS) -o bin/main main.o re_der.o -lm

re_der.o: re_der.c re_der.h
	gcc $(CFLAGS) -lm -c re_der.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
