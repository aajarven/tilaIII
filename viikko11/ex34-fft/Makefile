CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o fft.o
	gcc $(CFLAGS) -o bin/main main.o fft.o -lm

fft.o: fft.c fft.h
	gcc $(CFLAGS) -lm -c fft.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
