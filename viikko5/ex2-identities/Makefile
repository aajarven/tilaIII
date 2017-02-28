CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o myexp.o
	gcc $(CFLAGS) -o bin/main main.o myexp.o -lm -llapack -lblas

harmonic.o: myexp.c roots.h lapack.h
	gcc $(CFLAGS) -lm -llapack -lblas -c myexp.c 

main.o: main.c
	gcc $(CFLAGS) -c main.c
