CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm *.o bin/*

main: main.o jacobi.o lapack.h
	gcc $(CFLAGS) -o bin/main main.o jacobi.o -lm -llapack -lblas

jacobi.o: jacobi.c jacobi.h
	gcc $(CFLAGS) -lm -llapack -lblas -c jacobi.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
