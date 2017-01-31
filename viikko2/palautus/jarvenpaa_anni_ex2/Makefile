CFLAGS = -g -Wall -Wextra -pedantic -std=c89

all: main

clean: 
	rm *.o bin/*

main: main.o residual.o
	gcc $(CFLAGS) -o bin/main main.o residual.o -lm -llapack -lblas

residual.o: residual.c residual.h
	gcc $(CFLAGS) -lm -llapack -lblas -c residual.c 

main.o: main.c lapac.h
	gcc $(CFLAGS) -c main.c

