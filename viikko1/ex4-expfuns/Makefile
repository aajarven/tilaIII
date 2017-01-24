CFLAGS = -g -Wall -Wextra -pedantic -std=c99

all: main

clean: 
	rm bin/* *.o

main: main.o exp_funs.o
	gcc $(CFLAGS) -o bin/main main.o exp_funs.o -lm

exp_funs.o: exp_funs.c exp_funs.h
	gcc $(CFLAGS) -c exp_funs.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
