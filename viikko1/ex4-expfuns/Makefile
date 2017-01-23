CFLAGS = -g -Wall -Wextra -pedantic

all: main

clean: 
	rm main *.o

main: main.o exp_funs.o
	gcc $(CFLAGS) -o main main.o harmonic.o

exp_funs.o: exp_funs.c exp_funs.h
	gcc $(CFLAGS) -c exp_funs.c

main.o: main.c
	gcc $(CFLAGS) -c main.c
