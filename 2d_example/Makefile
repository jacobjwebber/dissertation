#Makefile for 2d example
CC = gcc
CFLAGS = -lm

square : rotating_square.o
	${CC} ${CFLAGS} rotating_square.o -o square

rotating_square.o : rotating_square.c
	${CC} ${CFLAGS} -c rotating_square.c -o rotating_square.o

clean:
	rm *.o square
