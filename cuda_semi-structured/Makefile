#Makefile for 2d example
CC = nvcc
CFLAGS = -lm

sphere : sphere.o
	${CC} ${CFLAGS} sphere.o -o sphere

sphere.o : sphere.cu
	${CC} ${CFLAGS} -c sphere.cu -o sphere.o

clean:
	rm *.o sphere
