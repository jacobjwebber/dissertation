#Makefile for 2d example
CC = nvcc
CFLAGS = -lm -I/opt/cuda7/samples/common/inc -g -G


sphere : sphere.o
	${CC} ${CFLAGS} sphere.o -o sphere

sphere.o : sphere.cu
	${CC} ${CFLAGS} -c sphere.cu -o sphere.o

clean:
	rm *.o sphere
