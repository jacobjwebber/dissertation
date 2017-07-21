#Makefile for 2d example
CC = nvcc
CFLAGS = -lm -arch=sm_35 -g -G -I.
DEPS = semi-structured_lib.h

sphere : sphere.o semi-structured_lib.o
	${CC} ${CFLAGS} semi-structured_lib.o sphere.o -o sphere

%.o: %.cu $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm *.o sphere
