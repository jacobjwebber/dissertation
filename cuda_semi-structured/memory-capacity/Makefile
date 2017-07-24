#Makefile for 2d example

IDIR =../lib
CC=nvcc
CFLAGS=-I$(IDIR) -arch=sm_35 -g -G 

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = semi-structured_lib.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = memory-capacity.o semi-structured_lib.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cu $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: ../lib/%.cu $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)


memory-capacity : $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ sphere
