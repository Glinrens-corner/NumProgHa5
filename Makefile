CC = gcc
RM = rm -f
# RM = del

CFLAGS = -Wall \
	 -O3 -funroll-loops -march=native -mtune=native \
         -std=c99 -g
LIBS   = -lm
# LIBS   = -lm -lwinmm
         
CFLAGS += -DUSE_OPENMP -fopenmp

#CFLAGS += -DUSE_BLAS
#LIBS   += -lblas -lgfortran -llapack

PROGRAMS = p5_newton

all: $(PROGRAMS)

basic.o: basic.c basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
miniblas.o: miniblas.c miniblas.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
matrix.o: matrix.c matrix.h miniblas.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
factorizations.o: factorizations.c factorizations.h matrix.h miniblas.h basic.h
	$(CC) -c $(CFLAGS) $< -o $@ 
	
$(PROGRAMS): %: %.c basic.o miniblas.o matrix.o factorizations.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	$(RM) *.o
	$(RM) $(PROGRAMS)


