CC=gcc

all:	tcat

tcat:	nth_catlan.o
	$(CC) -O3 nth_catlan.o -o tcat -fopenmp

nth_catlan.o:	nth_catlan.c
	$(CC) -c nth_catlan.c -O3 -fopenmp

clean:
	rm -f *.o tmcc *~
