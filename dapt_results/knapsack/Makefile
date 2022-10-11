CC=gcc
LD=g++
CPP=g++

all:	tknapsack

tknapsack:	knapsack.o
	$(LD) -std=c++17 -O3 knapsack.o -o tknapsack -fopenmp

knapsack.o:	knapsack.cc
	$(CPP) -c knapsack.cc -O3 -fopenmp

clean:
	rm -f *.o tknapsack *~
