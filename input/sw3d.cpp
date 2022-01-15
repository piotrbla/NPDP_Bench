#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>


#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))
#define CHECK_VALID 0


int ***H;
int ***H1;

int ***H0;
int ***H2;
int ***H3;
int ***H4;


int ***tmp_H;
int ***m1;
int ***m2;
int ***m3;
int ***m4;
int ***m5;
int ***m6;

int *W;
unsigned char *a;
unsigned char *b;
unsigned char *c;


int N = 100, DIM = 102;

//Similarity score of the elements that constituted the two sequences
int s(unsigned char x, unsigned char z){
   return (x == z) ? 3 : -3;
}




#include "mem3d.h"
#include "sw3d_oryg.h"
#include "sw3d_traco16.h"
#include "sw3d_traco128.h"
#include "sw3d_tstile.h"
#include "sw3d_pluto.h"



int main(int argc, char *argv[]){

    int i,j,k,l;


    int num_proc=1;
    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+2;


    if(argc > 3)
        kind = atoi(argv[3]);

    // H is the scoring matrix
    H = mem();
    H1 = mem();
    
    H0 = mem();
    H2 = mem();
    H3 = mem();
    H4 = mem();    
    
    
    m1 = mem();
    m2 = mem();
    m3 = mem();
    m4 = mem();
    m5 = mem();
    m6 = mem();

    W = (int*)malloc(DIM * sizeof(int));
    a = (unsigned char *)malloc(DIM * sizeof(unsigned char ));
    b = (unsigned char *)malloc(DIM * sizeof(unsigned char ));
    c = (unsigned char *)malloc(DIM * sizeof(unsigned char ));



    for(i=0; i<=N; i++){
        H[i][0][0] = 0;
        H1[i][0][0]= 0;
        H[0][i][0] = 0;
        H1[0][i][0] = 0;
        H[0][0][i] = 0;
        H1[0][0][i] = 0;

    }


    // W is the gap alignment
    W[0] = 2;
    for(i=1; i<=N; i++)
        W[i] = i*W[0];

    rand_seq(a, N);
    rand_seq(b, N);
    rand_seq(c, N);



    omp_set_num_threads(num_proc);

    double start = omp_get_wtime();

    if(kind == 1)
        sw_seq();



    if(CHECK_VALID)
    {    tmp_H = H;
        H = H1;
        sw_seq();
        H = tmp_H;


     for(i=0; i<N; i++)
      for(j=0; j<N; j++)
        for(l=0; l<N; l++)
       if(H[i][j][l] != H1[i][j][l]){
          printf("error!\n");
          printf("%d\n", H[i][j][l]);
          printf("%d\n", H1[i][j][l]);
          printf("%i %i %i \n", i,j,l);
          exit(0);
      }

    }
    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);




return 0;
}
