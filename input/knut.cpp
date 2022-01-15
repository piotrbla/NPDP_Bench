// set check to 0 to measure time, set to 1 to check validity

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

int **c;
int **ck;
int **w;

int N;
int DIM;

#include "mem.h"

int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    srand(time(NULL));


    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;
    int n = N;
int check = 0;
    if(argc > 3)
        kind = atoi(argv[3]);


    omp_set_num_threads(num_proc);


    c = mem();
    ck = mem();
    w = mem();

    for(i=0; i<DIM; i++)
        for(j=0; j<DIM; j++){
          ck[i][j] = i+j;
          c[i][j] = ck[i][j];
          w[i][j] = i-j;
    }


    double start = omp_get_wtime();

        #pragma scop
        for(i=n-1; i>=1; i--)
           for(j=i+1; j<=n; j++)
               for(k=i+1; k<j; k++)
                  ck[i][j] = MIN(ck[i][j], w[i][j]+ck[i][k]+ck[k][j]);
        #pragma endscop
 

  
    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);




    return 0;
}
