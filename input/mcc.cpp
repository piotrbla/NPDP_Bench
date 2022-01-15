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

double ** Q;
double ** Q1;
double ** Qbp;
double ** Qbp1;
double ** Pbp;
double ** Pu;
double ** M;
int CHECK_VALID = 0;

int Ebp = 0; // Energy weight of base pair  -2, -1, 0, 1, 2
int RT = 1; // 'Normalized' temperature 1,2,3,4,5
float ERT;
int l = 0; //minimum loop length 0-5
int delta = 1;  // Base pair weighting  1-5

unsigned char * RNA;  //only ACGU

int N;
int DIM;

#include "mem.h"

int paired(int i, int j) {
   char nt1 = RNA[i];
   char nt2 = RNA[j];
         if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
             (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
             (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
         else
            return 0;
}


int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    ERT = exp((float)-Ebp/(float)RT);
//    printf("ERT: %.f", ERT);
 
    srand(time(NULL));



    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;


    if(argc > 3)
        kind = atoi(argv[3]);



    omp_set_num_threads(num_proc);
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    RNA =  (unsigned char*) malloc(DIM * sizeof(unsigned char*));  //read from FASTA file
    rand_seq(RNA, N);

  

//    printf("Sequence: ");
//    for(i=0; i<N; i++)
//       printf("%c", RNA[i]);
//    printf("\n\n");




    Q = memd();
    Q1 = memd();
    Qbp = memd();
    Qbp1 = memd();
    Pbp = memd();
    Pu = memd();
    M = memd();

    rna_array_init(Q, 0.4, 0.4);
    rna_array_init(Q1, 0.4, 0.4);
    rna_array_init(Qbp, 0.5, 0.5);
    rna_array_init(Qbp1, 0.5, 0.5);
    rna_array_init(Pbp, 0, 0);
    rna_array_init(Pu, 0, 0);
    rna_array_init(M, 0, 0);


//for(i=0; i<N; i++)
  //for(j=0; j<N; j++)
    // printf("%.6f %.6f -- %d %d\n", Q[i][j], Q1[i][j], i, j);



    double start = omp_get_wtime();

        #pragma scop
        if(N>=1 && l>=0 && l<=5)
        for(i=N-1; i>=0; i--){
         for(j=i+1; j<N; j++){
            Q1[i][j] =  Q1[i][j-1];
           for(k=0; k<j-i-l; k++){
             Qbp1[k+i][j] = Q1[k+i+1][j-1] * ERT * paired(k+i,j-1);
             Q1[i][j] +=  Q1[i][k+i] * Qbp1[k+i][j];
           }

         }
        }
       #pragma endscop

    


    double stop = omp_get_wtime();



    return 0;

}
