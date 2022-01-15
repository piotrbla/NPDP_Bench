// set check to 0 to measure time, set to 1 to check validity

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

int ** c;
int ** ck;

int zz = 2;

int ** F;  //only ACGU

int N;
int DIM;

char * RNA;  //only ACGU


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
    int i,j,k,ll,p,q,l=0;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6,t7,t8,t9,t10;
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


    if(argc > 3)
        kind = atoi(argv[3]);



    omp_set_num_threads(num_proc);
    //printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    F =  mem();
    c = mem();
    ck = mem();

   for(i=0; i<DIM; i++)
    for(i=0; i<DIM; i++){
     c[i][j] = i+j;
     ck[i][j] = i+j;
    }

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);

//for(i=0; i<DIM; i++)
 //printf("%c", RNA[i]);

printf("\n");
     int check=1;


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp

        #pragma scop
        for (i = N-2;  i>=1; i--){
           for ( j=i+2; j<= N; j++){
             for ( k = i; k<=j-l; k++){
              ck[i][j] +=  ck[i][j-1] + paired(k,j) ?  ck[i][k-1] + ck[k+1][j-1] : 0;
           }
          }
        }
       #pragma endscop

 
    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);

    if(check)
    for(i=0; i<DIM; i++)
    for(j=0; j<DIM; j++)
     if(c[i][j] != ck[i][j]){
        printf("err: %d %d %d %d\n", i, j,c[i][j], ck[i][j]);
        exit(0);
     }


    return 0;

}
