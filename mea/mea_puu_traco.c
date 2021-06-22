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
double ** Qbp;
double ** Pbp;
double ** Pu;
double ** M;


int Ebp = 0; // Energy weight of base pair  -2, -1, 0, 1, 2
int RT = 1; // 'Normalized' temperature 1,2,3,4,5
float ERT;
int l = 0; //minimum loop length 0-5
int delta = 1;  // Base pair weighting  1-5

char * RNA;  //only ACGU

int N;
int DIM;

#include "../mem.h"

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
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c15;

    int t1, t2, t3, t4, t5, t6,t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    ERT = exp((float)-Ebp/(float)RT);


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

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);


    //printf("Sequence: ");
    //for(i=0; i<N; i++)
    //   printf("%c", RNA[i]);
    //printf("\n\n");




    Q = memd();
    Qbp = memd();
    Pbp = memd();
    Pu = memd();
    M = memd();

    rna_array_init(Q, 1, 1);
    rna_array_init(Qbp, 0, 0);
    rna_array_init(Pbp, 0, 0);
    rna_array_init(Pu, 0, 0);
    rna_array_init(M, 0, 0);
    double * Puu = (double*)malloc(DIM * sizeof(double));


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1){
    #pragma scop
    for(i=0; i<=N; i++){
     Puu[i] = 1;
     for(j=i+1; j<N; j++){
       Puu[i] += -1 * Pbp[i][j+1];
     }
     for(k=0; k<i; k++){
       Puu[i] += -1 * Pbp[k][i+1];
     }
    }
    #pragma endscop

    }
    if(kind==2) // pluto
    {
        printf("pluto\n");

/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if (N >= 0) {
  lbp=0;
  ubp=floord(N,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
  for (t2=lbp;t2<=ubp;t2++) {
    lbv=16*t2;
    ubv=min(N,16*t2+15);
#pragma ivdep
#pragma vector always
    for (t3=lbv;t3<=ubv;t3++) {
      Puu[t3] = 1;;
    }
  }
  if (N >= 1) {
    lbp=0;
    ubp=floord(N,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
    for (t2=lbp;t2<=ubp;t2++) {
      if (t2 <= floord(N-2,16)) {
        for (t4=t2;t4<=floord(N-1,16);t4++) {
          for (t5=16*t2;t5<=min(min(N-2,16*t2+15),16*t4+14);t5++) {
            for (t7=max(16*t4,t5+1);t7<=min(N-1,16*t4+15);t7++) {
              Puu[t5] += -1 * Pbp[t5][t7+1];;
            }
          }
        }
      }
      for (t4=0;t4<=min(floord(N-1,16),t2);t4++) {
        for (t5=max(16*t2,16*t4+1);t5<=min(N,16*t2+15);t5++) {
          for (t7=16*t4;t7<=min(16*t4+15,t5-1);t7++) {
            Puu[t5] += -1 * Pbp[t7][t5+1];;
          }
        }
      }
    }
  }
}
/* End of CLooG code */





    }
    if(kind==3) // traco
    {
        printf("traco\n");


    }
    if(kind==4) // traco tstile
    {
        printf("traco cor\n");

        #pragma omp parallel for
for( c1 = 0; c1 <= N/16; c1 += 1) {
  for( c2 = 0; c2 <= min(1, -8 * c1 + N / 2); c2 += 1) {
    if (c2 == 1) {
      for( c3 = 0; c3 <= -c1 + (N - 2) / 16; c3 += 1)
        for( c5 = 16 * c1; c5 <= min(16 * c1 + 15, N - 16 * c3 - 2); c5 += 1)
          for( c7 = 16 * c3 + c5 + 1; c7 <= min(N - 1, 16 * c3 + c5 + 16); c7 += 1)
            Puu[c5] += -1 * Pbp[c5][c7+1];
    } else {
      for( c5 = 16 * c1; c5 <= min(N, 16 * c1 + 15); c5 += 1)
        Puu[c5] = 1;
    }
  }
  for( c3 = 0; c3 <= min(c1, floord(N - 1, 16)); c3 += 1)
    for( c5 = max(16 * c1, 16 * c3 + 1); c5 <= min(N, 16 * c1 + 15); c5 += 1)
      for( c7 = 16 * c3; c7 <= min(16 * c3 + 15, c5 - 1); c7 += 1)
        Puu[c5] += -1 * Pbp[k][c5+1];
}




    }

    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);

    exit(0);



    return 0;

}
