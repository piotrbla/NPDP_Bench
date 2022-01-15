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

int N = 8, DIM = 12;

#include "../mem.h"

int paired(i, j) {
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

    ERT = exp((float)-Ebp/(float)RT);



    srand(time(NULL));


    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    if(argc > 2)
        N = atoi(argv[2]);
    DIM = 2*N+2;


    if(argc > 3)
        kind = atoi(argv[3]);


    printf(" -exp(Ebp/RT) = %5.3f\n", ERT);

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);

    RNA = "GGUCCAC";


    printf("Sequence: ");
    for(i=0; i<N; i++)
       printf("%c", RNA[i]);
    printf("\n\n");


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



    //  compute the partition functions Q and Qbp

    #pragma scop
    if(N>=1 && l>=0 && l<=5)
    for(i=N-1; i>=0; i--){
     for(j=i+1; j<N; j++){
        Q[i][j] =  Q[i][j-1];
       for(k=0; k<j-i-l; k++){
         Qbp[k+i][j] = Q[k+i+1][j-1] * ERT * paired(k+i,j-1);
         Q[i][j] +=  Q[i][k+i] * Qbp[k+i][j];
       }

     }
    }
   #pragma endscop

    printf("Q\n");
    rna_array_print(Q);
    printf("Qbp\n");
    rna_array_print(Qbp);






    #pragma scop
    for(i=0; i<N; i++){
     for(j=i+1; j<N; j++){
       Pbp[i][j] = (Q[0][i]*Q[j][N-1]*Qbp[i][j])/Q[0][N-1];   //  Pbp[i][j] = (Q[1][i]*Q[j+1][N]*Qbp[i][j])/Q[0][N-1];
       for(p=0; p<i; p++){
        for(q=j+1; q<N; q++){
         Pbp[i][j] += (Pbp[p][q] * ERT * Q[p+1][i] * Qbp[i][j] * Q[j+1][q-1]) / (Qbp[p][q] ==0 ? 1 : Qbp[p][q]) ;

        }
      }
     }
    }
    #pragma endscop


   printf("Pbp\n");
    rna_array_print(Pbp);


    #pragma scop
      for(i=N-1; i>=0; i--){
     for(j=i+1; j<N; j++){
       Pu[i][j] = (Q[0][i]*Q[j][N-1]*1)/Q[0][N-1];
       for(p=0; p<i; p++){
        for(q=j+1; q<N; q++){
         Pu[i][j] += (Pbp[p][q] * ERT * Q[p+1][i] * 1 * Q[j+1][q-1]) /  (Qbp[p][q] ==0 ? 1 : Qbp[p][q]) ;
        }
      }
     }
    }
    #pragma endscop

    printf("Pu\n");
    rna_array_print(Pu);


    double * Puu = (double*)malloc(DIM * sizeof(double));


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

    printf("Puu\n");
    for(i=0; i<N-1; i++)
        printf("%3.3f ", Puu[i]);
    printf("\n");

    #pragma scop
    for(i=N-1; i>=0; i--){
     for(j=i+1; j<N; j++){
       for(k=0; k<j-i-l; k++){
        M[i][j] = MAX(M[i][j], M[i][k+i-1] + M[k+i+1][j-1] + delta*Pbp[k+i][j])*paired(k+i,j-1);
      }
      M[i][j] = MAX(M[i][j], M[i][j-1] + Puu[j-1]);
     }
    }
    #pragma endscop

    printf("M\n");
    rna_array_print(M);

    return 0;

}
