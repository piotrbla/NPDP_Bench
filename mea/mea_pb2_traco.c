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


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1){
     #pragma scop
     for(i=0; i<N; i++){
     for(j=i+1; j<N; j++){
       for(p=0; p<i; p++){
        for(q=j+1; q<N; q++){
         Pbp[i][j] += (Pbp[p][q] * ERT * Q[p+1][i] * Qbp[i][j] * Q[j+1][q-1]) / (Qbp[p][q] ==0 ? 1 : Qbp[p][q]);

        }
      }
      Pbp[i][j] += (Q[0][i]*Q[j][N-1]*Qbp[i][j])/Q[0][N-1];
     }
    }
    #pragma endscop

    }
    if(kind==2) // pluto
    {
        printf("pluto\n");
     /* Start of CLooG code   */
if(1==0)
if (N >= 2) {
  for (t1=0;t1<=N-2;t1++) {
    lbp=t1+1;
    ubp=N-1;
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      if (t2 <= N-2) {
        for (t3=0;t3<=floord(t1-1,16);t3++) {
          for (t5=16*t3;t5<=min(t1-1,16*t3+15);t5++) {
            for (t6=t2+1;t6<=N-1;t6++) {
              Pbp[t1][t2] += (Pbp[t5][t6] * ERT * Q[t5+1][t1] * Qbp[t1][t2] * Q[t2+1][t6-1]) / (Qbp[t5][t6] ==0 ? 1 : Qbp[t5][t6]);;
            }
          }
        }
      }
      t3 = floord(t1,16);
      Pbp[t1][t2] += (Q[0][t1]*Q[t2][N-1]*Qbp[t1][t2])/Q[0][N-1];;
    }
  }
}

/* Start of CLooG code  */
if(1==1)
 if (N >= 2) {
   for (t1=0;t1<=N-2;t1++) {
       lbp=ceild(t1-14,16);
           ubp=floord(N-1,16);
           #pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
               for (t2=lbp;t2<=ubp;t2++) {
                     for (t3=0;t3<=floord(t1,16);t3++) {
                             if ((t1 >= 16*t3+1) && (t1 <= 16*t3+15)) {
                                       for (t4=max(16*t2,t1+1);t4<=min(N-2,16*t2+15);t4++) {
                                                   for (t5=16*t3;t5<=t1-1;t5++) {
                                                                 for (t6=t4+1;t6<=N-1;t6++) {
                                                                                 Pbp[t1][t4] += (Pbp[t5][t6] * ERT * Q[t5+1][t1] * Qbp[t1][t4] * Q[t4+1][t6-1]) / (Qbp[t5][t6] ==0 ? 1 :  Qbp[t5][t6]);;
                                                                                               }
                                                                                                           }
                                                                                                                       Pbp[t1][t4] += (Q[0][t1]*Q[t4][N-1]*Qbp[t1][t4])/Q[0][N-1];;
                                                                                                                                 }
                                                                                                                                         }
                                                                                                                                                 if (t1 >= 16*t3+16) {
                                                                                                                                                           for (t4=max(16*t2,t1+1);t4<=min(N-2,16*t2+15);t4++) {
                                                                                                                                                                       for (t5=16*t3;t5<=16*t3+15;t5++) {
                                                                                                                                                                                     for (t6=t4+1;t6<=N-1;t6++) {
                                                                                                                                                                                                     Pbp[t1][t4] += (Pbp[t5][t6] * ERT * Q[t5+1][t1] * Qbp[t1][t4] * Q[t4+1][t6-1]) / (Qbp[t5][t6] ==0 ? 1 :  Qbp[t5][t6]);;
                                                                                                                                                                                                                   }
                                                                                                                                                                                                                               }
                                                                                                                                                                                                                                         }
                                                                                                                                                                                                                                                 }
                                                                                                                                                                                                                                                         if ((t1 >= 16*t3+1) && (t1 <= 16*t3+15) && (t2 >= ceild(N-16,16))) {
                                                                                                                                                                                                                                                                   Pbp[t1][(N-1)] += (Q[0][t1]*Q[(N-1)][N-1]*Qbp[t1][(N-1)])/Q[0][N-1];;
                                                                                                                                                                                                                                                                           }
                                                                                                                                                                                                                                                                                   if (t1 == 16*t3) {
                                                                                                                                                                                                                                                                                             for (t4=max(16*t2,t1+1);t4<=min(N-1,16*t2+15);t4++) {
                                                                                                                                                                                                                                                                                                         if (t1%16 == 0) {
                                                                                                                                                                                                                                                                                                                       Pbp[t1][t4] += (Q[0][t1]*Q[t4][N-1]*Qbp[t1][t4])/Q[0][N-1];;
                                                                                                                                                                                                                                                                                                                                   }
                                                                                                                                                                                                                                                                                                                                             }
                                                                                                                                                                                                                                                                                                                                                     }
                                                                                                                                                                                                                                                                                                                                                           }
                                                                                                                                                                                                                                                                                                                                                               }
                                                                                                                                                                                                                                                                                                                                                                 }
                                                                                                                                                                                                                                                                                                                                                                 }
                                                                                                                                                                                                                                                                                                                                                                  /*End of CLooG code */


    }
    if(kind==3) // traco
    {
        printf("traco\n");


  /*  for( c1 = 0; c1 < N - 1; c1 += 1)
  #pragma omp parallel for schedule(dynamic , 1)
  for( c3 = 0; c3 <= (N - c1 - 2) / 16; c3 += 1)
    for( c4 = max(max(0, -c1 + 1), -N + c1 + 16 * c3 + 3); c4 <= 1; c4 += 1) {
      if (c4 == 1) {
        for( c11 = c1 + 16 * c3 + 1; c11 <= min(N - 1, c1 + 16 * c3 + 16); c11 += 1)
          Pbp[c1][c11] += (Q[0][c1]*Q[c11][N-1]*Qbp[c1][c11])/Q[0][N-1];
      } else {
        for( c5 = 0; c5 <= (c1 - 1) / 16; c5 += 1)
          for( c7 = 0; c7 <= -c3 + (N - c1 - 3) / 16; c7 += 1)
            for( c11 = c1 + 16 * c3 + 1; c11 <= min(c1 + 16 * c3 + 16, N - 16 * c7 - 2); c11 += 1) {
              if (N >= 16 * c7 + c11 + 18) {
                for( c15 = 16 * c7 + c11 + 1; c15 <= 16 * c7 + c11 + 16; c15 += 1)
                  Pbp[c1][c11] += (Pbp[16*c5][c15] * ERT * Q[16*c5+1][c1] * Qbp[c1][c11] * Q[c11+1][c15-1]) / (Qbp[16*c5][c15] ==0 ? 1 : Qbp[16*c5][c15]);
              } else {
                for( c13 = 16 * c5; c13 <= min(c1 - 1, 16 * c5 + 15); c13 += 1) {
                  if (c13 >= 16 * c5 + 1)
                    for( c15 = c11 + 1; c15 <= 16 * c7 + c11; c15 += 1)
                      Pbp[c1][c11] += (Pbp[c13][c15] * ERT * Q[c13+1][c1] * Qbp[c1][c11] * Q[c11+1][c15-1]) / (Qbp[c13][c15] ==0 ? 1 :  Qbp[c13][c15]);
                  for( c15 = 16 * c7 + c11 + 1; c15 < N; c15 += 1)
                    Pbp[c1][c11] += (Pbp[c13][c15] * ERT * Q[c13+1][c1] * Qbp[c1][c11] * Q[c11+1][c15-1]) / (Qbp[c13][c15] ==0 ? 1 :  Qbp[c13][c15]);
                }
              }
            }
      }
    }
*/

for( c1 = 0; c1 < N - 1; c1 += 1)
  #pragma omp parallel for schedule(dynamic , 1)
  for( c3 = 0; c3 < N - c1 - 1; c3 += 1)
    for( c4 = max(max(0, -c1 + 1), -N + c1 + c3 + 3); c4 <= 1; c4 += 1) {
      if (c4 == 1) {
        Pbp[c1][(c1+c3+1)] += (Q[0][c1]*Q[(c1+c3+1)][N-1]*Qbp[c1][(c1+c3+1)])/Q[0][N-1];
      } else {
        for( c5 = 0; c5 <= (c1 - 1) / 128; c5 += 1)
          for( c7 = 0; c7 <= (N - c1 - c3 - 3) / 16; c7 += 1) {
            if (N >= c1 + c3 + 16 * c7 + 19) {
              for( c15 = c1 + c3 + 16 * c7 + 2; c15 <= c1 + c3 + 16 * c7 + 17; c15 += 1)
                Pbp[c1][(c1+c3+1)] += (Pbp[128*c5][c15] * ERT * Q[128*c5+1][c1] * Qbp[c1][(c1+c3+1)] * Q[(c1+c3+1)+1][c15-1]) / (Qbp[128*c5][c15] ==0 ? 1 : Qbp[128*c5][c15]);
            } else {
              for( c13 = 128 * c5; c13 <= min(c1 - 1, 128 * c5 + 127); c13 += 1) {
                if (c13 >= 128 * c5 + 1)
                  for( c15 = c1 + c3 + 2; c15 <= c1 + c3 + 16 * c7 + 1; c15 += 1)
                    Pbp[c1][(c1+c3+1)] += (Pbp[c13][c15] * ERT * Q[c13+1][c1] * Qbp[c1][(c1+c3+1)] * Q[(c1+c3+1)+1][c15-1]) / (Qbp[c13][c15] ==0 ? 1 :  Qbp[c13][c15]);
                for( c15 = c1 + c3 + 16 * c7 + 2; c15 < N; c15 += 1)
                  Pbp[c1][(c1+c3+1)] += (Pbp[c13][c15] * ERT * Q[c13+1][c1] * Qbp[c1][(c1+c3+1)] * Q[(c1+c3+1)+1][c15-1]) / (Qbp[c13][c15] ==0 ? 1 : Qbp[c13][c15]);
              }
            }
          }
      }
    }



    }
    if(kind==4) // traco tstile
    {



    }

    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);

    exit(0);


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
