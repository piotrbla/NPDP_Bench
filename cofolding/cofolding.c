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



int ** F;  //only ACGU

int N;
int DIM;

char * RNA;  //only ACGU


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
    int i,j,k,ll,p,q,l=0;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6,t7,t8,t9,t10;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    rand_seq(RNA, N);


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





    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1){
        #pragma scop
        for (i = N-2;  i>=1; i--){
           for ( j=i+2; j<= N; j++){
             for ( k = i; k<=j-l; k++){
              c[i][j] +=  c[i][j-1] + paired(k,j) ?  c[i][k-1] + c[k+1][j-1] : 0;
           }
          }
        }
       #pragma endscop
    }
    if(kind==3) // pluto
    {
        /* Start of CLooG code */
/* Start of CLooG code */
if ((N >= 3) && (N >= l+1)) {
  for (t1=max(3,l+1);t1<=N;t1++) {
    lbp=0;
    ubp=min(floord(t1-2,16),floord(t1-l,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=t2;t3<=floord(t1-l,16);t3++) {
        for (t4=max(1,16*t2);t4<=min(min(t1-2,t1-l),16*t2+15);t4++) {
          for (t5=max(16*t3,t4);t5<=min(t1-l,16*t3+15);t5++) {
            c[t4][t1] += c[t4][t1-1] + paired(t5,t1) - c[t4][t5-1] + c[t5+1][t1-1] * 0;;
          }
        }
      }
    }
  }
}
/* End of CLooG code */


        /* End of CLooG code */
    }
    if(kind==2) // traco
    {

    }

    if(kind==4) // traco
    {
    for( c0 = max(0, l + floord(l - 2, 16) - 2); c0 < N + floord(N - 3, 16) - 2; c0 += 1)
  #pragma omp parallel for
  for( c1 = c0 - (c0 + 17) / 17 + 1; c1 <= min(min(N - 3, c0), c0 + (-l + 1)/16 + 1); c1 += 1)
    for( c3 = max(l, 16 * c0 - 16 * c1 + 2); c3 <= min(c1 + 2, 16 * c0 - 16 * c1 + 17); c3 += 1)
      for( c4 = (N - c1 - 2) / 16; c4 <= (-l + N - c1 + c3 - 2) / 16; c4 += 1)
        for( c10 = max(N - c1 - 2, 16 * c4); c10 <= min(-l + N - c1 + c3 - 2, 16 * c4 + 15); c10 += 1)
          c[(N-c1-2)][(N-c1+c3-2)] += c[(N-c1-2)][(N-c1+c3-2)-1] + paired(c10,(N-c1+c3-2)) - c[(N-c1-2)][c10-1] + c[c10+1][(N-c1+c3-2)-1] * 0;
    }





    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);


    return 0;

}
