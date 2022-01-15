#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
//#include<cmath>

#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

#define CHECK_VALID 0

int **W;
int **V;
int **V1;
int **tmp_V;
int **W1;
int **tmp_W;
int **EFL;
int **EHF;
int **VMI;
int **VBI;
int **WZ;

int N = 100, DIM = 102;

void zuker_seq()
{
int i,j,k,m;

#pragma scop
for (i = N-1; i >= 0; i--){
 for (j = i+1; j < N; j++) {
   for (k = i+1; k < j; k++){
   for(m=k+1; m <j; m++){
    if(k-i + j - m > 2 && k-i + j - m < 30)
       V[i][j] = MIN(V[k][m] + EFL[i][j], V[i][j]);
   }
   W[i][j] += MIN ( MIN(W[i][k], W[k+1][j]), W[i][j]);
   if(k < j-1)
     V[i][j] = MIN(W[i+1][k] + W[k+1][j-1], V[i][j]);
  }
 V[i][j] = MIN( MIN (V[i+1][j-1], EHF[i][j]), V[i][j]);
 W[i][j] = MIN( MIN ( MIN ( W[i+1][j], W[i][j-1]), V[i][j]), W[i][j]);
 }
}
#pragma endscop

}




int main(int argc, char *argv[]){

    int i,j,k;
//

    int num_proc=1;
    num_proc = atoi(argv[1]);

    int kind;

    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+2;


    if(argc > 3)
        kind = atoi(argv[3]);

    W = mem();
    V = mem();
    V1 = mem();
    W1 = mem();
    EFL = mem();
    EHF = mem();
    VMI = mem();
    VBI = mem();
    WZ = mem();

    for(i=0; i<N; i++)
     for(j=0; j<N; j++){
     W[i][j] = i*j;
     V[i][j] = i+1;
     EHF[i][j] = i+1;
     EFL[i][j] = i+1;
     V1[i][j] = V[i][j];
     W1[i][j] = W[i][j];
    }
    omp_set_num_threads(num_proc);

    double start = omp_get_wtime();

    zuker_seq();

    double stop = omp_get_wtime();
    printf("%.2f seconds\n",stop - start);

return 0;
}
