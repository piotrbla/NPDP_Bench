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

#include "zuker_oryg.h"
#include "zuker_pluto.h"
#include "zuker_traco.h"
#include "zuker_traco3.h"
#include "mem.h"




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

    if(kind == 4){
        printf("traco\n");
        zuker_traco3();
    }


    if(kind == 3){
        printf("traco tilecorr\n");
        zuker_traco();
    }

    if(kind == 2)
   {
      printf("pluto\n");
      zuker_pluto();
   }

    if(kind == 1 || CHECK_VALID)
    {    tmp_V = V;
        V = V1;
      tmp_W = W;
        W = W1;

        zuker_seq();
        V = tmp_V;
        W = tmp_W;

    if(CHECK_VALID && kind > 1)
     for(i=0; i<N; i++)
      for(j=0; j<N; j++)
       if((V[i][j] != V1[i][j]) || (W[i][j] != W1[i][j])){
          printf("error! %i %i\n", V[i][j], V1[i][j] );
         exit(0);
      }

    }


    double stop = omp_get_wtime();
    printf("%.2f seconds\n",stop - start);






return 0;
}
