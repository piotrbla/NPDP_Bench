#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int ***mem()
{
int i,j;
int ***S;
S = (int ***) malloc(DIM * sizeof(int**));

for (i=0; i<DIM; i++)
    S[i] = (int**)malloc(DIM * sizeof(int*));

for (i=0; i<DIM; i++)
    for (j=0; j<DIM; j++)
      S[i][j] = (int*)malloc(DIM * sizeof(int));

return S;
}


float ****mem4()
{
int i,j,k;
float ****S;
S = (float ****) malloc(DIM * sizeof(float***));

for (i=0; i<DIM; i++)
    S[i] = (float***)malloc(DIM * sizeof(float**));

for (i=0; i<DIM; i++)
    for (j=0; j<DIM; j++)
      S[i][j] = (float**)malloc(DIM * sizeof(float*));

for (i=0; i<DIM; i++)
    for (j=0; j<DIM; j++)
      for (k=0; k<DIM; k++)
      S[i][j][k] = (float*)malloc(DIM * sizeof(float));

return S;
}



void rand_seq(unsigned char*a, int N){
  int i, tmp;
  srand(time(NULL));
  for(i=0; i<N; i++)
  {
      tmp = rand()%4;

      switch(tmp){
          case 0 : a[i] = 'A'; break;
          case 1 : a[i] = 'G'; break;
          case 2 : a[i] = 'C'; break;
          case 3 : a[i] = 'T'; break;
      }

  }

}
