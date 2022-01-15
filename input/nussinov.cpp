//nussinov

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <cstring>
#include <string>

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))



long double **S;
char *RNA;
int N;

void oryg(){

int i,j,k;



    for (i = N-1; i >= 0; i--) {
     if(i % 100==0)
       printf("%i \n", i);
     for (j = i+1; j < N; j++) {
      for (k = 0; k < j-i; k++) {
        S[i][j] = max(S[i][k+i] + S[k+i+1][j], S[i][j]);
      }
      for (k = 0; k < 1; k++) {
       S[i][j] = max(S[i][j], S[i+1][j-1]  + can_pair(RNA, i, j));

     }
    }
   }


}


int main(int argc, char *argv[]){

    int i,j,k;
    char *filename, *method;
    int num_proc=-1;

    method = "oryg";

    if(argc < 4 )
    {
        printf("./nuss [method: oryg,tstile,tilecorr,pluto] [number of threads] [N] \n");
        return -1;
     }
    
    method = argv[1];
    num_proc = atoi(argv[2]);
    N = atoi(argv[3]);

    omp_set_num_threads(num_proc);  // else default max number
    
    S = mem();
    RNA = new char[N+5];

    printf("\nmethod %s\n", method);
    printf("N %i\n", N);

    double start = omp_get_wtime();



       oryg();





    double stop = omp_get_wtime();
//    saveTable(method, num_proc, filename);
    //saveTable();



    printf("Time: %.2f\n", stop - start);

   // printf("Traceback:\n");
   //  char *wout = new char[256];
  //   strcpy(wout, filename);

   //  FILE *plik = fopen(strcat(wout, ".traceback.txt") ,"w");
   //  traceback(0, N-1, plik);
   //  fclose(plik);

    return 0;
}
