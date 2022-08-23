#define INT_MAX 2147483647
#define INT_MIN -2147483648
#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))
#define CHECK_VALID 0


int **F;
int **F1;
int **H;
int **H1;
int **tmp_H;
int **m1;
int **m2;

int *W;
unsigned char *a;
unsigned char *b;


int N = 100, DIM = 102;

//Similarity score of the elements that constituted the two sequences
int s(unsigned char x, unsigned char z){
   return (x == z) ? 1 : -1;
}


#include "mem.h"

void nw_seq()
{
int i,j,k;


#pragma scop
for (i=1; i <=N; i++)
    for (j=1; j <=N; j++){
    // Block S
        m1[i][j] = INT_MIN;
        for (k=1; k <=i; k++)
            m1[i][j] = MAX(m1[i][j] ,F[i-k][j] - W[k]);
        m2[i][j] = INT_MIN;
        for (k=1; k <=j; k++)
            m2[i][j] = MAX(m2[i][j] ,F[i][j-k] - W[k]);
        F[i][j] = MAX(0, MAX( F[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
    }
#pragma endscop


}


int main(int argc, char *argv[]){

    int i,j,k;


    int num_proc=1;

    int kind=1;

    N = 100;
    DIM = 2*N+2;


    kind = 1;

    // H is the scoring matrix
    F = mem();
    F1 = mem();
    m1 = mem();
    m2 = mem();

    //W = (int*)malloc(DIM * sizeof(int));
    //a = (unsigned char *)malloc(DIM * sizeof(unsigned char ));
    //b = (unsigned char *)malloc(DIM * sizeof(unsigned char ));



    for(i=0; i<=N; i++){
        F[i][0] = 0;
        F[0][i] = 0;
        F1[i][0] = 0;
        F1[0][i] = 0;

    }


    // W is the gap alignment
    W[0] = 2;
    for(i=1; i<=N; i++)
        W[i] = i*W[0];

    rand_seq(a, N);
    rand_seq(b, N);



    nw_seq();




return 0;
}
