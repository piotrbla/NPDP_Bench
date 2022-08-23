#define INT_MAX 2147483647
#define INT_MIN -2147483648
#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

#define sqrt(d) (double)(d)

int N = 1500, DIM = 1502;
#include "mem.h"

#define pluto 3
#define pluto2 6
#define traco 2
#define tstile 4

int **points;

// A utility function to find distance between two points in a plane
double dist(int * p1, int * p2)
{
    return sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) +
                (p1[1] - p2[1])*(p1[1] - p2[1]));
}

// A utility function to find cost of a triangle. The cost is considered
// as perimeter (sum of lengths of all edges) of the triangle
double cost(int i, int j, int k)
{
    int* p1 = points[i];
    int* p2 = points[j];
    int *p3 = points[k];
    return dist(p1, p2) + dist(p2, p3) + dist(p3, p1);
}




// Matrix Ai has dimension p[i-1] x p[i] for i = 1..n
double mcTDP(int kind)
{

    double** table = memd();


    int i, j, k, gap, q;

    /* m[i,j] = Minimum number of scalar multiplications needed
       to compute the matrix A[i]A[i+1]...A[j] = A[i..j] where
       dimension of A[i] is p[i-1] x p[i] */



   #pragma scop
   for (gap = 0; gap < N; gap++)
   {
      for (j = gap; j < N; j++)    // i = j - gap
      {
          if (gap < 2)
             table[j-gap][j] = 0.0;
          else
          {
              table[j-gap][j] = INT_MAX;
              for (k = j-gap+1; k < j; k++)
              {
                table[j-gap][j]  = MIN(table[j-gap][j], table[j-gap][k] + table[k][j] + cost(j-gap,j,k));

              }
          }
      }
   }
  #pragma endscop

     return  table[0][N-1];
}

int main(int argc, char *argv[]){

    int num_proc=1, i;
    int kind=1;

    DIM = N+2;
    //points = (int **) malloc(DIM * sizeof(int*));

    // for(i=0; i<DIM; i++)
    //    points[i] = (int *) malloc(2 * sizeof(int));



    mcTDP(kind);

    return 0;
}
