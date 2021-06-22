// See the Cormen book for details of the following algorithm
#include<stdio.h>
#include<limits.h>
#include <math.h>
#include <omp.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))


int N = 1500, DIM = 1502;
#include "mem.h"

#define pluto 3
#define traco 2
#define tstile 4


int **s;


int minsq(int a,int b,int i,int j,int k)
{
    if(a < b){
        s[i][j] = k;
        return a;
    }
    else
        return b;
}


// Matrix Ai has dimension p[i-1] x p[i] for i = 1..n
int MatrixChainOrder(int p[], int n, int kind)
{

    /* For simplicity of the program, one extra row and one
       extra column are allocated in m[][].  0th row and 0th
       column of m[][] are not used */
    int** m = mem();
    s = mem();

    int i, j, k, L, q;

    /* m[i,j] = Minimum number of scalar multiplications needed
       to compute the matrix A[i]A[i+1]...A[j] = A[i..j] where
       dimension of A[i] is p[i-1] x p[i] */

    double start = omp_get_wtime();
    // cost is zero when multiplying one matrix.
    for (i=1; i<n; i++)
        m[i][i] = 0;

    // L is chain length.
    if(kind==-1)
    for (L=2; L<n; L++)
    {
        for (i=1; i<n-L+1; i++)
        {
            j = i+L-1;
            m[i][j] = INT_MAX;
            for (k=i; k<=j-1; k++)
            {
                // q = cost/scalar multiplications
                q = m[i][k] + m[k+1][j] + p[i-1]*p[k]*p[j];
                if (q < m[i][j]){
                    m[i][j] = q;
                    s[i][j] = k;
                }
            }
        }
    }

    if(kind==1)
    for (L=2; L<n; L++)
    {
        for (i=1; i<n-L+1; i++)
        {
            m[i][i+L-1] = INT_MAX;
            for (k=i; k<=i+L-1-1; k++)
            {
                m[i][i+L-1] = minsq(m[i][k] + m[k+1][i+L-1] + p[i-1]*p[k]*p[i+L-1], m[i][i+L-1], i, i+L-1, k);
            }
        }
    }

    if(kind==pluto)
    {

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if (n >= 3) {
  lbp=0;
  ubp=floord(n-1,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=0;t3<=min(floord(n-2,16),floord(-16*t2+n,16));t3++) {
      for (t4=max(1,16*t3);t4<=min(min(n-2,-16*t2+n),16*t3+15);t4++) {
        lbv=max(2,16*t2);
        ubv=min(16*t2+15,-t4+n);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          m[t4][t4+t5-1] = INT_MAX;;
        }
      }
    }
  }
  for (t1=0;t1<=floord(n-1,8);t1++) {
    lbp=ceild(t1,2);
    ubp=min(floord(n,16),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t2) shared(t1)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(2,16*t1-16*t2);t3<=min(min(n-1,16*t2+14),16*t1-16*t2+15);t3++) {
        for (t4=max(16*t2,t3+1);t4<=min(n,16*t2+15);t4++) {
          for (t5=-t3+t4;t5<=t4-2;t5++) {
            m[(-t3+t4)][(-t3+t4)+t3-1] = minsq(m[(-t3+t4)][t5] + m[t5+1][(-t3+t4)+t3-1] + p[(-t3+t4)-1]*p[t5]*p[(-t3+t4)+t3-1], m[(-t3+t4)][(-t3+t4)+t3-1], (-t3+t4),(-t3+t4)+t3-1,t5);;
          }
        }
      }
    }
  }
}
}
/* End of CLooG code */


if(kind==traco)
{
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;

  lbp=0;
  ubp=floord(n-1,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=0;t3<=min(floord(n-2,16),floord(-16*t2+n,16));t3++) {
      for (t4=max(1,16*t3);t4<=min(min(n-2,-16*t2+n),16*t3+15);t4++) {
        lbv=max(2,16*t2);
        ubv=min(16*t2+15,-t4+n);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          m[t4][t4+t5-1] = INT_MAX;;
        }
      }
    }
  }

int c1,c3,c5,c9,c7,c11;


for( c1 = 0; c1 < n - 2; c1 += 1)
  #pragma omp parallel for shared(c1) private(c3,c5,c9,c11)
  for( c3 = 0; c3 <= (n - c1 - 3) / 128; c3 += 1)
    for( c5 = 0; c5 <= c1 / 16; c5 += 1)
      for( c9 = 128 * c3 + 1; c9 <= min(n - c1 - 2, 128 * c3 + 128); c9 += 1)
        for( c11 = 16 * c5 + c9; c11 <= min(c1 + c9, 16 * c5 + c9 + 15); c11 += 1)
          m[c9][c9+(c1+2)-1] = minsq(m[c9][c11] + m[c11+1][c9+(c1+2)-1] + p[c9-1]*p[c11]*p[c9+(c1+2)-1], m[c9][c9+(c1+2)-1], c9,c9+(c1+2)-1,c11);



if(1==0)
for( c1 = 0; c1 <= floord(n - 3, 16); c1 += 1)
  #pragma omp parallel for
  for( c3 = 0; c3 <= -c1 + (n - 3) / 16; c3 += 1)
    for( c5 = 0; c5 <= c1; c5 += 1) {
      if (c1 >= c5 + 1) {
        for( c9 = 16 * c1 + 16 * c3 + 3; c9 <= min(n, 16 * c1 + 16 * c3 + 18); c9 += 1)
          for( c11 = -16 * c1 + 16 * c5 + c9 - 2; c11 <= -16 * c1 + 16 * c5 + c9 + 13; c11 += 1)
            m[(-16*c1+c9-2)][(-16*c1+c9-2)+(16*c1+2)-1] = minsq(m[(-16*c1+c9-2)][c11] + m[c11+1][(-16*c1+c9-2)+(16*c1+2)-1] + p[(-16*c1+c9-2)-1]*p[c11]*p[(-16*c1+c9-2)+(16*c1+2)-1], m[(-16*c1+c9-2)][(-16*c1+c9-2)+(16*c1+2)-1], (-16*c1+c9-2),(-16*c1+c9-2)+(16*c1+2)-1,c11);
      } else
        for( c7 = 16 * c1 + 2; c7 <= min(16 * c1 + 17, n - c3 - 1); c7 += 1)
          for( c9 = max(16 * c1 + 16 * c3 + 3, c7 + 1); c9 <= min(n, 16 * c1 + 16 * c3 + 18); c9 += 1) {
            if (c9 >= 16 * c3 + c7 + 1) {
              if (c7 >= 16 * c1 + 3)
                for( c11 = -c7 + c9; c11 < 16 * c1 - c7 + c9; c11 += 1)
                  m[(-c7+c9)][(-c7+c9)+c7-1] = minsq(m[(-c7+c9)][c11] + m[c11+1][(-c7+c9)+c7-1] + p[(-c7+c9)-1]*p[c11]*p[(-c7+c9)+c7-1], m[(-c7+c9)][(-c7+c9)+c7-1], (-c7+c9),(-c7+c9)+c7-1,c11);
              for( c11 = 16 * c1 - c7 + c9; c11 < c9 - 1; c11 += 1)
                m[(-c7+c9)][(-c7+c9)+c7-1] = minsq(m[(-c7+c9)][c11] + m[c11+1][(-c7+c9)+c7-1] + p[(-c7+c9)-1]*p[c11]*p[(-c7+c9)+c7-1], m[(-c7+c9)][(-c7+c9)+c7-1], (-c7+c9),(-c7+c9)+c7-1,c11);
            } else
              for( c11 = -c7 + c9; c11 < c9 - 1; c11 += 1)
                m[(-c7+c9)][(-c7+c9)+c7-1] = minsq(m[(-c7+c9)][c11] + m[c11+1][(-c7+c9)+c7-1] + p[(-c7+c9)-1]*p[c11]*p[(-c7+c9)+c7-1], m[(-c7+c9)][(-c7+c9)+c7-1], (-c7+c9),(-c7+c9)+c7-1,c11);
          }
    }

}

if(kind == tstile){

  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
 int c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c0;

  lbp=0;
  ubp=floord(n-1,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
  for (t2=lbp;t2<=ubp;t2++) {
    for (t3=0;t3<=min(floord(n-2,16),floord(-16*t2+n,16));t3++) {
      for (t4=max(1,16*t3);t4<=min(min(n-2,-16*t2+n),16*t3+15);t4++) {
        lbv=max(2,16*t2);
        ubv=min(16*t2+15,-t4+n);
#pragma ivdep
#pragma vector always
        for (t5=lbv;t5<=ubv;t5++) {
          m[t4][t4+t5-1] = INT_MAX;;
        }
      }
    }
  }

for( c0 = 0; c0 <= floord(n - 3, 16); c0 += 1)
  #pragma omp parallel for private(c3,c1,c4,c10,c6) shared(c0)
  for( c1 = 0; c1 <= c0; c1 += 1)
    for( c3 = 16 * c0 + 16 * c1 + 5; c3 <= min(min(2 * n - 16 * c0 + 16 * c1 - 1, n + 16 * c1 + 17), 16 * c0 + 16 * c1 + 50); c3 += 1)
      for( c4 = max(c0 - c1, -2 * c1 + (c3 - 2) / 16 - 2); c4 <= min((n - 2) / 16, -c1 + (16 * c0 + 16 * c1 + c3 + 12) / 32); c4 += 1)
        for( c6 = max(max(max(16 * c1 + 2, -n + c3), -8 * c4 + c3 / 2 - 7), -8 * c0 + 8 * c1 + (c3 + 1) / 2 - 8); c6 <= min(min(16 * c1 + 17, c3 - 16 * c4 - 2), -8 * c0 + 8 * c1 + (c3 + 1) / 2 - 1); c6 += 1)
          for( c10 = max(16 * c4, c3 - 2 * c6); c10 <= min(16 * c4 + 15, c3 - c6 - 2); c10 += 1)
            m[(c3-2*c6)][(c3-2*c6)+c6-1] = minsq(m[(c3-2*c6)][c10] + m[c10+1][(c3-2*c6)+c6-1] + p[(c3-2*c6)-1]*p[c10]*p[(c3-2*c6)+c6-1], m[(c3-2*c6)][(c3-2*c6)+c6-1], (c3-2*c6), (c3-2*c6)+c6-1, c10);



if(1==0)
{

for( c0 = 0; c0 <= floord(n - 3, 16); c0 += 1)
  #pragma omp parallel for
  for( c1 = 0; c1 <= c0; c1 += 1) {
    for( c3 = 16 * c1 + 2; c3 <= min(n - 16 * c0 + 16 * c1 - 1, 16 * c1 + 17); c3 += 1)
      for( c4 = c0 - c1; c4 <= min(c0 - c1 + 1, (n - c3) / 16); c4 += 1) {
        if (c0 == 0 && c1 == 0 && c4 == 0)
          for( c6 = 2; c6 < (c3 + 1) / 2; c6 += 1)
            for( c10 = c3 - 2 * c6; c10 < c3 - c6 - 1; c10 += 1)
              m[(c3-2*c6)][(c3-2*c6)+c6-1] = minsq(m[(c3-2*c6)][c10] + m[c10+1][(c3-2*c6)+c6-1] + p[(c3-2*c6)-1]*p[c10]*p[(c3-2*c6)+c6-1], m[(c3-2*c6)][(c3-2*c6)+c6-1], (c3-2*c6), (c3-2*c6)+c6, c10);
        for( c8 = max(16 * c0 - 16 * c1 + 1, 16 * c4); c8 <= min(min(16 * c0 - 16 * c1 + 16, n - c3), 16 * c4 + 15); c8 += 1)
          m[c8][c8+c3-1] = INT_MAX;
      }
    for( c3 = max(5, n); c3 <= min(17, 2 * n - 1); c3 += 1)
      for( c6 = max(2, -n + c3); c6 < (c3 + 1) / 2; c6 += 1)
        for( c10 = c3 - 2 * c6; c10 < c3 - c6 - 1; c10 += 1)
          m[(c3-2*c6)][(c3-2*c6)+c6-1] = minsq(m[(c3-2*c6)][c10] + m[c10+1][(c3-2*c6)+c6-1] + p[(c3-2*c6)-1]*p[c10]*p[(c3-2*c6)+c6-1], m[(c3-2*c6)][(c3-2*c6)+c6-1], (c3-2*c6), (c3-2*c6)+c6, c10);
    for( c3 = max(16 * c0 + 16 * c1 + 5, 16 * c1 + 18); c3 <= min(min(2 * n - 16 * c0 + 16 * c1 - 1, n + 16 * c1 + 17), 16 * c0 + 16 * c1 + 50); c3 += 1)
      for( c4 = max(c0 - c1, -2 * c1 + (c3 - 2) / 16 - 2); c4 <= min((n - 2) / 16, -c1 + (16 * c0 + 16 * c1 + c3 + 12) / 32); c4 += 1)
        for( c6 = max(max(max(16 * c1 + 2, -n + c3), -8 * c4 + c3 / 2 - 7), -8 * c0 + 8 * c1 + (c3 + 1) / 2 - 8); c6 <= min(min(16 * c1 + 17, c3 - 16 * c4 - 2), -8 * c0 + 8 * c1 + (c3 + 1) / 2 - 1); c6 += 1)
          for( c10 = max(16 * c4, c3 - 2 * c6); c10 <= min(16 * c4 + 15, c3 - c6 - 2); c10 += 1)
            m[(c3-2*c6)][(c3-2*c6)+c6-1] = minsq(m[(c3-2*c6)][c10] + m[c10+1][(c3-2*c6)+c6-1] + p[(c3-2*c6)-1]*p[c10]*p[(c3-2*c6)+c6-1], m[(c3-2*c6)][(c3-2*c6)+c6-1], (c3-2*c6), (c3-2*c6)+c6, c10);
  }



}




}




    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

    return m[1][n-1];
}

int main(int argc, char *argv[]){

    int num_proc=1, i;
    if(argc > 1)
        num_proc = atoi(argv[1]);

    omp_set_num_threads(num_proc);

    int kind=1;

    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+2;


    if(argc > 3)
        kind = atoi(argv[3]);



    int *p = (int *) malloc(DIM * sizeof(int));
    for(i=0; i<N; i++)
        p[i] = (i % 20) + 1;

    //int p[] = {1, 2, 3, 4};
    //N = sizeof(p)/sizeof(p[0]);


    printf("Minimum number of multiplications is %d \n",
                       MatrixChainOrder(p, N, kind));

    return 0;
}
