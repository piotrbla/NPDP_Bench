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

int **c;
int **ck;
int **w;

int N;
int DIM;

#include "mem.h"

int main(int argc, char *argv[]){



    int num_proc=1;
    int i,j,k,ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    srand(time(NULL));


    if(argc > 1)
        num_proc = atoi(argv[1]);

    int kind=1;

    N = 8;
    DIM = 12;
    if(argc > 2)
        N = atoi(argv[2]);
    DIM = N+10;
    int n = N;

    if(argc > 3)
        kind = atoi(argv[3]);


    omp_set_num_threads(num_proc);

    c = mem();
    ck = mem();
    w = mem();

    double start = omp_get_wtime();
    if(kind==1){
        #pragma scop
        for(i=n-1; i>=1; i--)
           for(j=i+1; j<=n; j++)
               for(k=i+1; k<j; k++)
                  c[i][j] = MIN(c[i][j], w[i][j]+c[i][k]+c[k][j]);
        #pragma endscop
    }

    if(kind==2)
    {
        int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        int lb, ub, lbp, ubp, lb2, ub2;
        register int lbv, ubv;
        /* Start of CLooG code */
        if (n >= 3) {
        for (t2=-1;t2<=floord(n-16,16);t2++) {
            lbp=t2+1;
            ubp=min(floord(n,16),floord(16*t2+n+13,16));
            #pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10,t4) shared(t2)
            for (t4=lbp;t4<=ubp;t4++) {
            for (t5=max(max(-n+2,16*t2-16*t4),-16*t4-13);t5<=16*t2-16*t4+15;t5++) {
                for (t7=max(16*t4,-t5+2);t7<=min(n,16*t4+15);t7++) {
                for (t9=-t5+1;t9<=t7-1;t9++) {
                    c[-t5][t7] = MIN(c[-t5][t7], w[-t5][t7]+c[-t5][t9]+c[t9][t7]);;
                }
                }
            }
            }
        }
        }
        /* End of CLooG code */
    }


    if(kind == 3)  //traco
    {
        /*for( c1 = 1; c1 < n + floord(n - 2, 128) - 1; c1 += 1)
        #pragma omp parallel for shared(c1) private(c3,c5,c9,c11) schedule(dynamic, 1)
        for( c3 = max(0, -n + c1 + 2); c3 <= c1 / 129; c3 += 1)
            for( c5 = 0; c5 <= min(8 * c3 + 7, (c1 - c3 - 1) / 16); c5 += 1)
            for( c9 = max(n - c1 + 129 * c3, n - c1 + c3 + 16 * c5 + 1); c9 <= min(n, n - c1 + 129 * c3 + 127); c9 += 1)
                for( c11 = n - c1 + c3 + 16 * c5; c11 <= min(n - c1 + c3 + 16 * c5 + 15, c9 - 1); c11 += 1)
                c[(n-c1+c3-1)][c9] = MIN(c[(n-c1+c3-1)][c9], w[(n-c1+c3-1)][c9]+c[(n-c1+c3-1)][c11]+c[c11][c9]);
*/

/*
{
  for( c1 = 0; c1 <= floord(n - 3, 8); c1 += 1)
    #pragma omp parallel for shared(c1) private(c3,c5,c9,c11) schedule(dynamic, 1)
    for( c3 = max(0, c1 - (n + 14) / 16 + 1); c3 <= c1 / 2; c3 += 1)
      for( c5 = 0; c5 <= c3; c5 += 1)
        for( c7 = max(-n + 16 * c1 - 16 * c3 + 1, -n + 15 * c1 - 15 * c3 + c5 + 2); c7 <= min(-1, -n + 16 * c1 - 16 * c3 + 16); c7 += 1) {
          if (n + 16 * c3 + c7 >= 16 * c1 + 2) {
            if (c3 >= c5 + 1)
              for( c11 = 16 * c5 - c7 + 1; c11 <= 16 * c5 - c7 + 16; c11 += 1)
                c[(-c7)][(16*c3-c7+1)] = MIN(c[(-c7)][(16*c3-c7+1)], w[(-c7)][(16*c3-c7+1)]+c[(-c7)][c11]+c[c11][(16*c3-c7+1)]);
            for( c9 = 16 * c3 - c7 + 2; c9 <= min(min(n, 16 * c3 - c7 + 16), 16 * c5 - c7 + 18); c9 += 1) {
              if (c5 == c3 && n + 16 * c3 + c7 >= 16 * c1 + 3 && 16 * c1 + c9 == n + 32 * c3) {
                for( c11 = -c7 + 1; c11 <= 16 * c3 - c7; c11 += 1)
                  c[(-c7)][(n-16*c1+32*c3)] = MIN(c[(-c7)][(n-16*c1+32*c3)], w[(-c7)][(n-16*c1+32*c3)]+c[(-c7)][c11]+c[c11][(n-16*c1+32*c3)]);
              } else if (c7 + c9 >= 15 * c3 + c5 + 3 && n + 32 * c3 >= 16 * c1 + c9 + 1) {
                for( c11 = -c7 + 1; c11 <= min(16 * c3 - c7, -15 * c3 + 15 * c5 + c9 - 3); c11 += 1)
                  c[(-c7)][c9] = MIN(c[(-c7)][c9], w[(-c7)][c9]+c[(-c7)][c11]+c[c11][c9]);
              } else if (16 * c1 + c9 >= n + 16 * c3 + 16 * c5 + 1 && 16 * c1 + c9 >= n + 32 * c3) {
                for( c11 = -c7 + 1; c11 <= 16 * c5 - c7; c11 += 1)
                  c[(-c7)][c9] = MIN(c[(-c7)][c9], w[(-c7)][c9]+c[(-c7)][c11]+c[c11][c9]);
              }
              for( c11 = 16 * c3 - c7 + 1; c11 < c9 - 2; c11 += 1)
                c[(-c7)][c9] = MIN(c[(-c7)][c9], w[(-c7)][c9]+c[(-c7)][c11]+c[c11][c9]);
              for( c11 = max(16 * c5 - c7 + 1, -15 * c3 + 15 * c5 + c9 - 2); c11 <= min(16 * c5 - c7 + 16, c9 - 1); c11 += 1)
                c[(-c7)][c9] = MIN(c[(-c7)][c9], w[(-c7)][c9]+c[(-c7)][c11]+c[c11][c9]);
            }
          } else {
            for( c9 = max(n - 16 * c1 + 32 * c3, n - 16 * c1 + 16 * c3 + 16 * c5 + 1); c9 <= min(n, n - 16 * c1 + 32 * c3 + 15); c9 += 1)
              for( c11 = n - 16 * c1 + 16 * c3 + 16 * c5; c11 <= min(n - 16 * c1 + 16 * c3 + 16 * c5 + 15, c9 - 1); c11 += 1)
                c[(n-16*c1+16*c3-1)][c9] = MIN(c[(n-16*c1+16*c3-1)][c9], w[(n-16*c1+16*c3-1)][c9]+c[(n-16*c1+16*c3-1)][c11]+c[c11][c9]);
          }
        }
  if ((n - 2) % 16 == 0)
    for( c5 = 0; c5 < (n - 2) / 16; c5 += 1)
      for( c11 = 16 * c5 + 2; c11 <= 16 * c5 + 17; c11 += 1)
        c[1][n] = MIN(c[1][n], w[1][n]+c[1][c11]+c[c11][n]);
}
*/


{
  for( c1 = 0; c1 < n + floord(-3 * n - 3, 8); c1 += 1)
    #pragma omp parallel for
    for( c3 = max(0, c1 - (n + 6) / 8 + 1); c3 <= min(n / 2 - 1, c1 - (c1 + 6) / 5 + 1); c3 += 1)
      for( c5 = 0; c5 <= c3 / 128; c5 += 1)
        for( c7 = max(max(-n + 2 * c3 + 1, -n + 8 * c1 - 8 * c3 + 1), -n + c3 + 128 * c5 + 2); c7 <= min(-1, -n + 8 * c1 - 8 * c3 + 8); c7 += 1) {
          if (n + 8 * c3 + c7 >= 8 * c1 + 2) {
            for( c11 = 256 * c5 - c7 + 1; c11 <= min(2 * c3 - c7, 256 * c5 - c7 + 256); c11 += 1)
              c[(-c7)][(2*c3-c7+1)] = MIN(c[(-c7)][(2*c3-c7+1)], w[(-c7)][(2*c3-c7+1)]+c[(-c7)][c11]+c[c11][(2*c3-c7+1)]);
            if (128 * c5 + 128 >= c3 && n + c7 >= 2 * c3 + 2) {
              if (c3 >= 128 * c5 + 1)
                for( c11 = -c7 + 1; c11 <= 256 * c5 - c7; c11 += 1)
                  c[(-c7)][(2*c3-c7+2)] = MIN(c[(-c7)][(2*c3-c7+2)], w[(-c7)][(2*c3-c7+2)]+c[(-c7)][c11]+c[c11][(2*c3-c7+2)]);
              for( c11 = 256 * c5 - c7 + 1; c11 <= min(2 * c3 - c7 + 1, 256 * c5 - c7 + 256); c11 += 1)
                c[(-c7)][(2*c3-c7+2)] = MIN(c[(-c7)][(2*c3-c7+2)], w[(-c7)][(2*c3-c7+2)]+c[(-c7)][c11]+c[c11][(2*c3-c7+2)]);
            }
          } else {
            for( c9 = max(n - 8 * c1 + 10 * c3, n - 8 * c1 + 8 * c3 + 256 * c5 + 1); c9 <= min(n, n - 8 * c1 + 10 * c3 + 1); c9 += 1)
              for( c11 = n - 8 * c1 + 8 * c3 + 256 * c5; c11 <= min(n - 8 * c1 + 8 * c3 + 256 * c5 + 255, c9 - 1); c11 += 1)
                c[(n-8*c1+8*c3-1)][c9] = MIN(c[(n-8*c1+8*c3-1)][c9], w[(n-8*c1+8*c3-1)][c9]+c[(n-8*c1+8*c3-1)][c11]+c[c11][c9]);
          }
        }
  if ((n - 2) % 8 == 0)
    for( c5 = 0; c5 <= floord(n - 10, 256); c5 += 1)
      for( c11 = 256 * c5 + 2; c11 <= min(n - 1, 256 * c5 + 257); c11 += 1)
        c[1][n] = MIN(c[1][n], w[1][n]+c[1][c11]+c[c11][n]);
}


    }

    if(kind == 4)
    {
        for( c1 = 1; c1 < n + floord(n - 2, 128) - 1; c1 += 1)
            #pragma omp parallel for shared(c1) private(c3,c5,c9)
            for( c3 = max(0, -n + c1 + 2); c3 <= c1 / 129; c3 += 1){
                int lb1 = n - c1 + 129 * c3; 
                int id1 = n-c1+c3-1;      
                for( c5 = 0; c5 <= min(8 * c3 + 7, (c1 - c3 - 1) / 16); c5 += 1){
                int lb2 = n - c1 + c3 + 16 * c5;
                for( c9 = max(lb1, lb2 + 1); c9 <= min(n, lb1 + 127); c9 += 1)
                    for( c11 = lb2; c11 <= min(lb2 + 15, c9 - 1); c11 += 1)
                    c[(id1)][c9] = MIN(c[(id1)][c9], w[(id1)][c9]+c[(id1)][c11]+c[c11][c9]);
                }
            }
    }

    double stop = omp_get_wtime();
    printf("%.4f\n",stop - start);

     /* for(i=0; i<DIM; i++)
    for(j=0; j<DIM; j++)
     if(c[i][j] != c1[i][j]){

       printf("Bad! %.2f %.2f %i %i\n", Q[i][j], Q1[i][j],i,j);
       exit(1);
    }
*/

    return 0;
}
