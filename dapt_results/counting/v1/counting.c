#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int ** c;int ** ck; int zz = 2; int ** F; int N; char * RNA;  
#include "mem.h"
int paired(int i, int j) {return 0;}

int main(int argc, char *argv[]){
    int ll,p,q,l=0;
    int kind=1;
    N = 8;
    int check=1;

    if(kind==1 || check){
        #pragma scop
        for (int i = N-2;  i>=1; i--){
           for (int j=i+2; j<= N; j++){
             for (int k = i; k<=j-l; k++){
              ck[i][j] +=  ck[i][j-1] + paired(k,j) ?  ck[i][k-1] + ck[k+1][j-1] : 0;
           }
          }
        }
       #pragma endscop
    }
    return 0;
}

