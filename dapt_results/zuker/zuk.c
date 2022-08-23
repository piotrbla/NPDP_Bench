
#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
//define floord(n,d) floor(((double)(n))/((double)(d)))
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


void zuker_dapt0()
{
#pragma scop
	for (int i = N-1; i >= 0; i--){
	 for (int j = i+1; j < N; j++) {
	   for (int k = i+1; k < j; k++){
	   for(int m=k+1; m <j; m++){
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


  // print_matrix_org(W, N);
  //print_matrix_org(V, N);

}




int main(int argc, char *argv[]){

    int i,j,k;
//

    int num_proc=1;
    //num_proc = atoi(argv[1]);
    //num_proc = atoi("4");
    //if(argc > 2)
    //    N = atoi(argv[2]);
    N = 600; //atoi ("600");
    DIM = N+2;



    zuker_dapt0();



return 0;
}
