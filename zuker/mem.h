int **mem()
{
int i;
int **S;
S = (int **) malloc(DIM * sizeof(int*));

for (i=0; i<DIM; i++)
    S[i] = (int*)malloc(DIM * sizeof(int));

return S;
}
