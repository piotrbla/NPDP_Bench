

// We consider two sequences of the same length N.

void nw_seq()
{
int i,j,k;

printf("- oryginal code - \n\n");

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
