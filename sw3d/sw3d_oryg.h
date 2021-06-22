// We consider three sequences of the same length N.
void sw_seq()
{
int i,j,k,l;

printf("- oryginal code - \n\n");

#pragma scop
for (i=1; i <=N; i++)
    for (j=1; j <=N; j++){
        for (l=1; l<=N; l++){
          // Block S
            m1[i][j][l] = INT_MIN;
            for (k=1; k <=i; k++)
                m1[i][j][l] = MAX(m1[i][j][l] ,H[i-k][j][l] - 2*W[k]);
            m2[i][j][l] = INT_MIN;
            for (k=1; k <=j; k++)
                m2[i][j][l] = MAX(m2[i][j][l], H[i][j-k][l] - 2*W[k]);
            m3[i][j][l] = INT_MIN;
            for (k=1; k <=l; k++)
                m3[i][j][l] = MAX(m3[i][j][l], H[i][j][l-k] - 2*W[k]);
            m4[i][j][l] = INT_MIN;
            for (k=1; k <=min(i,j); k++)
                m4[i][j][l] = MAX(m4[i][j][l], H[i-k][j-k][l] - W[k] + s(a[i], b[j]));
            m5[i][j][l] = INT_MIN;
            for (k=1; k <=min(j,l); k++)
                m5[i][j][l] = MAX(m5[i][j][l], H[i][j-k][l-k] - W[k] + s(b[j], c[l]));
            m6[i][j][l] = INT_MIN;
            for (k=1; k <=min(i,l); k++)
                m6[i][j][l] = MAX(m6[i][j][l], H[i-k][j][l-k] - W[k] + s(a[i], c[l]));
            H[i][j][l] = MAX(0, MAX( H[i-1][j-1][l-1] + s(a[i], b[j]) + s(a[i], c[l]) + s(b[j], c[l]), MAX(m1[i][j][l], MAX(m2[i][j][l], MAX(m3[i][j][l], MAX(m4[i][j][l], MAX(m5[i][j][l], m6[i][j][l])))))));


        }
    }
#pragma endscop

}

