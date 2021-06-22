void sw_traco(){

printf("- traco [16x16x16 - \n\n");

int c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15;

#pragma omp parallel for
for( c1 = 0; c1 <= (N - 1)/16; c1 += 1)
  for( c3 = 0; c3 <= (N - 1) / 16; c3 += 1)
    for( c5 = 16 * c1 + 1; c5 <= min(N, 16 * c1 + 16); c5 += 1)
      for( c7 = 16 * c3 + 1; c7 <= min(N, 16 * c3 + 16); c7 += 1){
          m1[c5][c7] = INT_MIN;
          m2[c5][c7] = INT_MIN;
      }



  for( c1 = 0; c1 <= floord(N - 1, 8); c1 += 1)
  #pragma omp parallel for schedule(dynamic, 1) shared(c1) private(c3,c4,c7,c9,c10,c5,c11)
  for( c3 = max(0, c1 - (N + 15) / 16 + 1); c3 <= min(c1, (N - 1) / 16); c3 += 1)
    for( c4 = 0; c4 <= 2; c4 += 1) {
      if (c4 == 2) {
        for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
          for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1) {
            for( c10 = max(0, 16 * c1 - 16 * c3 - c7 + 2); c10 <= min(1, -16 * c3 + c9 - 1); c10 += 1) {
              if (c10 == 1) {
                for( c11 = 1; c11 <= c9; c11 += 1)
                  m2[c7][c9] = MAX(m2[c7][c9] + H[c7][c9-c11], W[c11]);
              } else
                for( c11 = 1; c11 <= c7; c11 += 1)
                  m1[c7][c9] = MAX(m1[c7][c9] + H[c7-c11][c9], W[c11]);
            }
            H[c7][c9] = MAX(0, MAX(H[c7-1][c9-1] + s(a[c7],b[c9]), MAX(m1[c7][c9], m2[c7][c9])));
          }
      } else if (c4 == 1) {
        for( c5 = 0; c5 <= c3; c5 += 1)
          for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
            for( c11 = 16 * c5 + 1; c11 <= min(16 * c3 + 1, 16 * c5 + 16); c11 += 1)
              m2[c7][(16*c3+1)] = MAX(m2[c7][(16*c3+1)], H[c7][(16*c3+1)-c11] + W[c11]);
      } else
        for( c5 = 0; c5 <= c1 - c3; c5 += 1)
          for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
            for( c11 = 16 * c5 + 1; c11 <= min(16 * c1 - 16 * c3 + 1, 16 * c5 + 16); c11 += 1)
              m1[(16*c1-16*c3+1)][c9] = MAX(m1[(16*c1-16*c3+1)][c9], H[(16*c1-16*c3+1)-c11][c9] + W[c11]);
    }



/*
for (i=1; i <=N; i++)
    for (j=1; j <=N; j++){
    // Block S
        m1[i][j] = INT_MIN;
        for (k=1; k <=i; k++)
            m1[i][j] = MAX(m1[i][j] ,H[i-k][j] + W[k]);
        m2[i][j] = INT_MIN;
        for (k=1; k <=j; k++)
            m2[i][j] = MAX(m2[i][j] ,H[i][j-k] + W[k]);
        H[i][j] = MAX(0, MAX( H[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
    }
*/


}
