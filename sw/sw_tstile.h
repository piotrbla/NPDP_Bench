void sw_tstile(){

printf("- tstile [16] -\n\n");

int c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15;

#pragma omp parallel for
for( c1 = 0; c1 <= (N - 1)/16; c1 += 1)
  for( c3 = 0; c3 <= (N - 1) / 16; c3 += 1)
    for( c5 = 16 * c1 + 1; c5 <= min(N, 16 * c1 + 16); c5 += 1)
      for( c7 = 16 * c3 + 1; c7 <= min(N, 16 * c3 + 16); c7 += 1){
          m1[c5][c7] = INT_MIN;
          m2[c5][c7] = INT_MIN;
      }


if(1==1)
for( c0 = 0; c0 <= floord(N - 1, 8); c0 += 1)
  #pragma omp parallel for schedule(dynamic, 1) shared(c0) private(c1,c3,c4,c6,c10)
  for( c1 = max(0, c0 - (N + 15) / 16 + 1); c1 <= min(c0, (N - 1) / 16); c1 += 1)
    for( c3 = 16 * c0 + 2; c3 <= min(min(min(2 * N, 16 * c0 + 32), N + 16 * c1 + 16), N + 16 * c0 - 16 * c1 + 16); c3 += 1) {
      for( c4 = max(max(-c0 + c1 - 1, -((N + 14) / 16)), c1 - (c3 + 13) / 16); c4 < c0 - c1 - (c3 + 13) / 16; c4 += 1)
        for( c6 = max(max(16 * c1 + 1, -16 * c0 + 16 * c1 + c3 - 16), -N + c3); c6 <= min(min(16 * c1 + 16, -16 * c0 + 16 * c1 + c3 - 1), c3 + 16 * c4 + 14); c6 += 1)
          for( c10 = max(1, c3 + 16 * c4 - c6); c10 <= c3 + 16 * c4 - c6 + 15; c10 += 1)
            m2[c6][(c3-c6)] = MAX(m2[c6][(c3-c6)] ,H[c6][(c3-c6)-c10] + W[c10]);
      if (c0 >= 2 * c1 + 1 && c3 >= 16 * c0 + 19)
        for( c6 = max(-16 * c0 + 16 * c1 + c3 - 16, -N + c3); c6 <= 16 * c1 + 16; c6 += 1)
          for( c10 = max(1, -16 * c1 + c3 - c6 - 32); c10 < -16 * c1 + c3 - c6 - 16; c10 += 1)
            m2[c6][(c3-c6)] = MAX(m2[c6][(c3-c6)] ,H[c6][(c3-c6)-c10] + W[c10]);
      for( c4 = max(max(-c1 - 1, -((N + 14) / 16)), c0 - c1 - (c3 + 13) / 16); c4 <= 0; c4 += 1) {
        if (N + 16 * c1 + 1 >= c3 && 16 * c0 + 17 >= c3 && c1 + c4 == -1)
          for( c10 = max(1, -32 * c1 + c3 - 17); c10 < -32 * c1 + c3 - 1; c10 += 1)
            m2[(16*c1+1)][(-16*c1+c3-1)] = MAX(m2[(16*c1+1)][(-16*c1+c3-1)] ,H[(16*c1+1)][(-16*c1+c3-1)-c10] + W[c10]);
        for( c6 = max(max(max(16 * c1 + 1, -16 * c0 + 16 * c1 + c3 - 16), -N + c3), -16 * c4 - 14); c6 <= min(min(N, 16 * c1 + 16), -16 * c0 + 16 * c1 + c3 - 1); c6 += 1) {
          for( c10 = max(1, 16 * c4 + c6); c10 <= min(c6, 16 * c4 + c6 + 15); c10 += 1)
            m1[c6][(c3-c6)] = MAX(m1[c6][(c3-c6)] ,H[c6-c10][(c3-c6)] + W[c10]);
          for( c10 = max(1, c3 + 16 * c4 - c6); c10 <= min(c3 - c6, c3 + 16 * c4 - c6 + 15); c10 += 1)
            m2[c6][(c3-c6)] = MAX(m2[c6][(c3-c6)] ,H[c6][(c3-c6)-c10] + W[c10]);
          if (c0 == 0 && c1 == 0 && c3 <= 15 && c4 == 0)
            H[c6][(c3-c6)] = MAX(0, MAX( H[c6-1][(c3-c6)-1] + s(a[c6], b[c6]), MAX(m1[c6][(c3-c6)], m2[c6][(c3-c6)])));
        }
      }
      if (c3 >= 16)
        for( c6 = max(max(16 * c1 + 1, -16 * c0 + 16 * c1 + c3 - 16), -N + c3); c6 <= min(min(N, 16 * c1 + 16), -16 * c0 + 16 * c1 + c3 - 1); c6 += 1)
          H[c6][(c3-c6)] = MAX(0, MAX( H[c6-1][(c3-c6)-1] + s(a[c6], b[c6]), MAX(m1[c6][(c3-c6)], m2[c6][(c3-c6)])));
    }



/*
for( c0 = 0; c0 < N + floord(N - 1, 16); c0 += 1)
  #pragma omp parallel for shared(c0) private(c1,c3,c4,c10) schedule(dynamic,1)
  for( c1 = max(0, c0 - (N + 15) / 16 + 1); c1 <= min(N - 1, c0); c1 += 1)
    for( c3 = 16 * c0 - 15 * c1 + 2; c3 <= min(N + c1 + 1, 16 * c0 - 15 * c1 + 17); c3 += 1) {
      for( c4 = -((-c1 + c3 + 13) / 16); c4 < -((c1 + 15) / 16); c4 += 1)
        for( c10 = max(1, -c1 + c3 + 16 * c4 - 1); c10 <= -c1 + c3 + 16 * c4 + 14; c10 += 1)
          m2[(c1+1)][(-c1+c3-1)] = MAX(m2[(c1+1)][(-c1+c3-1)] ,H[(c1+1)][(-c1+c3-1)-c10] + W[c10]);
      for( c4 = -((c1 + 15) / 16); c4 <= 0; c4 += 1) {
        for( c10 = max(1, c1 + 16 * c4 + 1); c10 <= min(c1 + 1, c1 + 16 * c4 + 16); c10 += 1)
          m1[(c1+1)][(-c1+c3-1)] = MAX(m1[(c1+1)][(-c1+c3-1)] ,H[(c1+1)-c10][(-c1+c3-1)] + W[c10]);
        for( c10 = max(1, -c1 + c3 + 16 * c4 - 1); c10 <= min(-c1 + c3 - 1, -c1 + c3 + 16 * c4 + 14); c10 += 1)
          m2[(c1+1)][(-c1+c3-1)] = MAX(m2[(c1+1)][(-c1+c3-1)] ,H[(c1+1)][(-c1+c3-1)-c10] + W[c10]);
      }
      H[(c1+1)][(-c1+c3-1)] = MAX(0, MAX( H[(c1+1)-1][(-c1+c3-1)-1] + s(a[(c1+1)], b[(c1+1)]), MAX(m1[(c1+1)][(-c1+c3-1)], m2[(c1+1)][(-c1+c3-1)])));
    }
*/


}
