void sw3d_pluto(){
printf("pluto\n");
  int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if (N >= 1) {
  for (t2=0;t2<=floord(N,8);t2++) {
    lbp=max(0,ceild(16*t2-N,16));
    ubp=min(floord(N,16),t2);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)
    for (t4=lbp;t4<=ubp;t4++) {
      for (t6=0;t6<=floord(N,16);t6++) {
        for (t7=max(1,16*t2-16*t4);t7<=min(N,16*t2-16*t4+15);t7++) {
          for (t9=max(1,16*t4);t9<=min(N,16*t4+15);t9++) {
            for (t11=max(1,16*t6);t11<=min(N,16*t6+15);t11++) {
              m1[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=t7;t13++) {
                m1[t7][t9][t11] = MAX(m1[t7][t9][t11] ,H[t7-t13][t9][t11] - 2*W[t13]);;
              }
              m2[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=t9;t13++) {
                m2[t7][t9][t11] = MAX(m2[t7][t9][t11], H[t7][t9-t13][t11] - 2*W[t13]);;
              }
              m3[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=t11;t13++) {
                m3[t7][t9][t11] = MAX(m3[t7][t9][t11], H[t7][t9][t11-t13] - 2*W[t13]);;
              }
              m4[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=min(t7,t9);t13++) {
                m4[t7][t9][t11] = MAX(m4[t7][t9][t11], H[t7-t13][t9-t13][t11] - W[t13] + s(a[t7], b[t9]));;
              }
              m5[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=min(t11,t9);t13++) {
                m5[t7][t9][t11] = MAX(m5[t7][t9][t11], H[t7][t9-t13][t11-t13] - W[t13] + s(b[t9], c[t11]));;
              }
              m6[t7][t9][t11] = INT_MIN;;
              for (t13=1;t13<=min(t11,t7);t13++) {
                m6[t7][t9][t11] = MAX(m6[t7][t9][t11], H[t7-t13][t9][t11-t13] - W[t13] + s(a[t7], c[t11]));;
              }
              H[t7][t9][t11] = MAX(0, MAX( H[t7-1][t9-1][t11-1] + s(a[t7], b[t9]) + s(a[t7], c[t11]) + s(b[t9], c[t11]), MAX(m1[t7][t9][t11], MAX(m2[t7][t9][t11], MAX(m3[t7][t9][t11], MAX(m4[t7][t9][t11], MAX(m5[t7][t9][t11], m6[t7][t9][t11])))))));;
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */


}