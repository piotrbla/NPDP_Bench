#include<stdio.h>
#include<stdlib.h>

#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

int _01_nCatlan_org(int n) {
    int *dp = (int *)calloc(n+1, sizeof(int));
    dp[0] = dp[1] = 1;
    for(int i=2;i<=n;i++){
        for(int j=0;j<i;j++){
            dp[i] += dp[j]*dp[i-j-1];
        }
    }
    return dp[n-1];
}

int _02_nCatlan_mod(int n) {
    int *dp = (int *)calloc(n+1, sizeof(int));
    dp[0] = dp[1] = 1;
for (int c0 = 1; c0 < n; c0 += 1)
  for (int c1 = c0 + 1; c1 <= min(n, 2 * c0 + 1); c1 += 1) {
    if (2 * c0 >= c1)
      dp[c1] += dp[-c0 + c1 - 1]*dp[c1-(-c0 + c1 - 1)-1];
    dp[c1] += dp[c0]*dp[c1-c0-1];
  }
  return dp[n-1];
}

int _03_nCatlan_dapt(int n) {
    int *dp = (int *)calloc(n+1, sizeof(int));
    dp[0] = dp[1] = 1;
for (int w0 = 0; w0 <= floord(5 * n - 3, 48); w0 += 1) {
  #pragma omp parallel for
  for (int h0 = max(w0 - (n + 24) / 24 + 1, (3 * w0 + 3) / 7); h0 <= min((n - 1) / 16, w0 - (2 * w0 + 2) / 5); h0 += 1) {
    for (int i0 = max(max(1, 16 * h0), 12 * w0 - 12 * h0); i0 <= min(min(n - 1, 16 * h0 + 15), 24 * w0 - 24 * h0 + 22); i0 += 1) {
      for (int i1 = max(24 * w0 - 24 * h0, i0 + 1); i1 <= min(min(n, 24 * w0 - 24 * h0 + 23), 2 * i0 + 1); i1 += 1) {
        {
          if (2 * i0 >= i1) {
            dp[i1] += (dp[-i0 + i1 - 1] * dp[i0]);
          }
          dp[i1] += (dp[i0] * dp[-i0 + i1 - 1]);
        }
      }
    }
  }
}
  return dp[n-1];
}

int main() {
  int n;
  for (n = 1; n <= 20; n++)
  {
    printf("n=%d\n", n);
    printf("01:%d\n", _01_nCatlan_org(n));
    printf("02:%d\n", _02_nCatlan_mod(n));
    printf("03:%d\n", _03_nCatlan_dapt(n));
  }
  return 0;
}