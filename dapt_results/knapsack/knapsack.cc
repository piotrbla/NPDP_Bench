// Unbounded Knapsack

#include <bits/stdc++.h>
using namespace std;
#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define floord(n, d) (((n) < 0) ? -((-(n) + (d)-1) / (d)) : (n) / (d))
#define MAX 1000


int knapsack(int size, int W, int *wt, int *val) {
    int dp[size + 1][W + 1];
    for (int i = 0; i <= size; i++) {
        for (int j = 0; j <= W; j++) {
            if (i == 0 || j == 0) {
                dp[i][j] = 0;
            }
            else if (wt[i - 1] <= j)
            {
                dp[i][j] = max(dp[i - 1][j], val[i - 1] + dp[i - 1][j - wt[i - 1]]);
            }
            else
            {
                dp[i][j] = dp[i - 1][j];
            }
        }
    }
    for (size_t i = 0; i < size; i++)
    {
        for (size_t j = 0; j < W; j++)
        {
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
    
    return dp[size][W];
}
//		#pragma endscop

int knapsack2(int size, int W, int *wt, int *val) {
    int dp[size + 1][W + 1];
    {
        if (W >= 0)
        {
            for (int c11 = 0; c11 <= size; c11 += 1)
            {
                if (c11 >= 1)
                {
                    dp[c11][0] = 0;
                }
                else
                {
                    for (int c12 = 0; c12 <= W; c12 += 1)
                    {
                        dp[0][c12] = 0;
                    }
                }
            }
        }
        for (int w0 = 0; w0 <= floord(size + W, 16); w0 += 1)
        {
#pragma omp parallel for
            for (int h0 = max(0, w0 - (W + 16) / 16 + 1); h0 <= min(w0, floord(size, 16)); h0 += 1)
            {
                for (int i1 = max(1, 16 * h0); i1 <= min(size, 16 * h0 + 15); i1 += 1)
                {
                    for (int i2 = max(1, 16 * w0 - 16 * h0); i2 <= min(W, 16 * w0 - 16 * h0 + 15); i2 += 1)
                    {
                        if (wt[i1 - 1] <= (i2))
                        {
                            dp[i1][i2] = ((dp[i1 - 1][i2] > (val[i1 - 1] + dp[i1 - 1][i2 - wt[i1 - 1]])) ? dp[i1 - 1][i2] : (val[i1 - 1] + dp[i1 - 1][i2 - wt[i1 - 1]]));
                        }
                        else
                        {
                            dp[i1][i2] = dp[i1 - 1][i2];
                        }
                    }
                }
            }
        }
    }
    for (size_t i = 0; i <= size; i++)
    {
        for (size_t j = 0; j <= W; j++)
        {
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
    
    return dp[size][W];
}

int main()
{
    int n = 100;int m=50;
    int *wt = (int *)malloc((n+1) * sizeof(int));
	int *val = (int *)malloc((n+1) * sizeof(int));
    //seed the random number generator
    time_t time_seed = time(NULL);
    srand(time_seed);
    printf("time_seed = %ld\n", time_seed);
    for (int i = 0; i < n; i++) {
        wt[i] = rand() % 100;
        val[i] = rand() % 100;
    }
    int kval = knapsack(n, m, wt, val);
    printf("Knapsack val: %d\n", kval);
    int kval2 = knapsack2(n, m, wt, val);
    printf("Knapsack2 val: %d\n", kval2);
    return 0;
}

