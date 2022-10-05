
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define MAX 1000

int findCost3(int arr[], int len)
{
    int dp[len][len];
    int q = 0;
    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            dp[i][j] = 0;
        }
    }
	#pragma scop
    {
        for (int c1 = 2; c1 < len; c1 += 1)
            for (int c2 = 0; c2 < len - c1; c2 += 1)
            {
                // S_0(c1, c2);
                dp[c2][c2 + c1] = 2147483647;
            }
        for (int c0 = 2; c0 < len; c0 += 1)
            for (int c1 = c0; c1 < min(len, 2 * c0 - 1); c1 += 1)
                for (int c2 = 0; c2 < len - c1; c2 += 1)
                {
                    if (2 * c0 >= c1 + 3)
                    {
                        //S_1(c1, c2, -c0 + c1 + c2 + 1);
                        //    l,  i,  k
                        if ((dp[c2][-c0 + c1 + c2 + 1] + dp[-c0 + c1 + c2 + 1][c2 + c1] + arr[c2] * arr[c2 + c1] * arr[-c0 + c1 + c2 + 1]) < dp[c2][c2 + c1])
                            dp[c2][c2 + c1] = dp[c2][-c0 + c1 + c2 + 1] + dp[-c0 + c1 + c2 + 1][c2 + c1] + arr[c2] * arr[c2 + c1] * arr[-c0 + c1 + c2 + 1];
                    }
                    //S_1(c1, c2, c0 + c2 - 1);
                    if ((dp[c2][c0 + c2 - 1] + dp[c0 + c2 - 1][c2 + c1] + arr[c2] * arr[c2 + c1] * arr[c0 + c2 - 1]) < dp[c2][c2 + c1])
                        dp[c2][c2 + c1] = dp[c2][c0 + c2 - 1] + dp[c0 + c2 - 1][c2 + c1] + arr[c2] * arr[c2 + c1] * arr[c0 + c2 - 1];
                }
    }
#pragma endscop
    return dp[0][len - 1];
}


int main()
{
    int n, a[MAX], b[MAX];
    //scanf("%d", &n);
    n=10;
 
 
    int res2 =  findCost3(a, n);
 
    return 0;
}