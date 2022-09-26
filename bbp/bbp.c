
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define MAX 1000

int lis2(int arr[], int n)
{
    int dp[n];
    for(int i=0; i<n; i++) dp[i] = 1;

#pragma scop
    for (int c0 = 0; c0 < n - 1; c0 += 1)
        for (int c1 = c0 + 1; c1 < n; c1 += 1)
        {
            //(c1, c0);
            if (arr[c1] >= arr[c0] && dp[c1] < dp[c0] + 1)
                dp[c1] = dp[c0] + 1;
        }
#pragma endscop
    int ans = dp[0];

    for (int i = 1; i < n; i++)
    {
        if(dp[i] > ans) ans = dp[i];
    }

    
    return ans;
}




int main()
{
    int n, a[MAX], b[MAX];
    //scanf("%d", &n);
    n=10;
 
 
    int res2 =  lis2(a, n);
 
    return 0;
}