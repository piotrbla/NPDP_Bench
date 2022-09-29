//rod cutting max product problem
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define MAX 1000

int diceproblem2(int dp[],int s)
{
    int temp=6;
    for (int i = 1; i <= s; i++)
    {
        dp[i]=0;
    }
    dp[s+1]=0;
    dp[0]=1;

#pragma scop
    for (int c0 = 0; c0 < s; c0 += 1)
        for (int c1 = c0 + 1; c1 <= min(s, temp + c0); c1 += 1)
        {
            //(c1, -c0 + c1);
            dp[c1] = dp[c1] + dp[c0];
        }
#pragma endscop
    return dp[s];

}


int main()
{
    int sum=20;
    int dp[sum+2];
    int x = diceproblem2(dp,sum);
    return 0;

}
