//rod cutting max product problem
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define MAX 1000



int findLis2(int *a, int n)
{
    int lis[100],i,j,max_val=0;

    for(i=0;i<n;i++)
        lis[i]=1;

#pragma scop
    for (int c0 = 0; c0 < n - 1; c0 += 1)
    {
        for (int c1 = c0 + 1; c1 < n; c1 += 1)
        { //(c1, c0);
            if (a[c1] > a[c0] && lis[c1] < (lis[c0] + 1))
                lis[c1]++;
        }
    }
#pragma endscop

    for(i=0;i<n;i++)
    {
        if(lis[i]>max_val)
            max_val=lis[i];
    }
    return max_val;
}

int main()
{
    int n = 100;
    int *arr;
    int max_val2 = findLis2(arr, n);
    return 0;
}

