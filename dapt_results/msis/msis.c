
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

int MSIS3(int a[],int size){
    int dp[size];
    for(int i=0;i<size;i++){
        dp[i] = a[i];
    }

#pragma scop

    for (int c0 = 0; c0 < size - 1; c0 += 1)
        for (int c1 = c0 + 1; c1 < size; c1 += 1)
        { //(c1, c0);
            if (a[c1] > a[c0] && dp[c1] < a[c1] + dp[c0])
                dp[c1] = a[c1] + dp[c0];
        }

#pragma endscop

    int max = dp[0];
    for(int i=1;i<size;i++){
        if(dp[i]>max){
            max = dp[i];
        }
    }
    return max;
}

int main(){
    int a[] = {1, 101, 2, 3, 15, 19, 93, 100, 4, 5, 7, 12, 123};
    int size = sizeof(a)/sizeof(a[0]);//13
    int x = MSIS3(a,size);
	return x;
}

