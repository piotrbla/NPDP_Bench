#define max(a,b) (((a)>(b))?(a):(b))

int LDS2(int a[],int size){
    int dp[size+1],res=1;
    for(int i=0;i<size;i++){
        dp[i] = 1;
    }
	
	#pragma scop
    for (int c0 = 1; c0 < size; c0 += 1)
        for (int c1 = size - c0; c1 < size; c1 += 1)
        {
            if(a[c1]>a[size - c0 - 1]){
                dp[c1] = max(dp[c1], 1+dp[size - c0 - 1]);
            }
            //(c1, size - c0 - 1);
        }
	#pragma endscop
	    
    return res;
}

int main(){
    int a[] = {10,9,2,5,3,7,101,18, 11, 15, 13};
    int size = sizeof(a)/sizeof(a[0]);
    int x = LDS2(a,size);
}
