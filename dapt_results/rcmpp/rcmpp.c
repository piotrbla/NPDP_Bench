//rod cutting max product problem
#define min(x,y)    ((x) < (y) ? (x) : (y))
#define max(x,y)    ((x) > (y) ? (x) : (y))
#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
#define MAX 1000






void rcmpp2(int result[], int n)
{
#pragma scop
    for (int c0 = 1; c0 < n; c0 += 1)
        for (int c1 = c0 + 1; c1 <= min(n, 2 * c0); c1 += 1)
        {
            //(c1, (-c0 + c1);
            result[c1]=max(result[c1],max((-c0 + c1)*result[c1-(-c0 + c1)],(-c0 + c1)*(c1-(-c0 + c1))));
        }
#pragma endscop
}


int main()
{
	int n=21;
    int result[n],i,j;
	for(i=0;i<n;i++)
        result[i]=0;

    rcmpp2(result,n);

	
    return 0;
}