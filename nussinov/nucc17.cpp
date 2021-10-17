#include <vector>
#include <algorithm>
#include <iostream>
#include <execution>
#include <chrono>
#include <cmath>
#include <tbb/parallel_for.h>
#include <omp.h>



#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))

int N=1000;


int can_pair(char* input, int a, int b) {

    return (((
        (input[a] == 'A' && input[b] == 'U') || (input[a] == 'U' && input[b] == 'A') ||
        (input[a] == 'G' && input[b] == 'C') || (input[a] == 'C' && input[b] == 'G') ||
        (input[a] == 'G' && input[b] == 'U') || (input[a] == 'U' && input[b] == 'G')
        )))/* && (a < b - 1))) */ ? 1 : 0;
}


long double** mem()
{
    int i;
    long double** S = new long double* [N + 5]();
    for (i = 0; i < N + 5; i++)
        S[i] = new long double[N + 5]();

    return S;

}


long double** S;
char* RNA;




int main(int argc, char *argv[])
{
    int method = atoi(argv[1]);  // 0 openmp, 1 c++ par stl, 2 tbb
    int num_proc = atoi(argv[2]);
    N = atoi(argv[3]);




    std::vector<int> nums{ };
    int c0, c1;

    
    nums.clear();

    S = mem();
    RNA = new char[N + 5];

    auto start = std::chrono::steady_clock::now();

if(method==0){
                  int c3, c6, c10, c4;
for(c0 = 0; c0 <= floord(N - 2, 8); c0 += 1)
#pragma omp parallel for schedule(dynamic, 1) shared(c0) private(c1,c3,c4,c6,c10)
   for( c1 = (c0 + 1) / 2; c1 <= min(c0, (N - 1) / 16); c1 += 1)
       for( c3 = 16 * c0 - 16 * c1 + 1; c3 <= min(min(N - 1, 16 * c1 + 15), 16 * c0 - 16 * c1 + 16); c3 += 1) {
         for( c4 = 0; c4 <= c0 - c1; c4 += 1)
           for( c6 = max(-N + 16 * c1 + 1, -N + c3 + 1); c6 <= min(0, -N + 16 * c1 + 16); c6 += 1) {
              for( c10 = 16 * c4; c10 <= min(c3 - 1, 16 * c4 + 15); c10 += 1)
                   S[(-c6)][(c3-c6)] = MAX(S[(-c6)][c10+(-c6)] + S[c10+(-c6)+1][(c3-c6)], S[(-c6)][(c3-c6)]);
                 if (c1 + c4 == c0 && 16 * c0 + c6 + 15 >= 16 * c1 + c3)
                   S[(-c6)][(c3-c6)] = MAX(S[(-c6)][(c3-c6)], S[(-c6)+1][(c3-c6)-1] + can_pair(RNA, (-c6), (c3-c6)));
         }
       for( c4 = max(c0 - c1 + 1, -c1 + (N + c3) / 16 - 1); c4 <= min((N - 1) / 16, -c1 + (N + c3 - 1) / 16); c4 += 1)
          for( c6 = max(max(-N + 16 * c1 + 1, -N + c3 + 1), c3 - 16 * c4 - 15); c6 <= min(-N + 16 * c1 + 16, c3 - 16 * c4); c6 += 1)
             S[(-c6)][(c3-c6)] = MAX(S[(-c6)][(c3-c6)], S[(-c6)+1][(c3-c6)-1] + can_pair(RNA, (-c6), (c3-c6)));
   }
}


if(method==1)
    for (c0 = 0; c0 <= floord(N - 2, 8); c0 += 1) {
        for (c1 = (c0 + 1) / 2; c1 <= min(c0, (N - 1) / 16); c1 += 1) {
            nums.push_back(c1);
	}
	std::for_each(std::execution::par_unseq, nums.cbegin(), nums.cend(), [&](const int& c1) {

                int c3, c6, c10, c4;
                for (c3 = 16 * c0 - 16 * c1 + 1; c3 <= min(min(N - 1, 16 * c1 + 15), 16 * c0 - 16 * c1 + 16); c3 += 1) {
                    for (c4 = 0; c4 <= c0 - c1; c4 += 1)
                        for (c6 = max(-N + 16 * c1 + 1, -N + c3 + 1); c6 <= min(0, -N + 16 * c1 + 16); c6 += 1) {
                            for (c10 = 16 * c4; c10 <= min(c3 - 1, 16 * c4 + 15); c10 += 1)
                                S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][c10 + (-c6)] + S[c10 + (-c6) + 1][(c3 - c6)], S[(-c6)][(c3 - c6)]);
                            if (c1 + c4 == c0 && 16 * c0 + c6 + 15 >= 16 * c1 + c3)
                                S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][(c3 - c6)], S[(-c6) + 1][(c3 - c6) - 1] + can_pair(RNA, (-c6), (c3 - c6)));
                        }
                    for (c4 = max(c0 - c1 + 1, -c1 + (N + c3) / 16 - 1); c4 <= min((N - 1) / 16, -c1 + (N + c3 - 1) / 16); c4 += 1)
                        for (c6 = max(max(-N + 16 * c1 + 1, -N + c3 + 1), c3 - 16 * c4 - 15); c6 <= min(-N + 16 * c1 + 16, c3 - 16 * c4); c6 += 1)
                            S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][(c3 - c6)], S[(-c6) + 1][(c3 - c6) - 1] + can_pair(RNA, (-c6), (c3 - c6)));
               }

        });
	
	nums.clear();
}


if(method==2)
for (c0 = 0; c0 <= floord(N - 2, 8); c0 += 1) {
tbb::parallel_for(tbb::blocked_range<int>((c0 + 1) / 2, min(c0, (N - 1) / 16)+1),
        [&](tbb::blocked_range<int> r)
{
                  int c3, c6, c10, c4;
          for(int c1 = r.begin(); c1 <r.end(); c1++){
                for (c3 = 16 * c0 - 16 * c1 + 1; c3 <= min(min(N - 1, 16 * c1 + 15), 16 * c0 - 16 * c1 + 16); c3 += 1) {
                    for (c4 = 0; c4 <= c0 - c1; c4 += 1)
                        for (c6 = max(-N + 16 * c1 + 1, -N + c3 + 1); c6 <= min(0, -N + 16 * c1 + 16); c6 += 1) {
                            for (c10 = 16 * c4; c10 <= min(c3 - 1, 16 * c4 + 15); c10 += 1)
                                S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][c10 + (-c6)] + S[c10 + (-c6) + 1][(c3 - c6)], S[(-c6)][(c3 - c6)]);
                            if (c1 + c4 == c0 && 16 * c0 + c6 + 15 >= 16 * c1 + c3)
                                S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][(c3 - c6)], S[(-c6) + 1][(c3 - c6) - 1] + can_pair(RNA, (-c6), (c3 - c6)));
                        }
                    for (c4 = max(c0 - c1 + 1, -c1 + (N + c3) / 16 - 1); c4 <= min((N - 1) / 16, -c1 + (N + c3 - 1) / 16); c4 += 1)
                        for (c6 = max(max(-N + 16 * c1 + 1, -N + c3 + 1), c3 - 16 * c4 - 15); c6 <= min(-N + 16 * c1 + 16, c3 - 16 * c4); c6 += 1)
                            S[(-c6)][(c3 - c6)] = MAX(S[(-c6)][(c3 - c6)], S[(-c6) + 1][(c3 - c6) - 1] + can_pair(RNA, (-c6), (c3 - c6)));
                }
             }
        });

}

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


}




