#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)


int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("Usage: %s <steps>\n", argv[0]);
        printf("  <steps> controls the precision of the approximation.\n");
        return 0;
    }

    // Reads the test parameters from the command line
    int N = atoi(argv[1]);
    printf("- Input parameters\n");
    printf("steps\t= %i\n", N);

    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    double out_result;

    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double x = (i + 0.5) / N;
        sum += sqrt(1 - x * x);
    }

    out_result = 4.0 / N * sum;

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("result\t= %.8f\n", out_result);
    const double realPiValue = 3.141592653589793238;
    printf("error\t= %.1e\n", fabs(out_result - realPiValue));

    return 0;
}
