#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)


double dotproduct(double *u, double *v, int n) {
    double result = 0;

    for (int i = 0; i < n; i++) {
        result += u[i] * v[i];
    }

    return result;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <n>\n", argv[0]);
        printf("  <n> is the desired test size.\n");
        return 1;
    }

    // Reads the test parameters from the command line
    int N = 0;
    sscanf(argv[1], "%d", &N);

    // Vector creation and initialization
    double *U = (double *)malloc(sizeof(double) * N);
    double *V = (double *)malloc(sizeof(double) * N);

    if (!U || !V) {
        free(U); free(V);
        printf("Error: not enough memory to run the test using n = %i\n", N);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        U[i] = rand() % 10;
        V[i] = rand() % 10;
    }

    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    double result = dotproduct(U, V, N);

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %i\n", N);
    printf("result\t= %.2f\n", result);

    free(U); free(V);

    return 0;
}
