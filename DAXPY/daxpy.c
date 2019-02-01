#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)


void daxpy(int n, double *D, double a, double *X, double *Y) {
    for (int i = 0; i < n; i++) {
        D[i] = a * X[i] + Y[i];
    }
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
    printf("- Input parameters\n");
    printf("n\t= %i\n", N);

    // Vector creation and initialization
    float a = 2;
    double *D = (double *)malloc(sizeof(double) * N);
    double *X = (double *)malloc(sizeof(double) * N);
    double *Y = (double *)malloc(sizeof(double) * N);

    if (!D || !X || !Y) {
        free(D); free(X); free(Y);
        printf("Error: not enough memory to run the test using n = %i\n", N);
        return 1;
    }

    for (int i = 0; i < N; i++) {
        D[i] = 0;
        X[i] = rand() % 10;
        Y[i] = rand() % 10;
    }

    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    daxpy(N, D, a, X, Y);

    // ================================================
    double time_finish = getClock();

    double checkSum = 0.0;
    for (int i = 0; i < N; i++)
        checkSum += D[i];

    // Prints an execution report
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("chksum\t= %.2f\n", checkSum);

    free(D); free(X); free(Y);

    return 0;
}
