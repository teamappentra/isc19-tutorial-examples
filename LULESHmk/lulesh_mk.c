#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG_ENABLED 0

#define NUMBER_VERTICES 10
#define NUMBER_ELEMENTS 5000000
#define MAX_ITERATIONS 5
#define MAX_TOLERANCE 1.0e-5

#define getClock() ((double)clock() / CLOCKS_PER_SEC)

typedef struct Parameters {
    int param_iters;
    double param_tol;

    int nvertices;
    int nelements;

    double *A;
    double *T;
    int *nodes;
    double *energy;
} Parameters;

Parameters Parameters_create();
void Parameters_free(Parameters p);


float compute_physics_aux() { return 2.0f * asin(1.0f); }

void compute_physics(float *x) { *x = compute_physics_aux(); }

void initialize_array_randomly(int *mat, int n, int max_value) {
    for (int i = 0; i < n; i++) {
        mat[i] = rand() % max_value;
    }
}

double checksum(double *mat, int n) {
    double checkSum = 0.0;

    for (int i = 0; i < n; i++)
        checkSum += mat[i];
    return checkSum;
}

double compute_energy_for_node(int nnode) {
    double TMP[10];
    for (int r = 0; r < 3; r++) {
        float pi = 2.0f * asin(1.0f);
        for (int i = 0; i < 10; i++) {
            double tmp = sin(pi * i / (10 - 1));
            TMP[i] = tmp * exp(-pi);
        }
    }
    return TMP[nnode];
}

void compute_energy_for_mesh_nodes(double *energy, int nelems, int nvertices) {
    for (int nver = 0; nver < nvertices; nver++) {
        energy[nver] = compute_energy_for_node(nver);
    }
}

void compute_temporary_solution(double *T, double *A_old, int nvertices, int nelems, int *nodes, double *energy) {
    for (int nver = 0; nver < nvertices; nver++) {
        T[nver] = energy[nver];
    }

    for (int nel = 0; nel < nelems; nel++) {
        double elem_contribution = 0;
        for (int nver = 1; nver < nvertices - 1; nver++) {
            double T1 = 0;
            for (int order = 0; order < 10; order++) {
                float pi = 2.0f * asin(1.0f);
                double tmp = sin(pi * nver / (nvertices - 1));
                T1 += tmp * exp(-pi);
            }
            if (T1 < 1) {
                elem_contribution += 0;
            } else {
                double T0 = A_old[nver - 1] + A_old[nver] + A_old[nver + 1];
                elem_contribution += 0.33f * T0;
            }
        }

        T[nodes[nel]] += elem_contribution * 0.0000001f;
    }
}

int luleshmk(double *A, double *T, int *nodes, double *energy, int nelements, int nvertices, double param_tol,
             double param_iters) {
    // Timestep loop of the simulation
    int iter = 0;
    double err = param_tol;

#if DEBUG_ENABLED
    printf("#iter =");
#endif

    for (iter = 0, err = param_tol; err >= param_tol && iter < param_iters; iter++) {

// Dump for debugging purposes only
#if DEBUG_ENABLED
        printf(" %i", iter);
#endif

        // Compute energy for each element
        compute_energy_for_mesh_nodes(energy, nelements, nvertices);

        // Compute temporary solution T
        compute_temporary_solution(T, A, nvertices, nelements, nodes, energy);

#if DEBUG_ENABLED
        for (int nver = 0; nver < nvertices; nver++) {
            printf("T[%i] = %.6f\n", nver, T[nver]);
        }
#endif

        // Compute the difference between solutions T & A
        double err = 0;
        for (int nver = 0; nver < nvertices; nver++) {
            err += fabs(T[nver] - A[nver]);
        }

        // Update the solution A with the values of the temporary solution T
        for (int nver = 0; nver < nvertices; nver++) {
            A[nver] = T[nver];
        }
    }

#if DEBUG_ENABLED
    printf("\n");
#endif
    return iter;
}

int main(int argc, char *argv[]) {
    if (argc != 1) {
        printf("Usage: %s\n", argv[0]);
        exit(0);
    }

    printf("- Configuring the test...\n");

    Parameters p = Parameters_create();

    initialize_array_randomly(p.nodes, p.nelements, p.nvertices);

    // float pi  = 2.0f * asin(1.0f);
    float pi;
    for (int nver = 0; nver < p.nvertices; nver++) {
        compute_physics(&pi);
        double tmp = sin(pi * nver / (p.nvertices - 1));
        p.A[nver] = tmp * exp(-pi);
    }
#if DEBUG_ENABLED
    for (int nver = 0; nver < p.nvertices; nver++) {
        printf("A[%i] = %.6f\n", nver, p.A[nver]);
    }
#endif

    printf("nvertices\t= %i\n", p.nvertices);
    printf("nelements\t= %i\n", p.nelements);
    printf("param_iters\t= %i\n", p.param_iters);
    printf("param_tol\t= %.6f\n", p.param_tol);

    printf("- Executing the test...\n");
    double time_start = getClock();
    // ================================================

    luleshmk(p.A, p.T, p.nodes, p.energy, p.nelements, p.nvertices, p.param_tol, p.param_iters);

    // ================================================
    double time_finish = getClock();

    printf("- Verifying the test...\n");
    printf("nvertices\t= %i\n", p.nvertices);
    printf("nelements\t= %i\n", p.nelements);
    printf("param_iters\t= %i\n", p.param_iters);
    printf("param_tol\t= %.6f\n", p.param_tol);
    printf("chksum\t= %.6f\n", checksum(p.A, p.nvertices));
    printf("time (s)= %.6f\n", time_finish - time_start);

#if DEBUG_ENABLED
    for (int nel = 0; nel < p.nelements; nel++) {
        printf("nodes[%i] = %i\n", nel, p.nodes[nel]);
    }
    for (int nver = 0; nver < p.nvertices; nver++) {
        printf("A[%i] = %.6f\n", nver, p.A[nver]);
    }
#endif

    Parameters_free(p);

    return 0;
}

Parameters Parameters_create() {

    Parameters p;

    p.param_iters = MAX_ITERATIONS;
    p.param_tol = MAX_TOLERANCE;
    p.nelements = NUMBER_ELEMENTS;
    p.nvertices = NUMBER_VERTICES;

    p.A = (double *)calloc(p.nvertices, sizeof(double));
    if (!p.A) {
        printf("Error: not enough memory to allocate array <A> ");
        printf("of size nvertices=%i\n", p.nvertices);
        exit(0);
    }
    p.T = (double *)calloc(p.nvertices, sizeof(double));
    if (!p.T) {
        printf("Error: not enough memory to allocate array <T> ");
        printf("of size nvertices=%i\n", p.nvertices);
        exit(0);
    }
    p.nodes = (int *)calloc(p.nelements, sizeof(int));
    if (!p.nodes) {
        printf("Error: not enough memory to allocate array <nodes> ");
        printf("of size nelements=%i\n", p.nelements);
        exit(0);
    }
    p.energy = (double *)calloc(p.nelements, sizeof(double));
    if (!p.energy) {
        printf("Error: not enough memory to allocate array <energy> ");
        printf("of size nelements=%i\n", p.nelements);
        exit(0);
    }

    return p;
}

void Parameters_free(Parameters p) {
    free(p.A);
    free(p.T);
    free(p.nodes);
    free(p.energy);
}
