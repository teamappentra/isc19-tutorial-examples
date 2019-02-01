#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)

typedef struct Matrix2D Matrix2D;

Matrix2D *Matrix2D_new(int rows, int cols);

void Matrix2D_delete(Matrix2D *mat);

Matrix2D *Matrix2D_zero(Matrix2D *mat);

double **Matrix2D_getData(Matrix2D *mat);

double Matrix2D_checksum(Matrix2D *mat);


void laplace_init(double **A, int n, int m) {
    float pi = 2.0f * asin(1.0f);
    for (int j = 0; j < n; j++) {
        double tmp = sin(pi * j / (n - 1));
        A[j][0] = tmp;
        A[j][m - 1] = tmp * exp(-pi);
    }
}

void compute(double *err, double **A, double **T, int n, int m) {
    double error = 0;
    for (int j = 1; j < n - 1; j++) {
        for (int i = 1; i < m - 1; i++) {
            T[j][i] = 0.25f * (A[j][i + 1] + A[j][i - 1] + A[j - 1][i] + A[j + 1][i]);
            double dif = fabs(T[j][i] - A[j][i]);
            error = fmaxf(error, dif);
        }
    }
    *err = error;
}

void backup(double **A, double **T, int n, int m) {
    for (int j = 1; j < m - 1; j++)
        for (int i = 1; i < n - 1; i++)
            A[j][i] = T[j][i];
}

int laplace(double **A, double **T, int n, int m, double tol, int iter_max) {
    int iter = 0;
    double err;
    for (iter = 0, err = tol; err >= tol && iter < iter_max; iter++) {
        compute(&err, A, T, n, m);
        backup(A, T, n, m);
    }
    return iter;
}

int main(int argc, char *argv[]) {
    int param_iters = 99999;

    if (argc != 2) {
        printf("Usage: %s <n> \n", argv[0]);
        printf("  <n> is the desired test size.\n");
        return 1;
    }

    // Reads the test parameters from the command line
    int param_n = 0;
    sscanf(argv[1], "%d", &param_n);
    printf("- Input parameters\n");
    printf("n\t= %i\n", param_n);

    // Allocates input/output resources
    Matrix2D *inoutA_mat = Matrix2D_new(param_n, param_n);
    Matrix2D *inT_mat = Matrix2D_new(param_n, param_n);
    if (!inoutA_mat || !inT_mat) {
        printf("Error: not enough memory to run the test using n = %i\n", param_n);
        return 1;
    }

    // Initializes data if needed
    double param_tol = 1.0e-5;
    Matrix2D_zero(inoutA_mat);
    laplace_init(Matrix2D_getData(inoutA_mat), param_n, param_n);

    // Calls to the corresponding function to perform the computation
    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    param_iters = laplace(
        Matrix2D_getData(inoutA_mat),
        Matrix2D_getData(inT_mat),
        param_n,
        param_n,
        param_tol,
        param_iters);

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    double checksum = Matrix2D_checksum(inoutA_mat);
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %i\n", param_n);
    printf("chksum\t= %.6f\n", checksum);
    if (param_iters > 1)
        printf("iters\t= %i\n", param_iters);

    // Release allocated resources
    Matrix2D_delete(inoutA_mat);
    Matrix2D_delete(inT_mat);

    return 0;
}

// Matrix structure
struct Matrix2D {
    int rows;
    int cols;
    long long size;
    double **data;
};

// Creates a new dense matrix with the specified rows and columns
Matrix2D *Matrix2D_new(int rows, int cols) {
    if (rows < 1 || cols < 1)
        return 0;
    Matrix2D *_this = (Matrix2D *)malloc(sizeof(Matrix2D));
    if (!_this)
        return 0;

    _this->rows = rows;
    _this->cols = cols;
    _this->size = rows * cols;
    // Last row pointer is always null
    _this->data = (double **)calloc(rows + 1, sizeof(double));
    if (_this->data == 0)
        return 0;

    size_t matBytes = cols * rows * sizeof(double);
    double *memPtr = (double *)malloc(matBytes);
    if (memPtr) {
        for (int i = 0; i < rows; i++)
            _this->data[i] = memPtr + i * cols;
        return _this;
    }

    // on memory allocation error
    if (_this->data)
        free(_this->data);
    if (_this)
        free(_this);
    return 0;
}

// Deletes the matrix and the resources allocated by it
void Matrix2D_delete(Matrix2D *mat) {
    if (mat && mat->data) {
        if (mat->data[0])
            free(mat->data[0]);
        free(mat->data);
    }
    if (mat)
        free(mat);
}

// Get raw pointer to matrix data (as list of pointers to rows)
double **Matrix2D_getData(Matrix2D *mat) { return mat->data; }

// Generates an empty matrix
Matrix2D *Matrix2D_zero(Matrix2D *mat) {
    if (!mat)
        return 0;
    size_t nBytes = sizeof(double) * mat->size;
    memset(mat->data[0], 0, nBytes);
    return mat;
}

// Generates a checksum based on the matrix data
double Matrix2D_checksum(Matrix2D *mat) {
    if (!mat)
        return 0;

    double checkSum = 0.0;
    for (int row = 0; row < mat->rows; row++)
        for (int col = 0; col < mat->cols; col++)
            checkSum += mat->data[row][col];
    return checkSum;
}
