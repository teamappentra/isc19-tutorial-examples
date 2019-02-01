#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)

typedef struct Matrix2D Matrix2D;

Matrix2D *Matrix2D_new(int rows, int cols);

void Matrix2D_delete(Matrix2D *mat);

Matrix2D *Matrix2D_rand(Matrix2D *mat);

double **Matrix2D_getData(Matrix2D *mat);

double Matrix2D_checksum(Matrix2D *mat);


void matmul(int m, int n, int p, double **A, double **B, double **C) {
    // Initialization
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
        }
    }

    // Accumulation
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < p; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    int param_iters = 1;

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
    Matrix2D *in1_mat = Matrix2D_new(param_n, param_n);
    Matrix2D *in2_mat = Matrix2D_new(param_n, param_n);
    Matrix2D *out_mat = Matrix2D_new(param_n, param_n);
    if (!in1_mat || !in2_mat || !out_mat) {
        printf("Error: not enough memory to run the test using n = %i\n", param_n);
        return 1;
    }

    // Initializes data
    Matrix2D_rand(in1_mat);
    Matrix2D_rand(in2_mat);

    // Calls to the corresponding function to perform the computation
    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    for (int iters = 0; iters < param_iters; iters++) {
        matmul(
            param_n,
            param_n,
            param_n,
            Matrix2D_getData(in1_mat),
            Matrix2D_getData(in2_mat),
            Matrix2D_getData(out_mat));
    }

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    double checksum = Matrix2D_checksum(out_mat);
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %i\n", param_n);
    printf("chksum\t= %.0f\n", checksum);
    if (param_iters > 1)
        printf("iters\t= %i\n", param_iters);

    // Release allocated resources
    Matrix2D_delete(in1_mat);
    Matrix2D_delete(in2_mat);
    Matrix2D_delete(out_mat);

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

// Generates a random dense matrix
Matrix2D *Matrix2D_rand(Matrix2D *mat) {
    if (!mat)
        return 0;
    for (int row = 0; row < mat->rows; row++)
        for (int col = 0; col < mat->cols; col++)
            mat->data[row][col] = rand() % 10;
    return mat;
}

// Get raw pointer to matrix data (as list of pointers to rows)
double **Matrix2D_getData(Matrix2D *mat) { return mat->data; }

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
