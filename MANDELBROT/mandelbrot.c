#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define getClock() ((double)clock() / CLOCKS_PER_SEC)

typedef struct Matrix2D Matrix2D;

Matrix2D *Matrix2D_new(int rows, int cols);

void Matrix2D_delete(Matrix2D *mat);

double **Matrix2D_getData(Matrix2D *mat);

double Matrix2D_checksum(Matrix2D *mat);


int mandelbrot(int max_iter, int height, int width, double **output, double real_min, double real_max, double imag_min,
               double imag_max) {
    double scale_real = (real_max - real_min) / width;
    double scale_imag = (imag_max - imag_min) / height;

    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {

            double x0 = real_min + col * scale_real;
            double y0 = imag_min + row * scale_imag;

            double y = 0, x = 0;
            int iter = 0;
            while (x * x + y * y < 4 && iter < max_iter) {
                double xtemp = x * x - y * y + x0;
                y = 2 * x * y + y0;
                x = xtemp;
                iter++;
            }
            output[row][col] = iter;
        }
    }
    return 0;
}

int main(int argc, char *argv[]) {
    int param_iters = 100;
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
    Matrix2D *out_mat = Matrix2D_new(param_n, param_n);
    if (!out_mat) {
        printf("Error: not enough memory to run the test using n = %i\n", param_n);
        return 0;
    }

    // Calls to the corresponding function to perform the computation
    printf("- Executing test...\n");
    double time_start = getClock();
    // ================================================

    mandelbrot(
        param_iters,
        param_n,
        param_n,
        Matrix2D_getData(out_mat),
        -1.4, 0.5, -0.95, 0.95);

    // ================================================
    double time_finish = getClock();

    // Prints an execution report
    double checksum = Matrix2D_checksum(out_mat);
    printf("time (s)= %.6f\n", time_finish - time_start);
    printf("size\t= %i\n", param_n);
    printf("chksum\t= %.0f\n", checksum);
    printf("iters\t= %i\n", param_iters);

    // Release allocated resources
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
