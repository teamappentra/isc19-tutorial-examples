ISC19 tutorial examples
=======================

These examples are well-known microbenchmarks that are representative of algorithms frequently used in scientific applications 
and suitable for parallelizing.

| Example                   | Description                                             |
|---------------------------|---------------------------------------------------------|
| [PI](PI/)                 | Approximation of the number PI                          |
| [DOTPRODUCT](DOTPRODUCT/) | Scalar product of two vectors                           |
| [DAXPY](DAXPY/)           | Scalar-Vector product plus Vector (a*X+Y)               |
| [MATMUL](MATMUL/)         | Matrix-Matrix Multiplication using dense linear algebra |
| [MANDELBROT](MANDELBROT/) | Computation of the Mandelbrot sets                      |
| [ATMUX](ATMUX/)           | Sparse Transposed Matrix-Vector Multiplication          |
| [LAPLACE2D](LAPLACE2D/)   | Two-dimensional Laplace equation                        |
| [LULESHmk](LULESHmk/)     | A micro-kernel of the CORAL LULESH benchmark            |

All examples contain a Makefile. Please note that this may need editing (e.g., to change the compiler or compilation flags).
