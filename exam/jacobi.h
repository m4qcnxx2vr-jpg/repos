#pragma once

#include "matrix.h"
#include <cmath>

namespace pp {

struct Jacobi {
    vector w;   // eigenvalues
    matrix V;   // eigenvectors stored as columns

    explicit Jacobi(matrix A, double tolerance = 1e-12);

    static void timesJ(matrix& A, int p, int q, double theta);
    static void Jtimes(matrix& A, int p, int q, double theta);

    static int cyclic(matrix& A, vector& w, matrix& V, double tolerance = 1e-12);
};

} // namespace pp