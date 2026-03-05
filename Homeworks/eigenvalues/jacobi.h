#pragma once
#include "matrix.h"
#include <cmath>

namespace pp {

struct Jacobi {

    vector w;   // eigenvalues
    matrix V;   // eigenvectors (stored as columns)

    // constructor: computes eigenvalue decomposition of A
    explicit Jacobi(const matrix& A, double tol = 1e-12);

    // A <- A * J(p,q,theta)  (column rotation)
    static void timesJ(matrix& A, int p, int q, double theta);

    // A <- J(p,q,theta)^T * A  (row rotation)
    static void Jtimes(matrix& A, int p, int q, double theta);

    // Jacobi diagonalization
    static int cyclic(matrix& A, vector& w, matrix& V, double tol = 1e-12);
};

} // namespace pp