#pragma once

#include "matrix.h"
#include <cmath>

namespace pp {

struct SVD {
    matrix U;
    matrix D;
    matrix V;
    vector s; // singular values

    explicit SVD(matrix A, double tolerance = 1e-12);

    static void timesJ(matrix& A, int p, int q, double theta);
    static void Jtimes(matrix& A, int p, int q, double theta);

    static int two_sided(matrix& A, matrix& U, matrix& V, vector& s,
                         double tolerance = 1e-12);
};

} 