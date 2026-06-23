#pragma once
#include "matrix.h"

namespace pp {

struct qr {
    pp::matrix Q;   // n x m
    pp::matrix R;   // m x m (upper triangular)

    explicit qr(const pp::matrix& A);
    pp::vector solve(const pp::vector b) const; // solve Q R x = b for x

    double det() const; // determinant of R (product of diagonal entries)
    pp::matrix inverse() const; // inverse of A (if square and non-singular), computed as A^{-1} = R^{-1} Q^T

};

} // namespace pp