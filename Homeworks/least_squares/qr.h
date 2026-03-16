#pragma once
#include "matrix.h"

namespace pp {

struct qr {
    pp::matrix Q;   // n x m
    pp::matrix R;   // m x m

    explicit qr(const pp::matrix& A);
    pp::vector solve(const pp::vector& b) const;

    double det() const;
    pp::matrix inverse() const;
};

} // namespace pp