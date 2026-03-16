#include "qr.h"
#include <stdexcept>

namespace pp {

qr::qr(const pp::matrix& A) {
    int n = A.size1();
    int m = A.size2();
    if (n < m) throw std::invalid_argument("qr: require n >= m");

    Q = A;
    R = pp::matrix(m, m, 0.0);

    for (int i = 0; i < m; i++) {
        double rii = Q[i].norm();
        if (rii == 0.0)
            throw std::runtime_error("qr: dependent columns (zero norm)");

        R(i, i) = rii;
        Q[i] /= rii;

        for (int j = i + 1; j < m; j++) {
            double rij = Q[i].dot(Q[j]);
            R(i, j) = rij;
            Q[j] -= Q[i] * rij;
        }
    }
}

pp::vector qr::solve(const pp::vector& b) const {
    int n = Q.size1();
    int m = Q.size2();

    if (b.size() != n)
        throw std::invalid_argument("qr::solve: b has wrong size");

    pp::vector y(m, 0.0);
    for (int i = 0; i < m; i++)
        y[i] = Q[i].dot(b);

    pp::vector x(m, 0.0);
    for (int i = m - 1; i >= 0; i--) {
        double s = y[i];
        for (int j = i + 1; j < m; j++)
            s -= R(i, j) * x[j];

        double rii = R(i, i);
        if (rii == 0.0)
            throw std::runtime_error("qr::solve: singular R (zero diagonal)");

        x[i] = s / rii;
    }

    return x;
}

double qr::det() const {
    int m = R.size1();
    if (m != R.size2())
        throw std::runtime_error("qr::det: R is not square");

    double d = 1.0;
    for (int i = 0; i < m; i++)
        d *= R(i, i);
    return d;
}

pp::matrix qr::inverse() const {
    int n = Q.size1();
    int m = Q.size2();

    if (n != m)
        throw std::invalid_argument("qr::inverse: A must be square (n==m)");

    pp::matrix Ainv(n, n, 0.0);

    for (int j = 0; j < n; j++) {
        pp::vector e(n, 0.0);
        e[j] = 1.0;

        pp::vector x = solve(e);
        Ainv[j] = x;
    }

    return Ainv;
}

} // namespace pp