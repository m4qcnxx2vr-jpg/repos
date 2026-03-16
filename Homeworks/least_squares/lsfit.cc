#include "lsfit.h"
#include <stdexcept>
#include <tuple>

namespace pp {

static pp::vector backsub(const pp::matrix& R, const pp::vector& b){
    int n = R.size1();

    if(R.size2() != n || b.size() != n)
        throw std::invalid_argument("backsub: size mismatch");

    pp::vector x(n, 0.0);

    for(int i = n - 1; i >= 0; i--){
        double s = b[i];
        for(int j = i + 1; j < n; j++)
            s -= R(i,j) * x[j];

        if(R(i,i) == 0.0)
            throw std::runtime_error("backsub: zero diagonal");

        x[i] = s / R(i,i);
    }

    return x;
}

std::tuple<pp::vector, pp::matrix> lsfit(
    const std::vector<Func>& fs,
    const pp::vector& x,
    const pp::vector& y,
    const pp::vector& dy
){
    int n = x.size();
    int m = (int)fs.size();

    if(y.size() != n || dy.size() != n)
        throw std::invalid_argument("lsfit: x, y, dy must have same size");

    if(m == 0)
        throw std::invalid_argument("lsfit: no functions given");

    if(n < m)
        throw std::invalid_argument("lsfit: need n >= m");

    pp::matrix A(n, m, 0.0);
    pp::vector b(n, 0.0);

    for(int i = 0; i < n; i++){
        if(dy[i] == 0.0)
            throw std::invalid_argument("lsfit: dy must be nonzero");

        b[i] = y[i] / dy[i];

        for(int k = 0; k < m; k++)
            A(i,k) = fs[k](x[i]) / dy[i];
    }

    pp::qr decomp(A);
    pp::vector c = decomp.solve(b);

    pp::matrix Rinv(m, m, 0.0);

    for(int j = 0; j < m; j++){
        pp::vector e(m, 0.0);
        e[j] = 1.0;
        pp::vector col = backsub(decomp.R, e);
        Rinv[j] = col;
    }

    pp::matrix cov = Rinv * Rinv.T();

    return {c, cov};
}

} // namespace pp