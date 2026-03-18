#include "fqinterp.h"
#include <stdexcept>

std::function<double(double)> make_qspline(
    std::vector<double> x,
    std::vector<double> y
){
    int n = (int)x.size();

    if((int)y.size() != n)
        throw std::invalid_argument("make_qspline: x and y must have same size");

    if(n < 2)
        throw std::invalid_argument("make_qspline: need at least 2 points");

    for(int i = 1; i < n; i++){
        if(x[i] <= x[i-1])
            throw std::invalid_argument("make_qspline: x must be strictly increasing");
    }

    std::vector<double> b(n-1), c(n-1), h(n-1), p(n-1);

    for(int i = 0; i < n-1; i++){
        h[i] = x[i+1] - x[i];
        p[i] = (y[i+1] - y[i]) / h[i];
    }

    c[0] = 0;
    for(int i = 0; i < n-2; i++){
        c[i+1] = (p[i+1] - p[i] - c[i]*h[i]) / h[i+1];
    }

    for(int i = 0; i < n-1; i++){
        b[i] = p[i] - c[i]*h[i];
    }

    return [x = std::move(x), y = std::move(y), b = std::move(b), c = std::move(c)]
           (double z){
        int n = (int)x.size();

        if(z < x[0] || z > x[n-1])
            throw std::invalid_argument("qspline function: z out of range");

        int i = 0;
        int j = n - 1;
        while(j - i > 1){
            int mid = (i + j) / 2;
            if(z > x[mid]) i = mid;
            else j = mid;
        }

        double dx = z - x[i];
        return y[i] + b[i]*dx + c[i]*dx*dx;
    };
}