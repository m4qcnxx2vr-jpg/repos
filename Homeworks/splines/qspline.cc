#include "qspline.h"
#include <stdexcept>

qspline::qspline(const vec& xs, const vec& ys)
    : n((int)xs.size()), x(xs), y(ys), b(n-1), c(n-1)
{
    if (xs.size() != ys.size())
        throw std::invalid_argument("qspline: x and y must have same size");

    if (n < 2)
        throw std::invalid_argument("qspline: need at least 2 points");

    for (int i = 1; i < n; i++) {
        if (x[i] <= x[i-1])
            throw std::invalid_argument("qspline: x must be strictly increasing");
    }

    vec h(n-1), p(n-1);
    for(int i = 0; i < n-1; i++){
        h[i] = x[i+1] - x[i];
        p[i] = (y[i+1] - y[i]) / h[i];
    }
    // computing c coefficients
    c[0] = 0;
    for(int i = 0; i < n-2; i++){
        c[i+1] = (p[i+1] - p[i] - c[i]*h[i]) / h[i+1];
    }
    //computing b coefficients
    for(int i = 0; i < n-1; i++){
        b[i] = p[i] - c[i]*h[i];
    }   

   
}

int qspline::binsearch(double z) const {
    if (z < x[0] || z > x[n-1])
        throw std::invalid_argument("binsearch: z out of range");

    int i = 0;
    int j = n - 1;

    while (j - i > 1) {
        int mid = (i + j) / 2;
        if (z >= x[mid]) i = mid;
        else j = mid;
    }

    return i;
}

double qspline::eval(double z) const {
    int i = binsearch(z);
    double dx = z - x[i];
    return y[i] + b[i]*dx + c[i]*dx*dx;
}

double qspline::deriv(double z) const {
    int i = binsearch(z);
    double dx = z - x[i];
    return b[i] + 2*c[i]*dx;
}

double qspline::integ(double z) const {
    int i = binsearch(z);
    double sum = 0.0;

    for(int k = 0; k < i; k++){
        double dx = x[k+1] - x[k];
        sum += y[k]*dx + b[k]*dx*dx/2 + c[k]*dx*dx*dx/3;
    }

    double dx = z - x[i];
    sum += y[i]*dx + b[i]*dx*dx/2 + c[i]*dx*dx*dx/3;

    return sum;
}