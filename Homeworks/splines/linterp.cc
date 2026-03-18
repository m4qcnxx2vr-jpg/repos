#include "linterp.h"
#include <stdexcept>
#include <cassert>

int binsearch(const std::vector<double>& x, double z){
    if(x.size() < 2)
        throw std::invalid_argument("binsearch: need at least 2 points");

    if(z < x[0] || z > x[x.size()-1])
        throw std::invalid_argument("binsearch: z out of range");

    int i = 0;
    int j = (int)x.size() - 1;

    while(j - i > 1){
        int mid = (i + j) / 2;
        if(z > x[mid]) i = mid;
        else j = mid;
    }

    return i;
}

double linterp(const std::vector<double>& x,
               const std::vector<double>& y,
               double z){
    if(x.size() != y.size())
        throw std::invalid_argument("linterp: x and y must have same size");

    int i = binsearch(x, z);
    double dx = x[i+1] - x[i];
    assert(dx > 0);

    double dy = y[i+1] - y[i];
    return y[i] + dy/dx * (z - x[i]);
}

double linterpInteg(const std::vector<double>& x,
                    const std::vector<double>& y,
                    double z){
    if(x.size() != y.size())
        throw std::invalid_argument("linterpInteg: x and y must have same size");

    int i = binsearch(x, z);
    double sum = 0.0;

    for(int k = 0; k < i; k++){
        double dx = x[k+1] - x[k];
        double dy = y[k+1] - y[k];
        sum += y[k]*dx + (dy/dx)*dx*dx/2.0;
    }

    double dx = z - x[i];
    double h = x[i+1] - x[i];
    double dy = y[i+1] - y[i];
    sum += y[i]*dx + (dy/h)*dx*dx/2.0;

    return sum;
}