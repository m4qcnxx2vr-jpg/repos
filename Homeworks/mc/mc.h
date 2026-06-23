#pragma once

#include <cmath>
#include <cstdint>
#include <functional>
#include <random>
#include <tuple>
#include <vector>

namespace pp {

using vector = std::vector<double>;
using tuple = std::tuple<double, double>;
using Func = std::function<double(vector)>;

struct lcg {
    uint32_t seed = 1;
    uint32_t a = 1664525;
    uint32_t c = 1013904223;
    uint64_t m = uint64_t{1} << 32;

    lcg() = default;

    lcg(uint32_t seed_in)
        : seed(seed_in)
    {}

    double uniform(double low = 0.0, double high = 1.0);
};

struct halton {
    int dim = 0;
    std::vector<int> bases;

    halton() = default;
    halton(int dim_in);

    vector point(int index);
    double van_der_corput(int index, int base);
    std::vector<int> first_primes(int n);
};

double randd(double low, double high);

double box_volume(const vector& a, const vector& b);

tuple plainmc(Func f, const vector& a, const vector& b, int N);

tuple lcgmc(Func f, const vector& a, const vector& b, int N, uint32_t seed = 42);

tuple quasimc(Func f, const vector& a, const vector& b, int N);

tuple stratifiedmc(Func f, const vector& a, const vector& b, int N, int nmin = 100);

}