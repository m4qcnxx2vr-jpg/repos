#include "mc.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace pp {

double lcg::uniform(double low, double high)
{
    seed = a * seed + c;

    double u = (static_cast<double>(seed) + 1.0)
             / (static_cast<double>(m) + 1.0);

    return low + u * (high - low);
}

halton::halton(int dim_in)
{
    dim = dim_in;
    bases = first_primes(dim);
}

double halton::van_der_corput(int index, int base)
{
    double q = 0.0;
    double factor = 1.0 / static_cast<double>(base);

    while (index > 0) {
        q += static_cast<double>(index % base) * factor;
        index /= base;
        factor /= static_cast<double>(base);
    }

    return q;
}

vector halton::point(int index)
{
    vector x(dim);

    for (int i = 0; i < dim; ++i) {
        x[i] = van_der_corput(index, bases[i]);
    }

    return x;
}

std::vector<int> halton::first_primes(int n)
{
    std::vector<int> primes;

    int candidate = 2;

    while (static_cast<int>(primes.size()) < n) {
        bool is_prime = true;

        for (int p : primes) {
            if (p * p > candidate) {
                break;
            }

            if (candidate % p == 0) {
                is_prime = false;
                break;
            }
        }

        if (is_prime) {
            primes.push_back(candidate);
        }

        ++candidate;
    }

    return primes;
}

double randd(double low, double high)
{
    static std::mt19937_64 gen(12345);
    std::uniform_real_distribution<double> dist(low, high);

    return dist(gen);
}

double box_volume(const vector& a, const vector& b)
{
    double V = 1.0;

    for (int i = 0; i < static_cast<int>(a.size()); ++i) {
        V *= b[i] - a[i];
    }

    return V;
}

tuple plainmc(Func f, const vector& a, const vector& b, int N)
{
    if (N <= 0) {
        return {0.0, 0.0};
    }

    int dim = static_cast<int>(a.size());
    double V = box_volume(a, b);

    double sum = 0.0;
    double sumsq = 0.0;

    for (int k = 0; k < N; ++k) {
        vector x(dim);

        for (int i = 0; i < dim; ++i) {
            x[i] = randd(a[i], b[i]);
        }

        double fx = f(x);

        sum += fx;
        sumsq += fx * fx;
    }

    double mean = sum / static_cast<double>(N);
    double variance = sumsq / static_cast<double>(N) - mean * mean;
    variance = std::max(variance, 0.0);

    double value = mean * V;
    double error = std::sqrt(variance / static_cast<double>(N)) * V;

    return {value, error};
}

tuple lcgmc(Func f, const vector& a, const vector& b, int N, uint32_t seed)
{
    if (N <= 0) {
        return {0.0, 0.0};
    }

    int dim = static_cast<int>(a.size());
    double V = box_volume(a, b);

    lcg gen(seed);

    double sum = 0.0;
    double sumsq = 0.0;

    for (int k = 0; k < N; ++k) {
        vector x(dim);

        for (int i = 0; i < dim; ++i) {
            x[i] = gen.uniform(a[i], b[i]);
        }

        double fx = f(x);

        sum += fx;
        sumsq += fx * fx;
    }

    double mean = sum / static_cast<double>(N);
    double variance = sumsq / static_cast<double>(N) - mean * mean;
    variance = std::max(variance, 0.0);

    double value = mean * V;
    double error = std::sqrt(variance / static_cast<double>(N)) * V;

    return {value, error};
}

tuple quasimc(Func f, const vector& a, const vector& b, int N)
{
    if (N <= 1) {
        return {0.0, 0.0};
    }

    int dim = static_cast<int>(a.size());
    double V = box_volume(a, b);

    halton seq(dim);

    int half = N / 2;

    double sum1 = 0.0;
    double sum2 = 0.0;

    for (int k = 1; k <= half; ++k) {
        vector u = seq.point(k);
        vector x(dim);

        for (int i = 0; i < dim; ++i) {
            x[i] = a[i] + u[i] * (b[i] - a[i]);
        }

        sum1 += f(x);
    }

    for (int k = half + 1; k <= N; ++k) {
        vector u = seq.point(k + 1000);
        vector x(dim);

        for (int i = 0; i < dim; ++i) {
            x[i] = a[i] + u[i] * (b[i] - a[i]);
        }

        sum2 += f(x);
    }

    double I1 = V * sum1 / static_cast<double>(half);
    double I2 = V * sum2 / static_cast<double>(N - half);

    double value = 0.5 * (I1 + I2);
    double error = std::abs(I1 - I2);

    return {value, error};
}

tuple stratifiedmc(Func f, const vector& a, const vector& b, int N, int nmin)
{
    if (N <= 0) {
        return {0.0, 0.0};
    }

    if (N <= 2 * nmin) {
        return plainmc(f, a, b, N);
    }

    int dim = static_cast<int>(a.size());

    std::vector<vector> points(nmin, vector(dim));
    vector values(nmin);

    for (int k = 0; k < nmin; ++k) {
        for (int i = 0; i < dim; ++i) {
            points[k][i] = randd(a[i], b[i]);
        }

        values[k] = f(points[k]);
    }

    int split_dim = 0;

    double largest_subvariance = -1.0;
    double left_variance_best = 0.0;
    double right_variance_best = 0.0;

    for (int d = 0; d < dim; ++d) {
        double midpoint = 0.5 * (a[d] + b[d]);

        double left_sum = 0.0;
        double left_sumsq = 0.0;
        double right_sum = 0.0;
        double right_sumsq = 0.0;

        int left_count = 0;
        int right_count = 0;

        for (int k = 0; k < nmin; ++k) {
            double fx = values[k];

            if (points[k][d] < midpoint) {
                left_sum += fx;
                left_sumsq += fx * fx;
                ++left_count;
            } else {
                right_sum += fx;
                right_sumsq += fx * fx;
                ++right_count;
            }
        }

        double left_variance = 0.0;
        double right_variance = 0.0;

        if (left_count > 1) {
            double left_mean = left_sum / static_cast<double>(left_count);
            left_variance =
                left_sumsq / static_cast<double>(left_count)
                - left_mean * left_mean;
        }

        if (right_count > 1) {
            double right_mean = right_sum / static_cast<double>(right_count);
            right_variance =
                right_sumsq / static_cast<double>(right_count)
                - right_mean * right_mean;
        }

        left_variance = std::max(left_variance, 0.0);
        right_variance = std::max(right_variance, 0.0);

        double subvariance = left_variance + right_variance;

        if (subvariance > largest_subvariance) {
            largest_subvariance = subvariance;
            split_dim = d;
            left_variance_best = left_variance;
            right_variance_best = right_variance;
        }
    }

    vector a_left = a;
    vector b_left = b;
    vector a_right = a;
    vector b_right = b;

    double midpoint = 0.5 * (a[split_dim] + b[split_dim]);

    b_left[split_dim] = midpoint;
    a_right[split_dim] = midpoint;

    int remaining = N - nmin;

    if (remaining <= 1) {
        return plainmc(f, a, b, N);
    }

    int N_left = 0;

    if (left_variance_best + right_variance_best == 0.0) {
        N_left = remaining / 2;
    } else {
        N_left = static_cast<int>(
            std::round(
                remaining * left_variance_best
                / (left_variance_best + right_variance_best)
            )
        );
    }

    N_left = std::max(1, std::min(remaining - 1, N_left));
    int N_right = remaining - N_left;

    auto [I_left, E_left] = stratifiedmc(f, a_left, b_left, N_left, nmin);
    auto [I_right, E_right] = stratifiedmc(f, a_right, b_right, N_right, nmin);

    double value = I_left + I_right;
    double error = std::sqrt(E_left * E_left + E_right * E_right);

    return {value, error};
}

}