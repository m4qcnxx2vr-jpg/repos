#include "Harmonic.h"

#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include <iomanip>
#include <functional>

int main(int argc, char* argv[]) {
    int nthreads = 1;
    int nterms = static_cast<int>(1e8);

    for (int i = 0; i < argc; i++) {
        std::string arg = argv[i];

        if ((arg == "-threads" || arg == "-nthreads") && i + 1 < argc) {
            nthreads = std::stoi(argv[i + 1]);
        }

        if ((arg == "-terms" || arg == "-nterms") && i + 1 < argc) {
            nterms = std::stoi(argv[i + 1]);
        }
    }

    if (nthreads < 1) {
        std::cerr << "Error: number of threads must be at least 1.\n";
        return 1;
    }

    if (nterms < 1) {
        std::cerr << "Error: number of terms must be at least 1.\n";
        return 1;
    }

    std::vector<Data> params(nthreads);

    for (int i = 0; i < nthreads; i++) {
        params[i].a = 1 + nterms / nthreads * i;
        params[i].b = 1 + nterms / nthreads * (i + 1);
        params[i].sum = 0.0;
    }

    params[params.size() - 1].b = nterms + 1;

    std::vector<std::thread> threads;
    threads.reserve(nthreads);

    for (int i = 0; i < nthreads; i++) {
        threads.emplace_back(harm, std::ref(params[i]));
    }

    for (auto& thread : threads) {
        thread.join();
    }

    double total = 0.0;

    for (const auto& p : params) {
        total += p.sum;
    }

    double expected = harmonic_asymptotic(nterms);
    double tolerance = 1e-8;

    std::cout << std::setprecision(16);
    std::cout << "Number of threads: " << nthreads << "\n";
    std::cout << "Number of terms:   " << nterms << "\n";
    std::cout << "Harmonic sum:      " << total << "\n";
    std::cout << "Asymptotic check:  " << expected << "\n";
    std::cout << "Difference:        " << std::abs(total - expected) << "\n";
    std::cout << "Approx OK:         "
              << std::boolalpha
              << approx(total, expected, tolerance)
              << "\n";

    return 0;
}