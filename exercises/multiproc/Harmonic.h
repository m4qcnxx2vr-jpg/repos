#pragma once

struct Data {
    int a;
    int b;
    double sum;
};

void harm(Data& arg);

double harmonic_asymptotic(int n);

bool approx(double x, double y, double tolerance);