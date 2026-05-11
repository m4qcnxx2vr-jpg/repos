#pragma once

#include <iostream>
#include <string>

class vec {
public:
    double x;
    double y;
    double z;

    vec();
    vec(double x, double y, double z);

    vec(const vec& other) = default;
    vec(vec&& other) = default;

    ~vec() = default;

    vec& operator=(const vec& other) = default;
    vec& operator=(vec&& other) = default;

    vec& operator+=(const vec& other);
    vec& operator-=(const vec& other);
    vec& operator*=(double c);
    vec& operator/=(double c);

    double dot(const vec& other) const;
    vec cross(const vec& other) const;
    double norm() const;

    void print(std::string s) const;
};

vec operator+(vec a, const vec& b);
vec operator-(vec a, const vec& b);

vec operator*(vec v, double c);
vec operator*(double c, vec v);
vec operator/(vec v, double c);

bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9);
bool approx(const vec& a, const vec& b, double acc = 1e-9, double eps = 1e-9);

std::ostream& operator<<(std::ostream& os, const vec& v);