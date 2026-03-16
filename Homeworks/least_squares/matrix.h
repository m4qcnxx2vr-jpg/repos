#pragma once
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

namespace pp {

// -------- vector (dynamic column vector) --------
struct vector {
    std::vector<double> data;

    vector() = default;
    explicit vector(int n, double value = 0.0);

    int size() const;
    void resize(int n);

    double& operator[](int i);
    const double& operator[](int i) const;

    vector& operator+=(const vector& other);
    vector& operator-=(const vector& other);
    vector& operator*=(double x);
    vector& operator/=(double x);

    vector& add(double x);
    vector& push_back(double x);

    double norm() const;
    double dot(const vector& other) const;
    vector map(std::function<double(double)> f) const;

    void print(std::string s = "") const;
};

// vector free operators
vector operator/(const vector& v, double x);
vector operator*(const vector& v, double x);
vector operator*(double x, const vector& a);
vector operator+(const vector& a, const vector& b);
vector operator-(const vector& a);
vector operator-(const vector& a, const vector& b);

// approx helpers
bool approx(double x, double y, double acc = 1e-6, double eps = 1e-6);
bool approx(const vector& u, const vector& v, double acc = 1e-6, double eps = 1e-6);

// -------- matrix (stored by columns) --------
struct matrix {
    std::vector<pp::vector> cols;

    matrix() = default;
    matrix(int n, int m, double value = 0.0);

    int size1() const; // rows
    int size2() const; // cols

    void resize(int n, int m);

    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;

    pp::vector& operator[](int j);
    const pp::vector& operator[](int j) const;

    matrix& operator+=(const matrix& other);
    matrix& operator-=(const matrix& other);
    matrix& operator*=(double x);
    matrix& operator/=(double x);

    void setid();
    matrix transpose() const;
    matrix T() const;

    void print(std::string s = "") const;
};

// matrix free operators
matrix operator/(const matrix& A, double x);
matrix operator*(const matrix& A, double x);
matrix operator*(double x, const matrix& A);
matrix operator+(const matrix& A, const matrix& B);
matrix operator-(const matrix& A, const matrix& B);

vector operator*(const matrix& M, const vector& v);
matrix operator*(const matrix& A, const matrix& B);

} // namespace pp