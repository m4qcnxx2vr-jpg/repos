#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdio>

#define SELF (*this)
#define FORV(i,v) for(int i=0;i<(v).size();i++)
#define FOR_COLS(j,A) for(int j=0;j<(A).size2();j++)

namespace pp {

// ---------------- approx ----------------
bool approx(double x, double y, double acc, double eps){
    if(std::fabs(x-y) < acc) return true;
    if(std::fabs(x-y) < eps*(std::fabs(x)+std::fabs(y))) return true;
    return false;
}

bool approx(const vector& u,const vector& v,double acc,double eps){
    if(u.size()!=v.size()) return false;
    for(int i=0;i<u.size();i++)
        if(!approx(u[i],v[i],acc,eps)) return false;
    return true;
}

// ---------------- vector ----------------
vector::vector(int n, double value) : data((std::size_t)n, value) {}

int vector::size() const { return (int)data.size(); }

void vector::resize(int n) { data.resize((std::size_t)n); }

double& vector::operator[](int i) { return data[(std::size_t)i]; }
const double& vector::operator[](int i) const { return data[(std::size_t)i]; }

vector& vector::operator+=(const vector& other) {
    if(size()!=other.size()) throw std::invalid_argument("vector += : size mismatch");
    FORV(i,SELF) data[(std::size_t)i] += other.data[(std::size_t)i];
    return SELF;
}

vector& vector::operator-=(const vector& other) {
    if(size()!=other.size()) throw std::invalid_argument("vector -= : size mismatch");
    FORV(i,SELF) data[(std::size_t)i] -= other.data[(std::size_t)i];
    return SELF;
}

vector& vector::operator*=(double x) {
    FORV(i,SELF) data[(std::size_t)i] *= x;
    return SELF;
}

vector& vector::operator/=(double x) {
    if(x==0.0) throw std::invalid_argument("vector /= 0");
    FORV(i,SELF) data[(std::size_t)i] /= x;
    return SELF;
}

vector& vector::add(double x){
    data.push_back(x);
    return SELF;
}

vector& vector::push_back(double x){
    data.push_back(x);
    return SELF;
}

double vector::norm() const {
    double s=0;
    FORV(i,SELF) s += SELF[i]*SELF[i];
    return std::sqrt(s);
}

double vector::dot(const vector& other) const {
    if(size() != other.size())
        throw std::invalid_argument("dot: size mismatch");
    double s = 0.0;
    for(int i = 0; i < size(); i++)
        s += (*this)[i] * other[i];
    return s;
}

vector vector::map(std::function<double(double)> f) const{
    vector r = SELF;
    for(int i=0;i<r.size();i++) r[i] = f(r[i]);
    return r;
}

void vector::print(std::string s) const {
    std::cout << s;
    FORV(i,SELF) std::printf("%9.3g ", (double)SELF[i]);
    std::printf("\n");
}

// vector free ops
vector operator/(const vector& v, double x){
    vector r=v; r/=x; return r;
}
vector operator*(const vector& v, double x){
    vector r=v; r*=x; return r;
}
vector operator*(double x,const vector& a){ return a*x; }
vector operator+(const vector& a, const vector& b){
    vector r=a; r+=b; return r;
}
vector operator-(const vector& a){
    vector r=a;
    for(int i=0;i<r.size();i++) r[i] = -r[i];
    return r;
}
vector operator-(const vector& a, const vector& b){
    vector r=a; r-=b; return r;
}



// ---------------- matrix ----------------
matrix::matrix(int n, int m, double value)
    : cols((std::size_t)m, pp::vector(n, value)) {}

int matrix::size1() const { return cols.empty() ? 0 : cols[0].size(); }
int matrix::size2() const { return (int)cols.size(); }

void matrix::resize(int n, int m){
    cols.resize((std::size_t)m);
    for(int j=0;j<m;++j) cols[(std::size_t)j].resize(n);
}

double& matrix::operator()(int i, int j){ return cols[(std::size_t)j][i]; }
const double& matrix::operator()(int i, int j) const { return cols[(std::size_t)j][i]; }

vector& matrix::operator[](int j){ return cols[(std::size_t)j]; }
const vector& matrix::operator[](int j) const { return cols[(std::size_t)j]; }

matrix& matrix::operator+=(const matrix& other) {
    if(size1()!=other.size1() || size2()!=other.size2())
        throw std::invalid_argument("matrix += : size mismatch");
    FOR_COLS(j,SELF) SELF[j] += other[j];
    return SELF;
}

matrix& matrix::operator-=(const matrix& other) {
    if(size1()!=other.size1() || size2()!=other.size2())
        throw std::invalid_argument("matrix -= : size mismatch");
    FOR_COLS(j,SELF) SELF[j] -= other[j];
    return SELF;
}

matrix& matrix::operator*=(double x) {
    FOR_COLS(j,SELF) SELF[j] *= x;
    return SELF;
}

matrix& matrix::operator/=(double x) {
    if(x==0.0) throw std::invalid_argument("matrix /= 0");
    FOR_COLS(j,SELF) SELF[j] /= x;
    return SELF;
}

matrix operator/(const matrix& A,double x){
    matrix R=A; R/=x; return R;
}
matrix operator*(const matrix& A,double x){
    matrix R=A; R*=x; return R;
}
matrix operator*(double x,const matrix& A){
    return A*x;
}
matrix operator+(const matrix& A, const matrix& B){
    matrix R=A; R+=B; return R;
}
matrix operator-(const matrix& A, const matrix& B){
    matrix R=A; R-=B; return R;
}

// matrix-vector (efficient for column storage)
vector operator*(const matrix& M, const vector& v){
    if(M.size2()!=v.size())
        throw std::invalid_argument("M*v: size mismatch");

    vector r(M.size1(), 0.0);

    for(int j=0;j<M.size2();j++){
        double vj = v[j];
        for(int i=0;i<M.size1();i++)
            r[i] += M(i,j) * vj;
    }
    return r;
}

matrix operator*(const matrix& A, const matrix& B){
    if(A.size2()!=B.size1())
        throw std::invalid_argument("A*B: size mismatch");

    matrix R(A.size1(), B.size2(), 0.0);

    for(int j=0;j<B.size2();j++){
        for(int k=0;k<A.size2();k++){
            double Bkj = B(k,j);
            for(int i=0;i<A.size1();i++)
                R(i,j) += A(i,k) * Bkj;
        }
    }
    return R;
}

void matrix::setid(){
    assert(size1()==size2());
    for(int i=0;i<size1();i++){
        SELF(i,i) = 1;
        for(int j=i+1;j<size1();j++) SELF(i,j) = SELF(j,i) = 0;
    }
}

matrix matrix::transpose() const {
    matrix R(size2(), size1(), 0.0);
    for(int i=0;i<size1();i++)
        for(int j=0;j<size2();j++)
            R(j,i) = SELF(i,j);
    return R;
}

matrix matrix::T() const { return SELF.transpose(); }

void matrix::print(std::string s) const {
    std::cout << s << std::endl;
    for(int i=0;i<size1();i++){
        for(int j=0;j<size2();j++)
            std::printf("%9.3g ", (double)SELF(i,j));
        std::printf("\n");
    }
    std::printf("\n");
}

} // namespace pp