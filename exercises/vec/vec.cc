#include "vec.h"

#include <cmath>

vec::vec()
    : x(0.0), y(0.0), z(0.0)
{
}

vec::vec(double x, double y, double z)
    : x(x), y(y), z(z)
{
}

vec& vec::operator+=(const vec& other)
{
    x += other.x;
    y += other.y;
    z += other.z;

    return *this;
}

vec& vec::operator-=(const vec& other)
{
    x -= other.x;
    y -= other.y;
    z -= other.z;

    return *this;
}

vec& vec::operator*=(double c)
{
    x *= c;
    y *= c;
    z *= c;

    return *this;
}

vec& vec::operator/=(double c)
{
    x /= c;
    y /= c;
    z /= c;

    return *this;
}

double vec::dot(const vec& other) const
{
    return x * other.x + y * other.y + z * other.z;
}

vec vec::cross(const vec& other) const
{
    return vec(
        y * other.z - z * other.y,
        z * other.x - x * other.z,
        x * other.y - y * other.x
    );
}

double vec::norm() const
{
    return std::sqrt(this->dot(*this));
}

void vec::print(std::string s) const
{
    std::cout << s << x << " " << y << " " << z << std::endl;
}

vec operator+(vec a, const vec& b)
{
    a += b;
    return a;
}

vec operator-(vec a, const vec& b)
{
    a -= b;
    return a;
}

vec operator*(vec v, double c)
{
    v *= c;
    return v;
}

vec operator*(double c, vec v)
{
    v *= c;
    return v;
}

vec operator/(vec v, double c)
{
    v /= c;
    return v;
}

bool approx(double a, double b, double acc, double eps)
{
    if (std::abs(a - b) < acc) {
        return true;
    }

    if (std::abs(a - b) < eps * (std::abs(a) + std::abs(b))) {
        return true;
    }

    return false;
}

bool approx(const vec& a, const vec& b, double acc, double eps)
{
    if (!approx(a.x, b.x, acc, eps)) return false;
    if (!approx(a.y, b.y, acc, eps)) return false;
    if (!approx(a.z, b.z, acc, eps)) return false;

    return true;
}

std::ostream& operator<<(std::ostream& os, const vec& v)
{
    os << "{ " << v.x << ", " << v.y << ", " << v.z << " }";
    return os;
}