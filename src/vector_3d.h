/*
* Templated 3D Vector Class
* http://cpp-wiki.wikidot.com/code:templated-3d-vector-class
*/

#pragma once

#include <ostream>
#include <cstddef>

template<typename T>
class Vector3D
{
public:
    T x;
    T y;
    T z;

public:
    Vector3D();

    template<typename R>
    explicit Vector3D(const Vector3D<R>& other);
    explicit Vector3D(const T& x, const T& y, const T& z);
    explicit Vector3D(const T coords[3]);

    const T& get_x() const;
    void set_x(const T& newX);

    const T& get_y() const;
    void set_y(const T& newY);

    const T& get_z() const;
    void set_z(const T& newZ);

    void getv(T buffer[3]) const;
    void setv(const T coords[3]);

    void get(T& x, T& y, T& z) const;
    void set(const T& x, const T& y, const T& z);

    // Interface for indexing

    const T& operator[] (size_t index) const;
    T& operator[] (size_t index);

    // Considering vectors as matrices with one row
    const T& operator() (size_t column) const;
    T& operator() (size_t column);

    // Standard operations

    // This does absolutely nothing, but it should be included for consistency
    const Vector3D operator+ () const;

    const Vector3D operator+ (const Vector3D& other) const;
    Vector3D& operator+= (const Vector3D& other);

    // The same as multiplying *this by -1
    const Vector3D operator- () const;

    const Vector3D operator- (const Vector3D& other) const;
    Vector3D& operator-= (const Vector3D& other);

    // Multiplying Vector3D by a scalar
    template<class U> friend const Vector3D<U> operator*(const U& scalar, const Vector3D<U>& v);
    template<class U> friend const Vector3D<U> operator*(const Vector3D<U>& v, const U& scalar);

    // Multiplying *this by a scalar
    Vector3D& operator*= (const T& scalar);

    // Same as multiplication by 1/scalar, maybe more accurate but also slower
    const Vector3D operator/ (const T& scalar) const;
    Vector3D& operator/= (const T& scalar);

    // Calculate the dot/inner/scalar product
    const T operator* (const Vector3D& other) const;

    // Calculate the cross/outer/vector product
    const Vector3D operator% (const Vector3D& other) const;
    Vector3D& operator%= (const Vector3D& other);

    // Auxiliary methods

    // Returns the squared length of *this
    const T get_sqr_len() const;

    // Returns the length of *this
    const T get_len() const;

    // Returns the length of *this
    const T abs() const;

    // Returns a vector with the same orientation, but with a length of 1
    const Vector3D get_unit() const;
};

template<class T>
std::ostream& operator<<(std::ostream& out, const Vector3D<T>& v)
{
    return (out << v.x << "\t" << v.y << "\t" << v.z);
}

typedef Vector3D<double> DblVector;

#include "vector_3d.inl"
