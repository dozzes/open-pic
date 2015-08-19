#if !defined (VECTOR_3D_H)
#define VECTOR_3D_H

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
    //! Sets all members to zero
    Vector3D();

    //! Explicitly converts from one type to another
    template<typename R>
    explicit Vector3D(const Vector3D<R>& other);

    explicit Vector3D(const T& x, const T& y, const T& z);

    explicit Vector3D(const T coords[3]);

    // Accessors

    const T& x() const;
    void set_x(const T& newX);

    const T& y() const;
    void set_y(const T& newY);

    const T& z() const;
    void set_z(const T& newZ);

    // Interface for indexing

    const T& operator[] (size_t index) const;
    T& operator[] (size_t index);

    //! Considering vectors as matrices with one row
    const T& operator() (size_t column) const;
    T& operator() (size_t column);

    // Standard operations

    //! This does absolutely nothing, but it should be included for consistency
    const Vector3D operator+ () const;

    const Vector3D operator+ (const Vector3D& other) const;
    Vector3D& operator+= (const Vector3D& other);

    //! The same as multiplying *this by -1
    const Vector3D operator- () const;

    const Vector3D operator- (const Vector3D& other) const;
    Vector3D& operator-= (const Vector3D& other);

    //! Multiplying Vector3D by a scalar
    template<class U> friend const Vector3D<U> operator*(const U& scalar, const Vector3D<U>& v);
    template<class U> friend const Vector3D<U> operator*(const Vector3D<U>& v, const U& scalar);

    //! Multiplying *this by a scalar
    Vector3D& operator*= (const T& scalar);

    //! Same as multiplication by 1/scalar, maybe more accurate but also slower
    const Vector3D operator/ (const T& scalar) const;
    Vector3D& operator/= (const T& scalar);

    //! Calculate the dot/inner/scalar product
    const T operator* (const Vector3D& other) const;

    //! Calculate the cross/outer/vector product
    const Vector3D operator% (const Vector3D& other) const;
    Vector3D& operator%= (const Vector3D& other);

    // Auxiliary methods

    //! Returns the squared length of *this
    const T sqr_length() const;

    //! Returns the length of *this
    const T length() const;

    //! Returns a vector with the same orientation, but with a length of 1
    const Vector3D unit() const;
};

template<class T>
std::ostream& operator<<(std::ostream& out, const Vector3D<T>& v)
{
    return (out << v.x << "\t" << v.y << "\t" << v.z);
}

typedef Vector3D<double> DblVector;

#include "vector_3d.inl"

#endif // VECTOR_3D_H
