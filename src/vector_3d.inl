#include <cmath>

template<typename T>
inline Vector3D<T>::Vector3D()
: x(0), y(0), z(0)
{}

template<typename T>
template<typename R>
inline Vector3D<T>::Vector3D(const Vector3D<R>& other)
: x(other.x), y(other.y), z(other.z)
{}

template<typename T>
inline Vector3D<T>::Vector3D(const T& x_, const T& y_, const T& z_)
: x(x_), y(y_), z(z_)
{}

template<typename T>
inline Vector3D<T>::Vector3D(const T coords[3])
: x(coords[0]), y(coords[1]), z(coords[2])
{}

template<typename T>
inline const T& Vector3D<T>::get_x() const
{
    return x;
}

template<typename T>
inline void Vector3D<T>::set_x(const T& new_x)
{
    x = new_x;
}

template<typename T>
inline const T& Vector3D<T>::get_y() const
{
    return y;
}

template<typename T>
inline void Vector3D<T>::set_y(const T& new_y)
{
    y = new_y;
}

template<typename T>
inline const T& Vector3D<T>::get_z() const
{
    return z;
}

template<typename T>
inline void Vector3D<T>::set_z(const T& new_z)
{
    z = new_z;
}

template<typename T>
inline void Vector3D<T>::getv(T buffer[3]) const
{
    buffer[0] = x;
    buffer[1] = y;
    buffer[2] = z;
}

template<typename T>
inline void Vector3D<T>::setv(const T coords[3])
{
    x = coords[0];
    y = coords[1];
    z = coords[2];
}

template<typename T>
inline void Vector3D<T>::get(T& x, T& y, T& z) const
{
    x = this->x;
    y = this->y;
    z = this->z;
}

template<typename T>
inline void Vector3D<T>::set(const T& x, const T& y, const T& z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

template<typename T>
inline const T& Vector3D<T>::operator[] (size_t index) const
{
    switch (index)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    }

    return T();
}

template<typename T>
inline T& Vector3D<T>::operator[] (size_t index)
{
    switch (index)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    }

    return T();
}

template<typename T>
inline const T& Vector3D<T>::operator() (size_t column) const
{
    switch (column)
    {
    case 1:
        return x;
    case 2:
        return y;
    case 3:
        return z;
    }

    return T();
}

template<typename T>
inline T& Vector3D<T>::operator() (size_t column)
{
    switch (column)
    {
    case 1:
        return x;
    case 2:
        return y;
    case 3:
        return z;
    }

    return T();
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator+ () const
{
    return *this;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator+ (const Vector3D& other) const
{
    return Vector3D(x + other.x, y + other.y, z + other.z);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator+= (const Vector3D& other)
{
    return *this = *this + other;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator- () const
{
    return Vector3D(-x, -y, -z);
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator- (const Vector3D& other) const
{
    return Vector3D(x - other.x, y - other.y, z - other.z);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator-= (const Vector3D& other)
{
    return *this = *this - other;
}

template<typename U>
const Vector3D<U> operator*(const U& scalar, const Vector3D<U>& v)
{
    return Vector3D<U>(v.x*scalar, v.y*scalar, v.z*scalar);
}

template<typename U>
const Vector3D<U> operator*(const Vector3D<U>& v, const U& scalar)
{
    return Vector3D<U>(v.x*scalar, v.y*scalar, v.z*scalar);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator*= (const T& scalar)
{
    return *this = *this * scalar;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator/ (const T& scalar) const
{
    return Vector3D(x/scalar, y/scalar, z/scalar);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator/= (const T& scalar)
{
    return *this = *this / scalar;
}

template<typename T>
inline const T Vector3D<T>::operator* (const Vector3D& other) const
{
    return x*other.x + y*other.y + z*other.z;
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::operator% (const Vector3D& other) const
{
    return Vector3D(y*other.z - z*other.y,
                    z*other.x - x*other.z,
                    x*other.y - y*other.x);
}

template<typename T>
inline Vector3D<T>& Vector3D<T>::operator%= (const Vector3D& other)
{
    return *this = *this % other;
}

template<typename T>
inline const T Vector3D<T>::get_sqr_len() const
{
    return x*x + y*y + z*z;
}

template<typename T>
inline const T Vector3D<T>::get_len() const
{
    return std::sqrt(get_sqr_len());
}

template<typename T>
inline const T Vector3D<T>::abs() const
{
    return get_len();
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::get_unit() const
{
    if (get_sqr_len() != 0)
        return *this / get_len();

    return *this;
}
