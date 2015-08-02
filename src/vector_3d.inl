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
inline const T& Vector3D<T>::getX() const
{
    return x;
}
 
template<typename T>
inline void Vector3D<T>::setX(const T& newX)
{
    x = newX;
}
 
template<typename T>
inline const T& Vector3D<T>::getY() const
{
    return y;
}
 
template<typename T>
inline void Vector3D<T>::setY(const T& newY)
{
    y = newY;
}
 
template<typename T>
inline const T& Vector3D<T>::getZ() const
{
    return z;
}
 
template<typename T>
inline void Vector3D<T>::setZ(const T& newZ)
{
    z = newZ;
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
inline const T Vector3D<T>::getSqrLen() const
{
    return x*x + y*y + z*z;
}
 
template<typename T>
inline const T Vector3D<T>::getLen() const
{
    return std::sqrt(getSqrLen());
}

template<typename T>
inline const T Vector3D<T>::abs() const
{
    return getLen();
}

template<typename T>
inline const Vector3D<T> Vector3D<T>::getUnit() const
{
    if (getSqrLen() != 0)
        return *this / getLen();
 
    return *this;
}
 
/*template<typename T>
inline const Vector3D<T> Vector3D<T>::getInterpolation(const Vector3D& other, const T& ratio) const
{
    return *this + (other - *this)*ratio;
}
 
template<typename T>
inline const Vector3D<T> Vector3D<T>::getReflection(const Vector3D& surfaceNormal) const
{
    return *this - surfaceNormal*((*this * surfaceNormal)*2);
}
 
template<typename T>
inline const Vector3D<T> Vector3D<T>::getRotationEuler(const Vector3D& angles, const Vector3D& origin, bool degs, bool ccw) const
{
    T sin_phi, sin_theta, sin_psi;
    T cos_phi, cos_theta, cos_psi;
 
    if (degs)
    {
        const T degToRad = T(M_PI)/360;
 
        sin_phi = std::sin(angles.getX()*degToRad);
        sin_theta = std::sin(angles.getY()*degToRad);
        sin_psi = std::sin(angles.getZ()*degToRad);
 
        cos_phi = std::cos(angles.getX()*degToRad);
        cos_theta = std::cos(angles.getY()*degToRad);
        cos_psi = std::cos(angles.getZ()*degToRad);
    }
    else
    {
        sin_phi = std::sin(angles.getX());
        sin_theta = std::sin(angles.getY());
        sin_psi = std::sin(angles.getZ());
 
        cos_phi = std::cos(angles.getX());
        cos_theta = std::cos(angles.getY());
        cos_psi = std::cos(angles.getZ());
    }
 
    if (!ccw)
    {
        sin_phi = -sin_phi;
        sin_theta = -sin_theta;
        sin_psi = -sin_psi;
    }
 
    Vector3D temp = *this - origin;
    Vector3D result = temp;
 
    result.setY( temp.getY()*cos_phi + temp.getZ()*sin_phi);
    result.setZ(-temp.getY()*sin_phi + temp.getZ()*cos_phi);
 
    temp = result;
 
    result.setX(temp.getX()*cos_theta - temp.getZ()*sin_theta);
    result.setZ(temp.getX()*sin_theta + temp.getZ()*cos_theta);
 
    temp = result;
 
    result.setX( temp.getX()*cos_psi + temp.getY()*sin_psi);
    result.setY(-temp.getX()*sin_psi + temp.getY()*cos_psi);
 
    return result + origin;
}
 
template<typename T>
inline const Vector3D<T> Vector3D<T>::getRotationArbAxis(const Vector3D& axis, const T& amount, const Vector3D& origin, bool degs, bool ccw) const
{
    T cos_theta, sin_theta;
 
    if (degs)
    {
        const T degToRad = T(M_PI)/360;
 
        cos_theta = std::cos(amount*degToRad);
        sin_theta = std::sin(amount*degToRad);
    }
    else
    {
        cos_theta = std::cos(amount);
        sin_theta = std::sin(amount);
    }
 
    if (!ccw)
        sin_theta = -sin_theta;
 
    return *this*cos_theta + (axis % *this)*sin_theta + axis*((axis * *this)*(1 - cos_theta));
}
*/
