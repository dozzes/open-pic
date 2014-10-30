#if !defined (VECTOR_3D_H)
#define VECTOR_3D_H

#include <vector>
#include <ostream>


template <class T>
struct Vector3D
{
   Vector3D() : x(0.0), y(0.0), z(0.0) {}

   explicit
   Vector3D(T _x, T _y, T _z)
    : x(_x), y(_y), z(_z) {}

   Vector3D(const Vector3D<T>& other)
   : x(other.x), y(other.y), z(other.z) {}

   Vector3D& operator=(const Vector3D& rhs)
   {
      if ((void*)this == (void*)&rhs)
      {
         return *this;
      }

      x = rhs.x;
      y = rhs.y;
      z = rhs.z;

      return *this;
   }

   double abs2() const { return (x*x + y*y + z*z); }
   double abs() const { return sqrt(abs2()); }
   T x, y, z;
};

template<class T>
std::ostream& operator<<(std::ostream& out, const Vector3D<T>& v)
{
   return (out << v.x  << "\t" << v.y  << "\t" << v.z);
}

typedef Vector3D<double> DblVector;

#endif // VECTOR_3D_H
