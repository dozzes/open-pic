#if !defined (UTILITY_H)
#define UTILITY_H

#include <ctime>
#include <cstdlib>

#include "vector_3d.h"
#include "constants.h"


inline
double rnd(double r_min = 0.0, double r_max = 1.0)
{
   return (( (double)rand() / (double) RAND_MAX) * r_max + r_min);
}

typedef Vector3D<double> DblVector;

void rnd_ball(double R, DblVector& result)
{
   double r = R * pow(rnd(), 1.0/3.0);
   double t = acos(2.0 * rnd() - 1);
   double f = 2.0 * PIC::Constants::pi() * rnd();

   result.x = r * sin(t) * cos(f);
   result.y = r * sin(t) * sin(f);
   result.z = r * cos(t);
}

#endif // UTILITY_H
