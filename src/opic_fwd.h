#if !defined (OPIC_FWD_H)
#define OPIC_FWD_H


struct Cell;

template<class NodeT> class GridContainer;

typedef GridContainer<Cell> Grid;

template <class ScalarT> struct Vector3D;

typedef Vector3D<double> DblVector;

typedef double (DblVector::*VectorComp);

struct Particle;

class Particles;

typedef unsigned long index_t;

#endif // OPIC_FWD_H
