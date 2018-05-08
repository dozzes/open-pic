#pragma once

template<class NodeT> class GridContainer;

struct Cell;
typedef GridContainer<Cell> Grid;

template <class ScalarT> class Vector3D;

typedef Vector3D<double> DblVector;

typedef double (DblVector::*VectorComp);

struct Particle;

class Particles;

typedef unsigned long index_t;
