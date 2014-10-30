/*
Full-integer grids [F] at i*h, j*h, k*h (i = 1, 2, ...,Nx, j = 1, 2, ...,Ny, k = 1, 2, ...,Nz)
and half-integer grids [H] at (i + 1/2)*h, (j + 1/2)*h, (k + 1/2)*h.
(i,j,k) is labeled [FFF].
(i+1/2,j+1/2,k+1/2) is labeled [HHH].

Charge density NP defined at [FFF].
Ex, Ey, Ez are defined at [HFF], [FHF], [FFH]. Face centered
Jx, Jy, Jz are defined at [HFF], [FHF], [FFH]. Face centered
Bx, By, Bz are defined at [FHH], [HFH], [HHF]. Edge centered
*/

#if !defined (GATHER_SCATTER_H)
#define GATHER_SCATTER_H

#include <vector>
#include <string>

#include "grid.h"
#include "particles.h"


namespace PIC
{

/***************************************************************************
* 1D particle form-factor                                                  *
***************************************************************************/
double R(double x, double h);

/***************************************************************************
* 3D particle form-factor                                                  *
***************************************************************************/
double R(double x, double y, double z, double h);

/***************************************************************************
* Pointer to Cell value                                                    *
***************************************************************************/
typedef DblVector (Cell::*CellVectorValue);
typedef double (Cell::*CellScalarValue);
typedef double (DblVector::*VectorComp);

/*****************************************************************************
* Gather edge-centered values.                                               *
*    x, y, z - specified position                                            *
*    val_x, val_y, val_z - Cell required values                              *
*    ret_vec - result gathered value                                         *
* Used to interpolate Bx, By, Bz for specified position                      *
*****************************************************************************/
void gather_edge(const Grid& grid,
                 const DblVector& atPoint,
                 CellVectorValue val,
                 DblVector& ret_vec);

/*****************************************************************************
* Gather face-centered values: Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz      *
*    x, y, z - specified position                                            *
*    val_x, val_y, val_z - Pointers to cell members required                 *
*    ret_vec - result gathered value                                         *
* Used to interpolate Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz values        *
* for specified position.                                                    *
*****************************************************************************/
void gather_face(const Grid& grid,
                 const DblVector& atPoint,
                 CellVectorValue val,
                 DblVector& ret_vec);

/*****************************************************************************
* Gather cell-centered value.                                                *
*   x, y, z - specified position                                             *
*   val - Pointer to cell member required                                    *
*   ret_val - result gathered value                                          *
* Used to interpolate density NP values for specified position.              *
*                                                                            *
*****************************************************************************/
double gather_center(const Grid& grid,
                     const DblVector& atPoint,
                     CellScalarValue val);

/***************************************************************************
* interpolates grid values to specified point <vec_to_point>               *
***************************************************************************/
void from_grid_to_point(const Grid& grid,
                        const DblVector& vec_to_point,
                        Grid::NodeType& val_at_point);

/***************************************************************************
* interpolates grid values to specified point (x,y,z)                      *
***************************************************************************/
void from_grid_to_point(const Grid& grid,
                        double x, double y, double z,
                        Grid::NodeType& val_at_point);

struct FaceXCentering
{
   static double x() { return 0.0; }
   static double y() { return 0.5; }
   static double z() { return 0.5; }
};

struct FaceYCentering
{
   static double x() { return 0.5; }
   static double y() { return 0.0; }
   static double z() { return 0.5; }
};

struct FaceZCentering
{
   static double x() { return 0.5; }
   static double y() { return 0.5; }
   static double z() { return 0.0; }
};

struct EdgeXCentering
{
   static double x() { return 0.5; }
   static double y() { return 0.0; }
   static double z() { return 0.0; }
};

struct EdgeYCentering
{
   static double x() { return 0.0; }
   static double y() { return 0.5; }
   static double z() { return 0.0; }
};

struct EdgeZCentering
{
   static double x() { return 0.0; }
   static double y() { return 0.0; }
   static double z() { return 0.5; }
};

struct CellCentering
{
   static double x() { return 0.5; }
   static double y() { return 0.5; }
   static double z() { return 0.5; }
};

struct NodeCentering
{
   static double x() { return 0.0; }
   static double y() { return 0.0; }
   static double z() { return 0.0; }
};

struct Density
{
   Density() : UP(0.0, 0.0, 0.0), NP(0.0) { }
   DblVector UP; // ion velocity
   double NP; // ion charge density
};

template<class Centering, class GridT>
double gather_vector(const GridT& grid,
                       const DblVector& atPoint,
                       DblVector GridT::NodeType::* vec,
                       double DblVector::* comp)
{
    const double h = grid.step();

    // home cell node indexes
    const index_t ip = static_cast<index_t>(atPoint.x/h - Centering::x());
    const index_t jp = static_cast<index_t>(atPoint.y/h - Centering::y());
    const index_t kp = static_cast<index_t>(atPoint.z/h - Centering::z());

    // point shift from home cell node
    const double dx = atPoint.x/h - (ip + Centering::x());
    const double dy = atPoint.y/h - (jp + Centering::y());
    const double dz = atPoint.z/h - (kp + Centering::z());

    // field interpolation weights
    const double w_000 = (1.0 - dx)*(1.0 - dy)*(1.0 - dz);
    const double w_001 = (1.0 - dx)*(1.0 - dy)*(dz);
    const double w_010 = (1.0 - dx)*(dy)*(1.0 - dz);
    const double w_011 = (1.0 - dx)*(dy)*(dz);
    const double w_100 = (dx)*(1.0 - dy)*(1.0 - dz);
    const double w_101 = (dx)*(1.0 - dy)*(dz);
    const double w_110 = (dx)*(dy)*(1.0 - dz);
    const double w_111 = (dx)*(dy)*(dz);

    const double v_000 = grid(ip,jp,kp).*vec.*comp;
    const double v_001 = grid(ip,jp,kp+1).*vec.*comp;
    const double v_010 = grid(ip,jp+1,kp).*vec.*comp;
    const double v_011 = grid(ip,jp+1,kp+1).*vec.*comp;
    const double v_100 = grid(ip+1,jp,kp).*vec.*comp;
    const double v_101 = grid(ip+1,jp,kp+1).*vec.*comp;
    const double v_110 = grid(ip+1,jp+1,kp).*vec.*comp;
    const double v_111 = grid(ip+1,jp+1,kp+1).*vec.*comp;

    const double ret_val = w_000*v_000 + w_001*v_001 + w_010*v_010 +
                           w_011*v_011 + w_100*v_100 + w_101*v_101 +
                           w_110*v_110 + w_111*v_111;

    return ret_val;
}

template<class Centering, class GridT>
double gather_vector(const GridT& grid,
                     double x, double y, double z,
                     DblVector GridT::NodeType::* vec,
                     double DblVector::* comp)
{
    return gather_vector<Centering>(grid, DblVector(x, y, z), vec, comp);
}

template<class GridT>
void gather_face(const GridT& grid,
                   const DblVector& atPoint,
                   DblVector GridT::NodeType::* val,
                   DblVector& ret_vec)
{
    ret_vec.x = gather_vector<PIC::FaceXCentering>(grid, atPoint, val, &DblVector::x);
    ret_vec.y = gather_vector<PIC::FaceYCentering>(grid, atPoint, val, &DblVector::y);
    ret_vec.z = gather_vector<PIC::FaceZCentering>(grid, atPoint, val, &DblVector::z);
}

template<class Centering, class GridT>
double gather_scalar(const GridT& grid,
                     const DblVector& atPoint,
                     double GridT::NodeType::* val)
{
    const double h = grid.step();

    // home cell node indexes
    const index_t i = static_cast<index_t>(atPoint.x/h - Centering::x());
    const index_t j = static_cast<index_t>(atPoint.y/h - Centering::y());
    const index_t k = static_cast<index_t>(atPoint.z/h - Centering::z());

    // point shift from home cell node
    const double dx = atPoint.x/h - (i + Centering::x());
    const double dy = atPoint.y/h - (j + Centering::y());
    const double dz = atPoint.z/h - (k + Centering::z());

    // field interpolation weights
    const double w_000 = (1.0 - dx)*(1.0 - dy)*(1.0 - dz);
    const double w_001 = (1.0 - dx)*(1.0 - dy)*(dz);
    const double w_010 = (1.0 - dx)*(dy)*(1.0 - dz);
    const double w_011 = (1.0 - dx)*(dy)*(dz);
    const double w_100 = (dx)*(1.0 - dy)*(1.0 - dz);
    const double w_101 = (dx)*(1.0 - dy)*(dz);
    const double w_110 = (dx)*(dy)*(1.0 - dz);
    const double w_111 = (dx)*(dy)*(dz);

    const double v_000 = grid(i,j,k).*val;
    const double v_001 = grid(i,j,k+1).*val;
    const double v_010 = grid(i,j+1,k).*val;
    const double v_011 = grid(i,j+1,k+1).*val;
    const double v_100 = grid(i+1,j,k).*val;
    const double v_101 = grid(i+1,j,k+1).*val;
    const double v_110 = grid(i+1,j+1,k).*val;
    const double v_111 = grid(i+1,j+1,k+1).*val;

    const double ret_val = w_000*v_000 + w_001*v_001 + w_010*v_010 +
                           w_011*v_011 + w_100*v_100 + w_101*v_101 +
                           w_110*v_110 + w_111*v_111;
    return ret_val;
}

typedef GridContainer<Density> DensityGrid;

void scatter_particle_std(const Particle& p, DensityGrid& dg);

/***************************************************************************
* Set boundary conditions                                                  *
***************************************************************************/
void set_boundary_conditions(const std::string& lua_func_name);

void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg);

} // namespace PIC

#endif // GATHER_SCATTER_H
