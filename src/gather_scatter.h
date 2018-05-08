/*
Full-integer grids [F] at i*h, j*h, k*h (i = 1..Nx, j = 1..Ny, k = 1..Nz)
and half-integer grids [H] at (i + 1/2)*h, (j + 1/2)*h, (k + 1/2)*h.
(i,j,k) is labeled [FFF].
(i+1/2,j+1/2,k+1/2) is labeled [HHH].

Charge density NP defined at [FFF].
Ex, Ey, Ez are defined at [HFF], [FHF], [FFH]. Face centered
Jx, Jy, Jz are defined at [HFF], [FHF], [FFH]. Face centered
Bx, By, Bz are defined at [FHH], [HFH], [HHF]. Edge centered
*/

#pragma once

#include "grid.h"
#include "particles.h"

#include <vector>
#include <string>

namespace PIC {

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
*    at_point - specified position                                           *
*    val_x, val_y, val_z - Cell required values                              *
*    ret_vec - result gathered value                                         *
* Used to interpolate Bx, By, Bz for specified position                      *
*****************************************************************************/
void gather_edge(const Grid& grid,
                 const DblVector& at_point,
                 CellVectorValue val,
                 DblVector& ret_vec);

/*****************************************************************************
* Gather face-centered values: Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz      *
*    at_point - specified position                                           *
*    val_x, val_y, val_z - Pointers to cell members required                 *
*    ret_vec - result gathered value                                         *
* Used to interpolate Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz values        *
* for specified position.                                                    *
*****************************************************************************/
void gather_face(const Grid& grid,
                 const DblVector& at_point,
                 CellVectorValue val,
                 DblVector& ret_vec);

/*****************************************************************************
* Gather cell-centered value.                                                *
*   at_point - specified position                                            *
*   val - Pointer to cell member required                                    *
*   ret_val - result gathered value                                          *
* Used to interpolate density NP values for specified position.              *
*                                                                            *
*****************************************************************************/
double gather_center(const Grid& grid,
                     const DblVector& at_point,
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
                     const DblVector& at_point,
                     DblVector GridT::NodeType::* vec,
                     double DblVector::* comp)
{
    const double h = grid.step();

    // home cell node indexes
    const index_t pi = static_cast<index_t>(at_point.x/h - Centering::x());
    const index_t pj = static_cast<index_t>(at_point.y/h - Centering::y());
    const index_t pk = static_cast<index_t>(at_point.z/h - Centering::z());

    double ret_val = 0.0;

    for (index_t i = pi; i != pi+2; ++i)
    for (index_t j = pj; j != pj+2; ++j)
    for (index_t k = pk; k != pk+2; ++k)
    {
        ret_val += (grid(i,j,k).*vec.*comp)*R((i + Centering::x())*h - at_point.x,
                                              (j + Centering::y())*h - at_point.y,
                                              (k + Centering::z())*h - at_point.z, h);
    }

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
                 const DblVector& at_point,
                 DblVector GridT::NodeType::* val,
                 DblVector& ret_vec)
{
    ret_vec.x = gather_vector<PIC::FaceXCentering>(grid, at_point, val, &DblVector::x);
    ret_vec.y = gather_vector<PIC::FaceYCentering>(grid, at_point, val, &DblVector::y);
    ret_vec.z = gather_vector<PIC::FaceZCentering>(grid, at_point, val, &DblVector::z);
}

template<class Centering, class GridT>
double gather_scalar(const GridT& grid,
                     const DblVector& at_point,
                     double GridT::NodeType::* val)
{
    const double h = grid.step();

    // home cell node indexes
    const index_t pi = static_cast<index_t>(at_point.x/h - Centering::x());
    const index_t pj = static_cast<index_t>(at_point.y/h - Centering::y());
    const index_t pk = static_cast<index_t>(at_point.z/h - Centering::z());

    double ret_val = 0.0;

    for (index_t i = pi; i != pi+2; ++i)
    for (index_t j = pj; j != pj+2; ++j)
    for (index_t k = pk; k != pk+2; ++k)
    {
        ret_val += (grid(i,j,k).*val)*R((i + Centering::x())*h - at_point.x,
                                        (j + Centering::y())*h - at_point.y,
                                        (k + Centering::z())*h - at_point.z, h);
    }

    return ret_val;
}

typedef GridContainer<Density> DensityGrid;

void scatter_particle_std(const Particle& particlep, DensityGrid& grid, const DblVector& dr);
void scatter_particle_zigzag(const Particle& particle, DensityGrid& grid, const DblVector& dr);

/***************************************************************************
* Set boundary conditions                                                  *
***************************************************************************/
void set_boundary_conditions(const std::string& lua_func_name);

void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg);

} // namespace PIC
