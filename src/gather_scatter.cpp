#include "gather_scatter.h"
#include "config.h"
#include "use_lua.h"
#include "io_utilities.h"
#include "particle_groups.h"
#include "particles.h"
#include "check_particle.h"
#include "call_lua_function.h"

#include <boost/format.hpp>
#include <boost/ref.hpp>
#include <luabind/luabind.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace PIC {

double R(double x, double h)
{
    return ((fabs(x) < h) ? (1.0 - fabs(x)/h) : 0.0);
}

double R(double x, double y, double z, double h)
{
    return (R(x, h)*R(y, h)*R(z, h));
}

/*****************************************************************************
* Gather edge-centered values.                                               *
*    at_point - specified position                                           *
*    val_x, val_y, val_z - Cell required values                              *
*    ret_vec - result gathered value                                         *
* Used to interpolate Bx, By, Bz for specified position                      *
*****************************************************************************/
void gather_edge(const Grid& grid, const DblVector& at_point,
                 CellVectorValue val, DblVector& ret_vec)
{
    ret_vec.x = gather_vector<EdgeXCentering>(grid, at_point, val, &DblVector::x);
    ret_vec.y = gather_vector<EdgeYCentering>(grid, at_point, val, &DblVector::y);
    ret_vec.z = gather_vector<EdgeZCentering>(grid, at_point, val, &DblVector::z);
}

/*****************************************************************************
* Gather face-centered values: Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz      *
*    at_point - specified position                                           *
*    val_x, val_y, val_z - Pointers to cell members required                 *
*    ret_vec - result gathered value                                         *
* Used to interpolate Ex, Ey, Ez, UPx, UPy, UPz, UEx, UEy, UEz values        *
* for specified position.                                                    *
*****************************************************************************/
void gather_face(const Grid& grid, const DblVector& at_point,
                 CellVectorValue val, DblVector& ret_vec)

{
    ret_vec.x = gather_vector<FaceXCentering>(grid, at_point, val, &DblVector::x);
    ret_vec.y = gather_vector<FaceYCentering>(grid, at_point, val, &DblVector::y);
    ret_vec.z = gather_vector<FaceZCentering>(grid, at_point, val, &DblVector::z);
}

/*****************************************************************************
* Gather cell-centered value.                                                *
*   x, y, z - specified position                                             *
*   val - Pointer to cell member required                                    *
*   ret_val - result gathered value                                          *
* Used to interpolate density NP values for specified position.              *
*                                                                            *
*****************************************************************************/
double gather_center(const Grid& grid, const DblVector& at_point, CellScalarValue val)
{
    return gather_scalar<CellCentering>(grid, at_point, val);
}

/******************************************************
* interpolates grid values to specified point (x,y,z) *
******************************************************/
void from_grid_to_point(const Grid& grid, const DblVector& vec_to_point, Grid::NodeType& val_at_point)
{
    // interpolate NP
    val_at_point.NP = gather_center(grid, vec_to_point, &Cell::NP);

    // interpolate B
    gather_edge(grid, vec_to_point, &Cell::B, val_at_point.B);

    //interpolate E
    gather_face(grid, vec_to_point, &Cell::E, val_at_point.E);

    // interpolate UP
    gather_face(grid, vec_to_point, &Cell::UP, val_at_point.UP);

    // interpolate UE
    gather_face(grid, vec_to_point, &Cell::UE, val_at_point.UE);
}

void set_boundary_conditions(const std::string& lua_func_name)
{
    if (!call_lua_function(lua_func_name.c_str()))
    {
        std::string msg_fmt = "Please check function \"%s\" in \"%s\" and restart.";

        std::string msg = str(boost::format(msg_fmt) % lua_func_name
                                                     % PIC::Config::cfg_script_name());

        throw std::domain_error(msg);
    }
}

void set_boundary_conditions(const std::string& lua_func_name, DensityGrid& dg)
{
    UseLua lua;

    try
    {
        luabind::call_function<void>(lua, lua_func_name.c_str(), boost::ref(dg));
    }
    catch (...)
    {
        const std::string fmt = "Error: in function \"%s\" not defined or invalid in \"%s\".";
        std::string msg = str(boost::format(fmt)
                              % lua_func_name
                              % PIC::Config::cfg_script_name());
        throw std::domain_error(msg);
    }
}

void scatter_particle_std(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    const double h = Config::h();

    ParticleGroups part_groups;
    const double q = particle.ni*part_groups[particle.group_name].charge;

    const index_t pi = static_cast<index_t>((particle.r.x + dr.x)/h);
    const index_t pj = static_cast<index_t>((particle.r.y + dr.y)/h);
    const index_t pk = static_cast<index_t>((particle.r.z + dr.z)/h);

    typedef DensityGrid::NodeType DensityNode;

    for (index_t i = pi - 1; i != pi + 2; ++i)
    for (index_t j = pj - 1; j != pj + 2; ++j)
    for (index_t k = pk - 1; k != pk + 2; ++k)
    {
        DensityNode& cell = grid(i, j, k);

        cell.NP += q*R((i + 0.5)*h - particle.r.x,
                       (j + 0.5)*h - particle.r.y,
                       (k + 0.5)*h - particle.r.z, h);

        cell.UP.x += q*particle.v.x*R((i)*h       - particle.r.x,
                                      (j + 0.5)*h - particle.r.y,
                                      (k + 0.5)*h - particle.r.z, h);

        cell.UP.y += q*particle.v.y*R((i + 0.5)*h - particle.r.x,
                                      (j)*h       - particle.r.y,
                                      (k + 0.5)*h - particle.r.z, h);

        cell.UP.z += q*particle.v.z*R((i + 0.5)*h - particle.r.x,
                                      (j + 0.5)*h - particle.r.y,
                                      (k)*h       - particle.r.z, h);
    }
}

void scatter_particle_zigzag(const Particle& particle, DensityGrid& grid, const DblVector& dr)
{
    const double h = Config::h();

    const double tau2 = Config::tau_2();

    const double x1 = particle.r.x;
    const double x2 = x1 + dr.x;
    const index_t i[2] = { static_cast<index_t>(x1/h), static_cast<index_t>(x2/h) };

    const double y1 = particle.r.y;
    const double y2 = y1 + dr.y;
    const index_t j[2] = { static_cast<index_t>(y1/h), static_cast<index_t>(y2/h) };

    const double z1 = particle.r.z;
    const double z2 = z1 + dr.z;
    const index_t k[2] = { static_cast<index_t>(z1/h), static_cast<index_t>(z2/h) };
    ParticleGroups part_groups;
    const double q = particle.ni*part_groups[particle.group_name].charge;

    {
        // scatter charge Q and store it in NP
        const double wi[2] = { (i[1] + 1 - x2/h), (x2/h - i[1]) };
        const double wj[2] = { (j[1] + 1 - y2/h), (y2/h - j[1]) };
        const double wk[2] = { (k[1] + 1 - z2/h), (z2/h - k[1]) };

        for (int l = 0; l != 2; ++l)
        for (int m = 0; m != 2; ++m)
        for (int n = 0; n != 2; ++n)
        {
            grid(i[1] + l, j[1] + m, k[1] + n).NP += q*wi[l]*wj[m]*wk[n];
        }
    }

    // scatter currents Jx, Jy, Jz and store them  in UPx, UPy, UPz
    DblVector r(min(min(i[0]*h, i[1]*h) + h, max(max(i[0]*h, i[1]*h), 0.5*(x1 + x2))),
                min(min(j[0]*h, j[1]*h) + h, max(max(j[0]*h, j[1]*h), 0.5*(y1 + y2))),
                min(min(k[0]*h, k[1]*h) + h, max(max(k[0]*h, k[1]*h), 0.5*(z1 + z2))));

    DblVector F[2] = { DblVector(q*(r.x - x1)/tau2, q*(r.y - y1)/tau2, q*(r.z - z1)/tau2),
                       DblVector(q*(x2 - r.x)/tau2, q*(y2 - r.y)/tau2, q*(z2 - r.z)/tau2) };

    DblVector W[2] = { DblVector(0.5*(x1 + r.x)/h - i[0], 0.5*(y1 + r.y)/h - j[0], 0.5*(z1 + r.z)/h - k[0]),
                       DblVector(0.5*(r.x + x2)/h - i[1], 0.5*(r.y + y2)/h - j[1], 0.5*(r.z + z2)/h - k[1]) };

    for (int m = 0; m != 2; ++m)
    {
        grid(i[m], j[m], k[m]).UP.x += F[m].x*(1.0 - W[m].y)*(1.0 - W[m].z);
        grid(i[m], j[m] + 1, k[m]).UP.x += F[m].x*W[m].y*(1.0 - W[m].z);
        grid(i[m], j[m], k[m] + 1).UP.x += F[m].x*(1.0 - W[m].y)*W[m].z;
        grid(i[m], j[m] + 1, k[m] + 1).UP.x += F[m].x*W[m].y*W[m].z;

        grid(i[m], j[m], k[m]).UP.y += F[m].y*(1.0 - W[m].x)*(1.0 - W[m].z);
        grid(i[m] + 1, j[m], k[m]).UP.y += F[m].y*W[m].x*(1.0 - W[m].z);
        grid(i[m], j[m], k[m] + 1).UP.y += F[m].y*(1.0 - W[m].x)*W[m].z;
        grid(i[m] + 1, j[m], k[m] + 1).UP.y += F[m].y*W[m].x*W[m].z;

        grid(i[m], j[m], k[m]).UP.z += F[m].z*(1.0 - W[m].x)*(1.0 - W[m].y);
        grid(i[m] + 1, j[m], k[m]).UP.z += F[m].z*W[m].x*(1.0 - W[m].y);
        grid(i[m], j[m] + 1, k[m]).UP.z += F[m].z*(1.0 - W[m].x)*W[m].y;
        grid(i[m] + 1, j[m] + 1, k[m]).UP.z += F[m].z*W[m].x*W[m].y;
    }
}

} // end of namespace PIC
