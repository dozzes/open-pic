#include <cmath>
#include <iostream>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/ref.hpp>

#include <luabind/luabind.hpp>

#include "particle_groups.h"
#include "particles.h"

#include "check_particle.h"
#include "config.h"

#include "use_lua.h"
#include "io_utilities.h"
#include "call_lua_function.h"

#include "gather_scatter.h"


namespace PIC {

double R(double x, double h)
{
    return ((fabs(x) < h) ? (1.0 - fabs(x) / h) : 0.0);
}

double R(double x, double y, double z, double h)
{
    return (R(x, h) * R(y, h) * R(z, h));
}

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
                 DblVector& ret_vec)
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
void gather_face(const Grid& grid,
                 const DblVector& at_point,
                 CellVectorValue val,
                 DblVector& ret_vec)

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
    catch (std::exception& e)
    {
        std::string msg_fmt = "Error: \"%s\".\nFunction \"%s\" not defined or invalid in \"%s\".";
        std::string msg = str(boost::format(msg_fmt) % e.what()
            % lua_func_name
            % PIC::Config::cfg_script_name());
        throw std::domain_error(msg);
    }
}

void scatter_particle_std(const Particle& p, DensityGrid& dg)
{
    const double h = Config::h();

    ParticleGroups part_groups;
    const double q = p.ni*part_groups[p.group_name].charge;

    const size_t pi = static_cast<size_t>(p.r.x / h);
    const size_t pj = static_cast<size_t>(p.r.y / h);
    const size_t pk = static_cast<size_t>(p.r.z / h);

    typedef DensityGrid::NodeType DensityNode;

    for (size_t i = pi - 1; i != pi + 2; ++i)
    for (size_t j = pj - 1; j != pj + 2; ++j)
    for (size_t k = pk - 1; k != pk + 2; ++k)
    {
        DensityNode& cell = dg(i, j, k);

        cell.NP += q * R((i + 0.5)*h - p.r.x,
                         (j + 0.5)*h - p.r.y,
                         (k + 0.5)*h - p.r.z, h);

        cell.UP.x += q * p.v.x * R((i)*h - p.r.x,
                                   (j + 0.5)*h - p.r.y,
                                   (k + 0.5)*h - p.r.z, h);

        cell.UP.y += q * p.v.y * R((i + 0.5)*h - p.r.x,
                                   (j)*h - p.r.y,
                                   (k + 0.5)*h - p.r.z, h);

        cell.UP.z += q * p.v.z * R((i + 0.5)*h - p.r.x,
                                   (j + 0.5)*h - p.r.y,
                                   (k)*h - p.r.z, h);
    }
}

} // end of namespace PIC
