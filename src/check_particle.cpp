#include <stdexcept>
#include <omp.h>
#include <boost/format.hpp>

#include "opic_fwd.h"
#include "io_utilities.h"
#include "config.h"
#include "grid.h"
#include "particles.h"
#include "gather_scatter.h"
#include "call_lua_function.h"

#include "check_particle.h"


namespace {

 /*********************************************************
 * check if particle is active.                           *
 * equation of motion is solved for active particles only *
 * return true if particles is active, otherwise - false  *
 *********************************************************/
 bool is_particle_active(Particle& particle, const Grid& grid)
 {
     if (particle.is_absorbed)
         return false;

     const double h = grid.step();

     const index_t pi = static_cast<index_t>(floor(particle.r.x/h));
     const index_t pj = static_cast<index_t>(floor(particle.r.y/h));
     const index_t pk = static_cast<index_t>(floor(particle.r.z/h));

     if (pi < 1 || pi > grid.size_x() - 2 ||
         pj < 1 || pj > grid.size_y() - 2 ||
         pk < 1 || pk > grid.size_z() - 2)
     {
         return false;
     }

     PIC::CellState home_cell_state = grid(pi,pj,pk).state();

     if (home_cell_state == PIC::cs_active)
     {
         return true;
     }

     if (home_cell_state == PIC::cs_absorptive)
     {
         return false;
     }

     if (home_cell_state == PIC::cs_custom)
     {
         return lua_validate_particle(particle);
     }

     return true;
 }

} // anonymous namespace

bool validate_particle(Particle& particle, const Grid& grid)
{
    if (particle.is_absorbed)
        return false;

    particle.is_absorbed = !is_particle_active(particle, grid);

    return !particle.is_absorbed;
}

bool check_particle_move(Particle& particle, const Grid& grid, const DblVector& dr)
{
    const double h = PIC::Config::h();
    const double h_2 = 0.5*h;

    if (fabs(dr.x) > h_2 || fabs(dr.y) > h_2 || fabs(dr.z) > h_2)
    {
        const index_t step = PIC::Config::current_time_step();

        Grid::NodeType point_val;
        PIC::from_grid_to_point(grid, particle.r, point_val);

        const index_t pi = static_cast<size_t>(particle.r.x/h);
        const index_t pj = static_cast<size_t>(particle.r.y/h);
        const index_t pk = static_cast<size_t>(particle.r.z/h);

        const std::string msg = boost::str(boost::format("Step = %1%:\n\tError CFL:"
                                                  "\n\tparticle's home cell = (%2%, %3%, %4%);"
                                                  "\n\tdr = (%5%, %6%, %7%),\n\tdr.abs() = %8%; h = %9%;"
                                                  "\n\tparticle [%10%] at (%11%, %12%, %13%)"
                                                  "\n\tfield values at particle position:"
                                                  "\n\t\t B = %14%, Bx = %15%, By = %16%, Bz = %17%"
                                                  "\n\t\t E = %18%, Ex = %19%, Ey = %20%, Ez = %21%"
                                                  "\n\t\t cell state = %22% (0 - active, 1 - absorptive, 2 - custom)")
                                                  % step
                                                  % pi
                                                  % pj
                                                  % pk
                                                  % dr.x % dr.y % dr.z % (dr.abs())
                                                  % h
                                                  % particle.group_name
                                                  % particle.r.x % particle.r.y % particle.r.z
                                                  % point_val.B.abs() % point_val.B.x % point_val.B.y % point_val.B.z
                                                  % point_val.E.abs() % point_val.E.x % point_val.E.y % point_val.E.z
                                                  % grid(pi, pj, pk).state());

        const int tid = omp_get_thread_num();
        const std::string log_file_name = str(boost::format("opic_thread_%1%_error.log") % tid);
        std::ofstream ofs_log(log_file_name.c_str(), std::ios_base::app);
        ofs_log << msg << std::endl;

        if (PIC::Config::CFL_severity() == PIC::Absorb)
        {
            particle.is_absorbed = true;
            return true;
        }

        if (PIC::Config::CFL_severity() == PIC::Ignore)
        {
            return true;
        }

        // Here PIC::Config::CFL_severity() == PIC::Stop
        return false;
    }

    return true;
}
