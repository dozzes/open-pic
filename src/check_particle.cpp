#include "check_particle.h"
#include "opic_fwd.h"
#include "config.h"
#include "io_utilities.h"
#include "grid.h"
#include "particles.h"
#include "gather_scatter.h"
#include "call_lua_function.h"

#include <boost/format.hpp>
#include <omp.h>
#include <stdexcept>

 /*********************************************************
 * Check if particle is active.                           *
 * equation of motion is solved for active particles only *
 * return true if particles is active, otherwise - false  *
 *********************************************************/
bool is_particle_can_scatter(Particles& particles, index_t p, const Grid& grid, const DblVector& dr)
{
    if (particles.is_inactive(p))
        return false;

    const Particle& particle = particles[p];
    const double h = PIC::Config::h();
    const index_t pi = static_cast<index_t>(floor((particle.r.x + dr.x)/h));
    const index_t pj = static_cast<index_t>(floor((particle.r.y + dr.y)/h));
    const index_t pk = static_cast<index_t>(floor((particle.r.z + dr.z)/h));
    const PIC::CellState home_cell_state = grid(pi, pj, pk).state();

    if (pi < 1 || pi > grid.size_x() - 2 ||
        pj < 1 || pj > grid.size_y() - 2 ||
        pk < 1 || pk > grid.size_z() - 2 ||
        home_cell_state == PIC::cs_absorptive)
    {
        particles.remove_later(p);
        return false;
    }

    return true;
}

bool is_particle_can_move(Particles& particles, index_t p, const Grid& grid)
{
    if (!is_particle_can_scatter(particles, p, grid, DblVector()))
        return false;

    Particle& particle = particles[p];
    const double h = PIC::Config::h();
    const index_t pi = static_cast<index_t>(floor(particle.r.x/h));
    const index_t pj = static_cast<index_t>(floor(particle.r.y/h));
    const index_t pk = static_cast<index_t>(floor(particle.r.z/h));
    const PIC::CellState home_cell_state = grid(pi,pj,pk).state();

    if (home_cell_state == PIC::cs_active)
        return true;

    if (home_cell_state == PIC::cs_custom)
        return lua_validate_particle(particle);

    // never reach here
    return false;
}

bool check_particle_move(const Particle& particle, const Grid& grid, const DblVector& dr)
{
    if (PIC::Config::CFL_severity() == PIC::Ignore)
        return true;

    const double h = PIC::Config::h();
    const double h_2 = PIC::Config::h_2();

    if (fabs(dr.x) > h_2 || fabs(dr.y) > h_2 || fabs(dr.z) > h_2)
    {
        const index_t step = PIC::Config::current_time_step();

        Grid::NodeType point_val;
        PIC::from_grid_to_point(grid, particle.r, point_val);

        const index_t pi = static_cast<index_t>(particle.r.x/h);
        const index_t pj = static_cast<index_t>(particle.r.y/h);
        const index_t pk = static_cast<index_t>(particle.r.z/h);

        const std::string msg = boost::str(boost::format("Step = %1%:\n\tError CFL, particle shift dr should be  < h/2:"
                                                  "\n\tparticle's home cell = (%2%, %3%, %4%);"
                                                  "\n\tgrid size = (%5% x %6% x %7%);"
                                                  "\n\tdr = (%8%, %9%, %10%),\n\tdr.abs() = %11%; h/2 = %12%;"
                                                  "\n\tparticle [%13%] at (%14%, %15%, %16%)"
                                                  "\n\tgrid values at particle position:"
                                                  "\n\t\t B = %17%, Bx = %18%, By = %19%, Bz = %20%"
                                                  "\n\t\t E = %21%, Ex = %22%, Ey = %23%, Ez = %24%"
                                                  "\n\t\t NP = %25%, UPx = %26%, UPy = %27%, UPz = %28%"
                                                  "\n\t\t UEx = %29%, UEy = %30%, UEz = %31%"
                                                  "\n\t\t cell state = %32% (0 - active, 1 - absorptive, 2 - custom)")
                                                  % step
                                                  % pi % pj % pk
                                                  % grid.size_x() % grid.size_y() % grid.size_z()
                                                  % dr.x % dr.y % dr.z % (dr.abs())
                                                  % h_2
                                                  % particle.group_name
                                                  % particle.r.x % particle.r.y % particle.r.z
                                                  % point_val.B.abs() % point_val.B.x % point_val.B.y % point_val.B.z
                                                  % point_val.E.abs() % point_val.E.x % point_val.E.y % point_val.E.z
                                                  % point_val.NP % point_val.UP.x % point_val.UP.y % point_val.UP.z
                                                  % point_val.UE.x % point_val.UE.y % point_val.UE.z
                                                  % grid(pi,pj,pk).state());

        const int tid = omp_get_thread_num();
        const std::string log_file_name = str(boost::format("opic_thread_%1%_check_particle_move_err.log") % tid);
        std::ofstream ofs_log(log_file_name.c_str(), std::ios_base::app);
        ofs_log << msg << std::endl;
        return false;
    }

    return true;
}
