#include <stdexcept>
#include <omp.h>
#include <boost/format.hpp>

#include "config.h"
#include "call_lua_function.h"
#include "grid.h"
#include "particles.h"
#include "check_particle.h"
#include "particle_groups.h"
#include "gather_scatter.h"
#include "save_grid.h"
#include "save_particles.h"

#include "simulate.h"


namespace PIC {

/************************************************************
*                                                           *
* I step                                                    *
*                                                           *
* Calculate magnetic field in half time step  (m+1/2)I step *
*                                                           *
************************************************************/
void calc_magnetic_field_half_time(Grid& grid)
{
   const index_t kd = 1;

   const index_t to_kx = (grid.size_x() - kd);
   const index_t to_ky = (grid.size_y() - kd);
   const index_t to_kz = (grid.size_z() - kd);

   const double ctau_2h = Config::ctau_2h();
   for (index_t kx = kd; kx != to_kx; ++kx)
   for (index_t ky = kd; ky != to_ky; ++ky)
   for (index_t kz = kd; kz != to_kz; ++kz)
   {
      Cell& cell = grid(kx,ky,kz);

      if (cell.state() == PIC::cs_active)
      {
         cell.B.x -= ctau_2h * ((cell.E.z - grid(kx,ky - 1,kz).E.z) -
                                (cell.E.y - grid(kx,ky,kz - 1).E.y));

         cell.B.y -= ctau_2h * ((cell.E.x - grid(kx,ky,kz - 1).E.x) -
                                (cell.E.z - grid(kx - 1,ky,kz).E.z));

         cell.B.z -= ctau_2h * ((cell.E.y - grid(kx - 1,ky,kz).E.y) -
                                (cell.E.x - grid(kx,ky - 1,kz).E.x));
      }
   }
}

/*******************************************
*                                          *
* II step                                  *
*                                          *
* Move particles in half time step (m+1/2) *
*                                          *
*******************************************/

void push_particle_std(const Grid& grid, Particle& rp)
{
   DblVector Bp;
   gather_edge(grid, rp.r, &Cell::B, Bp);

   DblVector Ep;
   gather_face(grid, rp.r, &Cell::E, Ep);

   ParticleGroups part_groups;
   const ParticleGroups::ParticleGroup& part_group = part_groups[rp.group_name];

   const double qm = part_group.charge / part_group.mass * Constants::e_mp();

   const double tau_qm = Config::tau() * qm;
   const double tau_2qmc = 0.5 * tau_qm / Constants::c();

   const double a = tau_2qmc * Bp.x;
   const double b = tau_2qmc * Bp.y;
   const double c = tau_2qmc * Bp.z;

   const double A = rp.v.x + tau_qm * Ep.x + c*rp.v.y - b*rp.v.z;
   const double B = rp.v.y + tau_qm * Ep.y + a*rp.v.z - c*rp.v.x;
   const double C = rp.v.z + tau_qm * Ep.z + b*rp.v.x - a*rp.v.y;

   const double a2 = a * a;
   const double b2 = b * b;
   const double c2 = c * c;

   const double D = 1.0/(a2 + b2 + c2 + 1.0);

   const double ab = a * b;
   const double ac = a * c;
   const double bc = b * c;

   /* Mathematica solution:
   x -> (A + a^2 A + a b B + B c - b C + a c C)/(1 + a^2 + b^2 + c^2),
   y -> (a A b + B + b^2 B - A c + a C + b c C)/(1 + a^2 + b^2 + c^2),
   z -> (A b - a B + a A c + b B c + C + c^2 C)/(1 + a^2 + b^2 + c^2).
   */

   rp.v.x = D*(A + a2*A + ab*B + B*c - b*C + ac*C);
   rp.v.y = D*(ab*A + B + b2*B - A*c + a*C + bc*C);
   rp.v.z = D*(A*b - a*B + ac*A + bc*B + C + c2*C);
}

void push_particle_boris(const Grid& grid, Particle& rp)
{
   DblVector B;
   gather_edge(grid, rp.r, &Cell::B, B);

   DblVector E;
   gather_face(grid, rp.r, &Cell::E, E);

   ParticleGroups part_groups;
   const ParticleGroups::ParticleGroup& part_group = part_groups[rp.group_name];

   const double qm = part_group.charge / part_group.mass * Constants::e_mp();
   const double tau_2qm = 0.5*Config::tau()*qm;
   const double tau_2qmc = tau_2qm / Constants::c();

   const DblVector Vm(rp.v.x + tau_2qm*E.x,
                      rp.v.y + tau_2qm*E.y,
                      rp.v.z + tau_2qm*E.z);

   const DblVector V0(Vm.x + tau_2qmc*(Vm.y*B.z - Vm.z*B.y),
                      Vm.y + tau_2qmc*(Vm.z*B.x - Vm.x*B.z),
                      Vm.z + tau_2qmc*(Vm.x*B.y - Vm.y*B.x));

   const double d = 2.0*tau_2qmc/(1.0 + tau_2qmc*tau_2qmc*(B.x*B.x + B.y*B.y + B.z*B.z));

   const DblVector Vp(Vm.x + d*(V0.y*B.z - V0.z*B.y),
                      Vm.y + d*(V0.z*B.x - V0.x*B.z),
                      Vm.z + d*(V0.x*B.y - V0.y*B.x));

   rp.v.x = Vp.x + tau_2qm*E.x;
   rp.v.y = Vp.y + tau_2qm*E.y;
   rp.v.z = Vp.z + tau_2qm*E.z;
}

void (*push_particle)(const Grid&, Particle&);

vector<DensityGrid> densityGrids;

void (*scatter_particle)(const Particle&, DensityGrid&);


void move_particles_half_time(const Grid& grid,
                              Particles& particles,
                              const std::string& group_name)
{
   size_t removed_num = particles.remove_absorbed();

   if (removed_num > 0)
   {
      std::cout << removed_num << " particle(s) removed.\n";
   }

   const double tau_2 = Config::tau_2();

   long particles_num = particles.size();

   long err_count = 0;

   #pragma omp parallel for reduction(+: err_count)
   for (long p = 0; p < particles_num; ++p)
   {
      try
	  {
		  Particle& particle = particles[p];

		  if ((group_name == "all" || particle.group_name == group_name) &&
			  validate_particle(particle, grid) &&
			  (err_count == 0))
		  {
			  push_particle(grid, particle);

			  const DblVector dr(tau_2 * particle.v.x, tau_2 * particle.v.y, tau_2 * particle.v.z);

			  if (check_particle_move(particle, dr))
			  {
				  particle.r.x += dr.x;
				  particle.r.y += dr.y;
				  particle.r.z += dr.z;

				  if(validate_particle(particle, grid))
				  {
					  const int thread_num = omp_get_thread_num();
					  scatter_particle(particle, densityGrids[thread_num]);
				  }
			  }
			  else
			  {
				  ++err_count;
			  }
		  }
	  } 
	  catch (...)
	  {
		  ++err_count;
	  }
   } // for (long p = 0; p < particles_num; ++p)

   if (err_count != 0)
   {
      throw domain_error("CFL condition violated! See opic_thread_*_error.log files for details.");
   }
}

/***********************************
*                                  *
* IV step                          *
*                                  *
* Move particles in full time step *
*                                  *
***********************************/
void move_particles_full_time(const Grid& grid,
                              Particles& particles,
                              const std::string& group_name) // m+1
{
   const double tau_2 = Config::tau_2();

   long particles_num = particles.size();

   #pragma omp parallel for
   for (long p = 0; p < particles_num; ++p)
   {
      Particle& particle = particles[p];

      if ((group_name == "all" || particle.group_name == group_name) &&
          validate_particle(particle, grid))
      {
         particle.r.x += tau_2 * particle.v.x;
         particle.r.y += tau_2 * particle.v.y;
         particle.r.z += tau_2 * particle.v.z;
      }
   }
}

/***************************************************
*                                                  *
* V step ( repeat step I )                         *
*                                                  *
* Calculate magnetic field in full time step (m+1) *
*                                                  *
***************************************************/

/***************************************************
*                                                  *
* VI step                                          *
*                                                  *
* Calculate electron velocities                    *
*                                                  *
***************************************************/

void calc_electrons_velocity(Grid& grid)
{
   const index_t kd = 1;

   const index_t to_kx = (grid.size_x() - kd);
   const index_t to_ky = (grid.size_y() - kd);
   const index_t to_kz = (grid.size_z() - kd);

   const double c_4pi_e_h = Config::c_4pi_e_h();

   double NP = 0.0;
   for (index_t kx = kd; kx != to_kx; ++kx)
   for (index_t ky = kd; ky != to_ky; ++ky)
   for (index_t kz = kd; kz != to_kz; ++kz)
   {
      Cell& cell = grid(kx,ky,kz);

      if (cell.state() == PIC::cs_active)
      {
         // NP for UEx
         NP = 0.5 * (cell.NP + grid(kx-1,ky,kz).NP);
         cell.UE.x = cell.UP.x - c_4pi_e_h / NP * ((grid(kx,ky + 1,kz).B.z - cell.B.z) -
                                                   (grid(kx,ky,kz + 1).B.y - cell.B.y));
         // NP for UEy
         NP = 0.5 * (cell.NP + grid(kx,ky-1,kz).NP);
         cell.UE.y = cell.UP.y - c_4pi_e_h / NP * ((grid(kx,ky,kz + 1).B.x - cell.B.x) -
                                                   (grid(kx + 1,ky,kz).B.z - cell.B.z));
         //NP for UEz
         NP = 0.5 * (cell.NP + grid(kx,ky,kz-1).NP);
         cell.UE.z = cell.UP.z - c_4pi_e_h / NP * ((grid(kx + 1,ky,kz).B.y - cell.B.y) -
                                                   (grid(kx,ky + 1,kz).B.x - cell.B.x));
      }
   }
}

/***************************************************
*                                                  *
* VI step                                          *
*                                                  *
* Calculate electric field                         *
*                                                  *
***************************************************/
void calc_electric_field(Grid& grid)
{
   DblVector UE, B;

   const index_t kd = 1;

   const index_t to_kx = (grid.size_x() - kd);
   const index_t to_ky = (grid.size_y() - kd);
   const index_t to_kz = (grid.size_z() - kd);

   for (index_t kx = kd; kx != to_kx; ++kx)
   for (index_t ky = kd; ky != to_ky; ++ky)
   for (index_t kz = kd; kz != to_kz; ++kz)
   {
      Cell& cell = grid(kx,ky,kz);

      if (cell.state() == PIC::cs_active)
      {
         // calculate Ex
         B.y = 0.5 * (cell.B.y + grid(kx,ky,kz + 1).B.y);
         B.z = 0.5 * (cell.B.z + grid(kx,ky + 1,kz).B.z);

         UE.y = 0.25 * (cell.UE.y + grid(kx,ky + 1,kz).UE.y +
                         grid(kx-1,ky,kz).UE.y + grid(kx - 1,ky + 1,kz).UE.y);

         UE.z = 0.25 * (cell.UE.z + grid(kx,ky,kz + 1).UE.z+
                         grid(kx-1,ky,kz).UE.z + grid(kx - 1,ky,kz + 1).UE.z);

         cell.E.x = (UE.z * B.y - UE.y * B.z) / Constants::c();

         // calculate Ey
         B.z = 0.5 * (cell.B.z + grid(kx + 1,ky,kz).B.z);
         B.x = 0.5 * (cell.B.x + grid(kx,ky,kz + 1).B.x);

         UE.x = 0.25 * (cell.UE.x + grid(kx + 1,ky,kz).UE.x +
                         grid(kx,ky - 1,kz).UE.x + grid(kx + 1,ky - 1,kz).UE.x);

         UE.z = 0.25 * (cell.UE.z + grid(kx,ky,kz + 1).UE.z+
                         grid(kx,ky - 1,kz).UE.z + grid(kx,ky - 1,kz + 1).UE.z);

         cell.E.y = (UE.x * B.z - UE.z * B.x) / Constants::c();

         // calculate Ez
         B.x = 0.5 * (cell.B.x + grid(kx,ky + 1,kz).B.x);
         B.y = 0.5 * (cell.B.y + grid(kx + 1,ky,kz).B.y);

         UE.y = 0.25 * (cell.UE.y + grid(kx,ky + 1,kz).UE.y +
                         grid(kx,ky,kz - 1).UE.y + grid(kx,ky + 1,kz - 1).UE.y);

         UE.x = 0.25 * (cell.UE.x + grid(kx + 1,ky,kz).UE.x +
                         grid(kx,ky,kz - 1).UE.x + grid(kx + 1,ky,kz - 1).UE.x);

         cell.E.z = (UE.y * B.x - UE.x * B.y) / Constants::c();
      }
   }
}

template <class DensityGridType>
void set_grid_UP(const Grid& grid, DensityGridType& dg_group)
{
   double NP = 0.0;

   const index_t m = 1;

   for (index_t i = m; i != (grid.size_x() - m); ++i)
   for (index_t j = m; j != (grid.size_y() - m); ++j)
   for (index_t k = m; k != (grid.size_z() - m); ++k)
   {
      const CellState cell_state = grid(i,j,k).state();

      if (cell_state == PIC::cs_active)
      {
         typename DensityGridType::NodeType& dn = dg_group(i,j,k);

         // NP for UPx
         NP = 0.5 * (dn.NP + dg_group(i-1,j,k).NP);
         NP > 1.0e-6 ? dn.UP.x /= NP : dn.UP.x = 0.0;

         // NP for UPy
         NP = 0.5 * (dn.NP + dg_group(i,j-1,k).NP);
         NP > 1.0e-6 ? dn.UP.y /= NP : dn.UP.y = 0.0;

         //NP for UEz
         NP = 0.5 * (dn.NP + dg_group(i,j,k-1).NP);
         NP > 1.0e-6 ? dn.UP.z /= NP : dn.UP.z = 0.0;
      }
   }
}

template <class DensityGridType>
void normalize_NP(const Grid& grid, DensityGridType& dg_group)
{
    const double cell_volume = grid.cell_volume();

    for (index_t i = 0; i != grid.size_x(); ++i)
    for (index_t j = 0; j != grid.size_y(); ++j)
    for (index_t k = 0; k != grid.size_z(); ++k)
    {
        const CellState cell_state = grid(i,j,k).state();

        if (cell_state == cs_active)
        {
            dg_group(i,j,k).NP /= cell_volume;
        }
    }
}

void local_Alfven_CFL(Grid& grid, size_t i, size_t j, size_t k)
{
    Cell& cell = grid(i,j,k);

    const double h = grid.step();
    DblVector vec_to_point(i*h, j*h, k*h);
    DblVector vec_B;
    gather_edge(grid, vec_to_point, &Cell::B, vec_B);
    const double B = vec_B.abs();

    const double mp = ParticleGroups::max_mp() * PIC::Constants::mp();
    const double pi = Constants::pi();
    const double tau_2 = Config::tau_2();

    const double min_dens =  1.5*B*B*tau_2*tau_2/(4*pi*h*h*mp);
    const double cell_volume = grid.cell_volume();

    // At this stage grid stores NP*cell_volume
    if (cell.NP < min_dens*cell_volume)
    {
        PIC::Config::logger() << "\nlocal_Alfven_CFL:"
                              << "(i,j,k) = (" << i << "," << j << "," << k << ")"
                              << " cell.NP = " << cell.NP
                              << " min_dens = " << min_dens
                              << std::endl;

        cell.NP = min_dens*cell_volume;
    }
}

void backgr_fracture(Grid& grid, size_t i, size_t j, size_t k)
{
    Cell& cell = grid(i,j,k);

    const double dens_cutoff = PIC::Config::dens_cutoff();
    const double cell_volume = grid.cell_volume();

     // At this stage grid stores NP*cell_volume
    if (cell.NP < dens_cutoff*cell_volume)
    {
        PIC::Config::logger() << "\nbackgr_fracture:"
                              << "(i,j,k) = (" << i << "," << j << "," << k << ")"
                              << " cell.NP = " << cell.NP/grid.cell_volume()
                              << " dens_cutoff = " << dens_cutoff
                              << std::endl;

        cell.NP = dens_cutoff*cell_volume;
    }
}

void (*grid_threshold)(Grid& grid, size_t i, size_t j, size_t k);

void set_threshold(Grid& grid)
{
    // At this stage grid stores NP*grid.cell_volume()
    for (index_t i = 1; i != grid.size_x()-1; ++i)
    for (index_t j = 1; j != grid.size_y()-1; ++j)
    for (index_t k = 1; k != grid.size_z()-1; ++k)
    {
        const CellState cell_state = grid(i,j,k).state();

        if (cell_state == PIC::cs_active)
        {
            grid_threshold(grid, i, j, k);
        }
    }
}

template<class GridType>
void save_grid_node_data(const std::string& prefix, const GridType& grid)
{
    using namespace PIC;

    std::string f_name = create_out_file_name(prefix, "grd_nodes_z=4_", PIC::Config::current_time_step());
    ofstream ofsg(f_name.c_str());
    //ofsg <<"x\ty\tNP\tUP\tUPx\tUPy\tUPz\tB\tBx\tBy\tBz\tE\tEx\tEy\tEz\n";
    ofsg <<"x\ty\tNP\tUP\tUPx\tUPy\tUPz\n";

    const size_t k = 4;
    for(size_t i = 0; i != (grid.size_x()); ++i)  
    for(size_t j = 0; j != (grid.size_y()); ++j)
    {
        const GridType::NodeType& c = grid(i,j,k);
        const double UP = sqrt(c.UP.x*c.UP.x + c.UP.y*c.UP.y + c.UP.z*c.UP.z);
        //const double B = sqrt(c.Bx*c.Bx + c.By*c.By + c.Bz*c.Bz);
        //const double E = sqrt(c.Ex*c.Ex + c.Ey*c.Ey + c.Ez*c.Ez);
        ofsg << i << "\t" << j
             << "\t" << c.NP
             << "\t" << UP << "\t"<< c.UP.x << "\t"<< c.UP.y << "\t" << c.UP.z
            //<< "\t" << B  << "\t"<< c.Bx  << "\t"<< c.By  << "\t" << c.Bz
            //<< "\t" << E  << "\t"<< c.Ex  << "\t"<< c.Ey  << "\t" << c.Ez
            << endl;
    }
}

void simulate(Grid& grid, Particles& particles)
{
   std::cout << "\n\nPIC simulation started ...\n\n";

   push_particle = &push_particle_std;
   if (Config::particle_push_alg() == Boris)
   {
       push_particle = push_particle_boris;
   }

   std::cout << "\nParticle push algorithm: "
             << (Config::particle_push_alg() == Direct ? "Direct" :
                (Config::particle_push_alg() == Boris  ? "Boris"  : "Undefined"));

   scatter_particle = &scatter_particle_std;

   grid_threshold = &backgr_fracture;
   if (Config::grid_threshold() == Local_CFL)
   {
       grid_threshold = &local_Alfven_CFL;
   }

   const int max_threads = omp_get_max_threads();
   densityGrids = vector<DensityGrid>(max_threads,
                                      DensityGrid(grid.size_x(), grid.size_y(),grid.size_z(), grid.step()));;

   omp_set_dynamic(0); // Explicitly disable dynamic teams
   omp_set_num_threads(max_threads); // Use max threads for all consecutive parallel regions

   DensityGrid dg_group(grid.size_x(), grid.size_y(), grid.size_z(), grid.step());

   const index_t steps_num = Config::time_steps();

   for (index_t step = 0; step != steps_num; ++step)
   {
      Config::set_current_time_step(step);

      std::cout << "\nstep " << step << " of " << PIC::Config::time_steps()
                << " : " << particles.size() << " particle(s)" << std::endl;

      // reset charge and current densities on grid.
      grid.reset_current();

      print_tm("Call \"on iteration begin\" handler");
      call_lua_function("on_iteration_begin");

      print_tm("Calculate magnetic field with half time step");
      calc_magnetic_field_half_time(grid);

      // move particle on half time step.
      if ( Config::is_on_save_step() )
      {
         ParticleGroups part_groups;
         const std::vector<std::string>& group_names = part_groups.group_names();

         // move particle groups half time step separately
         // to allow saving their grid data in separate files.
         for (size_t g = 0; g != group_names.size(); ++g)
         {
            const std::string& group_name = group_names[g];

            print_tm("Call \"particles moved on half time\" handler");
            call_lua_function("on_particles_moved_half_time");

            print_tm("Move particles in half time step: ", group_name);
            move_particles_half_time(grid, particles, group_name);

            // reset particle group charge and current densities.
            dg_group.reset_current();

            // accumulate particle group charge and current densities
            // scattered in particles push threads.
            for(int t = 0; t != max_threads; ++t)
            {
               DensityGrid& dg_thread = densityGrids[t];
               dg_group.add_current(dg_thread);
               dg_thread.reset_current();
            }

            // add particle group charge and current densities to grid.
            grid.add_current(dg_group);

            // calculate UP of particle group.
            set_grid_UP(grid, dg_group);
            set_boundary_conditions("on_set_boundary_group_UP", dg_group);

            // calculate NP of particle group.
            normalize_NP(grid, dg_group);
            set_boundary_conditions("on_set_boundary_group_NP", dg_group);

            // save group grid data.
            print_tm("Save grid data: ", group_name);
            save_grid(grid, dg_group, group_name);

            //save_grid_node_data(group_name, dg_group);

         } // for (size_t g = 0; g!= group_names.size(); ++g)
      }
      else
      {

         print_tm("Call \"particles moved on half time\" handler");
         call_lua_function("on_particles_moved_half_time");

         // move all particle on half time step.
         print_tm("Move particles in half time step: ", "all");
         move_particles_half_time(grid, particles, "all");

         // accumulate all particles charge and current densities
         // scattered in particles push threads.
         for(int t = 0; t != max_threads; ++t)
         {
            DensityGrid& dg = densityGrids[t];
            grid.add_current(dg);
            dg.reset_current();
         }
      }

      // set velocity of all particles.
      print_tm("Set boundary conditions for UP");
      set_grid_UP(grid, grid);
      set_boundary_conditions("on_set_boundary_UP");

      // set charge density all particles.
      print_tm("Set boundary conditions for NP");
      set_threshold(grid);
      normalize_NP(grid, grid);
      set_boundary_conditions("on_set_boundary_NP");

      if ( Config::is_on_save_step() )
      {
         print_tm("Save particles data");
         save_particles(grid, particles);

         print_tm("Save grid data: ", "all");
         save_grid(grid, grid, "all");

         //save_grid_node_data("all", grid);
      }

      // move all particle on full time step
      print_tm("Move particles in full time step: ", "all");
      move_particles_full_time(grid, particles, "all");

      print_tm("Call \"particles moved on full time\" handler");
      call_lua_function("on_particles_moved_full_time");

      print_tm("Calculate magnetic field with full time step");
      calc_magnetic_field_half_time(grid);
      set_boundary_conditions("on_set_boundary_MF");

      print_tm("Calculate electron velocities");
      calc_electrons_velocity(grid);
      set_boundary_conditions("on_set_boundary_UE");

      print_tm("Calculate electric field");
      calc_electric_field(grid);
      set_boundary_conditions("on_set_boundary_EF");

      print_tm("Call \"on iteration end\" handler");
      call_lua_function("on_iteration_end");
   }
}

} // namespace PIC
