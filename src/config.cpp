#include "config.h"

#include <ostream>
#include <fstream>
#include <stdexcept>
#include <functional>

// static
PIC::Config::Parameters PIC::Config::params = Config::Parameters();

namespace {

static std::ofstream out("opic_trace.log");

void check_if_positive(const std::string& msg, double val)
{
    if (val < 0.0) throw std::invalid_argument(msg);
}

} // unnamed namespace

// static
std::ofstream& PIC::Config::ofs_log = out;

void PIC::Config::set_dens_cutoff(double dens_cutoff)
{
    check_if_positive("min NP must be > 0.0.", dens_cutoff);
    params.dens_cutoff = dens_cutoff;
}

void PIC::Config::set_tau(double tau)
{
    check_if_positive("tau must be > 0.0.", tau);
    params.tau = tau;
}

void PIC::Config::set_L_scale(double L_scale)
{
    check_if_positive("L_scale must be > 0.0.", L_scale);
    params.L_scale = L_scale;
}

void PIC::Config::set_T_scale(double T_scale)
{
    check_if_positive("T_scale must be > 0.0.", T_scale);
    params.T_scale = T_scale;
}

void PIC::Config::set_U_scale(double U_scale)
{
    check_if_positive("U_scale must be > 0.0.", U_scale);
    params.U_scale = U_scale;
}

void PIC::Config::set_N_scale(double N_scale)
{
    check_if_positive("N_scale must be > 0.0.", N_scale);
    params.N_scale = N_scale;
}

void PIC::Config::set_E_scale(double E_scale)
{
    check_if_positive("E_scale must be > 0.0.", E_scale);
    params.E_scale = E_scale;
}

void PIC::Config::set_B_scale(double B_scale)
{
    check_if_positive("B_scale must be > 0.0.", B_scale);
    params.B_scale = B_scale;
}

// static
void PIC::Config::to_stream(std::ostream& os)
{
    os << params;
}

PIC::Config::Parameters::Parameters()
: cfg_script_name(),
  h(0.0),
  tau(0.0),
  dens_cutoff(0.0),
  time_steps(0), save_time_steps(0),
  grid_size_x(0), grid_size_y(0), grid_size_z(0),
  total_particles_num(0),
  save_all_particles(false),
  save_whole_grid(false),
  save_grid_x_plains(false),
  save_grid_y_plains(false),
  save_grid_z_plains(false),
  os_name(""),
  process_num(1),
  process_rank(0),
  current_time_step(0),
  L_scale(1.0),
  T_scale(1.0),
  U_scale(1.0),
  N_scale(1.0),
  E_scale(1.0),
  B_scale(1.0),
  particle_push_alg(Direct),
  scatter_alg(Standard),
  grid_threshold(Min_Density),
  CFL_severity(Stop)
{
    set_os_name();
}

PIC::Config::Parameters::Parameters(const std::string& cfg_script_name_,
                                    double h_,
                                    double tau_,
                                    double dens_cutoff_,
                                    index_t time_steps_, index_t save_time_steps_,
                                    index_t grid_size_x_, index_t grid_size_y_, index_t grid_size_z_,
                                    index_t total_particles_num_,
                                    bool save_all_particles_,
                                    bool save_whole_grid_,
                                    bool save_grid_x_plains_, bool save_grid_y_plains_, bool save_grid_z_plains_,
                                    const std::string& os_name_,
                                    index_t process_num_,
                                    index_t process_rank_,
                                    index_t current_time_step_,
                                    double L_scale_,
                                    double T_scale_,
                                    double U_scale_,
                                    double N_scale_,
                                    double E_scale_,
                                    double B_scale_,
                                    ParticlePushAlg particle_push_alg_,
                                    ScatterAlg scatter_alg_,
                                    GridThreshold grid_threshold_,
                                    CFLSeverity CFL_severity_)
: cfg_script_name(cfg_script_name_),
  h(h_),
  tau(tau_),
  dens_cutoff(dens_cutoff_),
  time_steps(time_steps_), save_time_steps(save_time_steps_),
  grid_size_x(grid_size_x_), grid_size_y(grid_size_y_), grid_size_z(grid_size_z_),
  total_particles_num(total_particles_num_),
  save_all_particles(save_all_particles_),
  save_whole_grid(save_whole_grid_),
  save_grid_x_plains(save_grid_x_plains_),
  save_grid_y_plains(save_grid_y_plains_),
  save_grid_z_plains(save_grid_z_plains_),
  os_name(os_name_),
  process_num(process_num_),
  process_rank(process_rank_),
  current_time_step(current_time_step_),
  L_scale(L_scale_),
  T_scale(T_scale_),
  U_scale(U_scale_),
  N_scale(N_scale_),
  E_scale(E_scale_),
  B_scale(B_scale_),
  particle_push_alg(particle_push_alg_),
  scatter_alg(scatter_alg_),
  grid_threshold(grid_threshold_),
  CFL_severity(CFL_severity_)
{

}

bool PIC::Config::Parameters::is_valid() const
{
    const bool is_density_threshold_valid = (grid_threshold == Min_Density ) ? (dens_cutoff > 0.0) : true;

    return (h > 0.0 && tau > 0.0 && is_density_threshold_valid &&
            time_steps > 0 && save_time_steps > 0 &&
            grid_size_x > 0 && grid_size_y > 0 && grid_size_z > 0 &&
            L_scale > 0.0 &&
            T_scale > 0.0 &&
            U_scale > 0.0 &&
            N_scale > 0.0 &&
            E_scale > 0.0 &&
            B_scale > 0.0);
}

bool PIC::Config::is_on_save_step()
{
    return ((current_time_step() % save_time_steps() == 0) ||
            (current_time_step() == time_steps()));
}

void PIC::Config::set_os_name()
{
#ifdef _WIN64
    params.os_name = "Windows 64-bit";
#elif _WIN32
    params.os_name = "Windows 32-bit";
#elif __unix || __unix__
    params.os_name = "Unix";
#elif __APPLE__ || __MACH__
    params.os_name = "Mac OSX";
#elif __linux__
    params.os_name = "Linux";
#elif __FreeBSD__
    params.os_name = "FreeBSD";
#else
    params.os_name = "Unknown";
#endif
}

namespace PIC {

ostream& operator<<(ostream& os, const Config::Parameters& params)
{
    os << "\nConfiguaration:"
       << "\ncfg_script_name = " << params.cfg_script_name
       << "\nh = " << params.h
       << "\ntau = " << params.tau
       << "\ndens_cutoff = " << params.dens_cutoff
       << "\nL_scale = " << params.L_scale
       << "\nT_scale = " << params.T_scale
       << "\nU_scale = " << params.U_scale
       << "\nN_scale = " << params.N_scale
       << "\nE_scale = " << params.E_scale
       << "\nB_scale = " << params.B_scale
       << "\ntime steps number = " << params.time_steps
       << "\nsave time steps number = " << params.save_time_steps
       << "\ngrid size = " << params.grid_size_x << " x " << params.grid_size_y << " x " << params.grid_size_z
       << "\ntotal particles number = " << params.total_particles_num
       << boolalpha
       << "\nsave all particles = " << params.save_all_particles
       << "\nsave whole grid = " << params.save_whole_grid
       << "\nsave grid X-plains = " << params.save_grid_x_plains
       << "\nsave grid Y-plains = " << params.save_grid_y_plains
       << "\nsave grid Z-plains = " << params.save_grid_z_plains
       << "\nOS name = " << params.os_name
       << "\nprocess_num = " << params.process_num
       << "\nprocess_rank = " << params.process_rank
       << "\nalgorithm of solving particle motion equation = "
       << (params.particle_push_alg == Direct ? "Direct" :
           params.particle_push_alg == Boris ? "Boris" : "Undefined")
       << "\nalgorithm of scatter particles to grid = "
       << (params.scatter_alg == Standard ? "Standard" :
           params.scatter_alg == Zigzag ? "Zigzag" : "Undefined")
       << "\ngrid threshold type = "
       << (params.grid_threshold == Min_Density ? "Min_Density" :
           params.grid_threshold == Local_CFL ? "Local_CFL" : "Undefined")
       << "\nCFL_severity = " << (params.CFL_severity == Ignore ? "Ignore" :
                                  params.CFL_severity == Absorb ? "Absorb" :
                                  params.CFL_severity == Stop ? "Stop" : "Undefined")
       << "\nConstants:"
       << "\nc = " << Constants::c()
       << "\ne = " << Constants::e()
       << "\nmp = " << Constants::mp()
       << "\ne_mp = " << Constants::e_mp()
       << "\npi = " << Constants::pi()
       << "\nc_4pi_e = " << Constants::c_4pi_e()
       << endl;

    return os;
}

} // namespace PIC
