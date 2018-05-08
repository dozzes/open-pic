#pragma once

#include "opic_fwd.h"
#include "constants.h"

#include <string>
#include <iosfwd>

namespace PIC {

using namespace std;

enum ParticlePushAlg { Direct = 0, Boris};
enum ScatterAlg { Standard = 0, Zigzag };
enum GridThreshold { Min_Density = 0, Local_CFL};
enum CFLSeverity { Ignore, Absorb, Stop};

class Config
{
public:
    struct Parameters
    {
        Parameters();

        Parameters(const string& cfg_script_name_,
                   double h_,
                   double tau_,
                   double dens_cutoff_,
                   index_t time_steps_,
                   index_t save_time_steps_,
                   index_t grid_size_x_,
                   index_t grid_size_y_,
                   index_t grid_size_z_,
                   index_t total_particles_num_,
                   bool   save_all_particles_,
                   bool   save_whole_grid_,
                   bool   save_grid_x_plains,
                   bool   save_grid_y_plains,
                   bool   save_grid_z_plains,
                   const string& os_name_,
                   index_t process_num_,
                   index_t process_rank_,
                   index_t current_time_step_,
                   double length_scale_,
                   double time_scale_,
                   double U_scale_,
                   double N_scale_,
                   double E_scale_,
                   double B_scale_,
                   ParticlePushAlg particl_push_alg_,
                   ScatterAlg scatter_alg_,
                   GridThreshold grid_threshold_,
                   CFLSeverity CFL_severity_);

        bool is_valid() const;

        string cfg_script_name;
        double h;               // spatial step
        double tau;             // time step
        double dens_cutoff;          // minimum possible value of grid NP (density cutoff value)
        index_t time_steps ;     // time steps number
        index_t save_time_steps; // save time steps number

        index_t grid_size_x;
        index_t grid_size_y;
        index_t grid_size_z;

        index_t total_particles_num;

        bool save_all_particles;
        bool save_whole_grid;
        bool save_grid_x_plains;
        bool save_grid_y_plains;
        bool save_grid_z_plains;

        string os_name;
        index_t process_num;
        index_t process_rank;

        index_t current_time_step;

        double L_scale;
        double T_scale;
        double U_scale;
        double N_scale;
        double E_scale;
        double B_scale;

        ParticlePushAlg particle_push_alg;
        ScatterAlg scatter_alg;
        GridThreshold grid_threshold;
        CFLSeverity CFL_severity;
    };

    static const string& cfg_script_name() { return params.cfg_script_name; }

    static void set_tau(double tau);
    static double tau() { return params.tau; }
    static double h()   { return params.h; }
    static double h_2() { return 0.5*h(); }
    static void set_dens_cutoff(double dens_cutoff); // min density (NP)
    static double dens_cutoff() { return params.dens_cutoff; }

    static double tau_2()     { return (0.5*tau()); }
    static double tau_2h()    { return (tau_2()/h()); }
    static double ctau_2()    { return (0.5*Constants::c()*tau()); }
    static double ctau_h()    { return (Constants::c()*tau()/h()); }
    static double ctau_2h()   { return (0.5*ctau_h());}
    static double c_4pi_e_h() { return (Constants::c_4pi_e()/h()); }

    static index_t grid_size_x() { return params.grid_size_x; }
    static index_t grid_size_y() { return params.grid_size_y; }
    static index_t grid_size_z() { return params.grid_size_z; }

    static index_t total_particles_num() { return params.total_particles_num; }

    static bool save_all_particles() { return params.save_all_particles; }
    static bool save_whole_grid()    { return params.save_whole_grid; }
    static bool save_grid_x_plains() { return params.save_grid_x_plains; }
    static bool save_grid_y_plains() { return params.save_grid_y_plains; }
    static bool save_grid_z_plains() { return params.save_grid_z_plains; }

    static string os_name()       { return params.os_name; }
    static index_t process_num()  { return params.process_num; }
    static index_t process_rank() { return params.process_rank; }

    static index_t time_steps()                          { return params.time_steps; } // time steps number
    static index_t save_time_steps()                     { return params.save_time_steps; } // save time steps number
    static index_t current_time_step()                   { return params.current_time_step; }
    static void   set_current_time_step(index_t time_step) { params.current_time_step = time_step; }
    static bool   is_on_save_step();
    static void   set_os_name();

    static double L_scale() { return params.L_scale; }
    static double T_scale() { return params.T_scale; }
    static double U_scale() { return params.U_scale; }
    static double N_scale() { return params.N_scale; }
    static double E_scale() { return params.E_scale; }
    static double B_scale() { return params.B_scale; }

    static void set_L_scale(double L_scale);
    static void set_T_scale(double T_scale);
    static void set_U_scale(double U_scale);
    static void set_N_scale(double N_scale);
    static void set_E_scale(double E_scale);
    static void set_B_scale(double B_scale);

    static void set_particle_push_alg(ParticlePushAlg alg) { params.particle_push_alg = alg; }
    static ParticlePushAlg particle_push_alg() { return params.particle_push_alg; }

    static void set_scatter_alg(ScatterAlg alg) { params.scatter_alg = alg; }
    static ScatterAlg scatter_alg() { return params.scatter_alg; }

    static void set_grid_threshold(GridThreshold gt) { params.grid_threshold = gt; }
    static GridThreshold grid_threshold() { return params.grid_threshold; }

    static void set_CFL_severity (CFLSeverity CFL_severity) { params.CFL_severity = CFL_severity; }
    static CFLSeverity CFL_severity() { return params.CFL_severity; }

    static Parameters& parameters() { return params; }

    static void to_stream(ostream& os);

    static ofstream& logger() { return ofs_log; }

private:
    static Parameters params;
    static ofstream& ofs_log;
};

ostream& operator<<(ostream& out, const Config::Parameters& params);

} // namespace PIC
