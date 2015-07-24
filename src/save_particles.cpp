#include <cstdio>

#include "grid.h"
#include "particles.h"
#include "gather_scatter.h"
#include "particle_groups.h"
#include "marker_particles.h"
#include "gather_scatter.h"
#include "config.h"
#include "io_utilities.h"

#include "save_particles.h"


using namespace std;

namespace
{
    void save_particle(const Grid& grid, const Particle& particle, FILE* fout)
    {
        if (particle.is_absorbed)
        {
            return;
        }

        Grid::NodeType point_val;
        PIC::from_grid_to_point(grid, particle.r, point_val);

        const double L = PIC::Config::L_scale();
        const double U = PIC::Config::U_scale();
        const double E = PIC::Config::E_scale();
        const double B = PIC::Config::B_scale();

        fprintf(fout,
                "%e\t%e\t%e\t"   // (particle.r.x/L), particle.r.y/L), (particle.r.z/L)
                "%e\t"           // (particle.v.abs()/U)
                "%e\t%e\t%e\t"   // (particle.v.x/U), particle.v.y/U), (particle.r.z/U)
                "%e\t"           // (point_val.B.abs()/B)
                "%e\t%e\t%e\t"   // (point_val.B.x/B), (point_val.B.y/B), (point_val.B.z/B)
                "%e\t"           // (point_val.E.abs()/E)
                "%e\t%e\t%e\n",  // (point_val.B.x/B), (point_val.B.y/B), (point_val.B.z/B)
                (particle.r.x / L), (particle.r.y / L), (particle.r.z / L),
                (particle.v.abs() / U),
                (particle.v.x / U), (particle.v.y / U), (particle.v.z / U),
                (point_val.B.abs() / B),
                (point_val.B.x / B), (point_val.B.y / B), (point_val.B.z / B),
                (point_val.E.abs() / E),
                (point_val.E.x / E), (point_val.E.y / E), (point_val.E.z / E));
    }

    void save_particles_group(const Grid& grid, const Particles& particles, const string& group_name)
    {
        const string file_name = create_out_file_name(group_name, "parts", PIC::Config::current_time_step());

        FILE* fout = fopen(file_name.c_str(), "w");

        if (!fout)
        {
            cerr << "\nCan't save data: " << file_name << endl;
            return;
        }

        for (index_t p = 0; p != particles.size(); ++p)
        {
            const Particle& particle = particles[p];

            if ((group_name == "all" || (particle.group_name == group_name)) &&
                !particle.is_absorbed)
            {
                save_particle(grid, particle, fout);
            }
        }

        fclose(fout);
    }
} //  unnamed namespace

void save_particles(const Grid& grid, const Particles& particles)
{
    ParticleGroups part_groups;

    const vector<string>& group_names = part_groups.group_names();

    for (index_t i = 0; i != group_names.size(); ++i)
    {
        const string& group_name = group_names[i];
        const ParticleGroups::ParticleGroup& group = part_groups[group_name];

        if (group.diag & PIC::save_positions)
        {
            save_particles_group(grid, particles, group_name);
        }
    }

    MarkerParticles markers;
    MarkerParticles::MapNameToMarkers::const_iterator it = markers.begin();

    for (; it != markers.end(); ++it)
    {
        const deque<index_t>& vec_idx = (*it).second;

        if (vec_idx.empty()) continue;

        const string file_name = create_out_file_name((*it).first, "markers", PIC::Config::current_time_step());

        FILE* fout = fopen(file_name.c_str(), "w");

        if (!fout)
        {
            cerr << "\nCan't save data: " << file_name << endl;
        }
        else
        {
            const char* the_header = "X\tY\tZ"
                                     "\tV\tVx\tVy\tVz"
                                     "\tB\tBx\tBy\tBz"
                                     "\tE\tEx\tEy\tEz\n";
            fputs(the_header, fout);

            for (index_t i = 0; i != vec_idx.size(); ++i)
            {
                save_particle(grid, particles[vec_idx[i]], fout);
            }
        }

        fclose(fout);
    }

    if (PIC::Config::save_all_particles())
    {
        save_particles_group(grid, particles, "all");
    }
}

