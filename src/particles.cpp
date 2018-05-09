#include "particles.h"
#include "particle_groups.h"
#include "marker_particles.h"
#include "check_particle.h"
#include "utility.h"
#include "io_utilities.h"

#include <boost/format.hpp>
#include <algorithm>
#include <ostream>
#include <omp.h>

Particles::Particles()
: indexes_to_remove_(omp_get_max_threads())
{
}

index_t Particles::size() const
{
    return data_.size();
}

void Particles::resize(index_t new_size)
{
    data_.resize(new_size);
}

const Particle& Particles::operator[](index_t p) const
{
    return data_[p];
}

Particle& Particles::operator[](index_t p)
{
    return data_.at(p);
}

const Particle& Particles::at(index_t p) const
{
    return data_.at(p);
}

Particle& Particles::at(index_t p)
{
    return data_.at(p);
}

void Particles::set(index_t i, const Particle& part)
{
    data_.at(i) = part;
}

void Particles::set(index_t p,
                    const std::string& group_name,
                    const DblVector& vec_r,
                    const DblVector& vec_v,
                    double ni)
{
    Particle& particle = data_.at(p);

    particle.r = vec_r;
    particle.v = vec_v;
    particle.group_name = group_name;
    particle.ni = ni;
}

void Particles::add(const Particle& part)
{
    data_.push_back(part);
}

void Particles::add(const std::string& group_name,
                    const DblVector& vec_r,
                    const DblVector& vec_v,
                    double ni)
{
    data_.push_back(Particle(group_name, vec_r, vec_v, ni));
}

void Particles::erase(index_t p)
{
    data_.erase(data_.begin() + p);

    MarkerParticles markers;
    markers.remove_idx(p);
}

void Particles::remove_later(index_t p)
{
    const int tid = omp_get_thread_num();
    indexes_to_remove_[tid].insert(p);
}

bool Particles::is_inactive(index_t p) const
{
    for (size_t t = 0; t != indexes_to_remove_.size(); ++t)
    {
        const IndexSet& indexes_to_remove = indexes_to_remove_[t];
        if (indexes_to_remove.find(p) != indexes_to_remove.end())
            return true;
    }
    return false;
}

size_t Particles::remove_inactives()
{
    return 0;
    IndexSet& indexes_to_remove = indexes_to_remove_[0];
    for (size_t t = 1; t != indexes_to_remove_.size(); ++t)
    {
        IndexSet& tid_indexes = indexes_to_remove_[t];
        indexes_to_remove.insert(tid_indexes.begin(), tid_indexes.end());
        tid_indexes.clear();
    }

    const size_t remove_count = indexes_to_remove.size();

    if (remove_count > 0)
    {
        print_tm("Removed inactive particles: ", remove_count);

        IndexSet::const_iterator i = indexes_to_remove.begin();
        for (; i != indexes_to_remove.end(); ++i)
            this->erase(*i);

        indexes_to_remove.clear();
    }

    return remove_count;
}

void Particles::init_ball(index_t from_idx, index_t parts_num,
                          const DblVector& center, double R,
                          const std::string& group_name,
                          double v_min, double v_max, double ni)
{
    const index_t curr_size = data_.size();
    const index_t new_size = from_idx + parts_num;

    if (new_size > curr_size)
        data_.resize(new_size);

    const DblVector vec_v_min(v_min, v_min, v_min);
    for (index_t p = from_idx; p != new_size; ++p)
    {
        DblVector r;
        rnd_ball(R, r);

        Particle& part = data_[p];
        part.r = center + r;
        part.v = vec_v_min + v_max*r/R;
        part.ni = ni;
        part.group_name = group_name;
    }

    std::cout << "\nParticles::init_ball(..): cloud particles index range = ["
              << from_idx << "," << new_size << ")"
              << std::endl;
}

void Particle::to_string(std::string& s) const
{
    s = str(boost::format("x = %1%\ny = %2%\nz = %3%"
                          "\nv_x = %4%\nv_y = %5%\nv_z = %6%")
                           % r.x % r.y % r.z
                           % v.x % v.y % v.z);
}

std::ostream& operator<<(std::ostream& out, const Particle& particle)
{
    return (out << particle.r.x << "\t" << particle.r.y << "\t" << particle.r.z << "\t"
                << particle.v.x << "\t" << particle.v.y << "\t" << particle.v.z << "\t"
                << particle.ni << std::endl);
}
