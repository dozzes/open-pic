#include <algorithm>
#include <ostream>
#include <boost/format.hpp>

#include "particle_groups.h"
#include "marker_particles.h"
#include "check_particle.h"
#include "utility.h"
#include "io_utilities.h"

#include "particles.h"


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

void Particles::absorbe(index_t p)
{
    data_[p].is_absorbed = true;
    indexes_to_remove_.insert(p);
}

index_t Particles::remove_absorbed()
{
    const index_t absorbed_num = indexes_to_remove_.size();

    if (absorbed_num < data_.size() / 10)
    {
        return 0;
    }

    print_tm("Remove absorbed particles:");
    std::cout << absorbed_num << std::endl;

    IndexSet::const_iterator it = indexes_to_remove_.begin();
    for (; it != indexes_to_remove_.end(); ++it)
    {
        erase(*it);
    }

    indexes_to_remove_.clear();

    return absorbed_num;
}

void Particles::init_ball(index_t from_idx, index_t parts_num,
                          const DblVector& center, double R,
                          const std::string& group_name,
                          double v_min, double v_max, double ni)
{
    const index_t curr_size = data_.size();
    const index_t new_size = from_idx + parts_num;

    if (new_size > curr_size)
    {
        data_.resize(new_size);
    }

    for (index_t p = from_idx; p != new_size; ++p)
    {
        DblVector r;
        rnd_ball(R, r);

        Particle& part = data_[p];

        part.r.x = center.x + r.x;
        part.r.y = center.y + r.y;
        part.r.z = center.z + r.z;

        part.v.x = v_min + v_max*r.x / R;
        part.v.y = v_min + v_max*r.y / R;
        part.v.z = v_min + v_max*r.z / R;

        part.ni = ni;
        part.group_name = group_name;
        part.is_absorbed = false;
    }

    std::cout << "\nParticles::init_ball(..): cloud particles index range = ["
              << from_idx << "," << new_size << ")"
              << std::endl;
}

void Particle::to_string(std::string& s) const
{
   s = str(boost::format("x = %1%"
                         "\ny = %2%"
                         "\nz = %3%"
                         "\nv_x = %4%"
                         "\nv_y = %5%"
                         "\nv_z = %6%")
           % r.x % r.y % r.z % v.x % v.y % v.z);
}

std::ostream& operator<<(std::ostream& out, const Particle& particle)
{
    return (out << particle.r.x << "\t"
                << particle.r.y << "\t"
                << particle.r.z << "\t"
                << particle.v.x << "\t"
                << particle.v.y << "\t"
                << particle.v.z << "\t"
                << particle.ni << std::endl);
}
