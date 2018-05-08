#pragma once

#include "opic_fwd.h"
#include "vector_3d.h"

#include <iosfwd>
#include <string>
#include <deque>
#include <functional>
#include <set>
#include <vector>

struct Particle
{
    Particle()
    : r(0.0, 0.0, 0.0),
      v(0.0, 0.0, 0.0),
      ni(0.0)
    {
    }

    Particle(const std::string& group_name_,
             double x_, double y_, double z_,
             double v_x_, double v_y_, double v_z_,
             double ni_)
    : r(x_, y_, z_),
      v(v_x_, v_y_, v_z_),
      ni(ni_),
      group_name(group_name_)
    {
    }

    Particle(const std::string& group_name_,
             const DblVector& vec_r, const DblVector& vec_v,
             double ni_)
    : r(vec_r),
      v(vec_v),
      ni(ni_),
      group_name(group_name_)
    {
    }

    Particle(const Particle& other)
    : r(other.r),
      v(other.v),
      ni(other.ni),
      group_name(other.group_name)
    {
    }

    Particle& operator=(const Particle& rhs)
    {
        if ((void*)this == (void*)&rhs)
        {
            return *this;
        }

        r = rhs.r;
        v = rhs.v;
        ni = rhs.ni;
        group_name = rhs.group_name;

        return *this;
    }

    void to_string(std::string& s) const;

    DblVector r;              // particle position
    DblVector v;              // particle velocity
    double    ni;             // ions in simulation particle

    std::string group_name;   // particle group name
};

std::ostream& operator<<(std::ostream& out, const Particle& p);

class Particles
{
    typedef std::deque<Particle> ParticlesContainer;

public:
    Particles();

    index_t size() const;
    void resize(index_t new_size);

    const Particle& operator[](index_t p) const;
    Particle& operator[](index_t p);

    const Particle& at(index_t p) const;
    Particle& at(index_t p);

    size_t remove_inactives();
    void remove_later(index_t p);
    bool is_inactive(index_t p) const;

    // for using from Lua
    void set(index_t p, const Particle& part);

    void set(index_t p,
             const std::string& group_name,
             const DblVector& vec_r, const DblVector& vec_v,
             double ni);

    void add(const Particle& part);

    void add(const std::string& group_name,
             const DblVector& vec_r,
             const DblVector& vec_v,
             double ni);

    void erase(index_t p);

    ParticlesContainer::const_iterator begin() const { return data_.begin(); }
    ParticlesContainer::const_iterator end() const { return data_.end(); }

    ParticlesContainer::iterator begin() { return data_.begin(); }
    ParticlesContainer::iterator end() { return data_.end(); }

    void init_ball(index_t from_idx, index_t parts_num,
                   const DblVector& center, double R,
                   const std::string& grp_name,
                   double v_min, double v_max, double ni);
private:
    ParticlesContainer data_;
    typedef std::set<index_t, std::greater<index_t> > IndexSet;
    std::vector<IndexSet> indexes_to_remove_;
};
