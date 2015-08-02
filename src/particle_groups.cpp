#include <iostream>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <limits> 

#include <boost/format.hpp>

#include "opic_fwd.h"

#include "particle_groups.h"


using namespace PIC;
using namespace std;

// static
const char* ParticleGroups::all_particles_name = "all";

// static
ParticleGroups::MapNameToGroup ParticleGroups::groups_ = ParticleGroups::MapNameToGroup();

// static
double ParticleGroups::max_mp_ = std::numeric_limits<double>::min();

ParticleGroups::ParticleGroup::ParticleGroup(const string &name_, double charge_, double mass_,
                                             Diagnostics diag_)
: name(name_), charge(charge_), mass(mass_), diag(diag_)
{

}

bool ParticleGroups::create_group(const string& name, double charge, double mass,
                                  Diagnostics diag)
{
    if (groups_.find(name) == end())
    {
        ParticleGroup pg(name, charge, mass, diag);

        if (mass > max_mp_)
        {
            max_mp_ = mass;
        }

        return groups_.insert(make_pair(name, pg)).second;
    }

    return true;
}

void ParticleGroups::remove_index(index_t /*idx_to_remove*/)
{

}

const vector<string>&
ParticleGroups::group_names() const
{
    static vector<string> names;
    names.clear();

    MapNameToGroup::const_iterator it = begin();
    for (; it != end(); ++it)
    {
        names.push_back(it->first);
    }

    return names;
}

const ParticleGroups::ParticleGroup&
ParticleGroups::operator[](const string& name) const
{
    MapNameToGroup::const_iterator it = groups_.find(name);
    if (it == end())
    {
        string msg = str(boost::format("ERROR : Invalid name [%s] "
            "in ParticleGroups::operator[](const string& name)") % name);
        throw invalid_argument(msg);
    }

    return (*it).second;
}

void ParticleGroups::to_stream(std::ostream& os) const
{
    MapNameToGroup::const_iterator it = begin();
    for (; it != end(); ++it)
    {
        const ParticleGroup& group = (*it).second;

        os << "\n" << group.name << ":"
            << "\n\t" << "charge = " << group.charge
            << "\n\t" << "mass = " << group.mass
            << "\n\t" << "diag = " << group.diag
            << std::endl;
    }
}

std::ostream& operator<<(std::ostream& os, const ParticleGroups& pg)
{
    pg.to_stream(os);
    return os;
}
