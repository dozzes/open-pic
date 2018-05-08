#include "marker_particles.h"
#include "opic_fwd.h"

#include <algorithm>

using namespace PIC;
using namespace std;

// static
MarkerParticles::MapNameToMarkers MarkerParticles::markers_ = MarkerParticles::MapNameToMarkers();
Diagnostics MarkerParticles::diag_ = save_positions;

void MarkerParticles::insert_idx(const string& name, index_t particle_idx)
{
    deque<index_t>& seq_idx = markers_[name];
    if (find(seq_idx.begin(), seq_idx.end(), particle_idx) == seq_idx.end())
        seq_idx.push_back(particle_idx);
}

void MarkerParticles::remove_idx(index_t idx)
{
    MapNameToMarkers::iterator it = begin();
    for (; it != end(); ++it)
    {
        deque<index_t>& seq_idx = (*it).second;

        seq_idx.erase(remove(seq_idx.begin(), seq_idx.end(), idx),
                      seq_idx.end());

        for (index_t i = 0; i != seq_idx.size(); ++i)
        {
            if (seq_idx[i] > idx)
                --seq_idx[i];
        }
    }
}

void MarkerParticles::remove_group(const string& name)
{
    markers_.erase(name);
}

const deque<index_t>&
MarkerParticles::operator[](const string& name) const
{
    static deque<index_t> empty;
    MapNameToMarkers::const_iterator it = markers_.find(name);
    return (it != end() ? (*it).second : empty);
}
