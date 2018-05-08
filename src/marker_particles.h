#pragma once

#include "opic_fwd.h"
#include "diagnostics.h"

#include <string>
#include <map>
#include <deque>

class MarkerParticles
{
public:
    typedef std::map<std::string, std::deque<index_t> > MapNameToMarkers;

    void insert_idx(const std::string& name, index_t particle_idx);
    void remove_group(const std::string& name);
    void remove_idx(index_t idx);
    void clear() { markers_ = MapNameToMarkers(); }
    void set_diagnostics(PIC::Diagnostics diag) { diag_ = diag; }
    PIC::Diagnostics diagnostics() { return diag_; }

    // get group by name
    const std::deque<index_t>& operator[](const std::string& name) const;

    MapNameToMarkers::const_iterator begin() const { return markers_.begin(); }
    MapNameToMarkers::const_iterator end() const { return markers_.end(); }

    MapNameToMarkers::iterator begin() { return markers_.begin(); }
    MapNameToMarkers::iterator end() { return markers_.end(); }

    const MapNameToMarkers& markers() const { return markers_; }

private:
    static MapNameToMarkers markers_;
    static PIC::Diagnostics diag_;
};
