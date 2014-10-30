#if !defined (PARTICLE_GROUPS_H)
#define PARTICLE_GROUPS_H

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include "opic_fwd.h"
#include "diagnostics.h"


class ParticleGroups
{
public:
   struct ParticleGroup
   {
      ParticleGroup(const std::string &name_,
                    double charge_, double mass_,
                    PIC::Diagnostics diag_);

      std::string name;      // group name
      double charge;         // charge of one ion in particle (in electron charge units)
      double mass;           // particle mass in ion mass unit
      PIC::Diagnostics diag; // group diagnostics type
   };

   typedef std::map<std::string, ParticleGroup> MapNameToGroup;

   bool create_group(const std::string& name,
                     double charge,
                     double mass,
                     PIC::Diagnostics diag);

   void remove_index(index_t idx_to_remove);
   void clear() { groups_.clear(); }

   // get group by name
   const ParticleGroup& operator[](const std::string& name) const;

   MapNameToGroup::const_iterator begin() const { return groups_.begin(); }
   MapNameToGroup::const_iterator end() const { return groups_.end(); }

   MapNameToGroup::iterator begin() { return groups_.begin(); }
   MapNameToGroup::iterator end() { return groups_.end(); }

   const MapNameToGroup& groups() const { return groups_; }

   const std::vector<std::string>& group_names() const;

   static double max_mp() { return max_mp_; }

   void to_stream(std::ostream& os) const;

private:
   static double max_mp_;
   static MapNameToGroup groups_;
};

std::ostream& operator<<(std::ostream& out, const ParticleGroups& pg);

#endif // PARTICLE_GROUPS_H
