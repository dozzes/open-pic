#if !defined (CHECK_PARTICLE_H)
#define CHECK_PARTICLE_H

#include "opic_fwd.h"


/**********************************************************
* Validate particle.                                      *
* If particle is not active set it absorbed.              *
* Equation of motion is solved for active particles only. *
* Return true if particles is active, otherwise - false   *
*********************************************************/
bool validate_particle(Particle& particle, const Grid& grid);

bool check_particle_move(Particle& particle, const Grid& grid, const DblVector& dr);

#endif // CHECK_PARTICLE_H
