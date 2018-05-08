#pragma once

#include "opic_fwd.h"

/*********************************************************
* Check if particle is active.                           *
* equation of motion is solved for active particles only *
* return true if particles is active, otherwise - false  *
*********************************************************/
bool is_particle_can_scatter(Particles& particles, index_t p, const Grid& grid, const DblVector& dr);
bool is_particle_can_move(Particles& particles, index_t p, const Grid& grid);

// Check if particle will be move less thne h/2.
bool check_particle_move(const Particle& particle, const Grid& grid, const DblVector& dr);
