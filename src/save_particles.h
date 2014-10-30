#if !defined (SAVE_PARTICLES_H)
#define SAVE_PARTICLES_H


struct Cell;

typedef GridContainer<Cell> Grid;

class Particles;

void save_particles(const Grid& grid, const Particles& particles);

#endif // SAVE_PARTICLES_H
