#pragma once

#include "grid.h"

struct Cell;
typedef GridContainer<Cell> Grid;

class Particles;
void save_particles(const Grid& grid, const Particles& particles);
