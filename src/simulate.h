#pragma once

#include "opic_fwd.h"

#include <boost/function.hpp>

namespace PIC {

/************************************************************
*                                                           *
* I step                                                    *
*                                                           *
* Calculate magnetic field in half time step  (m+1/2)       *
*                                                           *
************************************************************/
void calc_magnetic_field_half_time(Grid& grid);

/*******************************************
*                                          *
* II step                                  *
*                                          *
* Move particles in half time step (m+1/2) *
*                                          *
*******************************************/
void move_particles_half_time(const Grid& grid, Particles& particles,
                              const std::string& group_name);

/************************************************************
*                                                           *
* III step                                                  *
*                                                           *
* Set average particle velocities and density in grid       *
* See function from_particles_to_grid() in gather_scatter.h *
*                                                           *
************************************************************/

/***********************************
*                                  *
* IV step                          *
*                                  *
* Move particles in full time step *
*                                  *
***********************************/
void move_particles_full_time(const Grid& grid, Particles& particles,
                              const std::string& group_name);

/***************************************************
*                                                  *
* V step (repeat step I)                           *
*                                                  *
* Calculate magnetic field in full time step (m+1) *
*                                                  *
*                                                  *
***************************************************/

/***************************************************
*                                                  *
* VI step                                          *
*                                                  *
* Calculate electrons velocity                     *
*                                                  *
***************************************************/
void calc_electrons_velocity(Grid& grid);

/***************************************************
*                                                  *
* VI step                                          *
*                                                  *
* Calculate electric field                         *
*                                                  *
***************************************************/
void calc_electric_field(Grid& grid);

/**************************************************
* Simulation cycle                                *
**************************************************/
void simulate(Grid& grid, Particles& particles);

} // namespace PIC
