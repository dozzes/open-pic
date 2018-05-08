#pragma once

bool call_lua_function(const char* func_name);

struct Particle;

bool lua_validate_particle(const Particle& particle);
