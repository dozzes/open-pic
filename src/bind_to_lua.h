#pragma once

#include "opic_fwd.h"

bool bind_to_lua(const char* lua_cfg_file_name, Grid& grid,
                 Particles& particles);
