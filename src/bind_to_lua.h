#if !defined (BIND_TO_LUA_H)
#define BIND_TO_LUA_H


#include "opic_fwd.h"

bool bind_to_lua(const char* lua_cfg_file_name,
                 Grid& grid,
                 Particles& particles);

#endif // BIND_TO_LUA_H
