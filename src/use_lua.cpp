#include <stdexcept>

#include "use_lua.h"

UseLua::LuaState UseLua::luaState_ = UseLua::LuaState();

UseLua::LuaState::LuaState() : lua_(nullptr)
{
    lua_ = luaL_newstate();
    if (!lua_)
        throw std::runtime_error("Can't open Lua.");

    luabind::open(lua_); // connect luabind to this Lua state
}

UseLua::LuaState::~LuaState()
{
    lua_close(lua_);
}
