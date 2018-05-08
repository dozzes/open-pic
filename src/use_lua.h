#pragma once

#include "inc_lua.h"

struct UseLua
{
    operator lua_State*() { return luaState_.lua_; }

private:
    struct LuaState
    {
        LuaState();
        ~LuaState();
        lua_State* lua_;
    };

    static LuaState luaState_;
};
