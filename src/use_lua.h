#if !defined (USE_LUA_H)
#define USE_LUA_H

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

#endif // USE_LUA_H
