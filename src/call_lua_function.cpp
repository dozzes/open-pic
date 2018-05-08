#include "call_lua_function.h"
#include "config.h"
#include "use_lua.h"
#include "io_utilities.h"
#include "particles.h"

#include <boost/format.hpp>

bool call_lua_function(const char* func_name)
{
    UseLua lua;

    bool ok = true;

    try
    {
        luabind::call_function<void>(lua, func_name);
    }
    catch(std::exception& e)
    {
        print(e.what());
        ok = false;
    }    
    catch(...)
    {
        std::string msg_fmt = "Function \"%s\" not defined or invalid in \"%s\".";
        std::string msg = str(boost::format(msg_fmt) % func_name % PIC::Config::cfg_script_name());
        print(msg);
        ok = false;
    }

    return ok;
}

bool lua_validate_particle(const Particle& particle)
{
    UseLua lua;

    bool ok = true;

    try
    {
        ok = luabind::call_function<bool>(lua, "validate_particle", particle);
    }
    catch (std::exception& /*e*/)
    {
        std::string msg_fmt = "Function <validate_particle()> not defined or invalid in \"%s\".";
        std::string msg = str(boost::format(msg_fmt) % PIC::Config::cfg_script_name());
        print(msg);
        ok = false;
    }

    return ok;
}
