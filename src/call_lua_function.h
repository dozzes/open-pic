#if !defined (CALL_LUA_FUNCTION_H)
#define CALL_LUA_FUNCTION_H


bool call_lua_function(const char* func_name);

struct Particle;
bool lua_validate_particle(Particle& particle);

#endif // CALL_LUA_FUNCTION_H
