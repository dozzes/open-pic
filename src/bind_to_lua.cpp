#include "bind_to_lua.h"
#include "config.h"
#include "use_lua.h"
#include "gather_scatter.h"
#include "grid.h"
#include "grid_filters.h"
#include "particles.h"
#include "particle_groups.h"
#include "marker_particles.h"
#include "save_grid.h"

#include <luabind/operator.hpp>
#include <boost/format.hpp>
#include <string>
#include <iostream>
#include <stdexcept>

namespace {

void lua_status_msg(int lua_status, const char* lua_cfg_file_name)
{
    using namespace std;

    cout << "\nLua status:\n"
         << (lua_status == 0 ? "succeed" :
             lua_status == LUA_YIELD ? "yield" :
             lua_status == LUA_ERRRUN ? "runtime error" :
             lua_status == LUA_ERRSYNTAX ? "syntax error during pre-compilation" :
             lua_status == LUA_ERRMEM ? "memory allocation error" :
             lua_status == LUA_ERRERR ? "error while running the error handler function" : "undefined")
        << endl;

    if (lua_status != 0)
    {
        string msg = str(boost::format("Error in \"%s\"") % lua_cfg_file_name);
        throw runtime_error(msg);
    }
}

int luabind_error_handler(lua_State* lua)
{
    // log the error message
    luabind::object msg(luabind::from_stack(lua, -1));
    std::ostringstream str;
    str << "lua> run-time error: " << msg;
    std::cout << str.str();

    // log the call stack
    std::string traceback = luabind::call_function<std::string>(luabind::globals(lua)["debug"]["traceback"]);
    traceback = std::string("lua> ") + traceback;
    std::cout << traceback.c_str();

    // return unmodified error object
    return 1;
}

} // unnamed namespace

// bind PIC objects to Lua
bool bind_to_lua(const char* lua_cfg_file_name, Grid& grid, Particles& particles)
{
    using namespace luabind;
    using namespace PIC;

    UseLua lua; // Lua interpreter

    bool ok = true; // result

    module(lua) // bind functions classes
    [
        // bind gather function from_grid_to_point as "pic_gather"
        def("pic_gather", (void(*)(const Grid&, const DblVector&, Grid::NodeType&)) &from_grid_to_point),

        // bind save function from_grid_to_point as "pic_save_grid_node"
        def("pic_save_grid_node", (void(*)(const std::string& prefix, const Grid& grid)) &save_grid_node),

        // bind DblVector
        class_<DblVector>("DblVector")
            .def(constructor<>())
            .def(constructor<const DblVector&>())
            .def(constructor<double /* x */, double /* y */, double /* z */>())
            .property("x", &DblVector::get_x, &DblVector::set_x)
            .property("y", &DblVector::get_y, &DblVector::set_y)
            .property("z", &DblVector::get_z, &DblVector::set_z)
            .def("abs", &DblVector::abs),

        // bind Cell
        class_<Cell>("Cell")
            .def(constructor<>())
            .def(constructor<double /* NP */,
                             const DblVector& /* vec_B */,
                             const DblVector& /* vec_E */,
                             const DblVector& /* vec_UE */,
                             const DblVector& /* vec_UP */>())
            .def(constructor<const Cell&>())
            .def_readwrite("NP", &Cell::NP)
            .def_readwrite("B", &Cell::B)
            .def_readwrite("E", &Cell::E)
            .def_readwrite("UE", &Cell::UE)
            .def_readwrite("UP", &Cell::UP)
            .enum_("CellState")
            [
                value("cs_active", cs_active),
                value("cs_absorptive", cs_absorptive),
                value("cs_custom", cs_custom)
            ]
            .property("state", &Cell::state, &Cell::set_state),

        // bind Grid
        class_<Grid>("Grid")
            .def("at", (Cell& (Grid::*)(index_t, index_t, index_t)) &Grid::at)
            .def("at", (const Cell& (Grid::*)(index_t, index_t, index_t) const) &Grid::at)
            .def("set", &Grid::set)
            .def("size_x", &Grid::size_x)
            .def("size_y", &Grid::size_y)
            .def("size_z", &Grid::size_z)
            .def("resize", &Grid::resize)
            .def("length_x", &Grid::length_x)
            .def("length_y", &Grid::length_y)
            .def("length_z", &Grid::length_z)
            .def("set_boundary_state", &Grid::set_boundary_state)
            .property("step", &Grid::step, &Grid::set_step),

        // bind Density
        class_<Density>("Density")
            .def(constructor<>())
            .def_readwrite("NP", &Density::NP)
            .def_readwrite("UP", &Density::UP),

        // bind DensityGrid
        class_<DensityGrid>("DensityGrid")
            .def("at", (Density& (DensityGrid::*)(index_t, index_t, index_t)) &DensityGrid::at)
            .def("at", (const Density& (DensityGrid::*)(index_t, index_t, index_t) const) &DensityGrid::at)
            .def("set", &DensityGrid::set)
            .def("size_x", &DensityGrid::size_x)
            .def("size_y", &DensityGrid::size_y)
            .def("size_z", &DensityGrid::size_z)
            .def("resize", &DensityGrid::resize)
            .def("length_x", &DensityGrid::length_x)
            .def("length_y", &DensityGrid::length_y)
            .def("length_z", &DensityGrid::length_z)
            .property("step", &DensityGrid::step, &DensityGrid::set_step),

        // bind UserGridFilters
        class_<UserGridFilters>("GridSaveFilters")
            .def("at", &UserGridFilters::at)
            .def("add", &UserGridFilters::append),

        // bind Particle
        class_<Particle>("Particle")
            .def(constructor<>())
            .def(constructor<const Particle&>())
            .def(constructor<const string& /* group_name */,
                             const DblVector& /* vec_r */,
                             const DblVector& /* vec_v */,
                             double /* ni */>())
            .def(constructor<const string& /* group_name */,
                             double /* x */, double /* y */, double/* z */,
                             double /* v_x */, double /* v_y */, double /* v_z */,
                             double/* ni_*/>())
            .def_readwrite("ni", &Particle::ni)
            .def_readwrite("r", &Particle::r)
            .def_readwrite("v", &Particle::v)
            .def_readonly("grp", &Particle::group_name),

        // bind Particles
        class_<Particles>("Particles")
            .def(constructor<>())
            .def("at", (Particle& (Particles::*)(index_t)) &Particles::at)
            .def("at", (const Particle& (Particles::*)(index_t) const) &Particles::at)
            .def("set", (void (Particles::*)(index_t /* i */, const Particle&)) &Particles::set)
            .def("set", (void (Particles::*)(index_t /* i */,
                                             const string& /* group_name */,
                                             const DblVector& /* vec_r */,
                                             const DblVector& /* vec_v */,
                                             double /* ni_*/))
                                             &Particles::set)
            .def("add", (void (Particles::*)(const Particle&)) &Particles::add)
            .def("add", (void (Particles::*)(const string& /* group_name */,
                                             const DblVector& /* vec_r */,
                                             const DblVector& /* vec_v */,
                                             double /* ni_*/))
                                             &Particles::add)
            .def("erase", &Particles::erase)
            .def("init_ball", &Particles::init_ball)

            .property("size", &Particles::size, &Particles::resize),

        // bind ParticleGroups
        class_<ParticleGroups>("ParticleGroups")
            .enum_("Diag")
            [
                value("no_diag", 0),
                value("save_positions", save_positions),
                value("save_grid_values", save_grid_values),
                value("save_all", save_all)
            ]
            .def("create_group", &ParticleGroups::create_group)
            .def(tostring(self)),

        // bind MarkerParticles
        class_<MarkerParticles>("MarkerParticles")
            .def("add", &MarkerParticles::insert_idx),

        // bind Config::Parameters
        class_<Config::Parameters>("Parameters")
            .def_readwrite("tau", &Config::Parameters::tau)
            .def_readwrite("time_steps", &Config::Parameters::time_steps)
            .def_readonly("current_time_step", &Config::Parameters::current_time_step)
            .def_readwrite("save_time_steps", &Config::Parameters::save_time_steps)
            .def_readwrite("dens_cutoff", &Config::Parameters::dens_cutoff)
            .def_readwrite("save_all_particles", &Config::Parameters::save_all_particles)
            .def_readwrite("save_whole_grid", &Config::Parameters::save_whole_grid)
            .def_readwrite("save_grid_x_plains", &Config::Parameters::save_grid_x_plains)
            .def_readwrite("save_grid_y_plains", &Config::Parameters::save_grid_y_plains)
            .def_readwrite("save_grid_z_plains", &Config::Parameters::save_grid_z_plains)
            .def_readonly("os_name", &Config::Parameters::os_name)
            .def_readonly("process_idx", &Config::Parameters::process_rank)
            .def_readwrite("L_scale", &Config::Parameters::L_scale)
            .def_readwrite("T_scale", &Config::Parameters::T_scale)
            .def_readwrite("U_scale", &Config::Parameters::U_scale)
            .def_readwrite("N_scale", &Config::Parameters::N_scale)
            .def_readwrite("E_scale", &Config::Parameters::E_scale)
            .def_readwrite("B_scale", &Config::Parameters::B_scale)

            .enum_("Push_Method")
            [
                value("Direct", Direct),
                value("Boris", Boris)
            ]
            .def_readwrite("push_method", &Config::Parameters::particle_push_alg)

            .enum_("Scatter_Method")
            [
                value("Standard", Standard),
                value("Zigzag", Zigzag)
            ]
            .def_readwrite("scatter_method", &Config::Parameters::scatter_alg)

            .enum_("Grid_Threshold")
            [
                value("Min_Density", Min_Density),
                value("Local_CFL", Local_CFL)
            ]

            .def_readwrite("grid_threshold", &Config::Parameters::grid_threshold)
            .def_readwrite("density_threshold", &Config::Parameters::grid_threshold)

            .enum_("CFL_Severity")
            [
                value("Ignore", Ignore),
                value("Absorb", Absorb),
                value("Stop", Stop)
            ]

            .def_readwrite("cfl_severity", &Config::Parameters::CFL_severity)

            .def(tostring(self))

    ]; // module(lua)

    try
    {
        // open Lua libraries
        luaL_openlibs(lua);

        // assign grid to a global in Lua
        globals(lua)["pic_grid"] = &grid;

        // assign particles to a global in Lua
        globals(lua)["pic_particles"] = &particles;

        // assign grid_save_filters to a global in Lua
        UserGridFilters grid_save_filters;
        globals(lua)["pic_grid_save_filters"] = &grid_save_filters;

        // assign part_groups to a global in Lua
        ParticleGroups part_groups;
        globals(lua)["pic_particle_groups"] = &part_groups;

        // assign markers to a global in Lua
        MarkerParticles markers;
        globals(lua)["pic_marker_particles"] = &markers;

        Config::Parameters& parameters = Config::parameters();

        parameters.cfg_script_name = lua_cfg_file_name;

        // assign parameters to a global in Lua
        globals(lua)["pic_parameters"] = &parameters;

        set_pcall_callback(luabind_error_handler);

        int lua_status = luaL_dofile(lua, lua_cfg_file_name);
        if (lua_status != 0)
            luabind_error_handler(lua);

        lua_status_msg(lua_status, lua_cfg_file_name);

        parameters.h = grid.step();
        parameters.grid_size_x = grid.size_x();
        parameters.grid_size_y = grid.size_y();
        parameters.grid_size_z = grid.size_z();
        parameters.total_particles_num = particles.size();

        if (!parameters.is_valid())
        {
            cerr << parameters;
            throw runtime_error("Invalid parameters");
        }
    }
    catch (error& e)
    {
        cerr << lua_tostring(lua, -1) << endl;
        cerr << "\n" << e.what() << endl;
        ok = false;
    }
    catch (exception& e)
    {
        cerr << "\n" << e.what() << endl;
        ok = false;
    }
    catch (...)
    {
        cerr << "\nUndefined error occurred."
                "\nPlease check \""
             << lua_cfg_file_name << "\"" << endl;
        ok = false;
    }

    return ok;
}
