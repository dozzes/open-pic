#include "config.h"
#include "use_lua.h"
#include "bind_to_lua.h"
#include "simulate.h"
#include "io_utilities.h"

#include "particles.h"
#include "grid.h"

#include <boost/format.hpp>


int main(int argc, char* argv[])
{
    print_header();

    // grid and particles objects
    Particles particles;
    Grid grid;

    try
    {
        if (argc != 2)
        {
            std::cerr << "Usage: pic3d <setup script in Lua>\n";
            std::cerr << "Press any key to exit ...";
            std::cin.get();
            return -1;
        }

        std::cout << "\nRun configuration script: \"" << argv[1] << "\" ...\n\n";
        if (!bind_to_lua(argv[1], grid, particles))
        {
            std::cerr << "\nLua script failed.";
            std::cerr << "\nPress any key to exit ...";
            std::cin.get();
            return -1;
        }

        PIC::Config::to_stream(std::cout);
        std::ofstream ofs_params("opic_config.txt");
        if (ofs_params)
        {
            ofs_params << std::boolalpha;
            PIC::Config::to_stream(ofs_params);
        }

        long start_time = elapsed_seconds();

        PIC::simulate(grid, particles);

        long run_time = (elapsed_seconds() - start_time);

        std::cout << "\nRun time: " << run_time << " sec.\nFinished!";

        PIC::Config::logger() << "\nRun time: " << run_time << " sec.\nFinished!";
    }
    catch (luabind::error& e)
    {
        UseLua lua; // Lua interpreter
        log(lua_tostring(lua, -1), true);
        log(e.what(), true);
    }
    catch (std::exception& e)
    {
        log(e.what(), true);
    }
    catch (...)
    {
        std::cerr << "\nUndefined error occurred."
                     "\nPlease check \""
                  << argv[1] << "\"" << std::endl;
    }

    return 0;
}
