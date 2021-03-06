#include "io_utilities.h"

#include <boost/format.hpp>
#include <ctime>

using namespace std;

void print_header()
{
    cout << "\n"
            "======================================================\n"
            " *                                                  * \n"
            "             _ \\               _ \\_ _|  __|         \n"
            "            (   |_ \\  -_)   \\  __/  |  (            \n"
            "           \\___/.__/\\___|_| _|_|  ___|\\___|        \n"
            "               _|                                     \n"
            " *                                                  * \n"
            "======================================================\n"
            "  3D hybrid particle-in-cell plasma simulation tool.  \n"
            "                                                      \n"
            "  Institute of Radiophysics and Electronics NAS RA.   \n"
            "                    www.irphe.am                      \n"
            "                                                      \n"
            "  version 2.0 05/08/2018         osipyan@gmail.com    \n"
            "======================================================\n"
            "\n";
}

clock_t elapsed_seconds()
{
    static time_t start = time(0);
    return (time(0) - start);
}

string create_out_file_name(const string& prefix, const string& type_name,
                            index_t step)
{
    using boost::format;
    using boost::io::group;

    static size_t time_steps_width = boost::lexical_cast<string>(PIC::Config::time_steps()).size();

    boost::format fmter("%1%_%2%_%3%.dat");
    fmter.modify_item(3, group(setw(time_steps_width), setfill('0')));
    fmter % prefix % type_name % step;

    return fmter.str();
}
