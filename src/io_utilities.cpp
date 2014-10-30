#include <ctime>
#include <boost/format.hpp>

#include "io_utilities.h"


using namespace std;

void print_header()
{
   cout << "\n******************************************************\n"
           "*                                                    *\n"
           "* 3D hybrid particle-in-cell plasma simulation tool. *\n"
           "*                                                    *\n"
           "* Institute of Radiophysics and Electronics NAS RA.  *\n"
           "*                                                    *\n"
           "*                   www.irphe.am                     *\n"
           "*                                                    *\n"
           "* version 1.2 08/22/2007.        osipyan@gmail.com   *\n"
           "*                                                    *\n"
           "******************************************************\n\n";
}

clock_t elapsed_seconds()
{
	static time_t start = time(0);
	return (time(0) - start);
}

string create_out_file_name(const string& prefix,
                            const string& type_name,
                            index_t step)
{
   using boost::format;
   using boost::io::group;

   static index_t time_steps_width = boost::lexical_cast<string>(PIC::Config::time_steps()).size();

   boost::format fmter("%1%_%2%_%3%.dat");
   fmter.modify_item(3, group(setw(time_steps_width), setfill('0')));
   fmter % prefix % type_name % step;

   return fmter.str();
}
