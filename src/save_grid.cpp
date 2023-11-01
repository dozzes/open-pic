#include "save_grid.h"

void save_grid_node(const std::string& prefix, const Grid& grid)
{
    using namespace PIC;

    std::string f_name = create_out_file_name(prefix, "grd_nodes_z=17_", PIC::Config::current_time_step());
    ofstream ofsg(f_name.c_str());
    ofsg << "x\ty\tz\tNP\tUP\tUPx\tUPy\tUPz\tB\tBx\tBy\tBz\tE\tEx\tEy\tEz\n";

    for (index_t i = 0; i != (grid.size_x()); ++i)
    for (index_t j = 0; j != (grid.size_y()); ++j)
    for (index_t k = 0; k != (grid.size_z()); ++k)
    {
        const Grid::NodeType& c = grid(i, j, k);
        const double UP = sqrt(c.UP.x*c.UP.x + c.UP.y*c.UP.y + c.UP.z*c.UP.z);
        const double B = sqrt(c.B.x*c.B.x + c.B.y*c.B.y + c.B.z*c.B.z);
        const double E = sqrt(c.E.x*c.E.x + c.E.y*c.E.y + c.E.z*c.E.z);
        ofsg << i << "\t" << j << "\t" << k
             << "\t" << c.NP
             << "\t" << UP << "\t" << c.UP.x << "\t" << c.UP.y << "\t" << c.UP.z
             << "\t" << B << "\t" << c.B.x << "\t" << c.B.y << "\t" << c.B.z
             << "\t" << E << "\t" << c.E.x << "\t" << c.E.y << "\t" << c.E.z
             << endl;
    }
}
