#if !defined (SAVE_GRID_H)
#define SAVE_GRID_H

#include <cstdio>
#include <vector>

#include "grid_filters.h"
#include "gather_scatter.h"
#include "particle_groups.h"
#include "constants.h"
#include "io_utilities.h"


template<class GridT, class DensGridT, class PointFilterFoT>
void save_subgrid(const GridT& grid, const DensGridT& dens_grid, const PointFilterFoT& point_filter)
{
    using namespace PIC;

    std::string fileName = create_out_file_name(point_filter.name(), "grd", PIC::Config::current_time_step());

    FILE* fout = fopen(fileName.c_str(), "w");
    if (!fout)
    {
        std::cerr << "can't save data: " << fileName << std::endl;
        return;
    }

    fprintf(fout, "X\tY\tZ"
        "\tNP"
        "\tB\tBx\tBy\tBz"
        "\tE\tEx\tEy\tEz\tdiv(E)"
        "\tUP\tUPx\tUPy\tUPz\tMa"
        "\tUE\tUEx\tUEy\tUEz\n");

    const double h = grid.step();

    const double L = Config::L_scale();
    const double U = Config::U_scale();
    const double N = Config::N_scale();
    const double E = Config::E_scale();
    const double B = Config::B_scale();

    const index_t m = 1;
    for (index_t i = m; i != (grid.size_x() - m); ++i)
    for (index_t j = m; j != (grid.size_y() - m); ++j)
    for (index_t k = m; k != (grid.size_z() - m); ++k)
    {
        DblVector node(i, j, k);
        DblVector point(i*h, j*h, k*h);

        if (point_filter(node, point))
        {
            typename GridT::NodeType point_val;
            from_grid_to_point(grid, point, point_val);

            typename DensGridT::NodeType dgNode;
            dgNode.NP = gather_scalar<PIC::CellCentering>(dens_grid, point, &DensGridT::NodeType::NP);

            gather_face(dens_grid, point, &DensGridT::NodeType::UP, dgNode.UP);

            double vA = point_val.B.abs() / sqrt(4 * Constants::pi() * point_val.NP * Constants::mp());

            const double d = 0.5*h;
            const double Ex_r = gather_vector<PIC::FaceXCentering>(grid, point.x + d, point.y, point.z, &GridT::NodeType::E, &DblVector::x);
            const double Ex_l = gather_vector<PIC::FaceXCentering>(grid, point.x - d, point.y, point.z, &GridT::NodeType::E, &DblVector::x);

            const double Ey_r = gather_vector<PIC::FaceYCentering>(grid, point.x, point.y + d, point.z, &GridT::NodeType::E, &DblVector::y);
            const double Ey_l = gather_vector<PIC::FaceYCentering>(grid, point.x, point.y - d, point.z, &GridT::NodeType::E, &DblVector::y);

            const double Ez_r = gather_vector<PIC::FaceZCentering>(grid, point.x, point.y, point.z + d, &GridT::NodeType::E, &DblVector::z);
            const double Ez_l = gather_vector<PIC::FaceZCentering>(grid, point.x, point.y, point.z - d, &GridT::NodeType::E, &DblVector::z);

            const double divE = (L / E)*((Ex_r - Ex_l) + (Ey_r - Ey_l) + (Ez_r - Ez_l)) / h;

            fprintf(fout,
                   "%e\t%e\t%e\t" /* point.x/L, point.y/L, point.z/L */
                   "%e\t" /* dgNode.NP/N */
                   "%e\t" /* point_val.B.abs()/B */
                   "%e\t%e\t%e\t" /* point_val.B.x/B, point_val.B.y/B, point_val.B.z/B */
                   "%e\t" /* point_val.E.abs()/E */
                   "%e\t%e\t%e\t" /* point_val.E.x/E, point_val.E.y/E, point_val.E.z/E */
                   "%e\t" /* div(E)/(E/L) */
                   "%e\t" /* dgNode.UP.abs()/U */
                   "%e\t%e\t%e\t" /* dgNode.UP.x/U, dgNode.UP.y/U, dgNode.UP.z/U */
                   "%e\t" /* dgNode.UP.abs()/vA */
                   "%e\t" /* point_val.UE.abs()/U */
                   "%e\t%e\t%e\n", /* dpoint_val.UE.x/U, point_val.UE.y/U, point_val.UE.z/U */
                   (point.x / L), (point.y / L), (point.z / L),
                   (dgNode.NP / N),
                   (point_val.B.abs() / B),
                   (point_val.B.x / B), (point_val.B.y / B), (point_val.B.z / B),
                   (point_val.E.abs() / E),
                   (point_val.E.x / E), (point_val.E.y / E), (point_val.E.z / E),
                   divE,
                   (dgNode.UP.abs() / U),
                   (dgNode.UP.x / U), (dgNode.UP.y / U), (dgNode.UP.z / U),
                   (dgNode.UP.abs() / vA),
                   (point_val.UE.abs() / U),
                   (point_val.UE.x / U), (point_val.UE.y / U), (point_val.UE.z / U));
        }
    }

    fclose(fout);
}

template<class GridT, class DensGridT>
void save_grid_levels(const GridT& grid,
                      const DensGridT& dens_grid,
                      PlainFilter& plain_filter,
                      index_t from_level, index_t to_level)
{
    for (index_t lv = from_level; lv != to_level; ++lv)
    {
        plain_filter.set_level(lv);
        save_subgrid(grid, dens_grid, plain_filter);
    }
}

template<class GridT, class DensGridT>
void save_grid(const GridT& grid, const DensGridT& dens_grid, const std::string& group_name)
{
    using namespace PIC;

    Diagnostics group_diag = save_grid_values;

    if (group_name != "all")
    {
        ParticleGroups part_groups;
        const ParticleGroups::ParticleGroup& group = part_groups[group_name];

        group_diag = group.diag;
    }

    if (group_diag & save_grid_values)
    {
        UserGridFilters user_grid_filters;
        for (index_t f = 0; f != user_grid_filters.size(); ++f)
        {
            UserGridFilter user_grid_filter(user_grid_filters.filter_names().at(f));
            save_subgrid(grid, dens_grid, user_grid_filter);
        }

        const size_t m = 1;

        if (Config::save_grid_x_plains())
        {
            PlainFilter grid_filter(group_name, PlainFilter::X, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_x() - m);
        }

        if (Config::save_grid_y_plains())
        {
            PlainFilter grid_filter(group_name, PlainFilter::Y, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_y() - m);
        }

        if (Config::save_grid_z_plains())
        {
            PlainFilter grid_filter(group_name, PlainFilter::Z, 0);
            save_grid_levels(grid, dens_grid, grid_filter, m, grid.size_z() - m);
        }

        if (Config::save_whole_grid())
        {
            SaveAllGrid grid_filter(group_name);
            save_subgrid(grid, dens_grid, grid_filter);
        }
    }
}

#endif // SAVE_GRID_H

