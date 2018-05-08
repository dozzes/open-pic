#include "grid.h"

#include <ostream>

std::ostream& operator<<(std::ostream& out, const Cell& cell)
{
    return (out << "\n" << cell.NP << "\n"
                        << cell.B  << "\n"
                        << cell.E  << "\n"
                        << cell.UE << "\n"
                        << cell.UP << "\n");
}
