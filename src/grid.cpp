#include <ostream>

#include "grid.h"

std::ostream& operator<<(std::ostream& out, const Cell& cell)
{
   return (out << "\n" << cell.NP << cell.B << "\n" << cell.E << "\n" << cell.UE << "\n" << cell.UP << "\n");
}
