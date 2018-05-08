#pragma once

#include "config.h"
#include "vector_3d.h"

#include <boost/format.hpp>
#include <vector>
#include <stdexcept>
#include <iosfwd>

namespace PIC {

enum CellState { cs_active = 0, cs_absorptive = 1, cs_custom = 2 };

} // namespace PIC

/************************************************************************/
/* Grid class                                                           */
/************************************************************************/

typedef Vector3D<double> DblVector;

struct Cell
{
    Cell() : NP(0.0), state_(PIC::cs_active) { }

    Cell(double NP_,
         const DblVector& vec_B, const DblVector& vec_E,
         const DblVector& vec_UE, const DblVector& vec_UP)
    : NP(NP_),
      B(vec_B), E(vec_E), UE(vec_UE),UP(vec_UP),
      state_(PIC::cs_active) { }

    Cell(const Cell& other)
    : NP(other.NP),
      B(other.B), E(other.E), UE(other.UE), UP(other.UP),
      state_(other.state()) { }

    Cell(double NP_,
         double Bx,  double By,  double Bz,
         double Ex,  double Ey,  double Ez,
         double UEx, double UEy, double UEz,
         double UPx, double UPy, double UPz)
    : NP(NP_),
      B(Bx,By, Bz),
      E(Ex, Ey, Ez),
      UE(UEx, UEy, UEz),
      UP(UPx, UPy, UPz),
      state_(PIC::cs_active) { }

    Cell& operator=(const Cell& rhs)
    {
       if ((void*)this == (void*)&rhs)
          return *this;

       NP = rhs.NP;
       B = rhs.B;
       E = rhs.E;
       UE = rhs.UE;
       UP = rhs.UP;
       state_ = rhs.state();

       return *this;
    }

    void set_state(PIC::CellState new_state) { state_ = new_state; }
    PIC::CellState state() const { return state_; }

    double NP;    // ions density (cm^-3)
    DblVector B;  // magnetic field
    DblVector E;  // electric field
    DblVector UE; // electron velocity
    DblVector UP; // ion velocity

private:
    PIC::CellState state_;
};

std::ostream& operator<<(std::ostream& out, const Cell& cell);

template<class NodeT>
class GridContainer
{
    template <class U> friend class GridContainer;

public:
    typedef NodeT NodeType;

    GridContainer() : size_x_(0), size_y_(0), size_z_(0), h_(0.0) { }

    GridContainer(index_t sx, index_t sy, index_t sz, double h)
    : size_x_(sx), size_y_(sy), size_z_(sz), h_(h)
    {
        resize(size_x_, size_y_, size_z_);
    }

    GridContainer(const GridContainer<NodeT>& other)
    : size_x_(other.size_x()), size_y_(other.size_y()), size_z_(other.size_z()),
      h_(other.step()), data_(other.data_) { }

    GridContainer<NodeT>& operator=(const GridContainer<NodeT>& rhs)
    {
        if ((void*)this == (void*)&rhs)
            return *this;

        size_x_ = rhs.size_x();
        size_y_ = rhs.size_y();
        size_z_ = rhs.size_z();
        h_ = rhs.step();
        data_ = rhs.data_;

        return *this;
    }

    void reset_current()
    {
        for (index_t i = 0; i != data_.size(); ++i)
        {
            NodeT& node = data_[i];

            node.NP = 0.0;
            node.UP.x = 0.0;
            node.UP.y = 0.0;
            node.UP.z = 0.0;
        }
    }

    template<class DenstyGridType>
    void add_current(const DenstyGridType& dg)
    {
        typedef typename DenstyGridType::NodeType DensityNode;

        for (index_t i = 0; i != data_.size(); ++i)
        {
            const DensityNode& dg_node = dg.data_[i];
            NodeT& node = this->data_[i];

            node.NP += dg_node.NP;
            node.UP += dg_node.UP;
        }
    }

    NodeT& operator()(index_t x, index_t y, index_t z)
    {
        return data_[(x*size_y_ + y)*size_z_ + z];
    }

    const NodeT& operator()(index_t x, index_t y, index_t z) const
    {
        return data_[(x*size_y_ + y)*size_z_ + z];
    }

    NodeT& at(index_t x, index_t y, index_t z)
    {
        return data_.at((x*size_y_ + y)*size_z_ + z);
    }

    const NodeT& at(index_t x, index_t y, index_t z) const
    {
        return data_.at((x*size_y_ + y)*size_z_ + z);
    }

    void set(index_t x, index_t y, index_t z, const NodeT& cell)
    {
        at(x, y, z) = cell;
    }

    void resize(index_t new_size_x, index_t new_size_y, index_t new_size_z)
    {
        size_x_ = new_size_x;
        size_y_ = new_size_y;
        size_z_ = new_size_z;

        data_.resize(size_x_*size_y_*size_z_);
    }

    void set_boundary_state(PIC::CellState new_state) //TODO: improve 3d matrix bound traversal
    {
        for (index_t i = 0; i != size_x(); ++i)
        for (index_t j = 0; j != size_y(); ++j)
        for (index_t k = 0; k != size_z(); ++k)
        {
            if (i == 0 || i == (size_x() - 1) ||
                j == 0 || j == (size_y() - 1) ||
                k == 0 || k == (size_z() - 1))
            {
                at(i, j, k).set_state(new_state);
            }
        }
    }

    void set_step(double h) { h_ = h; }
    double step() const     { return h_; }

    double cell_volume() const { return (h_*h_*h_);}

    index_t size_x() const { return size_x_; }
    index_t size_y() const { return size_y_; }
    index_t size_z() const { return size_z_; }

    double length_x() const { return ((size_x_ - 1)*h_); }
    double length_y() const { return ((size_y_ - 1)*h_); }
    double length_z() const { return ((size_z_ - 1)*h_); }

private:
    index_t size_x_;
    index_t size_y_;
    index_t size_z_;
    double h_;
    std::vector<NodeT> data_;
};

// Grid functions
template<class NodeT>
bool operator==(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2)
{
    return (g1.size_x() == g2.size_x() &&
            g1.size_y() == g2.size_y() &&
            g1.size_z() == g2.size_z() &&
            g1.step() == g2.step() &&
            g1.data_ == g2.data_);
}

template<class NodeT>
bool operator!=(const GridContainer<NodeT>& g1, const GridContainer<NodeT>& g2)
{
    return !(g1 == g2);
}

typedef GridContainer<Cell> Grid;
