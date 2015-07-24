#if !defined (GRID_FILTERS_H)
#define GRID_FILTERS_H

#include <string>

#include "opic_fwd.h"
#include "vector_3d.h"


typedef Vector3D<double> DblVector;

class GridFilter
{
public:
    GridFilter(const std::string& file_name) : name_(file_name) {}
    const std::string& name() const { return name_; }

protected:
    std::string name_;
};

class UserGridFilter : public GridFilter
{
public:
    UserGridFilter(const std::string& name);
    virtual bool operator()(const DblVector& node, const DblVector& pos) const;
};

class UserGridFilters
{
public:
    const std::string& operator[](index_t i) const  { return filters_[i]; }
    std::string& operator[](index_t i)              { return filters_[i]; }
    std::string& at(index_t i)                      { return filters_.at(i); }
    void append(const std::string& filter_name)     { filters_.push_back(filter_name); }
    void resize(index_t newsize)                    { filters_.resize(newsize); }
    index_t size()                                  { return filters_.size(); }
    const std::vector<std::string>& filter_names()  { return filters_; }

private:
    static std::vector<std::string> filters_;
};

class SaveAllGrid : public GridFilter
{
public:
    SaveAllGrid(const std::string& grid_group_name);
    bool operator()(const DblVector& /*node*/, const DblVector& /*pos*/) const { return true; }
};

class PlainFilter : public GridFilter
{
public:
    enum Plain { X = 0, Y, Z };

    PlainFilter(const std::string& grid_group_name, Plain plain, index_t level);
    void set_level(index_t level);
    bool operator()(const DblVector& node, const DblVector& pos) const;

private:
    Plain plain_;
    index_t level_;
    static const std::string plain_tags_;
};

#endif // GRID_FILTERS_H
