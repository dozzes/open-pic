#include "grid_filters.h"
#include "config.h"
#include "use_lua.h"
#include "io_utilities.h"

#include <boost/ref.hpp>
#include <boost/format.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace std;

// static
vector<string> UserGridFilters::filters_ = vector<string>();

UserGridFilter::UserGridFilter(const string& name) : GridFilter(name)
{
}

bool UserGridFilter::operator()(const DblVector& node, const DblVector& pos) const
{
    UseLua lua;

    bool do_save = false;

    try
    {
        do_save = luabind::call_function<bool>(lua, name_.c_str(), boost::ref(node), boost::ref(pos));
    }
    catch (exception& e)
    {
        cerr << "\n" << name_ << " : " << e.what() << endl;
    }

    return do_save;
}

SaveAllGrid::SaveAllGrid(const string& grid_group_name) : GridFilter(grid_group_name)
{
    name_ += "_xyz";
}

const string PlainFilter::plain_tags_ = "xyz";

PlainFilter::PlainFilter(const string& grid_group_name, Plain plain, index_t level)
: GridFilter(""), plain_(plain), level_(level)
{
    name_ = str(boost::format("%s_%s_%s") % grid_group_name
                                          % plain_tags_[plain_]
                                          % level_);
}

void PlainFilter::set_level(index_t level)
{
    level_ = level;
    index_t pos = name_.find_last_of('_');
    name_ = str(boost::format("%s%s") % name_.substr(0, pos + 1) % level_);
}

bool PlainFilter::operator()(const DblVector& node, const DblVector& /*pos*/) const
{
    return (level_ == (plain_ == X ? node.x : plain_ == Y ? node.y : node.z));
}
