// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2018 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "graph_properties.hh"
#include "graph.hh"

#include "graph_assortativity.hh"

using namespace std;
using namespace graph_tool;

pair<double,double>
assortativity_coefficient(GraphInterface& gi, GraphInterface::deg_t deg,
                          boost::any weight)
{
    typedef UnityPropertyMap<size_t,GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        weight_props_t;

    if (!weight.empty() && !belongs<edge_scalar_properties>()(weight))
        throw ValueException("weight edge property must have a scalar value type");

    if(weight.empty())
        weight = weight_map_t();

    double a = 0, a_err = 0;
    run_action<>()(gi,std::bind(get_assortativity_coefficient(),
                                std::placeholders::_1, std::placeholders::_2,
                                std::placeholders::_3,
                                std::ref(a), std::ref(a_err)),
                   scalar_selectors(), weight_props_t())
        (degree_selector(deg), weight);
    return make_pair(a, a_err);
}

pair<double,double>
scalar_assortativity_coefficient(GraphInterface& gi, GraphInterface::deg_t deg,
                                 boost::any weight)
{
    typedef UnityPropertyMap<size_t,GraphInterface::edge_t> weight_map_t;
    typedef boost::mpl::push_back<edge_scalar_properties, weight_map_t>::type
        weight_props_t;

    if (!weight.empty() && !belongs<edge_scalar_properties>()(weight))
        throw ValueException("weight edge property must have a scalar value type");

    if(weight.empty())
        weight = weight_map_t();

    double a = 0, a_err = 0;
    run_action<>()(gi, std::bind(get_scalar_assortativity_coefficient(),
                                 std::placeholders::_1, std::placeholders::_2,
                                 std::placeholders::_3,
                                 std::ref(a), std::ref(a_err)),
                   scalar_selectors(), weight_props_t())
        (degree_selector(deg), weight);
    return make_pair(a, a_err);
}

#include <boost/python.hpp>

using namespace boost::python;

void export_assortativity()
{
    def("assortativity_coefficient", &assortativity_coefficient);
    def("scalar_assortativity_coefficient", &scalar_assortativity_coefficient);
}
