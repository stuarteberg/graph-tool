// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"

#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

#include "graph_properties_copy.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

struct graph_views:
    boost::mpl::transform<graph_tool::detail::all_graph_views,
                          boost::mpl::quote1<std::add_pointer> >::type {};

void GraphInterface::copy_edge_property(const GraphInterface& src,
                                        boost::any prop_src,
                                        boost::any prop_tgt)
{
    run_action<>()
        (*this, std::bind(copy_property<edge_selector,edge_properties>(),
                          std::placeholders::_1, std::placeholders::_2,
                          std::placeholders::_3, prop_src),
         graph_views(), writable_edge_properties())
        (src.get_graph_view(), prop_tgt);
}
