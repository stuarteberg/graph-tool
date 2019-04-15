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

#include "graph_tool.hh"
#include "graph_blockmodel_uncertain_marginal.hh"

using namespace boost;
using namespace graph_tool;

class dummy_property
{};

template <class Key>
constexpr double get(dummy_property&, Key&&) { return 0.; }

template <class Key>
constexpr void put(dummy_property&, Key&&, double) {}

void collect_marginal_dispatch(GraphInterface& gi, GraphInterface& ui,
                               boost::any aecount)
{
    typedef eprop_map_t<int32_t>::type emap_t;
    auto ecount = any_cast<emap_t>(aecount);

    gt_dispatch<>()
        ([&](auto& g, auto& u) { collect_marginal(g, u, ecount,
                                                  dummy_property(),
                                                  dummy_property(),
                                                  dummy_property()); },
         all_graph_views(), all_graph_views())(gi.get_graph_view(),
                                               ui.get_graph_view());
}

void collect_xmarginal_dispatch(GraphInterface& gi, GraphInterface& ui,
                                boost::any ax, boost::any aecount,
                                boost::any axsum, boost::any ax2sum)
{
    typedef eprop_map_t<int32_t>::type emap_t;
    auto ecount = any_cast<emap_t>(aecount);

    typedef eprop_map_t<double>::type xmap_t;
    auto x = any_cast<xmap_t>(ax);
    auto xsum = any_cast<xmap_t>(axsum);
    auto x2sum = any_cast<xmap_t>(ax2sum);

    gt_dispatch<>()
        ([&](auto& g, auto& u) { collect_marginal(g, u, ecount,
                                                  x, xsum, x2sum); },
         all_graph_views(), all_graph_views())(gi.get_graph_view(),
                                               ui.get_graph_view());
}
