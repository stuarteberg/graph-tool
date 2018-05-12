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

#ifndef GRAPH_BLOCKMODEL_UNCERTAIN_MARGINAL_HH
#define GRAPH_BLOCKMODEL_UNCERTAIN_MARGINAL_HH

#include "config.h"

#include <vector>

#include "../support/graph_state.hh"
#include "../blockmodel/graph_blockmodel_util.hh"
#include "graph_blockmodel_sample_edge.hh"

namespace graph_tool
{
using namespace boost;
using namespace std;

template <class Graph, class UGraph, class Eprop>
void collect_marginal(Graph& g, UGraph& u, Eprop ecount)
{
    typedef typename graph_traits<Graph>::edge_descriptor edge_t;
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    gt_hash_map<std::tuple<vertex_t, vertex_t>, edge_t> emap;
    for (auto e : edges_range(g))
    {
        std::tuple<vertex_t, vertex_t> vs(source(e, g), target(e, g));
        if (!graph_tool::is_directed(g) && get<0>(vs) > get<1>(vs))
            std::swap(get<0>(vs), get<1>(vs));
        emap[vs] = e;
    }

    for (auto e : edges_range(u))
    {
        std::tuple<vertex_t, vertex_t> vs(source(e, u), target(e, u));
        if (!graph_tool::is_directed(g) && get<0>(vs) > get<1>(vs))
            std::swap(get<0>(vs), get<1>(vs));
        edge_t ge;
        auto iter = emap.find(vs);
        if (iter == emap.end())
        {
            ge = add_edge(get<0>(vs), get<1>(vs), g).first;
            emap[vs] = ge;
            ecount[ge] = 0;
        }
        else
        {
            ge = iter->second;
        }
        ecount[ge]++;
    }
}


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_UNCERTAIN_MARGINAL_HH
