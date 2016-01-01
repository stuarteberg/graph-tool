// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_KCORE_HH
#define GRAPH_KCORE_HH

namespace graph_tool
{
using namespace std;
using namespace boost;

struct kcore_decomposition
{
    template <class Graph, class VertexIndex, class CoreMap, class DegSelector>
    void operator()(Graph& g, VertexIndex vertex_index, CoreMap core_map,
                    DegSelector degS) const
    {
        unchecked_vector_property_map<size_t, VertexIndex>
            deg(vertex_index, num_vertices(g));
        unchecked_vector_property_map<size_t, VertexIndex>
            pos(vertex_index, num_vertices(g));

        typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

        vector<vector<vertex_t>> bins;

        for (auto v : vertices_range(g))
        {
            size_t k = degS(v, g);
            deg[v] = k;
            if (k >= bins.size())
                bins.resize(k + 1);
            bins[k].push_back(v);
            pos[v] = bins[k].size() - 1;
        }

        for (size_t k = 0; k < bins.size(); ++k)
        {
            auto& bins_k = bins[k];
            while (!bins_k.empty())
            {
                vertex_t v = bins_k.back();
                bins_k.pop_back();
                core_map[v] = k;
                for (auto e : out_edges_range(v, g))
                {
                    vertex_t u = target(e, g);
                    size_t ku = deg[u];
                    if (ku > deg[v])
                    {
                        auto& bins_ku = bins[ku];
                        vertex_t w = bins_ku.back();
                        auto pos_w = pos[w] = pos[u];
                        bins_ku[pos_w] = w;
                        bins_ku.pop_back();
                        auto& bins_ku_m = bins[ku - 1];
                        bins_ku_m.push_back(u);
                        pos[u] = bins_ku_m.size() - 1;
                        --deg[u];
                    }
                }
            }
        }
    }
};

} // graph_tool namespace

#endif // GRAPH_KCORE_HH
