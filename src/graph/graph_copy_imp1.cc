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
#include "graph_selectors.hh"

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/python/extract.hpp>

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class GraphSrc, class GraphTgt, class IndexMap, class SrcIndexMap>
struct copy_edge_property_dispatch
{
    copy_edge_property_dispatch(const GraphSrc& src, const GraphTgt& tgt,
                                boost::any& prop_src, boost::any& prop_tgt,
                                IndexMap& index_map,
                                SrcIndexMap& src_edge_index,
                                size_t max_src_edge_index,
                                bool& found)
        : src(src), tgt(tgt), prop_src(prop_src),
          prop_tgt(prop_tgt), index_map(index_map),
          src_edge_index(src_edge_index),
          max_src_edge_index(max_src_edge_index), found(found) {}


    const GraphSrc& src;
    const GraphTgt& tgt;
    boost::any& prop_src;
    boost::any& prop_tgt;
    IndexMap& index_map;
    SrcIndexMap& src_edge_index;
    size_t max_src_edge_index;
    bool& found;

    template <class PropertyMap>
    void operator()(PropertyMap) const
    {
        PropertyMap* psrc = any_cast<PropertyMap>(&prop_src);
        if (psrc == NULL)
            return;
        if (prop_tgt.empty())
            prop_tgt = PropertyMap(get(edge_index_t(), tgt));
        PropertyMap* ptgt = any_cast<PropertyMap>(&prop_tgt);
        if (ptgt == NULL)
            return;
        found = true;

        auto p_src = psrc->get_unchecked(max_src_edge_index + 1);
        auto p_tgt = ptgt->get_unchecked(num_edges(tgt));

        int i, N = num_vertices(src);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            auto v = vertex(i, src);
            if (v == graph_traits<GraphSrc>::null_vertex())
                continue;

            for (auto e : out_edges_range(v, src))
            {
                auto s = source(e, src);
                auto t = target(e, src);
                if (!is_directed::apply<GraphSrc>::type::value && s > t)
                    continue;
                size_t ei = src_edge_index[e];
                const auto& new_e = index_map[ei];
                p_tgt[new_e] = p_src[e];
            }
        }
    }
};

template <class PropertyMaps, class GraphSrc, class GraphTgt,
          class IndexMap, class SrcIndexMap>
void copy_edge_property(boost::any& prop_src, boost::any& prop_tgt,
                        const GraphSrc& src, const GraphTgt& tgt,
                        IndexMap& index_map, SrcIndexMap& src_vertex_index,
                        size_t max_src_edge_index)
{
    bool found = false;
    boost::mpl::for_each<PropertyMaps>(copy_edge_property_dispatch<GraphSrc, GraphTgt,
                                                                   IndexMap, SrcIndexMap>
            (src, tgt, prop_src, prop_tgt, index_map, src_vertex_index,
             max_src_edge_index, found));
    if (!found)
        throw ValueException("Cannot find property map type.");
}
