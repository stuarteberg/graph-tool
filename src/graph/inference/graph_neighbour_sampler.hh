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

#ifndef GRAPH_NEIGHBOUR_SAMPLER_HH
#define GRAPH_NEIGHBOUR_SAMPLER_HH

#include "config.h"

#include "graph_tool.hh"

// Sample neighbours efficiently
// =============================

namespace graph_tool
{

template <class Graph>
class UnweightedNeighbourSampler
{
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;

    template <class Eprop>
    UnweightedNeighbourSampler(Graph& g, Eprop&, bool self_loops=false)
        : _sampler_c(get(vertex_index_t(), g)),
          _sampler(_sampler_c.get_unchecked()),
          _sampler_pos_c(get(vertex_index_t(), g)),
          _sampler_pos(_sampler_pos_c.get_unchecked())

    {
        for (auto v : vertices_range(g))
        {
            auto& sampler = _sampler_c[v];
            auto& sampler_pos = _sampler_pos_c[v];
            sampler.clear();
            sampler_pos.clear();
            for (auto e : all_edges_range(v, g))
            {
                auto u = target(e, g);
                if (is_directed::apply<Graph>::type::value && u == v)
                    u = source(e, g);
                if (!self_loops && u == v)
                    continue;
                sampler_pos[u].push_back(sampler.size());
                sampler.push_back(u); // Self-loops will be added twice
            }
        }
    }

    template <class RNG>
    vertex_t sample(vertex_t v, RNG& rng)
    {
        return uniform_sample(_sampler[v], rng);
    }

    bool empty(vertex_t v)
    {
        return _sampler[v].empty();
    }

    template <class Weight>
    void remove(vertex_t v, vertex_t u, Weight)
    {
        auto& sampler = _sampler[v];
        auto& sampler_pos = _sampler_pos[v];

        auto& pos_u = sampler_pos[u];
        size_t i = pos_u.back();
        sampler[i] = sampler.back();
        for (auto& j : sampler_pos[sampler.back()])
        {
            if (j == sampler.size() - 1)
            {
                j = i;
                break;
            }
        }
        sampler.pop_back();
        pos_u.pop_back();
        if (pos_u.empty())
            sampler_pos.erase(u);
    }

    template <class Weight>
    void insert(vertex_t v, vertex_t u, Weight)
    {
        auto& sampler = _sampler[v];
        auto& sampler_pos = _sampler_pos[v];
        sampler_pos[u].push_back(sampler.size());
        sampler.push_back(u);
    }

private:
    typedef typename vprop_map_t<std::vector<vertex_t>>::type sampler_t;
    sampler_t _sampler_c;
    typename sampler_t::unchecked_t _sampler;

    typedef typename vprop_map_t<gt_hash_map<vertex_t, vector<size_t>>>::type sampler_pos_t;
    sampler_pos_t _sampler_pos_c;
    typename sampler_pos_t::unchecked_t _sampler_pos;
};

template <class Graph, template <class V> class Sampler>
class WeightedNeighbourSampler
{
public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor vertex_t;

    template <class Eprop>
    WeightedNeighbourSampler(Graph& g, Eprop& eweight, bool self_loops=false)
        : _sampler_c(get(vertex_index_t(), g)),
          _sampler(_sampler_c.get_unchecked()),
          _sampler_pos_c(get(vertex_index_t(), g)),
          _sampler_pos(_sampler_pos_c.get_unchecked())
    {
        for (auto v : vertices_range(g))
        {
            auto& sampler = _sampler_c[v];
            auto& sampler_pos = _sampler_pos_c[v];

            sampler.clear();
            sampler_pos.clear();

            for (auto e : all_edges_range(v, g))
            {
                auto u = target(e, g);
                if (is_directed::apply<Graph>::type::value && u == v)
                    u = source(e, g);
                if (!self_loops && u == v)
                    continue;
                auto w = eweight[e];
                if (w == 0)
                    continue;
                insert(v, u, w); // Self-loops will be added twice, and hence will
                                 // be sampled with probability 2 * eweight[e]
            }
        }
    }

    template <class RNG>
    vertex_t sample(vertex_t v, RNG& rng)
    {
        return _sampler[v].sample(rng);
    }

    bool empty(vertex_t v)
    {
        return _sampler[v].empty();
    }

    template <class Weight>
    void remove(vertex_t v, vertex_t u, Weight w)
    {
        auto& sampler = _sampler[v];
        auto& sampler_pos = _sampler_pos[v];

        auto& pos = sampler_pos[u];
        sampler.remove(pos.first);

        w -= pos.second;
        if (w > 0)
            insert(v, u, w);
        else
            sampler_pos.erase(u);
    }

    template <class Weight>
    void insert(vertex_t v, vertex_t u, Weight w)
    {
        auto& sampler = _sampler[v];
        auto& sampler_pos = _sampler_pos[v];
        auto pos = sampler_pos.find(u);
        if (pos != sampler_pos.end())
        {
            auto old_w = pos->second.second;
            remove(v, u, old_w);
            w += old_w;
        }
        sampler_pos[u] = std::make_pair(sampler.insert(u, w), w);
    }

private:
    typedef typename vprop_map_t<Sampler<vertex_t>>::type sampler_t;
    sampler_t _sampler_c;
    typename sampler_t::unchecked_t _sampler;

    typedef typename vprop_map_t<gt_hash_map<vertex_t, pair<size_t, double>>>::type sampler_pos_t;
    sampler_pos_t _sampler_pos_c;
    typename sampler_pos_t::unchecked_t _sampler_pos;
};

}

#endif // GRAPH_NEIGHBOUR_SAMPLER_HH
