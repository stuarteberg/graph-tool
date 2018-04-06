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

#ifndef GRAPH_SBM_SAMPLE_EDGE_HH
#define GRAPH_SBM_SAMPLE_EDGE_HH

#include <tuple>
#include <iostream>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"
#include "../../generation/sampler.hh"
#include "../../generation/urn_sampler.hh"

#include "random.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

template <class Graph>
class SBMEdgeSampler
{
public:
    SBMEdgeSampler()
        : _v_in_sampler((is_directed_::apply<Graph>::type::value) ?
                        __v_in_sampler : _v_out_sampler) {}

    template <class State>
    void sync(State& state)
    {
        vector<std::tuple<size_t, size_t>> rs;
        vector<size_t> rs_count;

        size_t E = 0;
        for (auto me : edges_range(state._bg))
        {
            size_t mrs = state._mrs[me];
            if (mrs == 0)
                continue;
            rs.emplace_back(source(me, state._bg),
                            target(me, state._bg));
            rs_count.push_back(mrs + 1);
            E += mrs;
        }

        rs.emplace_back(std::numeric_limits<size_t>::max(),
                        std::numeric_limits<size_t>::max());
        rs_count.push_back(E + 1);

        _rs_sampler = rs_sampler_t(rs, rs_count);

        vector<vector<size_t>> rvs;
        vector<vector<size_t>> v_in_probs, v_out_probs;

        bool deg_corr = state._deg_corr;
        for (auto v : vertices_range(state._g))
        {
            size_t r = state._b[v];
            if (r >= v_out_probs.size())
            {
                if (graph_tool::is_directed(state._g))
                    v_in_probs.resize(r+1);
                v_out_probs.resize(r+1);
                rvs.resize(r+1);
            }
            rvs[r].push_back(v);
            if (graph_tool::is_directed(state._g))
            {
                auto kin = (deg_corr) ?
                    in_degreeS()(v, state._g, state._eweight) : 0;
                v_in_probs[r].push_back(kin + 1);
            }
            auto kout = (deg_corr) ?
                out_degreeS()(v, state._g, state._eweight) : 0;
            v_out_probs[r].push_back(kout + 1);
        }

        __v_in_sampler.clear();
        _v_out_sampler.clear();
        _groups.clear();
        for (size_t r = 0; r < rvs.size(); ++r)
        {
            if (graph_tool::is_directed(state._g))
                __v_in_sampler.emplace_back(rvs[r], v_in_probs[r]);
            _v_out_sampler.emplace_back(rvs[r], v_out_probs[r]);
            if (!rvs[r].empty())
                _groups.push_back(r);
        }
    }

    template <class RNG>
    std::tuple<size_t, size_t> sample(RNG& rng)
    {
        auto rs = _rs_sampler.sample(rng);

        if (get<0>(rs) == std::numeric_limits<size_t>::max())
        {
            get<0>(rs) = uniform_sample(_groups, rng);
            get<1>(rs) = uniform_sample(_groups, rng);
        }

        auto& r_sampler = _v_out_sampler[get<0>(rs)];
        auto& s_sampler = _v_in_sampler[get<1>(rs)];
        return std::make_tuple(r_sampler.sample(rng),
                               s_sampler.sample(rng));
    }

private:
    typedef UrnSampler<std::tuple<size_t, size_t>, true> rs_sampler_t;
    rs_sampler_t _rs_sampler;

    typedef UrnSampler<size_t, true> vsampler_t;
    vector<vsampler_t> __v_in_sampler, _v_out_sampler;
    vector<vsampler_t>& _v_in_sampler;

    std::vector<size_t> _groups;
};


} // graph_tool namespace

#endif // GRAPH_SBM_SAMPLE_EDGE_HH
