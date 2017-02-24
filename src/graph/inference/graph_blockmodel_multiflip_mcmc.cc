// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2017 Tiago de Paula Peixoto <tiago@skewed.de>
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
#include "random.hh"

#include <boost/python.hpp>

#include "graph_blockmodel_util.hh"
#include "graph_blockmodel.hh"
#include "graph_blockmodel_multiflip_mcmc.hh"
#include "mcmc_loop.hh"

using namespace boost;
using namespace graph_tool;

GEN_DISPATCH(block_state, BlockState, BLOCK_STATE_params)

template <class State>
GEN_DISPATCH(mcmc_block_state, MCMC<State>::template MCMCBlockState,
             MCMC_BLOCK_STATE_params(State))

python::object do_multiflip_mcmc_sweep(python::object omcmc_state,
                                       python::object oblock_state,
                                       rng_t& rng)
{
    python::object ret;
    auto dispatch = [&](auto& block_state)
    {
        typedef typename std::remove_reference<decltype(block_state)>::type
            state_t;

        mcmc_block_state<state_t>::make_dispatch
           (omcmc_state,
            [&](auto& s)
            {
                auto ret_ = mcmc_sweep(s, rng);
                ret = python::make_tuple(ret_.first, ret_.second);
            });
    };
    block_state::dispatch(oblock_state, dispatch);
    return ret;
}


void get_worder(GraphInterface& gi, boost::any abs, boost::any aorder)
{
    typedef vprop_map_t<int64_t>::type vmap_t;

    vmap_t& order = any_cast<vmap_t&>(aorder);

    run_action<>()(gi,
                   [&](auto& g, auto& bs)
                   {
                       std::vector<size_t> vs;
                       vs.reserve(num_vertices(g));
                       for (auto v : vertices_range(g))
                           vs.push_back(v);
                       std::sort(vs.begin(), vs.end(),
                                 [&](auto u, auto v)
                                 { return bs[u] < bs[v]; });
                       for (size_t i = 0; i < vs.size(); ++i)
                           order[vs[i]] = i;
                   },
                   vertex_scalar_vector_properties())(abs);
}


void shuffle_labels(GraphInterface& gi, boost::any ab, rng_t& rng)
{
    run_action<>()(gi,
                   [&](auto& g, auto& b)
                   {
                       typedef typename property_traits
                           <typename std::remove_reference<decltype(b)>::type>
                           ::value_type val_t;

                       std::unordered_set<val_t> vals;
                       for (auto v : vertices_range(g))
                           vals.insert(b[v]);

                       std::vector<val_t> nvals(vals.begin(), vals.end());
                       std::shuffle(nvals.begin(), nvals.end(), rng);

                       std::unordered_map<val_t, val_t> val_map;
                       auto iter = nvals.begin();
                       for (auto r : vals)
                           val_map[r] = *iter++;

                       for (auto v : vertices_range(g))
                           b[v] = val_map[b[v]];
                   },
                   writable_vertex_scalar_properties())(ab);
}


void export_blockmodel_multiflip_mcmc()
{
    using namespace boost::python;
    def("multiflip_mcmc_sweep", &do_multiflip_mcmc_sweep);
    def("get_worder", &get_worder);
    def("shuffle_labels", &shuffle_labels);
}
