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

#ifndef GRAPH_BLOCKMODEL_UNCERTAIN_MCMC_HH
#define GRAPH_BLOCKMODEL_UNCERTAIN_MCMC_HH

#include "config.h"

#include <vector>

#include "graph_tool.hh"
#include "../support/graph_state.hh"
#include "graph_blockmodel_uncertain.hh"

namespace graph_tool
{
using namespace boost;
using namespace std;

typedef std::vector<size_t> vlist_t;

#define MCMC_UNCERTAIN_STATE_params(State)                                     \
    ((__class__,&, mpl::vector<python::object>, 1))                            \
    ((state, &, State&, 0))                                                    \
    ((slist, &, vlist_t&, 0))                                                  \
    ((tlist, &, vlist_t&, 0))                                                  \
    ((beta,, double, 0))                                                       \
    ((entropy_args,, entropy_args_t, 0))                                       \
    ((verbose,, bool, 0))                                                      \
    ((niter,, size_t, 0))


template <class State>
struct MCMC
{
    GEN_STATE_BASE(MCMCUncertainStateBase, MCMC_UNCERTAIN_STATE_params(State))

    template <class... Ts>
    class MCMCUncertainState
        : public MCMCUncertainStateBase<Ts...>
    {
    public:
        GET_PARAMS_USING(MCMCUncertainStateBase<Ts...>,
                         MCMC_UNCERTAIN_STATE_params(State))
        GET_PARAMS_TYPEDEF(Ts, MCMC_UNCERTAIN_STATE_params(State))

        template <class... ATs,
                  typename std::enable_if_t<sizeof...(ATs) ==
                                            sizeof...(Ts)>* = nullptr>
        MCMCUncertainState(ATs&&... as)
            : MCMCUncertainStateBase<Ts...>(as...),
              _vlist(_slist.size())
        {
            for (size_t i = 0; i < _vlist.size(); ++i)
                _vlist[i] = i;
        }

        std::pair<size_t, size_t> _e;
        std::vector<size_t> _vlist;
        int _null_move = 0;

        std::pair<size_t, size_t> get_edge(size_t ei)
        {
            size_t u = _slist[ei];
            size_t v = _tlist[ei];
            if (u > num_vertices(_state._g))
                std::tie(u, v) = _e;
            return {u, v};
        }

        size_t node_state(size_t u, size_t v)
        {
            auto&& m = _state.get_u_edge(u, v);
            if (m == _state._null_edge)
                return 0;
            return _state._eweight[m];
        }

        size_t node_state(size_t ei)
        {
            size_t u, v;
            std::tie(u, v) = get_edge(ei);
            return node_state(u, v);
        }

        bool skip_node(auto&)
        {
            return false;
        }

        size_t node_weight(auto&)
        {
            return 1;
        }

        template <class RNG>
        int move_proposal(size_t ei, RNG& rng)
        {
            if (_slist[ei] >= num_vertices(_state._g))
            {
                std::uniform_int_distribution<size_t>
                    sample(0, num_vertices(_state._g) - 1);
                _e = {sample(rng), sample(rng)};
            }

            std::bernoulli_distribution coin(.5);
            if (coin(rng))
            {
                size_t m = node_state(ei);
                if (m > 0)
                    return -1;
                else
                    return _null_move;
            }
            else
            {
                return 1;
            }
        }

        std::tuple<double, double>
        virtual_move_dS(size_t ei, int dm)
        {
            if (dm == 0)
                return {0., 0.};

            size_t u, v;
            std::tie(u, v) = get_edge(ei);

            double dS = 0;
            if (dm < 0)
                dS = _state.remove_edge_dS(u, v, _entropy_args);
            else
                dS = _state.add_edge_dS(u, v, _entropy_args);

            return {dS, 0.};
        }

        void perform_move(size_t ei, int dm)
        {
            if (dm == 0)
                return;
            size_t u, v;
            std::tie(u, v) = get_edge(ei);
            if (dm < 0)
                _state.remove_edge(u, v);
            else
                _state.add_edge(u, v);
        }

        bool is_deterministic()
        {
            return false;
        }

        bool is_sequential()
        {
            return false;
        }

        auto& get_vlist()
        {
            return _vlist;
        }

        double get_beta()
        {
            return _beta;
        }

        size_t get_niter()
        {
            return _niter;
        }

        void step(auto&, int)
        {
        }
    };
};


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_UNCERTAIN_MCMC_HH
