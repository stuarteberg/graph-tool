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

#ifndef GRAPH_BLOCKMODEL_UNCERTAIN_HH
#define GRAPH_BLOCKMODEL_UNCERTAIN_HH

#include "config.h"

#include <vector>

#include "../support/graph_state.hh"
#include "../blockmodel/graph_blockmodel_util.hh"
#include "graph_blockmodel_sample_edge.hh"

namespace graph_tool
{
using namespace boost;
using namespace std;

#define UNCERTAIN_STATE_params                                                 \
    ((__class__,&, mpl::vector<python::object>, 1))                            \
    ((g, &, all_graph_views, 1))                                               \
    ((q,, eprop_map_t<double>::type, 0))                                       \
    ((q_default,, double, 0))                                                  \
    ((S_const,, double, 0))                                                    \
    ((aE,, double, 0))                                                         \
    ((E_prior,, bool, 0))                                                      \
    ((self_loops,, bool, 0))

template <class BlockState>
struct Uncertain
{
    GEN_STATE_BASE(UncertainStateBase, UNCERTAIN_STATE_params)

    template <class... Ts>
    class UncertainState
        : public UncertainStateBase<Ts...>
    {
    public:
        GET_PARAMS_USING(UncertainStateBase<Ts...>,
                         UNCERTAIN_STATE_params)
        GET_PARAMS_TYPEDEF(Ts, UNCERTAIN_STATE_params)

        template <class... ATs,
                  typename std::enable_if_t<sizeof...(ATs) == sizeof...(Ts)>* = nullptr>
        UncertainState(BlockState& block_state, ATs&&... args)
            : UncertainStateBase<Ts...>(std::forward<ATs>(args)...),
              _block_state(block_state)
        {
            _u_edges.resize(num_vertices(_u));
            for (auto e : edges_range(_u))
            {
                get_u_edge<true>(source(e, _u), target(e, _u)) = e;
                _E += _eweight[e];
            }
            _edges.resize(num_vertices(_g));
            for (auto e : edges_range(_g))
                get_edge<true>(source(e, _g), target(e, _g)) = e;
            _block_state.enable_partition_stats();
        }

        UncertainState(const UncertainState& other)
            : UncertainStateBase<Ts...>(static_cast<const UncertainStateBase<Ts...>&>(other)),
              _block_state(other._block_state),
              _u_edges(other._u_edges),
              _edges(other._edges),
              _pe(other._pe),
              _E(other._E)
        {
            _block_state.enable_partition_stats();
        }

        BlockState& _block_state;
        typename BlockState::g_t& _u = _block_state._g;
        typename BlockState::eweight_t& _eweight = _block_state._eweight;
        GraphInterface::edge_t _null_edge;
        std::vector<double> _recs;

        std::vector<gt_hash_map<size_t, GraphInterface::edge_t>> _u_edges;
        std::vector<gt_hash_map<size_t, GraphInterface::edge_t>> _edges;

        SBMEdgeSampler<typename BlockState::g_t> _edge_sampler;

        double _pe = log(_aE);
        size_t _E = 0;

        void set_q_default(double q_default)
        {
            _q_default = q_default;
        }

        void set_S_const(double S_const)
        {
            _S_const = S_const;
        }

        template <bool insert, class Graph, class Elist>
        auto& _get_edge(size_t u, size_t v, Graph& g, Elist& edges)
        {
            if (!graph_tool::is_directed(g) && u > v)
                std::swap(u, v);
            auto& qe = edges[u];
            if (insert)
                return qe[v];
            auto iter = qe.find(v);
            if (iter != qe.end())
                return iter->second;
            return _null_edge;
        }

        template <bool insert=false>
        auto& get_u_edge(size_t u, size_t v)
        {
            return _get_edge<insert>(u, v, _u, _u_edges);
        }

        template <bool insert=false>
        auto& get_edge(size_t u, size_t v)
        {
            return _get_edge<insert>(u, v, _g, _edges);
        }

        double entropy()
        {
            double S = 0;
            for (auto m : edges_range(_g))
            {
                double q_e = _q[m];
                if (q_e == std::numeric_limits<double>::infinity())
                    continue;
                auto& e = get_u_edge<false>(source(m, _g), target(m, _g));
                if (e != _null_edge && _eweight[e] > 0 &&
                    (_self_loops || (source(e, _u) != target(e, _u))))
                    S += q_e;
            }

            for (auto e : edges_range(_u))
            {
                auto& m = get_edge<false>(source(e, _u), target(e, _u));
                if (m != _null_edge || _eweight[e] == 0 ||
                    (!_self_loops && source(m, _g) == target(m, _g)))
                    continue;
                if (_q_default == std::numeric_limits<double>::infinity())
                    continue;
                S += _q_default;
            }

            if (_E_prior)
                S += _E * _pe - lgamma_fast(_E + 1) - exp(_pe);
            S += _S_const;

            return -S;
        }

        double remove_edge_dS(size_t u, size_t v, entropy_args_t ea)
        {
            auto& e = get_u_edge(u, v);
            double dS = _block_state.template modify_edge_dS<false>(source(e, _u),
                                                                    target(e, _u),
                                                                    e, _recs, ea);
            if (_E_prior)
            {
                dS += _pe;
                dS += lgamma_fast(_E) - lgamma_fast(_E + 1);
            }
            if (_eweight[e] == 1 && (_self_loops || u != v))
            {
                auto& m = get_edge<false>(u, v);
                double q_e = (m == _null_edge) ? _q_default : _q[m];
                dS += q_e;
            }

            return dS;
        }

        double add_edge_dS(size_t u, size_t v, entropy_args_t ea)
        {
            auto& e = get_u_edge(u, v);
            double dS = _block_state.template modify_edge_dS<true>(u, v, e,
                                                                   _recs, ea);
            if (_E_prior)
            {
                dS -= _pe;
                dS += lgamma_fast(_E + 2) - lgamma_fast(_E + 1);
            }
            if ((e == _null_edge || _eweight[e] == 0) && (_self_loops || u != v))
            {
                auto& m = get_edge<false>(u, v);
                double q_e = (m == _null_edge) ? _q_default : _q[m];
                dS -= q_e;
            }
            return dS;
        }

        void remove_edge(size_t u, size_t v)
        {
            auto& e = get_u_edge(u, v);
            _block_state.template modify_edge<false>(u, v, e,
                                                     _recs);
            _E--;
        }

        void add_edge(size_t u, size_t v)
        {
            auto& e = get_u_edge<true>(u, v);
            _block_state.template modify_edge<true>(u, v, e,
                                                    _recs);
            _E++;
        }

    };
};

template <class State>
double get_edge_prob(State& state, size_t u, size_t v, entropy_args_t ea,
                     double epsilon)
{
    auto e = state.get_u_edge(u, v);
    size_t ew = 0;
    if (e != state._null_edge)
        ew = state._eweight[e];

    for (size_t i = 0; i < ew; ++i)
        state.remove_edge(u, v);

    double S = 0, Z = 1;
    double delta = 1. + epsilon;
    size_t ne = 0;

    while (delta > epsilon || ne < 2)
    {
        double dS = state.add_edge_dS(u, v, ea);
        state.add_edge(u, v);
        S += dS;
        ne++;
        double dZ = exp(-S);
        delta = dZ/Z;
        Z += dZ;
    }

    double L = log1p(-1./Z);

    for (int i = 0; i < int(ne - ew); ++i)
        state.remove_edge(u, v);
    for (int i = 0; i < int(ew - ne); ++i)
        state.add_edge(u, v);

    return L;
}

template <class State>
void get_edges_prob(State& state, python::object edges, python::object probs,
                    entropy_args_t ea, double epsilon)
{
    multi_array_ref<uint64_t,2> es = get_array<uint64_t,2>(edges);
    multi_array_ref<double,1> eprobs = get_array<double,1>(probs);
    for (size_t i = 0; i < eprobs.shape()[0]; ++i)
        eprobs[i] = get_edge_prob(state, es[i][0], es[i][1], ea, epsilon);
}

} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_UNCERTAIN_HH
