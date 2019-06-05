// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2019 Tiago de Paula Peixoto <tiago@skewed.de>
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

#ifndef GRAPH_BLOCKMODEL_MULTILEVEL_MCMC_HH
#define GRAPH_BLOCKMODEL_MULTILEVEL_MCMC_HH

#include "config.h"

#include <vector>
#include <algorithm>

#include "graph_tool.hh"
#include "../support/graph_state.hh"
#include "graph_blockmodel_util.hh"
#include <boost/mpl/vector.hpp>

namespace graph_tool
{
using namespace boost;
using namespace std;

#define MCMC_BLOCK_STATE_params(State)                                         \
    ((__class__,&, mpl::vector<python::object>, 1))                            \
    ((state, &, State&, 0))                                                    \
    ((beta,, double, 0))                                                       \
    ((c,, double, 0))                                                          \
    ((a,, double, 0))                                                          \
    ((entropy_args,, entropy_args_t, 0))                                       \
    ((allow_vacate,, bool, 0))                                                 \
    ((verbose,, bool, 0))                                                      \
    ((niter,, size_t, 0))


template <class State>
struct MCMC
{
    GEN_STATE_BASE(MCMCBlockStateBase, MCMC_BLOCK_STATE_params(State))

    template <class... Ts>
    class MCMCBlockState
        : public MCMCBlockStateBase<Ts...>
    {
    public:
        GET_PARAMS_USING(MCMCBlockStateBase<Ts...>,
                         MCMC_BLOCK_STATE_params(State))
        GET_PARAMS_TYPEDEF(Ts, MCMC_BLOCK_STATE_params(State))

        template <class... ATs,
                  typename std::enable_if_t<sizeof...(ATs) ==
                                            sizeof...(Ts)>* = nullptr>
        MCMCBlockState(ATs&&... as)
           : MCMCBlockStateBase<Ts...>(as...),
            _g(_state._g),
            _groups(num_vertices(_state._bg)),
            _vpos(get(vertex_index_t(), _state._g),
                  num_vertices(_state._g)),
            _rpos(get(vertex_index_t(), _state._bg),
                  num_vertices(_state._bg)),
            _clabel(get(vertex_index_t(), _state._g),
                    num_vertices(_state._g))
        {
            _state.init_mcmc(_c,
                             (_entropy_args.partition_dl ||
                              _entropy_args.degree_dl ||
                              _entropy_args.edges_dl));
            for (auto v : vertices_range(_state._g))
            {
                if (_state._vweight[v] > 0)
                {
                    auto r = _state._b[v];
                    if (_groups[r].empty())
                        add_element(_rlist, _rpos, r);
                    add_element(_groups[r], _vpos, v);
                }
            }
        }

        typename state_t::g_t& _g;

        std::vector<std::vector<size_t>> _groups;
        typename vprop_map_t<size_t>::type::unchecked_t _vpos;
        typename vprop_map_t<size_t>::type::unchecked_t _rpos;

        std::vector<size_t> _rlist;

        typename vprop_map_t<size_t>::type::unchecked_t _blabel;

        std::vector<gt_hash_map<size_t, size_t>> _mrpc;
        size_t get_mrpc(size_t r, size_t c)
        {
            auto& mrp = _mrpc[r];
            auto iter = mrp.find(c);
            if (iter != mrp.end())
                return iter->second;
            return 0;
        }

        std::vector<size_t> seed_rs;

        typename vprop_map_t<size_t>::type::unchecked_t _rtree;
        std::vector<std::pair<size_t, size_t>> _merges;

        template <class RNG>
        size_t sample_merge(size_t r, RNG& rng)
        {
            if (_state._mrp[r] + _state._mrm[r] == 0)
                return uniform_sample(_rclist[_blabel[r]], rng);

            const auto& e = _egroups.sample_edge(r, rng);
            auto t = _state._b[target(e, _g)];
            if (t == r)
                t = _b[source(e, _g)];

            double p_rand = 0;
            if (c > 0)
            {
                size_t B = _rlist.size();
                auto etr = _state.get_mrs(t, r);
                if (graph_tool::is_directed(_g))
                    etr += _state.get_mrs(r, t);
                p_rand = c * B / double(get_mrpc(t, _blabel[r]) - etr + c * B);
            }

            size_t s;

            std::uniform_real_distribution<> rdist;
            if (prand > 0 && rdist(rng) < p_rand)
            {
                s = uniform_sample(_rclist[_blabel[r]], rng);
            }
            else
            {
                do
                {
                    const auto& e = _egroups.sample_edge(t, rng);
                    auto s = _state._b[target(e, _g)];
                    if (s == t)
                        s = _b[source(e, _g)];
                }
                while (s == r || _blabel[s] != _blabel[r]);
            }

            return s;
        }

        template <class RNG>
        void build_tree(RNG& rng)
        {
            size_t m = sample_m(_groups[r].size(), rng);

            if (!_allow_vacate && _groups[r].size() == m)
                return null_group;

            assert(m <= _groups[r].size());

            _nproposals += m;
            _vs.clear();

            while (_vs.size() < m)
            {
                size_t v = uniform_sample(_groups[r], rng);
                _vs.push_back(v);
                remove_element(_groups[r], _vpos, v);
            }

            for (auto v : _vs)
                add_element(_groups[r], _vpos, v);

            size_t v = uniform_sample(_vs, rng);
            size_t s = _state.sample_block(v, _c, _d, rng);

            if (s >= _groups.size())
            {
                _groups.resize(s + 1);
                _rpos.resize(s + 1);
            }

            if (!_state.allow_move(r, s) || s == r)
                return null_group;

            if (!_groups[s].empty() || _groups[r].size() > m)
                _mproposals[m]++;
            return s;
        }

        std::tuple<double, double>
        virtual_move_dS(size_t r, size_t nr)
        {
            double dS = 0, a = 0;
            size_t m = _vs.size();

            a -= log_pm(m, _groups[r].size());
            a += log_pm(m, _groups[nr].size() + m);
            a -= -lbinom(_groups[r].size(), m);
            a += -lbinom(_groups[nr].size() + m, m);

            int B = _vlist.size();
            a -= -log(B);
            if (_groups[r].size() == m)
                B--;
            if (_groups[nr].empty())
                B++;
            a += -log(B);

            if (m == 1)
            {
                auto v = _vs.front();
                dS = _state.virtual_move(v, r, nr, _entropy_args);
                double pf = _state.get_move_prob(v, r, nr, _c, _d, false);
                double pb = _state.get_move_prob(v, nr, r, _c, _d, true);
                a += log(pb) - log(pf);
            }
            else
            {
                _state._egroups_enabled = false;
                double pf = 0, pb = 0;
                for (auto v : _vs)
                    pf += _state.get_move_prob(v, r, nr, _c, _d, false);
                pf /= m;
                for (auto v : _vs)
                {
                    dS += _state.virtual_move(v, r, nr, _entropy_args);
                    _state.move_vertex(v, nr);
                }
                for (auto v : _vs)
                    pb += _state.get_move_prob(v, nr, r, _c, _d, false);
                pb /= m;
                a += log(pb) - log(pf);
                for (auto v : _vs)
                    _state.move_vertex(v, r);
                _state._egroups_enabled = true;
            }
            return std::make_tuple(dS, a);
        }

        void perform_move(size_t r, size_t nr)
        {
            if (!_groups[nr].empty() || _groups[r].size() > _vs.size())
                _maccept[_vs.size()]++;
            if (_state._wr[nr] == 0)
                add_element(_vlist, _rpos, nr);
            for (auto v : _vs)
            {
                _state.move_vertex(v, nr);
                remove_element(_groups[r], _vpos, v);
                add_element(_groups[nr], _vpos, v);
            }
            if (_state._wr[r] == 0)
                remove_element(_vlist, _rpos, r);
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
            if (_nproposals < _niter * _N)
                return std::numeric_limits<size_t>::max();
            return 0;
        }

        void step(size_t, size_t)
        {
        }
    };
};


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_MCMC_HH
