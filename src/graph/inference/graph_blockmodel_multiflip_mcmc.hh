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

#ifndef GRAPH_BLOCKMODEL_MULTIFLIP_MCMC_HH
#define GRAPH_BLOCKMODEL_MULTIFLIP_MCMC_HH

#include "config.h"

#include <vector>
#include <algorithm>

#include "graph_tool.hh"
#include "graph_state.hh"
#include "graph_blockmodel_util.hh"
#include <boost/mpl/vector.hpp>

namespace graph_tool
{
using namespace boost;
using namespace std;

typedef vprop_map_t<int64_t>::type wmap_t;

#define MCMC_BLOCK_STATE_params(State)                                         \
    ((__class__,&, mpl::vector<python::object>, 1))                            \
    ((state, &, State&, 0))                                                    \
    ((E,, size_t, 0))                                                          \
    ((beta,, double, 0))                                                       \
    ((c,, double, 0))                                                          \
    ((d,, double, 0))                                                          \
    ((a,, double, 0))                                                          \
    ((w,, wmap_t, 0))                                                          \
    ((entropy_args,, entropy_args_t, 0))                                       \
    ((mproposals, & ,std::vector<size_t>&, 0))                                 \
    ((maccept, & ,std::vector<size_t>&, 0))                                    \
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
            _vpos(num_vertices(_g)),
            _mprobs(num_vertices(_state._g) + 1),
            _sequential(false),
            _deterministic(false)
        {
            _state.init_mcmc(_c,
                             (_entropy_args.partition_dl ||
                              _entropy_args.degree_dl ||
                              _entropy_args.edges_dl));

            for (auto v : vertices_range(_state._g))
            {
                if (_state._vweight[v] > 0)
                    _vlist.push_back(v);
            }

            std::sort(_vlist.begin(), _vlist.end(),
                      [&](auto u, auto v){ return this->_w[u] < this->_w[v]; });

            size_t i = 0;
            for (auto v : _vlist)
                _vpos[v] = i++;

            _v = _vlist.front();
            _m = 1;

            std::vector<size_t> ms;
            std::vector<double> probs;
            for (size_t m = 1; m < _mprobs.size(); ++m)
            {
                ms.push_back(m);
                probs.push_back(std::pow(m, -_a));
                _mprobs[m] = _mprobs[m - 1] + probs.back();
            }
            _msampler = Sampler<size_t, mpl::false_>(ms, probs);
        }

        typename state_t::g_t& _g;

        std::vector<size_t> _vlist;
        std::vector<size_t> _vpos;

        std::vector<double> _mprobs;
        Sampler<size_t, mpl::false_> _msampler;

        size_t _v;
        size_t _m;

        bool _sequential;
        bool _deterministic;

        size_t node_state(size_t v)
        {
            return _state._b[v];
        }

        bool skip_node(size_t)
        {
            _m--;
            return _m > 0;
        }

        size_t node_weight(size_t)
        {
            return _m;
        }


        template <class F>
        void iter_vs(size_t v, size_t m, F&& f)
        {
            auto iter = _vlist.begin() + _vpos[v];
            auto end = iter + m;
            for (; iter != end; ++iter)
            {
                if (!f(*iter))
                    break;
            }
        }

        template <class RNG>
        size_t move_proposal(size_t v, RNG& rng)
        {
            auto& b = _state._b;
            size_t r = b[v];

            auto m = _msampler.sample(rng);

            size_t x1 = _vpos[v];
            size_t x2 = x1;
            for (; x2 < _vlist.size() - 1; ++x2)
            {
                if (size_t(b[_vlist[x2 + 1]]) != r)
                    break;
                if (x2 - x1 + 1 == m)
                    break;
            }
            for (; x1 > 0; --x1)
            {
                if (size_t(b[_vlist[x1 - 1]]) != r)
                    break;
                if (x2 - x1 + 1 == m)
                    break;
            }

            _v = _vlist[x1];
            _m = x2 - x1 + 1;

            std::uniform_int_distribution<size_t> dis(0, _m - 1);
            size_t u = _vlist[_vpos[_v] + dis(rng)];
            size_t s = _state.sample_block(u, _c, _d, rng);

            if (r == s)
                return null_group;

            if (!_state.allow_move(r, s))
                return null_group;

            if (!_allow_vacate)
            {
                size_t vw = 0;
                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            vw += this->_state._vweight[u];
                            return true;
                        });
                if (vw == size_t(_state._wr[r]))
                    return null_group;
            }

            if (_state._wr[s] == 0)
                _state._bclabel[s] = _state._bclabel[r];
            else
                _mproposals[_m]++;

            return s;
        }

        double log_pm(size_t m)
        {
            return - _a * log(m) - log(_mprobs.back());
        }

        double log_pm_cum(size_t m)
        {
            assert(m >= 1);
            return log(_mprobs.back() - _mprobs[m - 1]) - log(_mprobs.back());
        }


        double proposal_prob(size_t v, size_t m, size_t s)
        {
            if (std::isinf(_beta))
                return 0;

            double L = 0;

            auto begin = _vpos[v];
            auto end = begin + m;

            auto& b = _state._b;

            if (end == _vlist.size() || b[v] != b[_vlist[end]])
            {
                if (begin == 0 || b[v] != b[_vlist[begin - 1]])
                    L = log_pm_cum(m) + log(m);
                else
                    L = log_pm(m) + log(m);
            }
            else
            {
                L = log_pm(m);
            }

            double p = 0;
            iter_vs(v, m,
                    [&](auto u)
                    {
                        p += _state.get_move_prob(u, this->_state._b[u], s,
                                                  this->_c, this->_d, false);
                        return true;
                    });
            L += log(p) - log(m);
            return L;
        }

        std::tuple<double, double>
        virtual_move_dS(size_t, size_t nr)
        {
            double dS = 0, a = 0;
            double pf = 0, pb = 0;

            size_t r = _state._b[_v];

            if (_m == 1)
            {
                dS = _state.virtual_move(_v, r, nr, _entropy_args);
                pf += log(_state.get_move_prob(_v, r, nr, _c, _d, false));
                pb += log(_state.get_move_prob(_v, nr, r, _c, _d, true));
            }
            else
            {
                _state._egroups_enabled = false;

                pf += proposal_prob(_v, _m, nr);

                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            dS += this->_state.virtual_move(u, r, nr,
                                                            this->_entropy_args);
                            this->_state.move_vertex(u, nr);
                            return true;
                        });

                pb += proposal_prob(_v, _m, r);

                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            this->_state.move_vertex(u, r);
                            return true;
                        });

                _state._egroups_enabled = true;
            }

            a -= pf;
            a += pb;

            return std::make_tuple(dS, a);
        }

        void perform_move(size_t, size_t nr)
        {
            if (_state._wr[nr] > 0)
                _maccept[_m]++;

            iter_vs(_v, _m,
                    [&](auto u)
                    {
                        this->_state.move_vertex(u, nr);
                        return true;
                    });
        }
    };
};


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_MCMC_HH
