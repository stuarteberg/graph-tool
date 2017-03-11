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
    ((q,, double, 0))                                                          \
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
            _rpos(num_vertices(_state._bg)),
            _empty_pos(num_vertices(_state._bg)),
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
        std::vector<size_t> _rs;

        std::vector<size_t> _rlist;
        std::vector<size_t> _rpos;
        std::vector<size_t> _nr;
        std::vector<size_t> _old_r;
        std::vector<size_t> _empty_list;
        std::vector<size_t> _empty_pos;

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
        void iter_vs_break(size_t v, size_t m, F&& f)
        {
            auto iter = _vlist.begin() + _vpos[v];
            auto end = iter + m;
            for (; iter != end; ++iter) { if (!f(*iter)) break; }
        }

        template <class F>
        void iter_vs(size_t v, size_t m, F&& f)
        {
            iter_vs_break(v, m, [&](auto u){ f(u); return true; });
        }

        template <class RNG>
        size_t move_proposal(size_t v, RNG& rng)
        {
            _v = v;
            _m = std::min(_msampler.sample(rng),
                          _vlist.size() - _vpos[_v]);

            _rs.clear();
            std::bernoulli_distribution rsample(_q);
            if (rsample(rng))
            {
                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            size_t s = this->_state.sample_block(u, _c, _d, rng);
                            this->_rs.push_back(s);
                        });
            }
            else
            {
                std::uniform_int_distribution<size_t> dis(0, _m - 1);
                size_t u = _vlist[_vpos[_v] + dis(rng)];
                size_t s = _state.sample_block(u, _c, _d, rng);
                iter_vs(_v, _m, [&](auto) {_rs.push_back(s);});
            }

            bool noop = true;
            auto s_iter = _rs.begin();
            auto t = _state._bclabel[_state._b[_v]];
            iter_vs_break(_v, _m,
                          [&](auto u)
                          {
                              size_t r = this->_state._b[u];
                              auto s = *s_iter++;
                              if (r != s)
                                  noop = false;
                              if (!this->_state.allow_move(r, s) ||
                                  this->_state._bclabel[r] != t)
                              {
                                  noop = true;
                                  return false;
                              }
                              return true;
                          });
            if (noop)
                return null_group;

            if (!_allow_vacate)
            {
                _nr.clear();
                _rlist.clear();
                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            size_t r = this->_state._b[u];
                            if (!has_element(this->_rlist, this->_rpos, r))
                            {
                                add_element(this->_rlist, this->_rpos, r);
                                this->_nr.push_back(0);
                            }
                            this->_nr[this->_rpos[r]] +=
                                this->_state._vweight[u];
                        });
                for (size_t i = 0; i < _nr.size(); ++i)
                {
                    if (_nr[i] == size_t(_state._wr[_rlist[i]]))
                        return null_group;
                }
            }

            s_iter = _rs.begin();
            iter_vs(_v, _m,
                    [&](auto u)
                    {
                        size_t r = this->_state._b[u];
                        size_t s = *s_iter++;
                        if (this->_state._wr[s] == 0)
                            this->_state._bclabel[s] = this->_state._bclabel[r];
                    });

            _mproposals[_m]++;

            return _rs.front();
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


        double proposal_prob(size_t v, size_t m, const std::vector<size_t>& ss)
        {
            if (std::isinf(_beta))
                return 0;

            double L = 0;

            auto begin = _vpos[v];
            auto end = begin + m;

            if (end == _vlist.size())
                L = log_pm_cum(m);
            else
                L = log_pm(m);

            double Lp = 0, p = 0;
            auto s_iter = ss.begin();
            size_t r = *s_iter;
            bool uniform = true;
            size_t nempty = 0;
            _empty_list.clear();
            iter_vs(v, m,
                    [&](auto u)
                    {
                        auto s = *s_iter++;
                        double p_u = _state.get_move_prob(u, this->_state._b[u],
                                                          s, this->_c, this->_d,
                                                          false);
                        if (this->_state._wr[s] == 0)
                        {
                            nempty++;
                            if (has_element(this->_empty_list,
                                            this->_empty_pos,
                                            s))
                            {
                                add_element(this->_empty_list,
                                            this->_empty_pos,
                                            s);
                            }
                        }
                        p += p_u;
                        Lp += log(p_u);
                        if (s != r)
                            uniform = false;
                    });

            if (nempty > 0)
                Lp += lbinom_fast(_state._empty_blocks.size(),
                                  _empty_list.size())
                    - nempty * log(_state._empty_blocks.size());

            double l1 = log(_q) + Lp;
            if (uniform)
            {
                double l2 = log1p(-_q) + log(p) - log(m);
                double x = std::max(l1, l2);
                double y = std::min(l1, l2);
                L += x + log1p(exp(y - x));
            }
            else
            {
                L += l1;
            }
            return L;
        }

        std::tuple<double, double>
        virtual_move_dS(size_t, size_t)
        {
            double dS = 0, a = 0;
            double pf = 0, pb = 0;

            if (_m == 1)
            {
                size_t r = _state._b[_v];
                size_t nr = _rs.front();
                dS = _state.virtual_move(_v, r, nr, _entropy_args);
                pf += log(_state.get_move_prob(_v, r, nr, _c, _d, false));
                pb += log(_state.get_move_prob(_v, nr, r, _c, _d, true));
            }
            else
            {
                _state._egroups_enabled = false;

                pf += proposal_prob(_v, _m, _rs);

                _old_r.clear();
                auto s_iter = _rs.begin();
                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            size_t r = this->_state._b[u];
                            size_t nr = *s_iter++;
                            this->_old_r.push_back(r);
                            dS += this->_state.virtual_move(u, r, nr,
                                                            this->_entropy_args);
                            this->_state.move_vertex(u, nr);
                        });

                pb += proposal_prob(_v, _m, _old_r);

                auto r_iter = _old_r.begin();
                iter_vs(_v, _m,
                        [&](auto u)
                        {
                            this->_state.move_vertex(u, *r_iter++);
                        });

                _state._egroups_enabled = true;
            }

            a -= pf;
            a += pb;

            return std::make_tuple(dS, a);
        }

        void perform_move(size_t, size_t)
        {
            auto r_iter = _rs.begin();
            iter_vs(_v, _m,
                    [&](auto u)
                    {
                        this->_state.move_vertex(u, *r_iter++);
                    });
            _maccept[_m]++;
        }
    };
};


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_MCMC_HH
