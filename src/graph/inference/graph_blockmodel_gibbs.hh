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

#ifndef GRAPH_BLOCKMODEL_GIBBS_HH
#define GRAPH_BLOCKMODEL_GIBBS_HH

#include "config.h"

#include <vector>

#include "graph_tool.hh"
#include "graph_state.hh"
#include "graph_blockmodel_util.hh"
#include <boost/mpl/vector.hpp>

namespace graph_tool
{
using namespace boost;
using namespace std;

#define GIBBS_BLOCK_STATE_params(State)                                        \
    ((__class__,&, mpl::vector<python::object>, 1))                            \
    ((state, &, State&, 0))                                                    \
    ((E,, size_t, 0))                                                          \
    ((vlist,&, std::vector<size_t>&, 0))                                       \
    ((block_list,&, std::vector<size_t>&, 0))                                  \
    ((beta,, double, 0))                                                       \
    ((multigraph,, bool, 0))                                                   \
    ((dense,, bool, 0))                                                        \
    ((partition_dl,, bool, 0))                                                 \
    ((degree_dl,, bool, 0))                                                    \
    ((edges_dl,, bool, 0))                                                     \
    ((allow_empty,, bool, 0))                                                  \
    ((parallel,, bool, 0))                                                     \
    ((sequential,, bool, 0))                                                   \
    ((verbose,, bool, 0))                                                      \
    ((niter,, size_t, 0))


template <class State, template <class Graph> class MEntries = EntrySet>
struct Gibbs
{
    GEN_STATE_BASE(GibbsBlockStateBase, GIBBS_BLOCK_STATE_params(State))

    template <class... Ts>
    class GibbsBlockState
        : public GibbsBlockStateBase<Ts...>
    {
    public:
        GET_PARAMS_USING(GibbsBlockStateBase<Ts...>,
                         GIBBS_BLOCK_STATE_params(State))
        GET_PARAMS_TYPEDEF(Ts, GIBBS_BLOCK_STATE_params(State))

        template <class... ATs,
                  typename std::enable_if_t<sizeof...(ATs) ==
                                            sizeof...(Ts)>* = nullptr>
        GibbsBlockState(ATs&&... as)
           : GibbsBlockStateBase<Ts...>(as...),
            _g(_state._g),
            _m_entries(num_vertices(_state._bg)),
            _B(num_vertices(_state._bg)),
            _tglobal(std::make_shared<tglobal_t>())
        {
            _state.init_mcmc(numeric_limits<double>::infinity(),
                             _partition_dl || _degree_dl || _edges_dl);

            size_t C = 0;
            for (auto v : vertices_range(_g))
                C = std::max(C, _state._bclabel[_state._b[v]]);
            C++;

            auto& moves = _tglobal->moves;
            auto& weights = _tglobal->weights;
            auto& empty = _tglobal->empty;

            empty.resize(C);

            for(size_t c = 0; c < C; ++c)
            {
                moves.push_back(_B + c);
                weights.push_back(0);
            }

            for (auto r : _block_list)
            {
                if (_state._wr[r] > 0)
                {
                    moves.push_back(r);
                    weights.push_back(1);
                }
                else
                {
                    auto c = _state._bclabel[r];
                    empty[c].push_back(r);
                }
            }

            for (size_t c = 0; c < C; ++c)
                weights[c] = empty[c].size();
        }

        typename state_t::g_t& _g;
        MEntries<typename state_t::g_t> _m_entries;

        size_t _B;

        struct tglobal_t
        {
            vector<size_t> moves;
            vector<size_t> weights;
            vector<vector<size_t>> empty;
        };

        std::shared_ptr<tglobal_t> _tglobal;

        auto& get_moves(size_t) { return _tglobal->moves; }
        auto& get_weights(size_t) { return _tglobal->weights; }

        size_t node_state(size_t v)
        {
            return _state._b[v];
        }

        double virtual_move_dS(size_t v, size_t nr)
        {
            if (nr >= _B)
            {
                if (_allow_empty)
                {
                    auto& empty = _tglobal->empty;
                    auto c = nr - _B;
                    if (empty[c].empty())
                        return numeric_limits<double>::infinity();
                    nr = empty[c].back();
                }
                else
                {
                    return numeric_limits<double>::infinity();
                }
            }
            size_t r = _state._b[v];
            if (_state._bclabel[r] != _state._bclabel[nr])
                return numeric_limits<double>::infinity();
            return _state.virtual_move(v, nr, _dense, _multigraph,
                                       _partition_dl, _degree_dl, _edges_dl,
                                       _m_entries);
        }

        void perform_move(size_t v, size_t nr)
        {
            size_t r = _state._b[v];

            if (r == nr)
                return;

            auto& moves = _tglobal->moves;
            auto& weights = _tglobal->weights;
            auto& empty = _tglobal->empty;

            if (nr >= _B)
            {
                auto c = nr - _B;
                assert(!empty[c].empty());
                nr = empty[c].back();
                empty[c].pop_back();
                weights[c]--;
                moves.push_back(nr);
                weights.push_back(1);
            }

            _state.move_vertex(v, nr);

            if (_state._wr[r] == 0)
            {
                auto c = _state._bclabel[r];
                empty[c].push_back(r);
                weights[c]++;

                auto iter = find(moves.begin(), moves.end(), r);
                assert(iter != moves.end());
                size_t pos = iter - moves.begin();
                std::swap(moves[pos], moves.back());
                moves.pop_back();
                weights.pop_back();
            }
        }
    };
};


} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_GIBBS_HH
