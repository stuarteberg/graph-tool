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

#ifndef GRAPH_BLOCKMODEL_ENTROPY_HH
#define GRAPH_BLOCKMODEL_ENTROPY_HH

#include "../support/util.hh"

namespace graph_tool
{

// ====================
// Entropy calculation
// ====================

enum deg_dl_kind
{
    ENT,
    UNIFORM,
    DIST
};

struct entropy_args_t
{
    bool dense;
    bool multigraph;
    bool exact;
    bool adjacency;
    bool recs;
    bool deg_entropy;
    bool partition_dl;
    bool degree_dl;
    deg_dl_kind  degree_dl_kind;
    bool edges_dl;
    bool recs_dl;
    double beta_dl;
};

// Sparse entropy terms
// ====================

// exact microcanonical deg-corr entropy
template <class Graph>
inline double eterm_exact(size_t r, size_t s, size_t mrs, const Graph& g)
{
    double val = lgamma_fast(mrs + 1);

    if (graph_tool::is_directed(g) || r != s)
    {
        return -val;
    }
    else
    {
#ifndef __clang__
        constexpr double log_2 = log(2);
#else
        const double log_2 = log(2);
#endif
        return -val - mrs * log_2;
    }
}

template <class Graph>
inline double vterm_exact(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                          const Graph& g)
{
    if (deg_corr)
    {
        if (graph_tool::is_directed(g))
            return lgamma_fast(mrp + 1) + lgamma_fast(mrm + 1);
        else
            return lgamma_fast(mrp + 1);
    }
    else
    {
        if (graph_tool::is_directed(g))
            return (mrp + mrm) * safelog(wr);
        else
            return mrp * safelog(wr);
    }
}

// "edge" term of the entropy
template <class Graph>
inline double eterm(size_t r, size_t s, size_t mrs, const Graph& g)
{
    if (!graph_tool::is_directed(g) && r == s)
        mrs *= 2;

    double val = xlogx(mrs);

    if (graph_tool::is_directed(g) || r != s)
        return -val;
    else
        return -val / 2;
}

// "vertex" term of the entropy
template <class Graph>
inline double vterm(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                    Graph& g)
{
    double one = 0.5;

    if (graph_tool::is_directed(g))
        one = 1;

    if (deg_corr)
        return one * (xlogx(mrm) + xlogx(mrp));
    else
        return one * (mrm * safelog(wr) + mrp * safelog(wr));
}



// Dense entropy
// =============

// "edge" term of the entropy
template <class Graph>
inline double eterm_dense(size_t r, size_t s, int ers, double wr_r,
                          double wr_s, bool multigraph, const Graph& g)
{
    // we should not use integers here, since they may overflow
    double nrns;

    if (ers == 0)
        return 0.;

    assert(wr_r + wr_s > 0);

    if (r != s || graph_tool::is_directed(g))
    {
        nrns = wr_r * wr_s;
    }
    else
    {
        if (multigraph)
            nrns = (wr_r * (wr_r + 1)) / 2;
        else
            nrns = (wr_r * (wr_r - 1)) / 2;
    }

    double S;
    if (multigraph)
        S = lbinom(nrns + ers - 1, ers); // do not use lbinom_fast!
    else
        S = lbinom(nrns, ers);
    return S;
}

// Edges description length
template <class Graph>
double get_edges_dl(size_t B, size_t E, Graph& g)
{
    size_t NB = (graph_tool::is_directed(g)) ? B * B : (B * (B + 1)) / 2;
    return lbinom(NB + E - 1, E);
}

} // namespace graph_tool

#endif // GRAPH_BLOCKMODEL_ENTROPY_HH
