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

#ifndef GRAPH_BLOCKMODEL_UTIL_HH
#define GRAPH_BLOCKMODEL_UTIL_HH

#include "config.h"

#include <tuple>

#include "hash_map_wrap.hh"

#include "../generation/sampler.hh"
#include "../generation/dynamic_sampler.hh"

#include "util.hh"
#include "cache.hh"

#include <boost/multi_array.hpp>

#ifdef USING_OPENMP
#include <omp.h>
#endif

namespace graph_tool
{

// ====================
// Entropy calculation
// ====================


// Sparse entropy terms
// ====================

// "edge" term of the entropy
template <class Graph>
inline double eterm(size_t r, size_t s, size_t mrs, const Graph&)
{
    if (!is_directed::apply<Graph>::type::value && r == s)
        mrs *= 2;

    double val = xlogx(mrs);

    if (is_directed::apply<Graph>::type::value || r != s)
        return -val;
    else
        return -val / 2;
}

// "vertex" term of the entropy
template <class Graph>
inline double vterm(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                    Graph&)
{
    double one = 0.5;

    if (is_directed::apply<Graph>::type::value)
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
                          double wr_s, bool multigraph, const Graph&)
{
    // we should not use integers here, since they may overflow
    double nrns;

    if (ers == 0)
        return 0.;

    if (r != s || is_directed::apply<Graph>::type::value)
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

// ===============
// Partition stats
// ===============

typedef vprop_map_t<std::vector<std::tuple<size_t, size_t, size_t>>>::type
    degs_map_t;

struct simple_degs_t {};

template <class Graph, class Vprop, class Eprop>
auto get_degs(size_t v, Vprop& vweight, Eprop& eweight, const simple_degs_t&,
              Graph& g)
{
    std::array<std::tuple<size_t, size_t, size_t>, 1>
        degs {{make_tuple(in_degreeS()(v, g, eweight),
                          out_degreeS()(v, g, eweight),
                          vweight[v])}};
    return degs;
}

template <class Graph, class Vprop, class Eprop>
auto& get_degs(size_t v, Vprop&, Eprop&, typename degs_map_t::unchecked_t& degs,
               Graph&)
{
    return degs[v];
}

class partition_stats_t
{
public:

    typedef gt_hash_map<pair<size_t,size_t>, int> map_t;

    template <class Graph, class Vprop, class VWprop, class Eprop, class Degs,
              class Mprop, class Vlist>
    partition_stats_t(Graph& g, Vprop& b, Vlist& vlist, size_t E, size_t B,
                      VWprop& vweight, Eprop& eweight, Degs& degs,
                      const Mprop& ignore_degree, std::vector<size_t>& bmap)
        : _bmap(bmap)
    {
        _N = 0;
        _E = E;
        _total_B = B;

        for (auto v : vlist)
        {
            auto r = get_r(b[v]);

            auto&& ks = get_degs(v, vweight, eweight, degs, g);
            if (v >= _ignore_degree.size())
                _ignore_degree.resize(v + 1, 0);
            _ignore_degree[v] = ignore_degree[v];
            for (auto& k : ks)
            {
                size_t kin = get<0>(k);
                size_t kout = get<1>(k);
                size_t n = get<2>(k);
                if (_ignore_degree[v] == 2)
                    kout = 0;
                if (_ignore_degree[v] != 1)
                {
                    _hist[r][make_pair(kin, kout)] += n;
                    _em[r] += kin * n;
                    _ep[r] += kout * n;
                }
                _total[r] += n;
                _N += n;
            }
        }

        _actual_B = 0;
        for (auto n : _total)
        {
            if (n > 0)
                _actual_B++;
        }
    }

    size_t get_r(size_t r)
    {
        constexpr size_t null =
            std::numeric_limits<size_t>::max();
        if (r >= _bmap.size())
            _bmap.resize(r + 1, null);
        size_t nr = _bmap[r];
        if (nr == null)
            nr = _bmap[r] = _hist.size();
        if (nr >= _hist.size())
        {
            _hist.resize(nr + 1);
            _total.resize(nr + 1);
            _ep.resize(nr + 1);
            _em.resize(nr + 1);
        }
        return nr;
    }

    double get_partition_dl()
    {
        double S = 0;
        S += lbinom(_actual_B + _N - 1, _N);
        S += lgamma(_N + 1);
        for (auto nr : _total)
            S -= lgamma(nr + 1);
        return S;
    }

    double get_deg_dl(bool ent, bool dl_alt, bool xi_fast)
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            if (_ep[r] + _em[r] == 0)
                continue;

            if (ent)
            {
                for (auto& k_c : _hist[r])
                {
                    double p = k_c.second / double(_total[r]);
                    S -= p * log(p) * _total[r];
                }
            }
            else
            {
                double S1 = 0;

                if (xi_fast)
                {
                    S1 += get_xi_fast(_total[r], _ep[r]);
                    S1 += get_xi_fast(_total[r], _em[r]);
                }
                else
                {
                    S1 += get_xi(_total[r], _ep[r]);
                    S1 += get_xi(_total[r], _em[r]);
                }

                S1 += lgamma(_total[r] + 1);
                for (auto& k_c : _hist[r])
                    S1 -= lgamma(k_c.second + 1);

                if (dl_alt)
                {
                    double S2 = 0;
                    S2 += lbinom(_total[r] + _ep[r] - 1, _ep[r]);
                    S2 += lbinom(_total[r] + _em[r] - 1, _em[r]);
                    S += min(S1, S2);
                }
                else
                {
                    S += S1;
                }
            }
        }
        return S;
    }

    template <class VProp>
    double get_delta_dl(size_t v, size_t r, size_t nr, VProp& vweight)
    {
        if (r == nr)
            return 0;
        r = get_r(r);
        nr = get_r(nr);
        int n = vweight[v];

        if (n == 0)
            return 0;

        double S_b = 0, S_a = 0;

        S_b += -lgamma_fast(_total[r] + 1);
        S_a += -lgamma_fast(_total[r] - n + 1);
        S_b += -lgamma_fast(_total[nr] + 1);
        S_a += -lgamma_fast(_total[nr] + n + 1);

        int dB = 0;
        if (_total[r] == n)
            dB--;
        if (_total[nr] == 0)
            dB++;

        if (dB != 0)
        {
            S_b += lbinom(_actual_B + _N - 1, _N);
            S_a += lbinom(_actual_B + dB + _N - 1, _N);
        }

        return S_a - S_b;
    }

    template <class Graph, class VProp>
    double get_delta_edges_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                              Graph&)
    {
        if (r == nr)
            return 0;

        r = get_r(r);
        nr = get_r(nr);

        double S_b = 0, S_a = 0;

        int n = vweight[v];

        if (n == 0)
            return 0;

        int dB = 0;
        if (_total[r] == n)
            dB--;
        if (_total[nr] == 0)
            dB++;

        if (dB != 0)
        {
            auto get_x = [](size_t B)
                {
                    if (is_directed::apply<Graph>::type::value)
                        return B * B;
                    else
                        return (B * (B + 1)) / 2;
                };

            S_b += lbinom(get_x(_total_B) + _E - 1, _E);
            S_a += lbinom(get_x(_total_B + dB) + _E - 1, _E);
        }

        return S_a - S_b;
    }

    template <class Graph, class VProp, class EProp, class Degs>
    double get_delta_deg_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                            EProp& eweight, Degs& degs, Graph& g)
    {
        if (r == nr || _ignore_degree[v] == 1)
            return 0;
        r = get_r(r);
        nr = get_r(nr);
        auto&& ks = get_degs(v, vweight, eweight, degs, g);
        auto* _ks = &ks;
        typename std::remove_reference<decltype(ks)>::type nks;
        if (_ignore_degree[v] == 2)
        {
            nks = ks;
            for (auto& k : nks)
                get<1>(k) = 0;
            _ks = &nks;
        }
        double dS = 0;
        dS += get_delta_deg_dl_change(v, r,  *_ks, -1);
        dS += get_delta_deg_dl_change(v, nr, *_ks, +1);
        return dS;
    }

    double get_Se (size_t s, int delta, int kin, int kout)
    {
        double S = 0;
        assert(_total[s] + delta >= 0);
        assert(_em[s] + kin >= 0);
        assert(_ep[s] + kout >= 0);
        S += get_xi_fast(_total[s] + delta, _em[s] + kin);
        S += get_xi_fast(_total[s] + delta, _ep[s] + kout);
        return S;
    };

    double get_Sr(size_t s, int delta)
    {
        assert(_total[s] + delta + 1 >= 0);
        return lgamma_fast(_total[s] + delta + 1);
    };

    double get_Sk(size_t s, pair<size_t, size_t>& deg, int delta)
    {
        int nd = 0;
        auto iter = _hist[s].find(deg);
        if (iter != _hist[s].end())
            nd = iter->second;
        assert(nd + delta >= 0);
        return -lgamma_fast(nd + delta + 1);
    };


    template <class Ks>
    double get_delta_deg_dl_change(size_t v, size_t r, Ks& ks, int diff)
    {
        double S_b = 0, S_a = 0;
        int tkin = 0, tkout = 0, n = 0;
        for (auto& k : ks)
        {
            size_t kin = get<0>(k);
            size_t kout = get<1>(k);
            int nk = get<2>(k);
            tkin += kin * nk;
            if (_ignore_degree[v] != 2)
                tkout += kout * nk;
            n += nk;

            auto deg = make_pair(kin, kout);
            S_b += get_Sk(r, deg,         0);
            S_a += get_Sk(r, deg, diff * nk);
        }

        S_b += get_Se(r,        0,           0,            0);
        S_a += get_Se(r, diff * n, diff * tkin, diff * tkout);

        S_b += get_Sr(r,        0);
        S_a += get_Sr(r, diff * n);

        return S_a - S_b;
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void change_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs,
                       int diff)
    {
        auto&& ks = get_degs(v, vweight, eweight, degs, g);
        for (auto& k : ks)
        {
            auto kin = get<0>(k);
            auto kout = get<1>(k);
            auto n = get<2>(k);
            change_k(v, r, deg_corr, n, kin, kout, diff);
        }
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void remove_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        r = get_r(r);
        change_vertex(v, r, deg_corr, g, vweight, eweight, degs, -1);
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void add_vertex(size_t v, size_t nr, bool deg_corr, Graph& g,
                    VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        nr = get_r(nr);
        change_vertex(v, nr, deg_corr, g, vweight, eweight, degs, 1);
    }

    void change_k(size_t v, size_t r, bool deg_corr, int vweight,
                  int kin, int kout, int diff)
    {
        if (_total[r] == 0 && diff * vweight > 0)
        {
            _actual_B++;
            _total_B++;
        }

        if (_total[r] == vweight && diff * vweight < 0)
        {
            _actual_B--;
            _total_B--;
        }

        _total[r] += diff * vweight;
        _N += diff * vweight;

        if (deg_corr && _ignore_degree[v] != 1)
        {
            if (_ignore_degree[v] == 2)
                kout = 0;
            auto deg = make_pair(kin, kout);
            _hist[r][deg] += diff * vweight;
            _em[r] += diff * deg.first * vweight;
            _ep[r] += diff * deg.second * vweight;
        }
    }

private:
    vector<size_t>& _bmap;
    size_t _N;
    size_t _E;
    size_t _actual_B;
    size_t _total_B;
    vector<map_t> _hist;
    vector<int> _total;
    vector<int> _ep;
    vector<int> _em;
    vector<uint8_t> _ignore_degree;
};

// ===============================
// Block moves
// ===============================

// this structure speeds up the access to the edges between given blocks, since
// we're using an adjacency list to store the block structure (it is simply an
// adjacency matrix)

template <class Graph, class BGraph>
class EMat
{
public:
    template <class Vprop, class RNG>
    EMat(Graph& g, Vprop b, BGraph& bg, RNG&)
        : _bedge(get(edge_index_t(), g), 0)
    {
        size_t B = num_vertices(bg);
        _mat.resize(boost::extents[B][B]);
        for (auto e : edges_range(bg))
        {
            _mat[source(e, bg)][target(e, bg)] = e;
            if (!is_directed::apply<BGraph>::type::value)
                _mat[target(e, bg)][source(e, bg)] = e;
        }

        auto bedge_c = _bedge.get_checked();
        for (auto e : edges_range(g))
        {
            auto r = b[source(e, g)];
            auto s = b[target(e, g)];
            bedge_c[e] = _mat[r][s];
        }
    }

    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<BGraph>::edge_descriptor edge_t;

    const auto& get_me(vertex_t r, vertex_t s) const
    {
        return _mat[r][s];
    }

    void put_me(vertex_t r, vertex_t s, const edge_t& e)
    {
        _mat[r][s] = e;
        if (!is_directed::apply<BGraph>::type::value && r != s)
            _mat[s][r] = e;
    }

    void remove_me(vertex_t r, vertex_t s, const edge_t& me, BGraph& bg,
                   bool delete_edge=false)
    {
        if (delete_edge)
        {
            _mat[r][s] = edge_t();
            if (!is_directed::apply<Graph>::type::value)
                _mat[s][r] = edge_t();
            remove_edge(me, bg);
        }
    }

    template <class Edge>
    auto& get_bedge(const Edge& e) { return _bedge[e]; }
    auto& get_bedge_map() { return _bedge; }
    const auto& get_bedge_map() const { return _bedge; }
    const auto& get_null_edge() const { return _null_edge; }

private:
    multi_array<edge_t, 2> _mat;
    typedef typename property_map_type::apply
        <edge_t,
         typename property_map<Graph,
                               edge_index_t>::type>::type bedge_t;
    typename bedge_t::unchecked_t _bedge;
    static const edge_t _null_edge;
};

template <class Graph, class BGraph>
const typename EMat<Graph, BGraph>::edge_t EMat<Graph, BGraph>::_null_edge;


template <class Key>
class perfect_hash_t
{
public:
    template <class RNG>
    perfect_hash_t(size_t N, RNG& rng)
        : _index(std::make_shared<std::vector<size_t>>())
    {
        auto& index = *_index;
        index.reserve(N);
        for (size_t i = 0; i < N; ++i)
            index.push_back(i);
        std::shuffle(index.begin(), index.end(), rng);
    }
    perfect_hash_t() {}
    size_t operator()(const Key& k) const {return (*_index)[k];}
private:
    std::shared_ptr<std::vector<size_t>> _index;
};

// this structure speeds up the access to the edges between given blocks, since
// we're using an adjacency list to store the block structure (this is like
// EMat above, but takes less space and is slower)

template <class Graph, class BGraph>
class EHash
{
public:

    template <class Vprop, class RNG>
    EHash(Graph& g, Vprop b, BGraph& bg, RNG& rng)
        : _hash_function(num_vertices(bg), rng),
          _hash(num_vertices(bg), ehash_t(0, _hash_function)),
          _bedge(get(edge_index_t(), g), 0)
    {

        for (auto e : edges_range(bg))
            put_me(source(e, bg), target(e, bg), e);

        auto bedge_c = _bedge.get_checked();
        for (auto e : edges_range(g))
        {
            auto r = b[source(e, g)];
            auto s = b[target(e, g)];
            bedge_c[e] = get_me(r, s);
        }
    }

    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<BGraph>::edge_descriptor edge_t;

    const auto& get_me(vertex_t r, vertex_t s) const
    {
        auto& map = _hash[r];
        const auto& iter = map.find(s);
        if (iter == map.end())
            return _null_edge;
        return iter->second;
    }

    void put_me(vertex_t r, vertex_t s, const edge_t& e)
    {
        assert(r < _hash.size());
        _hash[r][s] = e;
        if (!is_directed::apply<Graph>::type::value)
            _hash[s][r] = e;
    }

    void remove_me(vertex_t r, vertex_t s, const edge_t& me, BGraph& bg,
                   bool delete_edge = true)
    {
        if (delete_edge)
        {
            assert(r < _hash.size());
            _hash[r].erase(s);
            if (!is_directed::apply<Graph>::type::value)
                _hash[s].erase(r);
            remove_edge(me, bg);
        }
    }

    template <class Edge>
    auto& get_bedge(const Edge& e) { return _bedge[e]; }
    auto& get_bedge_map() { return _bedge; }
    const auto& get_bedge_map() const { return _bedge; }
    const auto& get_null_edge() const { return _null_edge; }

private:
    perfect_hash_t<vertex_t> _hash_function;
    typedef gt_hash_map<vertex_t, edge_t, perfect_hash_t<vertex_t>> ehash_t;
    std::vector<ehash_t> _hash;
    typedef typename property_map_type::apply
        <edge_t,
         typename property_map<Graph,
                               edge_index_t>::type>::type bedge_t;
    typename bedge_t::unchecked_t _bedge;
    static const edge_t _null_edge;
};

template <class Graph, class BGraph>
const typename EHash<Graph, BGraph>::edge_t EHash<Graph, BGraph>::_null_edge;

template <class Vertex, class Eprop, class Emat>
inline size_t get_mrs(Vertex r, Vertex s, const Eprop& mrs, const Emat& emat)
{
    const auto& me = emat.get_me(r, s);
    if (me != emat.get_null_edge())
        return mrs[me];
    return 0;
}


// Manage a set of block pairs and corresponding edge counts that will be
// updated

template <class Graph>
class EntrySet
{
public:
    EntrySet(size_t B)
    {
        _r_field_t.resize(B, _null);
        _nr_field_t.resize(B, _null);

        if (is_directed::apply<Graph>::type::value)
        {
            _r_field_s.resize(B, _null);
            _nr_field_s.resize(B, _null);
        }
    }

    void set_move(size_t r, size_t nr)
    {
        _rnr = make_pair(r, nr);
    }

    void insert_delta(size_t t, size_t s, int delta,
                      size_t mrs = numeric_limits<size_t>::max())
    {
        insert_delta(t, s, delta, mrs,
                     typename is_directed::apply<Graph>::type());
    }

    void insert_delta(size_t t, size_t s, int delta, size_t mrs, std::true_type)
    {
        bool src = false;
        if (t != _rnr.first && t != _rnr.second)
        {
            std::swap(t, s);
            src = true;
        }

        assert(t == _rnr.first || t == _rnr.second);

        auto& r_field = (src) ? _r_field_s : _r_field_t;
        auto& nr_field = (src) ? _nr_field_s : _nr_field_t;

        auto& field = (_rnr.first == t) ? r_field : nr_field;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            if (src)
                _entries.emplace_back(s, t);
            else
                _entries.emplace_back(t, s);
            _delta.push_back(delta);
            _mrs.push_back(mrs);
        }
        else
        {
            _delta[field[s]] += delta;
        }
    }

    void insert_delta(size_t t, size_t s, int delta, size_t mrs, std::false_type)
    {
        if (t > s)
            std::swap(t, s);
        if (t != _rnr.first && t != _rnr.second)
            std::swap(t, s);

        assert(t == _rnr.first || t == _rnr.second);

        auto& r_field = _r_field_t;
        auto& nr_field = _nr_field_t;

        auto& field = (_rnr.first == t) ? r_field : nr_field;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            _entries.emplace_back(t, s);
            _delta.push_back(delta);
            _mrs.push_back(mrs);
        }
        else
        {
            _delta[field[s]] += delta;
        }
    }

    int get_delta(size_t t, size_t s)
    {
        if (is_directed::apply<Graph>::type::value)
        {
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return get_delta_source(t, s);
            return 0;
        }
        else
        {
            if (t > s)
                std::swap(t, s);
            if (t != _rnr.first && t != _rnr.second)
                std::swap(t, s);
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            return 0;
        }
    }

    int get_delta_target(size_t r, size_t s)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_t : _nr_field_t;
        if (field[s] == _null)
            return 0;
        else
            return _delta[field[s]];
    }

    int get_delta_source(size_t s, size_t r)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_s : _nr_field_s;
        if (field[s] == _null)
            return 0;
        else
            return _delta[field[s]];
    }

    void clear()
    {
        for (const auto& rs : _entries)
        {
            size_t r = rs.first;
            size_t s = rs.second;
            _r_field_t[r] = _nr_field_t[r] = _null;
            _r_field_t[s] = _nr_field_t[s] = _null;
            if (is_directed::apply<Graph>::type::value)
            {
                _r_field_s[r] = _nr_field_s[r] = _null;
                _r_field_s[s] = _nr_field_s[s] = _null;
            }
        }
        _entries.clear();
        _delta.clear();
        _mrs.clear();
    }

    const vector<pair<size_t, size_t> >& get_entries() { return _entries; }
    const vector<int>& get_delta() { return _delta; }
    const vector<size_t>& get_mrs() { return _mrs; }

private:
    static constexpr size_t _null = numeric_limits<size_t>::max();

    pair<size_t, size_t> _rnr;
    vector<size_t> _r_field_t;
    vector<size_t> _nr_field_t;
    vector<size_t> _r_field_s;
    vector<size_t> _nr_field_s;
    vector<pair<size_t, size_t> > _entries;
    vector<int> _delta;
    vector<size_t> _mrs;
};

template <class Graph>
constexpr size_t EntrySet<Graph>::_null;

struct is_loop_nop
{
    bool operator()(size_t) const { return false; }
};


struct standard_neighbours_policy
{
    template <class Graph, class Vertex>
    IterRange<typename out_edge_iteratorS<Graph>::type>
    get_out_edges(Vertex v, Graph& g) const
    {
        return out_edges_range(v, g);
    }

    template <class Graph, class Vertex>
    IterRange<typename in_edge_iteratorS<Graph>::type>
    get_in_edges(Vertex v, Graph& g) const
    {
        return in_edges_range(v, g);
    }

    template <class Graph, class Vertex, class Weight>
    int get_out_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return out_degreeS()(v, g, eweight);
    }

    template <class Graph, class Vertex, class Weight>
    int get_in_degree(Vertex& v, Graph& g, Weight& eweight) const
    {
        return in_degreeS()(v, g, eweight);
    }
};

// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class CEprop, class EBedge, class MEntries,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
void move_entries(Vertex v, Vertex nr, Vprop& b, Eprop& eweights, CEprop& mrs,
                  EBedge& bedge, Graph& g, BGraph& bg, MEntries& m_entries,
                  const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    vertex_t r = b[v];

    m_entries.set_move(r, nr);

    remove_entries(v, r, b, eweights, mrs, bedge, g, bg, m_entries, npolicy,
                   is_loop);
    add_entries(v, nr, b, eweights, g, bg, m_entries, npolicy, is_loop);
}

template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class CEprop, class EBedge, class MEntries,
          class NPolicy = standard_neighbours_policy, class IL = is_loop_nop>
void remove_entries(Vertex v, Vertex r, Vprop& b, Eprop& eweights, CEprop& mrs,
                    EBedge& bedge, Graph& g, BGraph&, MEntries& m_entries,
                    const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    int self_weight = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];
        int ew = eweights[e];
        //assert(ew > 0);

        const auto& me = bedge[e];

        m_entries.insert_delta(r, s, -ew, mrs[me]);

        if (u == v || is_loop(v))
        {
            if (!is_directed::apply<Graph>::type::value)
                self_weight += ew;
        }
    }

    if (self_weight > 0 && self_weight % 2 == 0)
        m_entries.insert_delta(r, r, self_weight / 2, false);

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[e];

        const auto& me = bedge[e];

        m_entries.insert_delta(s, r, -ew, mrs[me]);
    }
}

template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class MEntries,  class NPolicy = standard_neighbours_policy,
          class IL = is_loop_nop>
void add_entries(Vertex v, Vertex nr, Vprop& b, Eprop& eweights, Graph& g,
                 BGraph&, MEntries& m_entries,
                 const NPolicy& npolicy = NPolicy(), IL is_loop = IL())
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    int self_weight = 0;
    for (auto e : npolicy.get_out_edges(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];
        int ew = eweights[e];
        //assert(ew > 0);

        if (u == v || is_loop(v))
        {
            s = nr;
            if (!is_directed::apply<Graph>::type::value)
                self_weight += ew;
        }
        m_entries.insert_delta(nr, s, +ew);
    }

    if (self_weight > 0 && self_weight % 2 == 0)
        m_entries.insert_delta(nr, nr, -self_weight / 2, false);

    for (auto e : npolicy.get_in_edges(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[e];
        m_entries.insert_delta(s, nr, +ew);
    }
}


// obtain the entropy difference given a set of entries in the e_rs matrix
template <class MEntries, class Eprop, class EMat, class BGraph>
double entries_dS(MEntries& m_entries, Eprop& mrs, EMat& emat, BGraph& bg)
{
    const auto& entries = m_entries.get_entries();
    const auto& delta = m_entries.get_delta();
    const auto& d_mrs = m_entries.get_mrs();

    double dS = 0;
    for (size_t i = 0; i < entries.size(); ++i)
    {
        auto& entry = entries[i];
        auto er = entry.first;
        auto es = entry.second;
        int d = delta[i];
        size_t ers = d_mrs[i];
        if (ers == numeric_limits<size_t>::max())
            ers = get_mrs(er, es, mrs, emat); // slower
        assert(int(ers) + d >= 0);
        dS += eterm(er, es, ers + d, bg) - eterm(er, es, ers, bg);
    }
    return dS;
}


// ====================================
// Construct and manage half-edge lists
// ====================================

//the following guarantees a stable (source, target) ordering even for
//undirected graphs
template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_source(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return source(e, g);
    return min(source(e, g), target(e, g));
}

template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_target(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return target(e, g);
    return max(source(e, g), target(e, g));
}


template <class Graph, class Weighted>
class EGroups
{
public:
    template <class Vprop, class Eprop, class BGraph>
    void init(Vprop b, Eprop eweight, Graph& g, BGraph& bg)
    {
        _egroups.clear();
        _egroups.resize(num_vertices(bg));

        for (auto e : edges_range(g))
        {
            size_t r = b[get_source(e, g)];
            auto& r_elist = _egroups[r];
            _epos[e].first = insert_edge(std::make_tuple(e, true), r_elist,
                                         eweight[e]);

            size_t s = b[get_target(e, g)];
            auto& s_elist = _egroups[s];
            _epos[e].second = insert_edge(std::make_tuple(e, false), s_elist,
                                          eweight[e]);
        }
    }

    void clear()
    {
        _egroups.clear();
    }

    bool empty()
    {
        return _egroups.empty();
    }


    template <class Edge, class EV>
    size_t insert_edge(const Edge& e, EV& elist, size_t)
    {
        elist.push_back(e);
        return elist.size() - 1;
    }

    template <class Edge>
    size_t insert_edge(const Edge& e, DynamicSampler<Edge>& elist,
                       size_t weight)
    {
        return elist.insert(e, weight);
    }


    template <class Edge>
    void remove_edge(size_t pos, vector<Edge>& elist)
    {
        if (get<1>(elist.back()))
            _epos[get<0>(elist.back())].first = pos;
        else
            _epos[get<0>(elist.back())].second = pos;
        if (get<1>(elist[pos]))
            _epos[get<0>(elist[pos])].first = numeric_limits<size_t>::max();
        else
            _epos[get<0>(elist[pos])].second = numeric_limits<size_t>::max();
        elist[pos] = elist.back();
        elist.pop_back();
    }

    template <class Edge>
    void remove_edge(size_t pos, DynamicSampler<Edge>& elist)
    {
        if (get<1>(elist[pos]))
            _epos[get<0>(elist[pos])].first = numeric_limits<size_t>::max();
        else
            _epos[get<0>(elist[pos])].second = numeric_limits<size_t>::max();
        elist.remove(pos);
    }

    template <class Vertex>
    void remove_vertex(Vertex v, Vertex r, Graph& g)
    {
        if (_egroups.empty())
            return;

        typedef Vertex vertex_t;

        auto& elist = _egroups[r];

        // update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = _epos[e].first;
                is_src = (pos < elist.size() && get<0>(elist[pos]) == e);
            }

            size_t pos = (is_src) ? _epos[e].first : _epos[e].second;
            assert(pos < elist.size());
            assert(get<0>(elist[pos]) == e);
            remove_edge(pos, elist);
        }
    }

    template <class Vertex, class Eprop>
    void add_vertex(Vertex v, Vertex s, Eprop& eweight, Graph& g)
    {
        if (_egroups.empty())
            return;

        typedef Vertex vertex_t;

        auto& elist = _egroups[s];

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = _epos[e].first;
                is_src = !(pos < elist.size() && get<0>(elist[pos]) == e);
            }

            typedef typename graph_traits<Graph>::edge_descriptor e_type;
            size_t pos = insert_edge(std::make_tuple(e_type(e), is_src),
                                     elist, size_t(eweight[e]));
            if (is_src)
                _epos[e].first = pos;
            else
                _epos[e].second = pos;
        }
    }

    template <class Vertex, class Eprop>
    void move_vertex(Vertex v, Vertex r, Vertex s, Eprop& eweight,
                     Graph& g)
    {
        if (r == s)
            return;
        remove_egroups(v, r, g);
        add_egroups(v, s, eweight, g);
    }

    template <class Edge, class RNG>
    const auto& sample_edge(const DynamicSampler<Edge>& elist, RNG& rng)
    {
        return get<0>(elist.sample(rng));
    }

    template <class Edge, class RNG>
    const auto& sample_edge(const vector<Edge>& elist, RNG& rng)
    {
        std::uniform_int_distribution<size_t> urand(0, elist.size() - 1);
        size_t ur = urand(rng);
        return get<0>(elist[ur]);
    }

    template <class Vertex, class RNG>
    const auto& sample_edge(Vertex r, RNG& rng)
    {
        return sample_edge(_egroups[r], rng);
    }

private:
    typedef typename std::conditional<Weighted::value,
                                      DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>,
                                      vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>>::type
        sampler_t;
    vector<sampler_t> _egroups;

    typedef typename eprop_map_t<pair<size_t, size_t>>::type epos_t;
    epos_t _epos;
};


// Sample neighbours efficiently
// =============================

template <class Vertex, class Graph, class Eprop, class SMap>
void build_neighbour_sampler(Vertex v, SMap& sampler, Eprop& eweight, Graph& g,
                             bool self_loops=true)
{
    vector<Vertex> neighbours;
    vector<double> probs;
    neighbours.reserve(total_degreeS()(v, g));
    probs.reserve(total_degreeS()(v, g));

    for (auto e : all_edges_range(v, g))
    {
        Vertex u = target(e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(e, g);
        if (!self_loops && u == v)
            continue;
        neighbours.push_back(u);
        probs.push_back(eweight[e]);  // Self-loops will be added twice, and
                                      // hence will be sampled with probability
                                      // 2 * eweight[e]
    }

    if (probs.empty())
    {
        neighbours.push_back(v);
        probs.push_back(1.);
    }

    sampler = Sampler<Vertex, mpl::false_>(neighbours, probs);
};

template <class Vertex, class Graph, class Eprop>
void build_neighbour_sampler(Vertex v, vector<size_t>& sampler, Eprop&, Graph& g,
                             bool self_loops=true)
{
    sampler.clear();
    for (auto e : all_edges_range(v, g))
    {
        Vertex u = target(e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(e, g);
        if (!self_loops && u == v)
            continue;
        sampler.push_back(u); // Self-loops will be added twice
    }
};

template <class Graph, class Eprop, class Sampler>
void init_neighbour_sampler(Graph& g, Eprop eweight, Sampler& sampler)
{
    for (auto v : vertices_range(g))
        build_neighbour_sampler(v, sampler[v], eweight, g);
}

template <class Sampler, class RNG>
auto sample_neighbour(Sampler& sampler, RNG& rng)
{
    return sampler.sample(rng);
}

template <class Vertex, class RNG>
auto sample_neighbour(vector<Vertex>& sampler, RNG& rng)
{
    std::uniform_int_distribution<Vertex> rand(0, sampler.size() - 1);
    return sampler[rand(rng)];
}


// Sampling marginal probabilities on the edges
template <class Graph, class Vprop, class MEprop>
void collect_edge_marginals(size_t B, Vprop b, MEprop p, Graph& g, Graph&)
{
    for (auto e : edges_range(g))
    {
        auto u = min(source(e, g), target(e, g));
        auto v = max(source(e, g), target(e, g));

        auto r = b[u];
        auto s = b[v];

        auto& pv = p[e];
        if (pv.size() < B * B)
            pv.resize(B * B);
        size_t j = r + B * s;
        pv[j]++;
    }
}

template <class Graph, class Vprop, class VVprop>
void collect_vertex_marginals(Vprop b, VVprop p, Graph& g)
{
    for (auto v : vertices_range(g))
    {
        auto r = b[v];
        auto& pv = p[v];
        if (pv.size() <= size_t(r))
            pv.resize(r + 1);
        pv[r]++;
    }
}

} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_UTIL_HH
