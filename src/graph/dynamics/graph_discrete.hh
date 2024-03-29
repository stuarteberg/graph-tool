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

#ifndef GRAPH_DISCRETE_HH
#define GRAPH_DISCRETE_HH

#include <boost/algorithm/minmax_element.hpp>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_util.hh"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "random.hh"
#include "parallel_rng.hh"
#include "../generation/sampler.hh"
#include "idx_map.hh"

namespace graph_tool
{
using namespace boost;

template <class Value = int32_t>
class discrete_state_base
{
public:
    typedef Value s_t;
    typedef typename vprop_map_t<s_t>::type::unchecked_t smap_t;

    discrete_state_base(smap_t s, smap_t s_temp)
        : _s(s), _s_temp(s_temp),
          _active(std::make_shared<std::vector<size_t>>()) {}

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph&, size_t, smap_t&, RNG&) { return 0; }

    template <class Graph>
    bool is_absorbing(Graph&, size_t) { return false; }

    constexpr bool has_absorbing() { return false; }

    template <class Graph>
    void update_sync(Graph&)
    {
        _s.swap(_s_temp);
    }

    smap_t _s;
    smap_t _s_temp;
    std::shared_ptr<std::vector<size_t>> _active;
};

template <bool exposed>
class SI_state: public discrete_state_base<>
{
public:

    enum State { S, I, R, E};

    template <class Graph, class RNG>
    SI_state(Graph& g, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _beta(python::extract<double>(params["beta"])),
          _epsilon(python::extract<double>(params["epsilon"])),
          _r(python::extract<double>(params["r"])),
          _m(num_vertices(g)),
          _m_temp(num_vertices(g))
    {
        size_t M = 0;
        for (auto v : vertices_range(g))
        {
            size_t k = 0;
            for (auto w : in_or_out_neighbors_range(v, g))
            {
                _m[v] += (_s[w] == State::I);
                ++k;
            }
            _m_temp[v] = _m[v];
            M = std::max(M, k);
        }
        for (size_t m = 0; m < M + 1; ++m)
            _prob.push_back(1-std::pow(1-_beta, m));
    };

    template <class Graph>
    void expose(Graph&, size_t v, smap_t& s_out)
    {
        s_out[v] = State::E;
    }

    template <bool sync, class Graph>
    void infect(Graph& g, size_t v, smap_t& s_out)
    {
        s_out[v] = State::I;
        if (sync)
        {
            for (auto w : out_neighbors_range(v, g))
            {
                auto& m = _m_temp[w];
                #pragma omp atomic
                m++;
            }
        }
        else
        {
            for (auto w : out_neighbors_range(v, g))
                _m[w]++;
        }
    }

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        if (_s[v] == State::I)
            return 0;

        if (exposed && _s[v] == State::E)
        {
            std::bernoulli_distribution einfect(_epsilon);
            if (_epsilon > 0 && einfect(rng))
            {
                infect<sync>(g, v, s_out);
                return 1;
            }
            return 0;
        }

        std::bernoulli_distribution spontaneous(_r);
        if (_r > 0 && spontaneous(rng))
        {
            if constexpr (exposed)
                expose(g, v, s_out);
            else
                infect<sync>(g, v, s_out);
            return 1;
        }

        auto m = _m[v];

        std::bernoulli_distribution minfect(_prob[m]);
        if (m > 0 && minfect(rng))
        {
            if constexpr (exposed)
                expose(g, v, s_out);
            else
                infect<sync>(g, v, s_out);
            return 1;
        }
        return 0;
    }

    template <class Graph>
    void update_sync(Graph& g)
    {
        parallel_vertex_loop(g, [&](auto v) { _m[v] = _m_temp[v]; });
        discrete_state_base::update_sync(g);
    }

    template <class Graph>
    bool is_absorbing(Graph&, size_t v) { return _s[v] == State::I; }

    constexpr bool has_absorbing() { return true; }

protected:
    double _beta;
    double _epsilon;
    double _r;
    discrete_state_base<>::smap_t _m, _m_temp;
    std::vector<double> _prob;
};

template <bool exposed, bool recovered>
class SIS_state: public SI_state<exposed>
{
public:

    typedef typename SI_state<exposed>::smap_t smap_t;
    typedef typename SI_state<exposed>::State State;
    using SI_state<exposed>::_s;
    using SI_state<exposed>::_m;
    using SI_state<exposed>::_m_temp;

    template <class Graph, class RNG>
    SIS_state(Graph& g, smap_t s, smap_t s_temp, python::dict params, RNG& rng)
        : SI_state<exposed>(g, s, s_temp, params, rng),
          _gamma(python::extract<double>(params["gamma"]))
    {};

    template <bool sync, class Graph>
    void recover(Graph& g, size_t v, smap_t& s_out)
    {
        s_out[v] = recovered ? State::R : State::S;
        if constexpr (sync)
        {
            for (auto w : out_neighbors_range(v, g))
            {
                auto& m = _m_temp[w];
                #pragma omp atomic
                m--;
            }
        }
        else
        {
            for (auto w : out_neighbors_range(v, g))
                _m[w]--;
        }
    }

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        if (_s[v] == State::I)
        {
            std::bernoulli_distribution srecover(_gamma);
            if (_gamma > 0 && srecover(rng))
            {
                recover<sync>(g, v, s_out);
                return 1;
            }
            return 0;
        }
        return SI_state<exposed>::template update_node<sync>(g, v, s_out, rng);
    }

    template <class Graph>
    bool is_absorbing(Graph&, size_t v) { return recovered && _s[v] == State::R; }

    constexpr bool has_absorbing() { return recovered; }


protected:
    double _gamma;
};

template <bool exposed>
class SIRS_state: public SIS_state<exposed, true>
{
public:

    typedef typename SI_state<exposed>::smap_t smap_t;
    typedef typename SI_state<exposed>::State State;
    using SI_state<exposed>::_s;

    template <class Graph, class RNG>
    SIRS_state(Graph& g, smap_t s, smap_t s_temp, python::dict params, RNG& rng)
        : SIS_state<exposed,true>(g, s, s_temp, params, rng),
          _mu(python::extract<double>(params["mu"]))
    {};


    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        if (_s[v] == State::R)
        {
            std::bernoulli_distribution srecover(_mu);
            if (_mu > 0 && srecover(rng))
            {
                s_out[v] = State::S;
                return 1;
            }
            return 0;
        }
        return SIS_state<exposed,true>::template update_node<sync>(g, v, s_out, rng);
    }

    template <class Graph>
    bool is_absorbing(Graph&, size_t) { return false; }

    constexpr bool has_absorbing() { return false; }

private:
    double _mu;
};


class voter_state: public discrete_state_base<>
{
public:
    template <class Graph, class RNG>
    voter_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _q(python::extract<int32_t>(params["q"])),
          _r(python::extract<double>(params["r"])){};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];
        std::bernoulli_distribution random(_r);
        if (_r > 0 && random(rng))
        {
            std::uniform_int_distribution<int32_t> sample(0, _q - 1);
            auto t = sample(rng);
            s_out[v] = t;
            return s != t;
        }

        size_t w;
        if (graph_tool::is_directed(g))
        {
            if (in_degreeS()(v, g) == 0)
                return 0;
            w = random_in_neighbor(v, g, rng);
        }
        else
        {
            if (out_degree(v, g) == 0)
                return 0;
            w = random_out_neighbor(v, g, rng);
        }
        auto t = _s[w];
        s_out[v] = t;
        return s != t;
    }

private:

    size_t _q;
    double _r;
};

class majority_voter_state: public discrete_state_base<>
{
public:
    template <class Graph, class RNG>
    majority_voter_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _q(python::extract<size_t>(params["q"])),
          _r(python::extract<double>(params["r"])) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];
        std::bernoulli_distribution random(_r);
        if (_r > 0 && random(rng))
        {
            std::uniform_int_distribution<int32_t> sample(0, _q - 1);
            auto t = sample(rng);
            s_out[v] = t;
            return (t != s);
        }

        for (auto w : in_or_out_neighbors_range(v, g))
            _m[_s[w]]++;

        if (_m.empty())
            return 0;

        auto qmax = boost::first_max_element(_m.begin(), _m.end(),
                                             [&](auto& a, auto& b)
                                             { return a.second < b.second; })->second;

        for (auto qc : _m)
        {
            if (qc.second == qmax)
                _qs.push_back(qc.first);
        }

        auto t = uniform_sample(_qs, rng);
        s_out[v] = t;

        _qs.clear();
        _m.clear();

        return (t != s);
    }

private:

    int32_t _q;
    double _r;
    idx_map<int32_t, size_t> _m;
    std::vector<int32_t> _qs;
};

class binary_threshold_state: public discrete_state_base<>
{
public:
    typedef vprop_map_t<double>::type::unchecked_t hmap_t;
    typedef eprop_map_t<double>::type::unchecked_t wmap_t;

    template <class Graph, class RNG>
    binary_threshold_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _r(python::extract<double>(params["r"])) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        double m = 0;
        size_t k = 0;
        std::bernoulli_distribution flip(_r);
        for (auto e : in_or_out_edges_range(v, g))
        {
            auto u = source(e, g);
            auto t = _s[u];

            if (_r > 0 && flip(rng))
                t ^= 1;

            m += _w[e] * t;
            k++;
        }

        auto s = _s[v];
        auto t = (m > _h[v] * k);
        s_out[v] = t;

        return s != t;
    }

private:

    hmap_t _h;
    wmap_t _w;
    double _r;
};

class ising_glauber_state: public discrete_state_base<>
{
public:

    typedef eprop_map_t<double>::type::unchecked_t wmap_t;
    typedef vprop_map_t<double>::type::unchecked_t hmap_t;

    template <class Graph, class RNG>
    ising_glauber_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _beta(python::extract<double>(params["beta"])) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];

        double m = 0;
        for (auto e : in_or_out_edges_range(v, g))
            m += _w[e] * _s[source(e, g)];

        std::bernoulli_distribution up(1./(1. + exp(-2 * (_h[v] + _beta * m))));

        auto t = up(rng) ? 1 : -1;
        s_out[v] = t;

        return s != t;
    }

private:

    wmap_t _w;
    hmap_t _h;
    double _beta;
};

class cising_glauber_state: public discrete_state_base<double>
{
public:

    typedef eprop_map_t<double>::type::unchecked_t wmap_t;
    typedef vprop_map_t<double>::type::unchecked_t hmap_t;

    template <class Graph, class RNG>
    cising_glauber_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _beta(python::extract<double>(params["beta"])) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];

        double m = 0;
        for (auto e : in_or_out_edges_range(v, g))
            m += _w[e] * _s[source(e, g)];
        double a = _beta * m + _h[v];

        std::uniform_real_distribution<> U(0, 1);
        double u = U(rng);

        double t;
        if (abs(a) > 1e-8)
        {
            //s_out[u] = t = log(exp(a + log(u)) + exp(-a + log1p(-u))) / a;
            if (a + log(u) > -a + log1p(-u))
                s_out[v] = t = 1 + (log(u) + log1p(exp(-2*a + log1p(-u) - log(u)))) / a;
            else
                s_out[v] = t = -1 + (log1p(-u) + log1p(exp(2 * a + log(u) - log1p(-u)))) / a;
        }
        else
        {
            s_out[v] = t = 2 * u - 1;
        }

        return s != t;
    }

private:

    wmap_t _w;
    hmap_t _h;
    double _beta;
};

class ising_metropolis_state: public discrete_state_base<>
{
public:

    typedef eprop_map_t<double>::type::unchecked_t wmap_t;
    typedef vprop_map_t<double>::type::unchecked_t hmap_t;

    template <class Graph, class RNG>
    ising_metropolis_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _beta(python::extract<double>(params["beta"])) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];

        double m = 0;
        for (auto e : in_or_out_edges_range(v, g))
            m += _w[e] * _s[source(e, g)];

        auto t = -s;

        std::uniform_real_distribution<double> sample;
        double r = exp(2 * t * (_h[v] + _beta * m));
        if (r > 1 || sample(rng) < r)
            s_out[v] = t;
        else
            t = s;

        return s != t;
    }

private:

    wmap_t _w;
    hmap_t _h;
    double _beta;
};

class kirman_state: public discrete_state_base<>
{
public:
    template <class Graph, class RNG>
    kirman_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _d(python::extract<double>(params["d"])),
          _c1(python::extract<double>(params["c1"])),
          _c2(python::extract<double>(params["c2"])){};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        bool active = _s[v] != 0;
        if (active)
        {
            std::bernoulli_distribution spontaneous(_c2);
            if (_c2 > 0 && spontaneous(rng))
            {
                s_out[v] = 0;
                return 1;
            }
        }
        else
        {
            std::bernoulli_distribution spontaneous(_c1);
            if (_c1 > 0 && spontaneous(rng))
            {
                s_out[v] = 1;
                return 1;
            }
        }

        size_t m = 0, k = 0;
        for (auto w : in_or_out_neighbors_range(v, g))
        {
            m += _s[w];
            k++;
        }

        if (active)
            m = k - m;

        std::bernoulli_distribution infect(1-std::pow(1-_d, m));

        if (infect(rng))
        {
            s_out[v] = not active;
            return 1;
        }
        return 0;
    }

private:

    double _d;
    double _c1;
    double _c2;
};

class potts_glauber_state: public discrete_state_base<>
{
public:

    typedef eprop_map_t<double>::type::unchecked_t wmap_t;
    typedef vprop_map_t<std::vector<double>>::type::unchecked_t hmap_t;

    template <class Graph, class RNG>
    potts_glauber_state(Graph&, smap_t s, smap_t s_temp, python::dict params, RNG&)
        : discrete_state_base(s, s_temp),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _f(get_array<double, 2>(params["f"])), _q(_f.shape()[0]), _m(_q) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        for (int r = 0; r < _q; ++r)
            _m[r] = _h[v][r];

        for (auto e : in_or_out_edges_range(v, g))
        {
            auto u = source(e, g);
            auto s = _s[u];
            for (int r = 0; r < _q; ++r)
                _m[r] += _w[e] * _f[r][s];
        }

        double mmax = *boost::first_max_element(_m.begin(), _m.end());

        for (int r = 0; r < _q; ++r)
        {
            _m[r] = exp(_m[r] - mmax);
            if (r > 0)
                _m[r] += _m[r-1];
        }

        std::uniform_real_distribution<double> sample(0, _m[_q-1]);

        auto u = sample(rng);

        int32_t t = 0;
        for (; t < _q; ++t)
        {
            if (_m[t] >= u)
                break;
        }

        int32_t s = _s[v];
        s_out[v] = t;

        return s != t;
    }

private:

    wmap_t _w;
    hmap_t _h;
    multi_array_ref<double, 2> _f;
    int32_t _q;
    std::vector<double> _m;
};

class potts_metropolis_state: public discrete_state_base<>
{
public:

    typedef eprop_map_t<double>::type::unchecked_t wmap_t;
    typedef vprop_map_t<std::vector<double>>::type::unchecked_t hmap_t;

    template <class Graph, class RNG>
    potts_metropolis_state(Graph&, smap_t s, smap_t s_temp, python::dict params,
                           RNG&)
        : discrete_state_base(s, s_temp),
          _w(any_cast<wmap_t::checked_t>(python::extract<any>(params["w"].attr("_get_any")())()).get_unchecked()),
          _h(any_cast<hmap_t::checked_t>(python::extract<any>(params["h"].attr("_get_any")())()).get_unchecked()),
          _f(get_array<double, 2>(params["f"])), _q(_f.shape()[0]), _m(_q) {};

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        auto s = _s[v];

        std::uniform_int_distribution<int> sample_s(0, _q - 1);

        auto ns = sample_s(rng);

        if (s == ns)
            return 0;

        double dH = _h[v][ns] - _h[v][s];

        for (auto e : in_or_out_edges_range(v, g))
        {
            auto u = source(e, g);
            auto r = _s[u];
            dH += _w[e] * (_f[ns][r] - _f[s][r]);
        }

        std::uniform_real_distribution<double> usample(0, 1);
        if (dH < 0 || usample(rng) < exp(-dH))
        {
            s_out[v] = ns;
            return 1;
        }
        return 0;
    }

private:

    wmap_t _w;
    hmap_t _h;
    multi_array_ref<double, 2> _f;
    int32_t _q;
    std::vector<double> _m;
};


class axelrod_state: public discrete_state_base<std::vector<int32_t>>
{
public:

    template <class Graph, class RNG>
    axelrod_state(Graph& g, smap_t s, smap_t s_temp, python::dict params, RNG& rng)
        : discrete_state_base(s, s_temp),
          _q(python::extract<int32_t>(params["q"])),
          _f(python::extract<int32_t>(params["f"])),
          _r(python::extract<double>(params["r"]))
    {
        std::uniform_int_distribution<int32_t> sample(0, _q - 1);
        for (auto v : vertices_range(g))
        {
            auto& s_v = _s[v];
            for (size_t i = s_v.size(); i < _f; ++i)
                s_v.push_back(sample(rng));
        }
    };

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        std::bernoulli_distribution random(_r);
        if (_r > 0 && random(rng))
        {
            std::uniform_int_distribution<int32_t> f_sample(0, _f - 1);
            std::uniform_int_distribution<int32_t> sample(0, _q - 1);
            auto i = f_sample(rng);
            auto t = sample(rng);
            auto s = _s[v][i];
            s_out[v][i] = t;
            return s != t;
        }

        size_t w;
        if (graph_tool::is_directed(g))
        {
            if (in_degreeS()(v, g) == 0)
                return 0;
            w = random_in_neighbor(v, g, rng);
        }
        else
        {
            if (out_degree(v, g) == 0)
                return 0;
            w = random_out_neighbor(v, g, rng);
        }

        auto& s_v = _s[v];
        auto& s_w = _s[w];

        _temp.clear();
        size_t d = 0;
        for (size_t i = 0; i < _f; ++i)
        {
            if (s_v[i] == s_w[i])
                d++;
            else
                _temp.push_back(i);
        }

        std::bernoulli_distribution copy(d / double(_f));
        if (!_temp.empty() && copy(rng))
        {
            size_t i = uniform_sample(_temp, rng);
            s_out[v][i] = _s[w][i];
            return 1;
        }
        return 0;
    }

private:

    size_t _q;
    size_t _f;
    double _r;

    std::vector<size_t> _temp;
};

class boolean_state: public discrete_state_base<uint8_t>
{
public:
    typedef vprop_map_t<std::vector<uint8_t>>::type::unchecked_t fmap_t;

    template <class Graph, class RNG>
    boolean_state(Graph& g, smap_t s, smap_t s_temp, python::dict params, RNG& rng)
        : discrete_state_base(s, s_temp),
          _f(any_cast<fmap_t::checked_t>(python::extract<any>(params["f"].attr("_get_any")())()).get_unchecked()),
          _r(python::extract<double>(params["r"]))
    {
        double p = python::extract<double>(params["p"]);
        std::bernoulli_distribution random(p);
        for (auto v : vertices_range(g))
        {
            auto& f_v = _f[v];
            size_t k = (graph_tool::is_directed(g)) ? in_degreeS()(v, g) : out_degree(v, g);
            size_t M = 1 << k;
            for (size_t pos = f_v.size(); pos < M; ++pos)
                f_v.push_back(random(rng));
        }
    };

    template <bool sync, class Graph, class RNG>
    bool update_node(Graph& g, size_t v, smap_t& s_out, RNG& rng)
    {
        std::bernoulli_distribution flip(_r);

        uint64_t idx = 0;
        size_t pos = 0;
        for (auto w : in_or_out_neighbors_range(v, g))
        {
            bool s_w = _s[w];
            if (_r > 0 && flip(rng))
                s_w ^= 1;
            if (s_w)
                idx += 1 << pos;
            ++pos;
        }

        auto s = s_out[v];
        auto& f_v = _f[v];
        s_out[v] = f_v[idx];
        return s != s_out[v];
    }

private:

    fmap_t _f;
    double _r;
};

template <class Graph, class State, class RNG>
size_t discrete_iter_sync(Graph& g, State state, size_t niter, RNG& rng_)
{
    size_t nflips = 0;

    parallel_rng<rng_t>::init(rng_);

    auto& active = *state._active;

    for (size_t i = 0; i < niter; ++i)
    {
        if (active.empty())
            break;

        #pragma omp parallel if (active.size() > OPENMP_MIN_THRESH)     \
             firstprivate(state) reduction(+:nflips)
        parallel_loop_no_spawn
            (active,
             [&] (auto, auto v)
             {
                 auto& rng = parallel_rng<rng_t>::get(rng_);
                 state._s_temp[v] = state._s[v];
                 nflips += state.template update_node<true>(g, v,
                                                            state._s_temp,
                                                            rng);
             });

        state.update_sync(g);

        if constexpr (state.has_absorbing())
        {
            auto iter = std::remove_if(active.begin(),
                                       active.end(),
                                       [&](auto v)
                                       {
                                           state._s_temp[v] = state._s[v];
                                           return state.is_absorbing(g, v);
                                       });
            active.erase(iter, active.end());
        }
    }
    return nflips;
}

template <class Graph, class State, class RNG>
size_t discrete_iter_async(Graph& g, State state, size_t niter, RNG& rng)
{
    size_t nflips = 0;

    auto& active = *state._active;
    for (size_t i = 0; i < niter; ++i)
    {
        if (active.empty())
            break;

        auto iter = uniform_sample_iter(active, rng);
        nflips += state.template update_node<false>(g, *iter, state._s, rng);

        if constexpr (state.has_absorbing())
        {
            if (state.is_absorbing(g, *iter))
            {
                std::swap(*iter, active.back());
                active.pop_back();
            }
        }
    }
    return nflips;
}


} // namespace graph_tool

#endif // GRAPH_DISCRETE_HH
