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

#ifndef GRAPH_BLOCKMODEL_ENTRIES_HH
#define GRAPH_BLOCKMODEL_ENTRIES_HH


namespace std
{
template <class T, class V>
void operator+=(vector<T>& ret, const V& v)
{
    ret.resize(max(ret.size(), v.size()));
    for (size_t i = 0; i < v.size(); ++i)
        ret[i] += v[i];
}

template <class T, class V>
void operator-=(vector<T>& ret, const V& v)
{
    ret.resize(max(ret.size(), v.size()));
    for (size_t i = 0; i < v.size(); ++i)
        ret[i] -= v[i];
}

template <class T1, class T2>
void operator*=(vector<T1>& ret, const T2& v)
{
    for (auto& x : ret)
        x *= v;
}

template <class T1, class T2>
void operator/=(vector<T1>& ret, const T2& v)
{
    for (auto& x : ret)
        x /= v;
}
}

namespace graph_tool
{

template <class PMap>
class VAdapter
{
public:
    typedef typename boost::property_traits<PMap>::key_type key_t;
    typedef typename boost::property_traits<PMap>::value_type val_t;
    VAdapter(std::vector<PMap>& v, const key_t& e)
        : _v(v), _e(e) {}

    size_t size() const { return _v.size(); }
    val_t& operator[](size_t i) { return _v[i][_e]; }
    const val_t& operator[](size_t i) const { return _v[i][_e]; }

    std::vector<PMap>& _v;
    const key_t& _e;
};

template <class PMap>
VAdapter<PMap> make_vadapter(std::vector<PMap>& v,
                             const typename VAdapter<PMap>::key_t& e)
{
    return VAdapter<PMap>(v, e);
}

// tuple utils
template <size_t i, class T, class OP>
void tuple_op_imp(T&, OP&&)
{
}

template <size_t i, class T, class OP, class Ti, class... Ts>
void tuple_op_imp(T& tuple, OP&& op, Ti&& v, Ts&&... vals)
{
    op(get<i>(tuple), v);
    tuple_op_imp<i+1>(tuple, op, vals...);
}

template <class OP, class T, class... Ts>
__attribute__((flatten))
void tuple_op(T& tuple, OP&& op, Ts&&... vals)
{
    tuple_op_imp<0>(tuple, op, vals...);
}

namespace detail {
   template <class F, class Tuple, std::size_t... I>
   constexpr decltype(auto) tuple_apply_impl(F&& f, Tuple&& t,
                                             std::index_sequence<I...>)
   {
       return f(std::get<I>(std::forward<Tuple>(t))...);
   }
} // namespace detail

template <class F, class Tuple>
constexpr decltype(auto) tuple_apply(F&& f, Tuple&& t)
{
    return detail::tuple_apply_impl
        (std::forward<F>(f), std::forward<Tuple>(t),
         std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{}>{});
}

// Manage a set of block pairs and corresponding edge counts that will be
// updated

template <class Graph, class BGraph, class... EVals>
class EntrySet
{
public:
    typedef typename graph_traits<BGraph>::edge_descriptor bedge_t;

    EntrySet(size_t B)
        : _dummy(_null)
    {
        _r_out_field.resize(B, _null);
        _nr_out_field.resize(B, _null);
        if (is_directed::apply<Graph>::type::value)
        {
            _r_in_field.resize(B, _null);
            _nr_in_field.resize(B, _null);
        }
    }

    void set_move(size_t r, size_t nr, size_t B)
    {
        clear();
        _rnr = make_pair(r, nr);
        if (B > _r_out_field.size())
        {
            _r_out_field.resize(B, _null);
            _nr_out_field.resize(B, _null);
            if (is_directed::apply<Graph>::type::value)
            {
                _r_in_field.resize(B, _null);
                _nr_in_field.resize(B, _null);
            }
        }
    }

    const pair<size_t, size_t>& get_move() { return _rnr; }

    template <bool First, bool Source>
    size_t& get_field_rnr(size_t s, size_t t)
    {
        auto& out_field = First ? _r_out_field : _nr_out_field;
        if (is_directed::apply<Graph>::type::value)
        {
            auto& in_field =  (First ? _r_in_field : _nr_in_field);
            return (Source || s == t) ? out_field[t] : in_field[s];
        }
        else
        {
            return (Source) ? out_field[t] : out_field[s];
        }
    }

    size_t& get_field(size_t s, size_t t)
    {
        if (s == _rnr.first)
            return get_field_rnr<true, true>(s, t);
        if (t == _rnr.first)
            return get_field_rnr<true, false>(s, t);
        if (s == _rnr.second)
            return get_field_rnr<false, true>(s, t);
        if (t == _rnr.second)
            return get_field_rnr<false, false>(s, t);
        return _dummy;
    }

    template <bool Add, class... DVals>
    void insert_delta_dispatch(size_t s, size_t t, size_t& f, int d, DVals&&... delta)
    {
        if (f == _null)
        {
            f = _entries.size();
            _entries.emplace_back(s, t);
            _delta.emplace_back();
            if (sizeof...(delta) > 0)
                _edelta.emplace_back();
        }

        if (Add)
        {
            _delta[f] += d;
            tuple_op(_edelta[f], [&](auto& r, auto& v){ r += v; },
                     delta...);
        }
        else
        {
            _delta[f] -= d;
            tuple_op(_edelta[f], [&](auto& r, auto& v){ r -= v; },
                     delta...);
        }
    }

    template <bool First, bool Source, bool Add, class... DVals>
    __attribute__((flatten))
    void insert_delta_rnr(size_t s, size_t t, int d, DVals&&... delta)
    {
        auto& f = get_field_rnr<First, Source>(s, t);
        insert_delta_dispatch<Add>(s, t, f, d, std::forward<DVals>(delta)...);
    }

    template <bool Add, class... DVals>
    __attribute__((flatten))
    void insert_delta(size_t s, size_t t, int d, DVals&&... delta)
    {
        auto& f = get_field(s, t);
        insert_delta_dispatch<Add>(s, t, f, d, std::forward<DVals>(delta)...);
    }

    int get_delta(size_t r, size_t s)
    {
        size_t f = get_field(r, s);
        if (f == _null)
            return 0;
        return _delta[f];
    }

    void clear()
    {
        for (const auto& rs : _entries)
        {
            size_t r = rs.first;
            size_t s = rs.second;
            auto& f = get_field(r, s);
            f = _null;
        }
        _entries.clear();
        _delta.clear();
        _edelta.clear();
        _mes.clear();
        _recs_entries.clear();
    }

    const vector<pair<size_t, size_t>>& get_entries() { return _entries; }
    const vector<int>&                  get_delta()   { return _delta;   }
    const vector<std::tuple<EVals...>>& get_edelta()  { return _edelta;  }

    template <class Emat>
    vector<bedge_t>& get_mes(Emat& emat)
    {
        _mes.reserve(_entries.size());
        for (size_t i = _mes.size(); i < _entries.size(); ++i)
        {
            auto& rs = _entries[i];
            _mes.push_back(emat.get_me(rs.first, rs.second));
            assert(_mes.back() != emat.get_null_edge() || get<0>(_delta[i]) >= 0);
        }
        return _mes;
    }

    template <class Emat>
    const bedge_t& get_me(size_t r, size_t s, Emat& emat)
    {
        size_t field = get_field(r, s);
        if (field >= _mes.size())
            return emat.get_me(r, s);
        return _mes[field];
    }

    std::tuple<EVals...> _self_eweight;

    std::vector<std::tuple<size_t, size_t,
                           GraphInterface::edge_t, int, std::vector<double>>>
        _recs_entries;

private:
    static constexpr size_t _null = numeric_limits<size_t>::max();

    pair<size_t, size_t> _rnr;
    vector<size_t> _r_out_field;
    vector<size_t> _r_in_field;
    vector<size_t> _nr_out_field;
    vector<size_t> _nr_in_field;
    vector<pair<size_t, size_t>> _entries;
    vector<int> _delta;
    vector<std::tuple<EVals...>> _edelta;
    vector<bedge_t> _mes;
    size_t _dummy;
};

template <class Graph, class BGraph, class... EVals>
constexpr size_t EntrySet<Graph, BGraph, EVals...>::_null;

struct is_loop_nop
{
    bool operator()(size_t) const { return false; }
};

template <bool Remove, bool Add, class Vertex, class Graph, class Vprop,
          class Eprop, class MEntries, class Efilt, class IL, class... Eprops>
__attribute__((flatten))
void modify_entries(Vertex v, Vertex r, Vertex nr, Vprop& _b, Graph& g,
                    Eprop& eweights, MEntries& m_entries, Efilt&& efilt,
                    IL&& is_loop, Eprops&... eprops)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    auto& eself_weight = m_entries._self_eweight;
    int self_weight = 0;
    if (!graph_tool::is_directed(g) && sizeof...(Eprops) > 0)
    {
        tuple_apply([&](auto&&... vals)
                    {
                        auto op = [](auto& x) -> auto& { x *= 0; return x; };
                        auto f = [](auto&...) {};
                        f(op(vals)...);
                    }, eself_weight);
    }

    for (auto e : out_edges_range(v, g))
    {
        if (efilt(e))
            continue;
        vertex_t u = target(e, g);
        vertex_t s = _b[u];
        int ew = eweights[e];

        if (Remove)
            m_entries.template insert_delta_rnr<true, true, false>
                (r, s, ew, make_vadapter(eprops, e)...);

        if (Add)
        {
            if (u == v)
                s = nr;

            if (s != r)
                m_entries.template insert_delta_rnr<false, true, true>
                    (nr, s, ew, make_vadapter(eprops, e)...);
            else
                m_entries.template insert_delta_rnr<true, false, true>
                    (nr, s, ew, make_vadapter(eprops, e)...);
        }

        if ((u == v || is_loop(v)) && !graph_tool::is_directed(g))
        {
            self_weight += ew;
            tuple_op(eself_weight, [&](auto& x, auto& val){ x += val; },
                     make_vadapter(eprops, e)...);
        }
    }

    if (self_weight > 0 && self_weight % 2 == 0 && !graph_tool::is_directed(g))
    {
        if (sizeof...(Eprops) > 0)
        {
            tuple_apply([&](auto&&... vals)
                        {
                            auto op = [](auto& x) -> auto& { x /= 2; return x; };
                            auto f = [](auto&...) {};
                            f(op(vals)...);

                            if (Add)
                                m_entries.template insert_delta_rnr<false, true, false>
                                    (nr, nr, self_weight / 2, vals...);
                            if (Remove)
                                m_entries.template insert_delta_rnr<true, true, true>
                                    (r, r, self_weight / 2, vals...);
                        }, eself_weight);
        }
        else
        {
            if (Add)
                m_entries.template insert_delta_rnr<false, true, false>
                    (nr, nr, self_weight / 2);
            if (Remove)
                m_entries.template insert_delta_rnr<true, true, true>
                    (r, r, self_weight / 2);
        }
    }

    if (graph_tool::is_directed(g))
    {
        for (auto e : in_edges_range(v, g))
        {
            if (efilt(e))
                continue;
            vertex_t u = source(e, g);
            if (u == v)
                continue;
            vertex_t s = _b[u];
            int ew = eweights[e];

            if (Remove)
                m_entries.template insert_delta_rnr<true, false, false>
                    (s, r, ew, make_vadapter(eprops, e)...);
            if (Add)
            {
                if (s != r)
                    m_entries.template insert_delta_rnr<false, false, true>
                        (s, nr, ew, make_vadapter(eprops, e)...);
                else
                    m_entries.template insert_delta_rnr<true, true, true>
                        (s, nr, ew, make_vadapter(eprops, e)...);
            }
        }
    }
}


// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class Vertex, class VProp, class Eprop,
          class MEntries, class EFilt, class IL, class... Eprops>
void move_entries(Vertex v, size_t r, size_t nr, VProp& _b, Graph& g,
                  Eprop& eweights, size_t B, MEntries& m_entries,
                  EFilt&& efilt, IL&& is_loop, Eprops&... eprops)
{
    m_entries.set_move(r, nr, B);

    if (r != null_group)
    {
        if (nr != null_group)
            modify_entries<true, true>(v, r, nr, _b, g, eweights, m_entries,
                                       std::forward<EFilt>(efilt),
                                       std::forward<IL>(is_loop), eprops...);
        else
            modify_entries<true, false>(v, r, nr, _b, g, eweights, m_entries,
                                        std::forward<EFilt>(efilt),
                                        std::forward<IL>(is_loop), eprops...);
    }
    else
    {
        modify_entries<false, true>(v, r, nr, _b, g, eweights, m_entries,
                                    std::forward<EFilt>(efilt),
                                    std::forward<IL>(is_loop), eprops...);
    }
}


// operation on a set of entries
template <class MEntries, class EMat, class OP>
void entries_op(MEntries& m_entries, EMat& emat, OP&& op)
{
    const auto& entries = m_entries.get_entries();
    const auto& delta = m_entries.get_delta();
    const auto& edelta = m_entries.get_edelta();
    auto& mes = m_entries.get_mes(emat);

    for (size_t i = 0; i < entries.size(); ++i)
    {
        auto& entry = entries[i];
        auto er = entry.first;
        auto es = entry.second;
        op(er, es, mes[i], delta[i], edelta[i]); // warning: edelta may be
                                                 // empty, but should be
                                                 // eliminated after compilation
    }
}

// obtain the entropy difference given a set of entries in the e_rs matrix
template <bool exact, class MEntries, class Eprop, class EMat, class BGraph>
double entries_dS(MEntries& m_entries, Eprop& mrs, EMat& emat, BGraph& bg)
{
    double dS = 0;
    entries_op(m_entries, emat,
               [&](auto r, auto s, auto& me, auto d, auto&)
               {
                   size_t ers = 0;
                   if (me != emat.get_null_edge())
                       ers = mrs[me];
                   assert(int(ers) + d >= 0);
                   if (exact)
                       dS += eterm_exact(r, s, ers + d, bg) - eterm_exact(r, s, ers, bg);
                   else
                       dS += eterm(r, s, ers + d, bg) - eterm(r, s, ers, bg);
               });
    return dS;
}

} // namespace graph_tool

#endif // GRAPH_BLOCKMODEL_ENTRIES_HH
