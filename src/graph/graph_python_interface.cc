// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
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

#include "graph_filtering.hh"
#include "graph.hh"
#include "graph_util.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <set>


using namespace std;
using namespace boost;
using namespace graph_tool;

namespace graph_tool
{

struct get_vertex_iterator
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi,
                    python::object& iter) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        iter = python::object(PythonIterator<Graph, PythonVertex<Graph>,
                                             vertex_iterator>(gp, vertices(g)));
    }
};

python::object get_vertices(GraphInterface& gi)
{
    python::object iter;
    run_action<>()(gi, std::bind(get_vertex_iterator(),
                                 placeholders::_1,
                                 std::ref(gi),
                                 std::ref(iter)))();
    return iter;
}

struct get_vertex_soft
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t i, python::object& v) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        if (i < num_vertices(g))
            v = python::object(PythonVertex<Graph>(gp, vertex(i, g)));
        else
            v = python::object(PythonVertex<Graph>(gp,
                                                   graph_traits<Graph>::null_vertex()));
    }
};

struct get_vertex_hard
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t i, python::object& v) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        size_t c = 0;
        for (auto vi : vertices_range(g))
        {
            if (c == i)
            {
                v = python::object(PythonVertex<Graph>(gp, vi));
                return;
            }
            ++c;
        }
        v = python::object(PythonVertex<Graph>(gp,
                                               graph_traits<Graph>::null_vertex()));
    }
};

python::object get_vertex(GraphInterface& gi, size_t i, bool use_index)
{
    python::object v;
    if (!use_index)
        run_action<>()(gi,
                       std::bind(get_vertex_hard(), placeholders::_1,
                                 std::ref(gi), i, std::ref(v)))();
    else
        run_action<>()(gi,
                       std::bind(get_vertex_soft(), placeholders::_1,
                                 std::ref(gi), i, std::ref(v)))();
    return v;
}

struct get_edge_iterator
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, python::object& iter)
        const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        iter = python::object(PythonIterator<Graph, PythonEdge<Graph>,
                                             edge_iterator>(gp, edges(g)));
    }
};

python::object get_edges(GraphInterface& gi)
{
    python::object iter;
    run_action<>()(gi, std::bind(get_edge_iterator(), placeholders::_1,
                                 std::ref(gi), std::ref(iter)))();
    return iter;
}

struct add_new_vertex
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t n,
                    python::object& new_v) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        if (n != 1)
        {
            for (size_t i = 0; i < n; ++i)
                add_vertex(g);
            new_v = python::object();
        }
        else
        {
            new_v = python::object(PythonVertex<Graph>(gp, add_vertex(g)));
        }
    }
};


python::object add_vertex(GraphInterface& gi, size_t n)
{
    python::object v;
    run_action<>()(gi, std::bind(add_new_vertex(), placeholders::_1,
                                 std::ref(gi), n, std::ref(v)))();
    return v;
}


void remove_vertex_array(GraphInterface& gi, const python::object& oindex, bool fast)
{
    boost::multi_array_ref<int64_t,1> index = get_array<int64_t,1>(oindex);
    auto& g = gi.get_graph();
    if (fast)
    {
        for (auto v : index)
            remove_vertex_fast(vertex(v, g), g);
    }
    else
    {
        for (auto v : index)
            remove_vertex(vertex(v, g), g);
    }
}

void remove_vertex(GraphInterface& gi, size_t v, bool fast)
{
    auto& g = gi.get_graph();
    if (fast)
    {
        remove_vertex_fast(vertex(v, g), g);
    }
    else
    {
        remove_vertex(vertex(v, g), g);
    }
}

struct do_clear_vertex
{
    template <class Graph>
    void operator()(Graph& g, size_t v) const
    {
        clear_vertex(vertex(v, g), g);
    }
};

void clear_vertex(GraphInterface& gi, size_t v)
{
    run_action<>()(gi, std::bind(do_clear_vertex(), placeholders::_1, v))();
}

struct add_new_edge
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t s, size_t t,
                    python::object& new_e) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        auto e = add_edge(vertex(s, g), vertex(t, g), g).first;
        new_e = python::object(PythonEdge<Graph>(gp, e));
    }
};

python::object add_edge(GraphInterface& gi, size_t s, size_t t)
{
    python::object new_e;
    run_action<>()(gi, std::bind(add_new_edge(), placeholders::_1, std::ref(gi),
                                 s, t, std::ref(new_e)))();
    return new_e;
}

struct get_edge_descriptor
{
    template <class Graph>
    void operator()(const Graph&, const python::object& e,
                    typename GraphInterface::edge_t& edge,
                    bool& found)  const
    {
        PythonEdge<Graph>& pe = python::extract<PythonEdge<Graph>&>(e);
        pe.check_valid();
        edge = pe.get_descriptor();
        found = true;
    }
};

void remove_edge(GraphInterface& gi, const python::object& e)
{
    GraphInterface::edge_t de;
    bool found = false;
    run_action<>()(gi, std::bind(get_edge_descriptor(), placeholders::_1,
                                 std::ref(e), std::ref(de), std::ref(found)))();
    remove_edge(de, gi.get_graph());
    if (!found)
        throw ValueException("invalid edge descriptor");
}

struct get_edge_dispatch
{
    template <class Graph>
    void operator()(Graph& g, GraphInterface& gi, size_t s, size_t t,
                    bool all_edges, boost::python::list& es) const
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(gi, g);
        size_t k_t = is_directed::apply<Graph>::type::value ?
            in_degreeS()(t, g) : out_degree(t, g);
        if (out_degree(s, g) <= k_t)
        {
            for (auto e : out_edges_range(vertex(s, g), g))
            {
                if (target(e, g) == vertex(t, g))
                {
                    es.append(PythonEdge<Graph>(gp, e));
                    if (!all_edges)
                        break;
                }
            }
        }
        else
        {
            for (auto e : in_or_out_edges_range(vertex(t, g), g))
            {
                auto w = is_directed::apply<Graph>::type::value ?
                    source(e, g) : target(e, g);
                if (w == vertex(s, g))
                {
                    if (!is_directed::apply<Graph>::type::value)
                        e.inv ^= true;
                    es.append(PythonEdge<Graph>(gp, e));
                    if (!all_edges)
                        break;
                }
            }
        }
    }
};

python::object get_edge(GraphInterface& gi, size_t s, size_t t, bool all_edges)
{
    python::list es;
    run_action<>()(gi, std::bind(get_edge_dispatch(), placeholders::_1,
                                 std::ref(gi), s, t, all_edges,
                                 std::ref(es)))();
    return es;
}


struct get_degree_map
{
    template <class Graph, class DegS, class Weight>
    void operator()(const Graph& g, python::object& odeg_map, DegS deg, Weight weight) const
    {
        typedef typename detail::get_weight_type<Weight>::type weight_t;
        typedef typename mpl::if_<std::is_same<weight_t, size_t>, int32_t, weight_t>::type deg_t;

        typedef typename property_map_type::apply<deg_t,
                                                  GraphInterface::vertex_index_map_t>::type
            map_t;

        map_t cdeg_map(get(vertex_index, g));
        typename map_t::unchecked_t deg_map = cdeg_map.get_unchecked(num_vertices(g));

        int i, N = num_vertices(g);
        #pragma omp parallel for default(shared) private(i) schedule(runtime) if (N > 100)
        for (i = 0; i < N; ++i)
        {
            typename graph_traits<Graph>::vertex_descriptor v = vertex(i, g);
            if (v == graph_traits<Graph>::null_vertex())
                continue;
            deg_map[v] = deg(v, g, weight);
        }

        odeg_map = python::object(PythonPropertyMap<map_t>(cdeg_map));
    }
};

python::object GraphInterface::degree_map(string deg, boost::any weight) const
{

    python::object deg_map;

    typedef mpl::push_back<edge_scalar_properties,
                           detail::no_weightS>::type weight_t;
    if (weight.empty())
        weight = detail::no_weightS();

    if (deg == "in")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), in_degreeS(), placeholders::_2), weight_t())
            (weight);
    else if (deg == "out")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), out_degreeS(), placeholders::_2), weight_t())
            (weight);
    else if (deg == "total")
        run_action<>()(const_cast<GraphInterface&>(*this),
                       std::bind(get_degree_map(), placeholders::_1,
                                 std::ref(deg_map), total_degreeS(), placeholders::_2), weight_t())
            (weight);
    return deg_map;
}

//
// Below are the functions with will properly register all the types to python,
// for every filter, type, etc.
//

// this will register all the Vertex/Edge classes to python
struct export_python_interface
{
    template <class Graph, class GraphViews>
    void operator()(Graph* gp, python::list vclasses,
                    python::list eclasses, GraphViews) const
    {
        using namespace boost::python;

        class_<PythonVertex<Graph>, bases<VertexBase>> vclass("Vertex", no_init);
        vclass
            .def("__in_degree", &PythonVertex<Graph>::get_in_degree,
                 "Return the in-degree.")
            .def("__weighted_in_degree", &PythonVertex<Graph>::get_weighted_in_degree,
                 "Return the weighted in-degree.")
            .def("__out_degree", &PythonVertex<Graph>::get_out_degree,
                 "Return the out-degree.")
            .def("__weighted_out_degree", &PythonVertex<Graph>::get_weighted_out_degree,
                 "Return the weighted out-degree.")
            .def("in_edges", &PythonVertex<Graph>::in_edges,
                 "Return an iterator over the in-edges.")
            .def("out_edges", &PythonVertex<Graph>::out_edges,
                 "Return an iterator over the out-edges.")
            .def("is_valid", &PythonVertex<Graph>::is_valid,
                 "Return whether the vertex is valid.")
            .def("graph_ptr", &PythonVertex<Graph>::get_graph_ptr)
            .def("graph_type", &PythonVertex<Graph>::get_graph_type)
            .def("__str__", &PythonVertex<Graph>::get_string)
            .def("__int__", &PythonVertex<Graph>::get_index)
            .def("__hash__", &PythonVertex<Graph>::get_hash);

        vclasses.append(vclass);

        class_<PythonEdge<Graph>, bases<EdgeBase>> eclass("Edge", no_init);
        eclass
            .def("source", &PythonEdge<Graph>::get_source,
                 "Return the source vertex.")
            .def("target", &PythonEdge<Graph>::get_target,
                 "Return the target vertex.")
            .def("is_valid", &PythonEdge<Graph>::is_valid,
                 "Return whether the edge is valid.")
            .def("graph_ptr", &PythonVertex<Graph>::get_graph_ptr)
            .def("graph_type", &PythonVertex<Graph>::get_graph_type)
            .def("__str__", &PythonEdge<Graph>::get_string)
            .def("__hash__", &PythonEdge<Graph>::get_hash);

        boost::mpl::for_each<GraphViews>(std::bind(export_python_interface(),
                                                   gp, std::placeholders::_1,
                                                   std::ref(eclass)));

        eclasses.append(eclass);

        typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
        class_<PythonIterator<Graph, PythonVertex<Graph>, vertex_iterator> >
            ("VertexIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("__next__", &PythonIterator<Graph, PythonVertex<Graph>,
                                             vertex_iterator>::next)
            .def("next", &PythonIterator<Graph, PythonVertex<Graph>,
                                         vertex_iterator>::next);

        typedef typename graph_traits<Graph>::edge_iterator edge_iterator;
        class_<PythonIterator<Graph, PythonEdge<Graph>,
                              edge_iterator> >("EdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("__next__", &PythonIterator<Graph, PythonEdge<Graph>,
                                         edge_iterator>::next)
            .def("next", &PythonIterator<Graph, PythonEdge<Graph>,
                                         edge_iterator>::next);

        typedef typename graph_traits<Graph>::out_edge_iterator
            out_edge_iterator;
        class_<PythonIterator<Graph, PythonEdge<Graph>,
                              out_edge_iterator> >("OutEdgeIterator", no_init)
            .def("__iter__", objects::identity_function())
            .def("__next__", &PythonIterator<Graph, PythonEdge<Graph>,
                                             out_edge_iterator>::next)
            .def("next", &PythonIterator<Graph, PythonEdge<Graph>,
                                         out_edge_iterator>::next);

        typedef typename graph_traits<Graph>::directed_category
            directed_category;
        typedef typename std::is_convertible<directed_category,
                                             boost::directed_tag>::type is_directed;
        if (is_directed::value)
        {
            typedef typename in_edge_iteratorS<Graph>::type in_edge_iterator;
            class_<PythonIterator<Graph, PythonEdge<Graph>,
                                  in_edge_iterator> >("InEdgeIterator", no_init)
                .def("__iter__", objects::identity_function())
                .def("__next__", &PythonIterator<Graph, PythonEdge<Graph>,
                                                 in_edge_iterator>::next)
                .def("next", &PythonIterator<Graph, PythonEdge<Graph>,
                                             in_edge_iterator>::next);
        }
    }

    template <class Graph, class OGraph, class Eclass>
    void operator()(Graph*, OGraph*, Eclass& eclass) const
    {
        std::function<bool(const PythonEdge<Graph>&,
                           const PythonEdge<OGraph>&)> eq =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 == e2; };
        std::function<bool(const PythonEdge<Graph>& e1,
                           const PythonEdge<OGraph>&)> ne =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 != e2; };
        std::function<bool(const PythonEdge<Graph>&,
                           const PythonEdge<OGraph>&)> gt =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 > e2; };
        std::function<bool(const PythonEdge<Graph>&,
                           const PythonEdge<OGraph>&)> lt =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 < e2; };
        std::function<bool(const PythonEdge<Graph>&,
                           const PythonEdge<OGraph>&)> ge =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 >= e2; };
        std::function<bool(const PythonEdge<Graph>&,
                           const PythonEdge<OGraph>&)> le =
            [] (const PythonEdge<Graph>& e1,
                const PythonEdge<OGraph>& e2) -> bool { return e1 <= e2; };

        eclass
            .def("__eq__", eq)
            .def("__ne__", ne)
            .def("__lt__", lt)
            .def("__gt__", gt)
            .def("__le__", le)
            .def("__ge__", ge);
    }
};

PythonPropertyMap<GraphInterface::vertex_index_map_t>
get_vertex_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::vertex_index_map_t>
        (g.get_vertex_index());
}

PythonPropertyMap<GraphInterface::edge_index_map_t>
do_get_edge_index(GraphInterface& g)
{
    return PythonPropertyMap<GraphInterface::edge_index_map_t>
        (g.get_edge_index());
}

template <class ValueList>
struct add_edge_list
{
    template <class Graph>
    void operator()(Graph& g, python::object aedge_list,
                    python::object& eprops, bool& found) const
    {
        boost::mpl::for_each<ValueList>(std::bind(dispatch(), std::ref(g),
                                                  std::ref(aedge_list),
                                                  std::ref(eprops),
                                                  std::ref(found),
                                                  placeholders::_1));
    }

    struct dispatch
    {
        template <class Graph, class Value>
        void operator()(Graph& g, python::object& aedge_list,
                        python::object& oeprops, bool& found, Value) const
        {
            if (found)
                return;
            try
            {
                boost::multi_array_ref<Value, 2> edge_list = get_array<Value, 2>(aedge_list);

                if (edge_list.shape()[1] < 2)
                    throw GraphException("Second dimension in edge list must be of size (at least) two");

                typedef typename graph_traits<Graph>::edge_descriptor edge_t;
                vector<DynamicPropertyMapWrap<Value, edge_t>> eprops;
                python::stl_input_iterator<boost::any> iter(oeprops), end;
                for (; iter != end; ++iter)
                    eprops.emplace_back(*iter, writable_edge_properties());

                size_t n_props = std::min(eprops.size(), edge_list.shape()[1] - 2);

                for (const auto& e : edge_list)
                {
                    size_t s = e[0];
                    size_t t = e[1];
                    while (s >= num_vertices(g) || t >= num_vertices(g))
                        add_vertex(g);
                    auto ne = add_edge(vertex(s, g), vertex(t, g), g).first;
                    for (size_t i = 0; i < n_props; ++i)
                    {
                        try
                        {
                            put(eprops[i], ne, e[i + 2]);
                        }
                        catch(bad_lexical_cast&)
                        {
                            throw ValueException("Invalid edge property value: " +
                                                 lexical_cast<string>(e[i + 2]));
                        }
                    }
                }
                found = true;
            }
            catch (invalid_numpy_conversion& e) {}
        }
    };
};

void do_add_edge_list(GraphInterface& gi, python::object aedge_list,
                      python::object eprops)
{
    typedef mpl::vector<bool, char, uint8_t, uint16_t, uint32_t, uint64_t,
                        int8_t, int16_t, int32_t, int64_t, uint64_t, double,
                        long double> vals_t;
    bool found = false;
    run_action<>()(gi, std::bind(add_edge_list<vals_t>(), placeholders::_1,
                                 aedge_list, std::ref(eprops),
                                 std::ref(found)))();
    if (!found)
        throw GraphException("Invalid type for edge list; must be two-dimensional with a scalar type");
}

template <class ValueList>
struct add_edge_list_hash
{
    template <class Graph, class VProp>
    void operator()(Graph& g, python::object aedge_list, VProp vmap,
                    bool& found, bool use_str, python::object& eprops) const
    {
        boost::mpl::for_each<ValueList>(std::bind(dispatch(), std::ref(g),
                                                  std::ref(aedge_list), std::ref(vmap),
                                                  std::ref(found), std::ref(eprops),
                                                  placeholders::_1));
        if (!found)
        {
            if (use_str)
                dispatch()(g, aedge_list, vmap, found, eprops, std::string());
            else
                dispatch()(g, aedge_list, vmap, found, eprops, python::object());
        }
    }

    struct dispatch
    {
        template <class Graph, class VProp, class Value>
        void operator()(Graph& g, python::object& aedge_list, VProp& vmap,
                        bool& found, python::object& oeprops, Value) const
        {
            if (found)
                return;
            try
            {
                boost::multi_array_ref<Value, 2> edge_list = get_array<Value, 2>(aedge_list);
                unordered_map<Value, size_t> vertices;

                if (edge_list.shape()[1] < 2)
                    throw GraphException("Second dimension in edge list must be of size (at least) two");

                typedef typename graph_traits<Graph>::edge_descriptor edge_t;
                vector<DynamicPropertyMapWrap<Value, edge_t>> eprops;
                python::stl_input_iterator<boost::any> iter(oeprops), end;
                for (; iter != end; ++iter)
                    eprops.emplace_back(*iter, writable_edge_properties());

                auto get_vertex = [&] (const Value& r) -> size_t
                    {
                        auto iter = vertices.find(r);
                        if (iter == vertices.end())
                        {
                            auto v = add_vertex(g);
                            vertices[r] = v;
                            vmap[v] = lexical_cast<typename property_traits<VProp>::value_type>(r);
                            return v;
                        }
                        return iter->second;
                    };

                for (const auto& e : edge_list)
                {
                    size_t s = get_vertex(e[0]);
                    size_t t = get_vertex(e[1]);
                    auto ne = add_edge(vertex(s, g), vertex(t, g), g).first;
                    for (size_t i = 0; i < e.size() - 2; ++i)
                    {
                        try
                        {
                            put(eprops[i], ne, e[i + 2]);
                        }
                        catch(bad_lexical_cast&)
                        {
                            throw ValueException("Invalid edge property value: " +
                                                 lexical_cast<string>(e[i + 2]));
                        }
                    }
                }
                found = true;
            }
            catch (invalid_numpy_conversion& e) {}
        }

        template <class Graph, class VProp>
        void operator()(Graph& g, python::object& edge_list, VProp& vmap,
                        bool& found, python::object& oeprops, std::string) const
        {
            if (found)
                return;
            try
            {
                unordered_map<std::string, size_t> vertices;

                typedef typename graph_traits<Graph>::edge_descriptor edge_t;
                vector<DynamicPropertyMapWrap<python::object, edge_t>> eprops;
                python::stl_input_iterator<boost::any> piter(oeprops), pend;
                for (; piter != pend; ++piter)
                    eprops.emplace_back(*piter, writable_edge_properties());

                auto get_vertex = [&] (const std::string& r) -> size_t
                    {
                        auto iter = vertices.find(r);
                        if (iter == vertices.end())
                        {
                            auto v = add_vertex(g);
                            vertices[r] = v;
                            vmap[v] = lexical_cast<typename property_traits<VProp>::value_type>(r);
                            return v;
                        }
                        return iter->second;
                    };

                python::stl_input_iterator<python::object> iter(edge_list), end;
                for (; iter != end; ++iter)
                {
                    const auto& row = *iter;

                    python::stl_input_iterator<python::object> eiter(row), eend;

                    size_t s = 0;
                    size_t t = 0;

                    typename graph_traits<Graph>::edge_descriptor e;
                    size_t i = 0;
                    for(; eiter != eend; ++eiter)
                    {
                        if (i >= eprops.size() + 2)
                            break;
                        const auto& val = *eiter;
                        switch (i)
                        {
                        case 0:
                            s = get_vertex(python::extract<std::string>(val));
                            while (s >= num_vertices(g))
                                add_vertex(g);
                            break;
                        case 1:
                            t = get_vertex(python::extract<std::string>(val));
                            while (t >= num_vertices(g))
                                add_vertex(g);
                            e = add_edge(vertex(s, g), vertex(t, g), g).first;
                            break;
                        default:
                            try
                            {
                                put(eprops[i - 2], e, val);
                            }
                            catch(bad_lexical_cast&)
                            {
                                throw ValueException("Invalid edge property value: " +
                                                     python::extract<string>(python::str(val))());
                            }
                        }
                        i++;
                    }
                }
                found = true;
            }
            catch (invalid_numpy_conversion& e) {}
        }

        template <class Graph, class VProp>
        void operator()(Graph& g, python::object& edge_list, VProp& vmap,
                        bool& found, python::object& oeprops, python::object) const
        {
            if (found)
                return;
            try
            {
                unordered_map<python::object, size_t> vertices;

                typedef typename graph_traits<Graph>::edge_descriptor edge_t;
                vector<DynamicPropertyMapWrap<python::object, edge_t>> eprops;
                python::stl_input_iterator<boost::any> piter(oeprops), pend;
                for (; piter != pend; ++piter)
                    eprops.emplace_back(*piter, writable_edge_properties());

                auto get_vertex = [&] (const python::object& r) -> size_t
                    {
                        auto iter = vertices.find(r);
                        if (iter == vertices.end())
                        {
                            auto v = add_vertex(g);
                            vertices[r] = v;
                            vmap[v] = python::extract<typename property_traits<VProp>::value_type>(r);
                            return v;
                        }
                        return iter->second;
                    };

                python::stl_input_iterator<python::object> iter(edge_list), end;
                for (; iter != end; ++iter)
                {
                    const auto& row = *iter;

                    python::stl_input_iterator<python::object> eiter(row), eend;

                    size_t s = 0;
                    size_t t = 0;

                    typename graph_traits<Graph>::edge_descriptor e;
                    size_t i = 0;
                    for(; eiter != eend; ++eiter)
                    {
                        if (i >= eprops.size() + 2)
                            break;
                        const auto& val = *eiter;
                        switch (i)
                        {
                        case 0:
                            s = get_vertex(val);
                            while (s >= num_vertices(g))
                                add_vertex(g);
                            break;
                        case 1:
                            t = get_vertex(val);
                            while (t >= num_vertices(g))
                                add_vertex(g);
                            e = add_edge(vertex(s, g), vertex(t, g), g).first;
                            break;
                        default:
                            try
                            {
                                put(eprops[i - 2], e, val);
                            }
                            catch(bad_lexical_cast&)
                            {
                                throw ValueException("Invalid edge property value: " +
                                                     python::extract<string>(python::str(val))());
                            }
                        }
                        i++;
                    }
                }
                found = true;
            }
            catch (invalid_numpy_conversion& e) {}
        }
    };
};

void do_add_edge_list_hashed(GraphInterface& gi, python::object aedge_list,
                             boost::any& vertex_map, bool is_str,
                             python::object eprops)
{
    typedef mpl::vector<bool, char, uint8_t, uint16_t, uint32_t, uint64_t,
                        int8_t, int16_t, int32_t, int64_t, uint64_t, double,
                        long double> vals_t;
    bool found = false;
    run_action<graph_tool::detail::all_graph_views, boost::mpl::true_>()
        (gi, std::bind(add_edge_list_hash<vals_t>(), placeholders::_1,
                       aedge_list, placeholders::_2, std::ref(found),
                       is_str, std::ref(eprops)),
         writable_vertex_properties())(vertex_map);
}


struct add_edge_list_iter
{
    template <class Graph>
    void operator()(Graph& g, python::object& edge_list,
                    python::object& oeprops) const
    {
        typedef typename graph_traits<Graph>::edge_descriptor edge_t;
        vector<DynamicPropertyMapWrap<python::object, edge_t>> eprops;
        python::stl_input_iterator<boost::any> piter(oeprops), pend;
        for (; piter != pend; ++piter)
            eprops.emplace_back(*piter, writable_edge_properties());

        python::stl_input_iterator<python::object> iter(edge_list), end;
        for (; iter != end; ++iter)
        {
            const auto& row = *iter;
            python::stl_input_iterator<python::object> eiter(row), eend;

            size_t s = 0;
            size_t t = 0;

            typename graph_traits<Graph>::edge_descriptor e;
            size_t i = 0;
            for(; eiter != eend; ++eiter)
            {
                if (i >= eprops.size() + 2)
                    break;
                const auto& val = *eiter;
                switch (i)
                {
                case 0:
                    s = python::extract<size_t>(val);
                    while (s >= num_vertices(g))
                        add_vertex(g);
                    break;
                case 1:
                    t = python::extract<size_t>(val);
                    while (t >= num_vertices(g))
                        add_vertex(g);
                    e = add_edge(vertex(s, g), vertex(t, g), g).first;
                    break;
                default:
                    try
                    {
                        put(eprops[i - 2], e, val);
                    }
                    catch(bad_lexical_cast&)
                    {
                        throw ValueException("Invalid edge property value: " +
                                             python::extract<string>(python::str(val))());
                    }
                }
                i++;
            }
        }
    }
};

void do_add_edge_list_iter(GraphInterface& gi, python::object edge_list,
                           python::object eprops)
{
    run_action<>()
        (gi, std::bind(add_edge_list_iter(), placeholders::_1,
                       std::ref(edge_list), std::ref(eprops)))();
}


} // namespace graph_tool

// register everything

void export_python_properties();

python::list* _vlist(0);
python::list* _elist(0);

python::list get_vlist()
{
    if (_vlist == nullptr)
        _vlist = new python::list();
    return *_vlist;
}

python::list get_elist()
{
    if (_elist == nullptr)
        _elist = new python::list();
    return *_elist;
}

void export_python_interface()
{
    using namespace boost::python;

    class_<VertexBase>("VertexBase", no_init);
    class_<EdgeBase>("EdgeBase", no_init);

    typedef boost::mpl::transform<graph_tool::detail::all_graph_views,
                                  boost::mpl::quote1<std::add_const> >::type const_graph_views;
    typedef boost::mpl::transform<graph_tool::detail::all_graph_views,
                                  boost::mpl::quote1<std::add_pointer> >::type all_graph_views;
    typedef boost::mpl::transform<const_graph_views,
                                  boost::mpl::quote1<std::add_pointer> >::type all_const_graph_views;
    typedef boost::mpl::joint_view<all_graph_views, all_const_graph_views>::type graph_views;
    boost::mpl::for_each<graph_views>(std::bind(graph_tool::export_python_interface(),
                                                placeholders::_1, get_vlist(),
                                                get_elist(), graph_views()));
    export_python_properties();
    def("new_vertex_property",
        &new_property<GraphInterface::vertex_index_map_t>);
    def("new_edge_property", &new_property<GraphInterface::edge_index_map_t>);
    def("new_graph_property",
        &new_property<ConstantPropertyMap<size_t,graph_property_tag> >);

    def("get_vertex", get_vertex);
    def("get_vertices", get_vertices);
    def("get_edges", get_edges);
    def("add_vertex", graph_tool::add_vertex);
    def("add_edge", graph_tool::add_edge);
    def("remove_vertex", graph_tool::remove_vertex);
    def("remove_vertex_array", graph_tool::remove_vertex_array);
    def("clear_vertex", graph_tool::clear_vertex);
    def("remove_edge", graph_tool::remove_edge);
    def("add_edge_list", graph_tool::do_add_edge_list);
    def("add_edge_list_hashed", graph_tool::do_add_edge_list_hashed);
    def("add_edge_list_iter", graph_tool::do_add_edge_list_iter);
    def("get_edge", get_edge);

    def("get_vertex_index", get_vertex_index);
    def("get_edge_index", do_get_edge_index);

    def("get_vlist", get_vlist);
    def("get_elist", get_elist);
}
