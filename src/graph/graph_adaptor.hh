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

#ifndef GRAPH_ADAPTOR_HH
#define GRAPH_ADAPTOR_HH

#include <list>

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/range/join.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "transform_iterator.hh"

namespace boost {

//==============================================================================
// UndirectedAdaptor
// This class encapsulates a directed graph with parallel edges and provides a
// view of the graph as undirected with parallel edges.
// Encapsulated graph can be: VertexListGraph, EdgeListGraph, IncidenceGraph,
//                            AdjacencyGraph, VertexMutableGraph,
//                            EdgeMutableGraph, VertexMutablePropertyGraph,
//                            EdgeMutablePropertyGraph, BidirectionalGraph
// The undirected graph obeys the same concepts.
//==============================================================================
template <class Graph> class UndirectedAdaptor
{
public:
    UndirectedAdaptor(const Graph& g) : _g(const_cast<Graph&>(g)){}

    typedef typename vertex_property_type<Graph>::type vertex_property_type;
    typedef typename edge_property_type<Graph>::type edge_property_type;
    typedef typename Graph::graph_tag graph_tag;
    typedef Graph graph_type;

    typedef Graph original_graph_t;

    typedef typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
        vertex_descriptor;
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor
        edge_descriptor;

    typedef undirected_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef typename graph_traits<Graph>::traversal_category traversal_category;

#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
    // Bundled properties support
    template<typename Descriptor>
    typename graph::detail::bundled_result<Graph, Descriptor>::type&
        operator[](Descriptor x) { return this->m_g[x]; }

    template<typename Descriptor>
    typename graph::detail::bundled_result<Graph, Descriptor>::type const&
        operator[](Descriptor x) const { return this->m_g[x]; }
#endif // BOOST_GRAPH_NO_BUNDLED_PROPERTIES

    const Graph& original_graph() const {return _g;}
    Graph& original_graph() {return _g;}

    static vertex_descriptor null_vertex() {graph_traits<Graph>::null_vertex();}

private:
    Graph& _g;
};

#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
template<typename Graph>
struct vertex_bundle_type<UndirectedAdaptor<Graph> >:
        vertex_bundle_type<Graph> { };

template<typename Graph>
struct edge_bundle_type<UndirectedAdaptor<Graph> >:
        edge_bundle_type<Graph> { };
#endif // BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template <class Graph>
struct get_iterator_category
{
    typedef typename graph_traits<Graph>::out_edge_iterator iter_t;
    typedef typename std::iterator_traits<iter_t>::iterator_category type;
};


template <class Graph, class Inverted = mpl::false_>
class joined_edge_iterator
    : public boost::iterator_facade<joined_edge_iterator<Graph, Inverted>,
                                    typename graph_traits<Graph>::edge_descriptor,
                                    typename get_iterator_category<Graph>::type,
                                    typename graph_traits<Graph>::edge_descriptor>
{
 public:
    typedef typename graph_traits<Graph>::in_edge_iterator in_iter_t;
    typedef typename graph_traits<Graph>::out_edge_iterator out_iter_t;

    typedef typename mpl::if_<Inverted, out_iter_t, in_iter_t>::type iter1_t;
    typedef typename mpl::if_<Inverted, in_iter_t, out_iter_t>::type iter2_t;

    joined_edge_iterator() {}
    explicit joined_edge_iterator(std::pair<iter1_t, iter1_t>&& range1,
                                  std::pair<iter2_t, iter2_t>&& range2,
                                  bool begin)
        : _range1(std::forward<std::pair<iter1_t, iter1_t>>(range1)),
          _range2(std::forward<std::pair<iter2_t, iter2_t>>(range2))
    {
        if (!begin)
        {
            _range1.first = _range1.second;
            _range2.first = _range2.second;
        }
    }

 private:
    friend class boost::iterator_core_access;
    void increment()
    {
        if (_range1.first == _range1.second)
            ++_range2.first;
        else
            ++_range1.first;
    }

    typedef typename std::iterator_traits<iter1_t>::difference_type diff_t;
    void advance(diff_t n)
    {
        diff_t d1 = _range1.second - _range1.first;
        if (n < d1)
        {
            _range1.first += n;
        }
        else
        {
            _range1.first = _range1.second;
            _range2.first += n - d1;
        }
    }

    diff_t distance_to(joined_edge_iterator const& other)
    {
        return (other._range1.first - _range1.first) +
            (other._range2.first - _range2.first);
    }

    bool equal(joined_edge_iterator const& other) const
    {
        return (_range2.first == other._range2.first &&
                _range1.first == other._range1.first);
    }

    template <class Edge>
    Edge inv(Edge&& e) const
    {
        e.inv = true;
        return e;
    }

    typename graph_traits<Graph>::edge_descriptor dereference() const
    {
        if (_range1.first == _range1.second)
            return *_range2.first;
        else
            return inv(*_range1.first);
    }

    std::pair<iter1_t, iter1_t> _range1;
    std::pair<iter2_t, iter2_t> _range2;
};



//==============================================================================
// UndirectedAdaptorAdjacencyIterator
// just keeps an internal reference to out_edge_iterator and calls target() when
// referenced
//==============================================================================

template <class Graph>
struct get_undirected_adjacency_iterator
{
    typedef joined_edge_iterator<Graph> out_edge_iter_t;
    typedef typename boost::adjacency_iterator_generator<UndirectedAdaptor<Graph>,
                                                         typename graph_traits<Graph>::vertex_descriptor,
                                                         out_edge_iter_t>::type type;
};


//==============================================================================
// graph_traits<UndirectedAdaptor>
// this defines all the necessary types associated with UndirectedAdaptor
//==============================================================================
template <class Graph>
struct graph_traits<UndirectedAdaptor<Graph> > {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;

    typedef typename get_undirected_adjacency_iterator<Graph>::type adjacency_iterator;
    typedef joined_edge_iterator<Graph, mpl::false_> out_edge_iterator;
    typedef joined_edge_iterator<Graph, mpl::true_> in_edge_iterator;
    typedef typename graph_traits<Graph>::vertex_iterator vertex_iterator;
    typedef typename graph_traits<Graph>::edge_iterator edge_iterator;


    typedef undirected_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef typename graph_traits<Graph>::traversal_category traversal_category;
    typedef typename graph_traits<Graph>::vertices_size_type vertices_size_type;
    typedef typename graph_traits<Graph>::edges_size_type edges_size_type;
    typedef typename graph_traits<Graph>::degree_size_type degree_size_type;

    static vertex_descriptor null_vertex()
    {
        return graph_traits<Graph>::null_vertex();
    }

private:
    typedef is_convertible<typename std::iterator_traits<typename graph_traits<Graph>::out_edge_iterator>::iterator_category,
                           std::random_access_iterator_tag> is_orig_ra;
    typedef is_convertible<typename std::iterator_traits<out_edge_iterator>::iterator_category,
                           std::random_access_iterator_tag> is_ra;
    BOOST_STATIC_ASSERT((!is_orig_ra::value || is_ra::value));
};

template <class Graph>
struct graph_traits< const UndirectedAdaptor<Graph> >:
    public graph_traits<UndirectedAdaptor<Graph> > {};

//==============================================================================
// Nonmember functions
// these provide manipulation of the graph
//==============================================================================

//==============================================================================
// source(e,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
source(const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
       const UndirectedAdaptor<Graph>& g)
{
    if (e.inv)
        return target(e, g.original_graph());
    else
        return source(e, g.original_graph());
}

//==============================================================================
// target(e,g)
//==============================================================================
template <class Graph>
inline typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
target(const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
       const UndirectedAdaptor<Graph>& g)
{
    if (e.inv)
        return source(e, g.original_graph());
    else
        return target(e, g.original_graph());
}

//==============================================================================
// vertex(n,g)
//==============================================================================
template <class Graph>
inline
typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type n,
       const UndirectedAdaptor<Graph>& g)
{
    return vertex(n, g.original_graph());
}

//==============================================================================
// vertices(g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::vertex_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::vertex_iterator >
vertices(const UndirectedAdaptor<Graph>& g)
{
    return vertices(g.original_graph());
}

//==============================================================================
// edges(g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::edge_iterator >
edges(const UndirectedAdaptor<Graph>& g)
{
    return edges(g.original_graph());
}

//==============================================================================
// edge(u, v, g)
//==============================================================================
template <class Graph>
inline
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
     const UndirectedAdaptor<Graph>& g)
{
    auto res = edge(u, v, g.original_graph());

    if (!res.second)
    {
        res = edge(v, u, g.original_graph());
        res.first.inv = true;
    }

    return res;
}

//==============================================================================
// out_edges(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph>>::out_edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph>>::out_edge_iterator >
out_edges(typename graph_traits<UndirectedAdaptor<Graph>>::vertex_descriptor u,
          const UndirectedAdaptor<Graph>& g)
{
    typedef joined_edge_iterator<Graph> iter_t;
    return std::make_pair(iter_t(in_edges(u, g.original_graph()),
                                 out_edges(u, g.original_graph()), true),
                          iter_t(in_edges(u, g.original_graph()),
                                 out_edges(u, g.original_graph()), false));
}

//==============================================================================
// in_edges(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph>>::in_edge_iterator,
          typename graph_traits<UndirectedAdaptor<Graph>>::in_edge_iterator >
in_edges(typename graph_traits<UndirectedAdaptor<Graph>>::vertex_descriptor u,
         const UndirectedAdaptor<Graph>& g)
{
    typedef joined_edge_iterator<Graph, mpl::true_> iter_t;
    return std::make_pair(iter_t(out_edges(u, g.original_graph()),
                                 in_edges(u, g.original_graph()), true),
                          iter_t(out_edges(u, g.original_graph()),
                                 in_edges(u, g.original_graph()), false));
}

//==============================================================================
// adjacent_vertices(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator,
          typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator>
adjacent_vertices
    (typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
     const UndirectedAdaptor<Graph>& g)
{
    typedef typename graph_traits<UndirectedAdaptor<Graph> >::adjacency_iterator
        adjacency_iterator;
    auto e_range = out_edges(u, g);
    return std::make_pair(adjacency_iterator(e_range.first, &g),
                          adjacency_iterator(e_range.second, &g));
}

//==============================================================================
// num_vertices(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertices_size_type
num_vertices(const UndirectedAdaptor<Graph>& g)
{
    return num_vertices(g.original_graph());
}

//==============================================================================
// num_edges(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::edges_size_type
num_edges(const UndirectedAdaptor<Graph>& g)
{
    return num_edges(g.original_graph());
}

//==============================================================================
// out_degree(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
out_degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
           const UndirectedAdaptor<Graph>& g)
{
    return (out_degree(u, g.original_graph()) +
            in_degree(u, g.original_graph()));
}

//==============================================================================
// in_degree(u,g)
//==============================================================================
template <class Graph>
inline
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
in_degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
          const UndirectedAdaptor<Graph>& g)
{
    return out_degree(u, g);
}

//==============================================================================
// degree(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::degree_size_type
degree(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
       const UndirectedAdaptor<Graph>& g)
{
    return out_degree(u, g);
}


//==============================================================================
// add_vertex(g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(UndirectedAdaptor<Graph>& g)
{
    return add_vertex(g.original_graph());
}

//==============================================================================
// add_vertex(vp,g)
//==============================================================================
template <class Graph, class VertexProperties>
inline __attribute__((always_inline))
typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor
add_vertex(const VertexProperties& p, UndirectedAdaptor<Graph>& g)
{
    return add_vertex(p, g.original_graph());
}

//==============================================================================
// clear_vertex(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void clear_vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                  UndirectedAdaptor<Graph>& g)
{
    clear_vertex(u, g.original_graph());
}

//==============================================================================
// remove_vertex(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_vertex(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                   UndirectedAdaptor<Graph>& g)
{
    remove_vertex(u, g.original_graph());
}

//==============================================================================
// remove_vertex_fast(u,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_vertex_fast(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                        UndirectedAdaptor<Graph>& g)
{
    remove_vertex_fast(u, g.original_graph());
}

//==============================================================================
// add_edge(u,v,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         UndirectedAdaptor<Graph>& g)
{
    return add_edge(u, v, g.original_graph());
}

//==============================================================================
// add_edge(u,v,ep,g)
//==============================================================================
template <class Graph, class EdgeProperties>
inline __attribute__((always_inline))
std::pair<typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor,
          bool>
add_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const EdgeProperties& ep, UndirectedAdaptor<Graph>& g)
{
    return add_edge(u, v, ep, g.original_graph());
}

//==============================================================================
// remove_edge(u,v,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor u,
                 typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
                 UndirectedAdaptor<Graph>& g)
{
    auto e = edge(u, v, g);
    if (e.second)
        remove_edge(e.first, g);
}

//==============================================================================
// remove_edge(e,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge(const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
                 UndirectedAdaptor<Graph>& g)
{
    remove_edge(e, g.original_graph());
}

//==============================================================================
// remove_edge(e_iter,g)
//==============================================================================
template <class Graph>
inline __attribute__((always_inline))
void remove_edge(const typename graph_traits<UndirectedAdaptor<Graph> >::out_edge_iterator& iter,
                 UndirectedAdaptor<Graph>& g)
{
    remove_edge(*iter, g);
}

//==============================================================================
// remove_out_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
inline
void remove_out_edge_if(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
                        Predicate predicate, UndirectedAdaptor<Graph>& g)
{
    std::vector<typename UndirectedAdaptor<Graph>::EdgeDescriptor> removed_edges;
    auto edge_range = out_edges(v,g);
    for(auto iter = edge_range.first; iter != edge_range.second; ++iter)
        if (predicate(*iter))
            removed_edges.push_back(*iter);
    for(auto& e : removed_edges)
        remove_edge(e, g);
}

//==============================================================================
// remove_in_edge_if(v,predicate,g)
//==============================================================================
template <class Graph, class Predicate>
inline
void remove_in_edge_if(typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
                       Predicate predicate, UndirectedAdaptor<Graph>& g)
{
    std::vector<typename UndirectedAdaptor<Graph>::EdgeDescriptor> removed_edges;
    auto edge_range = in_edges(v,g);
    for(auto iter = edge_range.first; iter != edge_range.second; ++iter)
        if (predicate(*iter))
            removed_edges.push_back(*iter);
    for(auto& e : removed_edges)
        remove_edge(e, g);
}


//==============================================================================
// Property maps
//==============================================================================

//==============================================================================
// vertex_property<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class vertex_property<UndirectedAdaptor<Graph> >
{
public:
    typedef typename vertex_property<Graph>::type type;
};

//==============================================================================
// vertex_property_type<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class vertex_property_type<UndirectedAdaptor<Graph> >
{
public:
    typedef typename vertex_property_type<Graph>::type type;
};

//==============================================================================
// edge_property<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class edge_property<UndirectedAdaptor<Graph> >
{
public:
    typedef typename edge_property<Graph>::type type;
};

//==============================================================================
// edge_property_type<UndirectedAdaptor>
//==============================================================================
template <class Graph>
class edge_property_type<UndirectedAdaptor<Graph> >
{
public:
    typedef typename edge_property_type<Graph>::type type;
};

//==============================================================================
// property_map<UndirecterdAdaptor, PropertyTag>
//==============================================================================
template <class Graph, class PropertyTag>
class property_map<UndirectedAdaptor<Graph>, PropertyTag>
{
public:
    typedef typename property_map<Graph, PropertyTag>::type type;
    typedef typename property_map<Graph, PropertyTag>::const_type const_type;
};

//==============================================================================
// property_map<UndirectedAdaptor, T Bundle::*>
//==============================================================================
template <typename Graph, typename T, typename Bundle>
class property_map<UndirectedAdaptor<Graph>, T Bundle::*>
{
public:
    typedef typename property_map<Graph, T Bundle::*>::type type;
    typedef typename property_map<Graph, T Bundle::*>::const_type const_type;
};


//==============================================================================
// get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::type
get(PropertyTag tag, UndirectedAdaptor<Graph>& g)
{
    return get(tag, g.original_graph());
}

//==============================================================================
// const get(tag,g)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_map<UndirectedAdaptor<Graph>, PropertyTag>::const_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g)
{
    return get(tag, g.original_graph());
}

//==============================================================================
// get(tag,g,v)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
    typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v)
{
    return get(tag, g.original_graph(), v);
}

//==============================================================================
// get(tag,g,e)
//==============================================================================
template <class PropertyTag, class Graph>
inline
typename property_traits
    <typename property_map<UndirectedAdaptor<Graph>,
                           PropertyTag>::const_type >::value_type
get(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
    const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e)
{
    return get(tag, g.original_graph(), e);
}

//==============================================================================
// put(tag, g, v, value)
//==============================================================================
template <class Graph, class PropertyTag, class Value>
inline
void put(PropertyTag tag, UndirectedAdaptor<Graph>& g,
         typename graph_traits<UndirectedAdaptor<Graph> >::vertex_descriptor v,
         const Value& value)
{
    put(tag, g.original_graph(), v, value);
}

//==============================================================================
// put(tag, g, e, value)
//==============================================================================
template <class Graph, class PropertyTag, class X, class Value>
inline
void put(PropertyTag tag, const UndirectedAdaptor<Graph>& g,
         const typename graph_traits<UndirectedAdaptor<Graph> >::edge_descriptor& e,
         const Value &value)
{
    put(tag, g.original_graph(), e, value);
}

//==============================================================================
// get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline
typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(UndirectedAdaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.original_graph(), tag);
}

//==============================================================================
// const get_property(g,tag)
//==============================================================================
template <class Graph, class GraphProperties, class GraphPropertyTag>
inline
const typename property_value<GraphProperties, GraphPropertyTag>::type&
get_property(const UndirectedAdaptor<Graph>& g, GraphPropertyTag tag)
{
    get_property(g.original_graph(), tag);
}

} // namespace boost


#endif // GRAPH_ADAPTOR_HH
