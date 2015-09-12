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

#ifndef PYTHON_INTERFACE_HH
#define PYTHON_INTERFACE_HH

#include <boost/python.hpp>
#include <boost/python/type_id.hpp>

#include <functional>

namespace std
{
    template<>
    struct hash<boost::python::object>
    {
        size_t operator()(const boost::python::object& o) const
        {
            return std::hash<int64_t>()(boost::python::extract<int64_t>(o.attr("__hash__")()));
        }
    };
}

#include <boost/graph/graph_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <type_traits>

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_selectors.hh"
#include "numpy_bind.hh"


// this file includes a simple python interface for the internally kept
// graph. It defines a PythonVertex, PythonEdge and PythonIterator template
// classes, which contain the proper member functions for graph traversal. These
// types are then specialized for each version of the adapted graph (directed,
// undirected, filtered, reversed).

namespace graph_tool
{

// generic iterator adaptor which can be used to iterate vertices, edges,
// out_edges and in_edges through python
template <class Graph, class Descriptor, class Iterator>
class PythonIterator
{
public:
    PythonIterator(std::shared_ptr<Graph>& gp,
                   std::pair<Iterator,Iterator> e)
        : _g(gp), _e(e) {}
    Descriptor Next()
    {
        if (_e.first == _e.second)
            boost::python::objects::stop_iteration_error();
        Descriptor e(_g, *_e.first);
        ++_e.first;
        return e;
    }
private:
    std::shared_ptr<Graph> _g;
    std::pair<Iterator,Iterator> _e;
};


// forward declaration of PythonEdge
template <class Graph>
class PythonEdge;

class VertexBase {}; // useful to unite all vertex

// below are classes related to the PythonVertex type
template <class Graph>
class PythonVertex : public VertexBase
{
public:
    PythonVertex(std::shared_ptr<Graph> g, GraphInterface::vertex_t v):
        _g(g), _v(v) {}

    bool IsValid() const
    {
        std::shared_ptr<Graph> gp(_g);
        Graph* g = gp.get();
        if (g == nullptr)
            return false;
        return ((_v != boost::graph_traits<Graph>::null_vertex()) &&
                (_v < num_vertices(*g)));
    }

    void CheckValid() const
    {
        if (!IsValid())
            throw ValueException("invalid vertex descriptor: " +
                                 boost::lexical_cast<string>(_v));
    }

    GraphInterface::vertex_t GetDescriptor() const
    {
        return _v;
    }

    template <class DegSelector>
    struct get_degree
    {
        void operator()(const Graph& g,
                        typename boost::graph_traits<Graph>::vertex_descriptor v,
                        size_t& deg) const
        {
            deg = DegSelector()(v, g);
        }

        template<class PMap>
        void operator()(const Graph& g,
                        typename boost::graph_traits<Graph>::vertex_descriptor v,
                        const boost::any& aweight, boost::python::object& deg,
                        bool& found, PMap) const
        {
            try
            {
                const PMap& weight = boost::any_cast<const PMap&>(aweight);
                deg = boost::python::object(DegSelector()(v, g, weight));
                found = true;
            }
            catch (boost::bad_any_cast&) {}
        }
    };

    size_t GetInDegree() const
    {
        CheckValid();
        std::shared_ptr<Graph> gp(_g);
        Graph& g = *gp.get();
        size_t in_deg;
        get_degree<in_degreeS>()(g, _v, in_deg);
        return in_deg;
    }

    boost::python::object GetWeightedInDegree(boost::any pmap) const
    {
        std::shared_ptr<Graph> gp(_g);
        Graph& g = *gp.get();
        boost::python::object in_deg;
        bool found = false;
        boost::mpl::for_each<edge_scalar_properties>(std::bind(get_degree<in_degreeS>(),
                                                               std::ref(g), _v,
                                                               std::ref(pmap),
                                                               std::ref(in_deg),
                                                               std::ref(found),
                                                               std::placeholders::_1));
        if (!found)
            throw ValueException("edge weight property must be of scalar type");
        return in_deg;
    }

    size_t GetOutDegree() const
    {
        CheckValid();
        std::shared_ptr<Graph> gp(_g);
        Graph& g = *gp.get();
        size_t out_deg;
        get_degree<out_degreeS>()(g, _v, out_deg);
        return out_deg;
    }


    boost::python::object GetWeightedOutDegree(boost::any pmap) const
    {
        std::shared_ptr<Graph> gp(_g);
        Graph& g = *gp.get();
        boost::python::object out_deg;
        bool found = false;
        boost::mpl::for_each<edge_scalar_properties>(std::bind(get_degree<out_degreeS>(),
                                                               std::ref(g), _v,
                                                               std::ref(pmap),
                                                               std::ref(out_deg),
                                                               std::ref(found),
                                                               std::placeholders::_1));
        if (!found)
            throw ValueException("edge weight property must be of scalar type");
        return out_deg;
    }

    // provide iterator support for out_edges
    boost::python::object OutEdges() const
    {
        CheckValid();
        std::shared_ptr<Graph> pg(_g);
        Graph& g = *pg;
        typedef typename boost::graph_traits<Graph>::out_edge_iterator
            out_edge_iterator;
        return boost::python::object(PythonIterator<Graph,PythonEdge<Graph>,
                                                    out_edge_iterator>
                                     (pg, out_edges(_v, g)));
    }

    boost::python::object InEdges() const
    {
        CheckValid();
        std::shared_ptr<Graph> pg(_g);
        Graph& g = *pg;
        typedef typename in_edge_iteratorS<Graph>::type
            in_edge_iterator;
        return boost::python::object(PythonIterator<Graph, PythonEdge<Graph>,
                                                    in_edge_iterator>
                                     (pg, in_edge_iteratorS<Graph>::get_edges(_v, g)));
    }

    std::string GetString() const
    {
        CheckValid();
        return boost::lexical_cast<std::string>(_v);
    }

    size_t GetHash() const
    {
        return std::hash<size_t>()(_v);
    }

    size_t GetIndex() const
    {
        return _v;
    }

    size_t GetGraphPtr() const
    {
        std::shared_ptr<Graph> pg(_g);
        return size_t(pg.get());
    }

    std::string GetGraphType() const
    {
        using boost::python::detail::gcc_demangle;
        return gcc_demangle(typeid(Graph).name());
    }

    template <class OGraph>
    bool operator==(const PythonVertex<OGraph>& other) const { return _v == other._v; }
    template <class OGraph>
    bool operator!=(const PythonVertex<OGraph>& other) const { return _v != other._v; }
    template <class OGraph>
    bool operator<(const PythonVertex<OGraph>& other) const { return _v < other._v; }
    template <class OGraph>
    bool operator<=(const PythonVertex<OGraph>& other) const { return _v <= other._v; }
    template <class OGraph>
    bool operator>(const PythonVertex<OGraph>& other) const { return _v > other._v; }
    template <class OGraph>
    bool operator>=(const PythonVertex<OGraph>& other) const { return _v >= other._v; }

private:
    std::weak_ptr<Graph> _g;
    GraphInterface::vertex_t _v;
};

// below are classes related to the PythonEdge type

class EdgeBase {}; // useful to unite all edge types

template <class Graph>
class PythonEdge : public EdgeBase
{
public:
    typedef typename boost::graph_traits<Graph>::edge_descriptor edge_descriptor;
    PythonEdge(std::shared_ptr<Graph> g, edge_descriptor e)
        : _g(g), _e(e) {}

    bool IsValid() const
    {
        std::shared_ptr<Graph> gp(_g);
        Graph* g = gp.get();
        if (g == nullptr)
            return false;

        auto s = source(_e, *g);
        auto t = target(_e, *g);

        return ((s != boost::graph_traits<Graph>::null_vertex()) &&
                (s < num_vertices(*g)) &&
                (t != boost::graph_traits<Graph>::null_vertex()) &&
                (t < num_vertices(*g)));
    }

    void CheckValid() const
    {
        if (!IsValid())
            throw ValueException("invalid edge descriptor");
    }

    GraphInterface::edge_t GetDescriptor() const
    {
        return _e;
    }

    PythonVertex<Graph> GetSource() const
    {
        CheckValid();
        std::shared_ptr<Graph> pg(_g);
        Graph& g = *pg;
        return PythonVertex<Graph>(pg, source(_e, g));
    }

    PythonVertex<Graph> GetTarget() const
    {
        CheckValid();
        std::shared_ptr<Graph> pg(_g);
        Graph& g = *pg;
        return PythonVertex<Graph>(pg, target(_e, g));
    }

    std::string GetString() const
    {
        CheckValid();
        Graph& g = *std::shared_ptr<Graph>(_g);
        auto s = source(_e, g);
        auto t = target(_e, g);
        return "(" + boost::lexical_cast<std::string>(s) + ", "
            + boost::lexical_cast<std::string>(t) + ")";
    }

    size_t GetHash() const
    {
        CheckValid();
        Graph& g = *std::shared_ptr<Graph>(_g);
        auto eindex = get(boost::edge_index_t(), g);
        return std::hash<size_t>()(eindex[_e]);
    }

    size_t GetGraphPtr() const
    {
        std::shared_ptr<Graph> pg(_g);
        return site_t(pg.get());
    }

    std::string GetGraphType() const
    {
        using boost::python::detail::gcc_demangle;
        return gcc_demangle(typeid(Graph).name());
    }

    template <class OGraph>
    bool operator==(const PythonEdge<OGraph>& other) const { return _e == other._e; }
    template <class OGraph>
    bool operator!=(const PythonEdge<OGraph>& other) const { return !(*this == other); }
    template <class OGraph>
    bool operator<(const PythonEdge<OGraph>& other)  const
    {
        CheckValid();
        other.CheckValid();
        Graph& g = *std::shared_ptr<Graph>(_g);
        OGraph& og = *std::shared_ptr<OGraph>(other._g);
        auto s = source(_e, g);
        auto t = target(_e, g);
        auto os = source(other._e, og);
        auto ot = target(other._e, og);
        if (not is_directed::apply<Graph>::type::value ||
            not is_directed::apply<OGraph>::type::value)
        {
            if (t < s)
                std::swap(s, t);
            if (ot < os)
                std::swap(os, ot);
        }
        return s < os && t < ot;
    }
    template <class OGraph>
    bool operator<=(const PythonEdge<OGraph>& other) const {return *this < other || *this == other;}
    template <class OGraph>
    bool operator>(const PythonEdge<OGraph>& other) const {return !(*this < other || *this == other);}
    template <class OGraph>
    bool operator>=(const PythonEdge<OGraph>& other) const {return *this > other || *this == other;}

private:
    std::weak_ptr<Graph> _g;
    edge_descriptor _e;

    template <class OGraph>
    friend class PythonEdge;
};

// metafunction to determine wether or not to return copies or internal
// references to property types
struct return_reference
{
    template <class ValueType>
    struct apply
    {
        // return actual references only for non-string and non-python::object
        // classes

        typedef typename boost::mpl::if_<
            typename boost::mpl::and_<
                std::is_class<ValueType>,
                typename boost::mpl::and_<
                    typename boost::mpl::not_<std::is_same<ValueType,
                                                           string> >::type,
                    typename boost::mpl::not_<std::is_same<ValueType,
                                                           boost::python::object> >::type>::type
                >::type,
            boost::mpl::bool_<true>,
            boost::mpl::bool_<false> >::type type;
    };
};

template <class PropertyMap>
class PythonPropertyMap
{
public:
    PythonPropertyMap(const PropertyMap& pmap)
        : _pmap(pmap) {}

    typedef typename boost::property_traits<PropertyMap>::value_type value_type;

    typedef typename boost::mpl::if_<
        typename return_reference::apply<value_type>::type,
        value_type&,
        value_type>::type reference;

    template <class PythonDescriptor>
    reference GetValue(const PythonDescriptor& key)
    {
        key.CheckValid();
        return get(_pmap, key.GetDescriptor());
    }

    // in this case, val should be a copy, not a reference. This is to avoid a
    // problem with vector-valued property maps
    template <class PythonDescriptor>
    void SetValue(const PythonDescriptor& key, value_type val)
    {
        set_value(key, val,
                  std::is_convertible<typename boost::property_traits<PropertyMap>::category,
                                      boost::writable_property_map_tag>());
    }

    template <class PythonDescriptor>
    void set_value(const PythonDescriptor& key, const value_type& val,
                   std::true_type)
    {
        key.CheckValid();
        put(_pmap, key.GetDescriptor(), val);
    }

    template <class PythonDescriptor>
    void set_value(const PythonDescriptor&, const value_type&,
                   std::false_type)
    {
        throw ValueException("property is read-only");
    }

    size_t GetHash() const
    {
        return std::hash<size_t>()(size_t(this));
    }

    std::string GetType() const
    {
        using boost::python::detail::gcc_demangle;
        if (std::is_same<typename boost::mpl::find<value_types,value_type>::type,
                         typename boost::mpl::end<value_types>::type>::value)
            return gcc_demangle(typeid(value_type).name());
        else
            return type_names[boost::mpl::find<value_types,
                                               value_type>::type::pos::value];
    }

    boost::any GetMap() const
    {
        return _pmap;
    }

    boost::any GetDynamicMap() const
    {
        return (boost::dynamic_property_map*)
            (new boost::detail::dynamic_property_map_adaptor<PropertyMap>
             (_pmap));
    }

    boost::python::object GetArray(size_t size)
    {
        typedef typename boost::mpl::or_<
            typename boost::mpl::or_<
                std::is_same<PropertyMap,
                             GraphInterface::vertex_index_map_t>,
                std::is_same<PropertyMap,
                             GraphInterface::edge_index_map_t> >::type,
            typename boost::mpl::not_<
                typename boost::mpl::has_key<numpy_types, value_type>::type >
            ::type>::type isnt_vector_map;
        return get_array(size, isnt_vector_map());
    }

    boost::python::object get_array(size_t size, boost::mpl::bool_<false>)
    {
        _pmap.reserve(size);
        return wrap_vector_not_owned(_pmap.get_storage());
    }

    boost::python::object get_array(size_t, boost::mpl::bool_<true>)
    {
        return boost::python::object();
    }

    bool IsWritable() const
    {
        return std::is_convertible<typename boost::property_traits<PropertyMap>::category,
                                   boost::writable_property_map_tag>::value;
    }

private:
    PropertyMap _pmap; // hold an internal copy, since it's cheap
};


//
// Create new properties
//

struct new_property_map
{
    template <class ValueType, class IndexMap>
    void operator()(ValueType, IndexMap index, const string& type_name,
                     boost::any pmap, boost::python::object& new_prop, bool& found) const
    {
        size_t i = boost::mpl::find<value_types,ValueType>::type::pos::value;
        if (type_name == type_names[i])
        {
            typedef typename property_map_type::apply<ValueType, IndexMap>::type
                map_t;
            map_t prop;
            if (pmap.empty())
                prop = map_t(index);
            else
                prop = boost::any_cast<map_t>(pmap);

            new_prop = boost::python::object(PythonPropertyMap<map_t>(prop));
            found = true;
        }
    }
};

template <class IndexMap>
boost::python::object new_property(const string& type, IndexMap index_map,
                                   boost::any pmap)
{
    boost::python::object prop;
    bool found = false;
    boost::mpl::for_each<value_types>(std::bind(new_property_map(),
                                                std::placeholders::_1, index_map,
                                                std::ref(type), pmap, std::ref(prop),
                                                std::ref(found)));
    if (!found)
        throw ValueException("Invalid property type: " + type);
    return prop;
}

//
// Python IO streams (minimal access to c++ streams)
//

class OStream
{
public:
    OStream(std::ostream& s): _s(s) {}

    void Write(const std::string& s, size_t n)
    {
        _s.write(s.c_str(), long(n));
    }

    void Flush()
    {
        _s.flush();
    }

private:
    std::ostream& _s;
};

class IStream
{
public:
    IStream(std::istream& s): _s(s) {}

    boost::python::object Read(size_t n)
    {
        std::string buf;
        buf.resize(n);
        _s.read(&buf[0], n);
        buf.resize(_s.gcount());

#if (PY_MAJOR_VERSION >= 3)
        // in python 3 we need to construct a 'bytes' instance
        PyObject* bytes = PyBytes_FromStringAndSize(&buf[0], buf.size());
        boost::python::handle<> x(bytes);
        boost::python::object pbuf(x);
#else
        boost::python::str pbuf(buf);
#endif
        return pbuf;
    }

private:
    std::istream& _s;
};


} //graph_tool namespace

#endif
