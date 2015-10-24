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
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/coroutine/all.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

class BFSVisitorWrapper
{
public:
    BFSVisitorWrapper(GraphInterface& gi, python::object vis)
        : _gi(gi), _vis(vis) {}

    template <class Vertex, class Graph>
    void initialize_vertex(Vertex u, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("initialize_vertex")(PythonVertex<Graph>(gp, u));
    }

    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("discover_vertex")(PythonVertex<Graph>(gp, u));
    }

    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("examine_vertex")(PythonVertex<Graph>(gp, u));
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("examine_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void tree_edge(Edge e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("tree_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void non_tree_edge(Edge e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("non_tree_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void gray_target(Edge e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("gray_target")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void black_target(Edge e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("black_target")(PythonEdge<Graph>(gp, e));
    }

    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("finish_vertex")(PythonVertex<Graph>(gp, u));
    }

private:
    GraphInterface& _gi;
    boost::python::object _vis;
};

struct do_bfs
{
    template <class Graph, class Visitor>
    void operator()(Graph& g, size_t s, Visitor vis) const
    {
        breadth_first_search(g, vertex(s, g), visitor(vis));
    }
};

void bfs_search(GraphInterface& g, size_t s, python::object vis)
{
    run_action<graph_tool::detail::all_graph_views,mpl::true_>()
        (g, std::bind(do_bfs(), placeholders::_1, s,
                      BFSVisitorWrapper(g, vis)))();
}

typedef boost::coroutines::asymmetric_coroutine<boost::python::object> coro_t;

class BFSGeneratorVisitor : public bfs_visitor<>
{
public:
    BFSGeneratorVisitor(GraphInterface& gi,
                        coro_t::push_type& yield)
        : _gi(gi), _yield(yield) {}

    template <class Edge, class Graph>
    void tree_edge(const Edge& e, Graph& g)
    {
        std::shared_ptr<Graph> gp = retrieve_graph_view<Graph>(_gi, g);
        _yield(boost::python::object(PythonEdge<Graph>(gp, e)));
    }

private:
    GraphInterface& _gi;
    coro_t::push_type& _yield;
};

class BFSGenerator
{
public:
    template <class Dispatch>
    BFSGenerator(Dispatch& dispatch)
        : _coro(std::make_shared<coro_t::pull_type>(dispatch)),
          _iter(begin(*_coro)), _end(end(*_coro)) {}
    boost::python::object next()
    {
        if (_iter == _end)
            boost::python::objects::stop_iteration_error();
        boost::python::object oe = *_iter;
        ++_iter;
        return oe;
    }
private:
    std::shared_ptr<coro_t::pull_type> _coro;
    coro_t::pull_type::iterator _iter;
    coro_t::pull_type::iterator _end;
};

boost::python::object bfs_search_generator(GraphInterface& g, size_t s)
{
#ifdef HAVE_BOOST_COROUTINE
    auto dispatch = [&](auto& yield)
        {
            BFSGeneratorVisitor vis(g, yield);
            run_action<graph_tool::detail::all_graph_views,mpl::true_>()
                (g, std::bind(do_bfs(), placeholders::_1, s, vis))();
        };
    return boost::python::object(BFSGenerator(dispatch));
#else
    throw GraphException("This functionality is not available because boost::coroutine was not found at compile-time")
#endif
}

void export_bfs()
{
    using namespace boost::python;
    def("bfs_search", &bfs_search);
    def("bfs_search_generator", &bfs_search_generator);
    class_<BFSGenerator>("BFSGenerator", no_init)
        .def("__iter__", objects::identity_function())
        .def("next", &BFSGenerator::next)
        .def("__next__", &BFSGenerator::next);
}
