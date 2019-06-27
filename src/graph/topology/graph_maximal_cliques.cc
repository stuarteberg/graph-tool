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

#include "graph_tool.hh"
#include "numpy_bind.hh"
#include "coroutine.hh"
#include "graph_python_interface.hh"

#include "graph_maximal_cliques.hh"

using namespace std;
using namespace graph_tool;

boost::python::object get_max_cliques(GraphInterface& gi)
{
#ifdef HAVE_BOOST_COROUTINE
    auto dispatch = [&](auto& yield)
        {
            run_action<>()
                (gi,
                 [&](auto& g)
                 {
                     max_cliques(g,
                                 [&](auto& s)
                                 {
                                     std::vector<size_t> v(s.begin(), s.end());
                                     auto c = wrap_vector_owned(v);
                                     yield(c);
                                 });
                 })();
        };
    return boost::python::object(CoroGenerator(dispatch));
#else
    throw GraphException("This functionality is not available because boost::coroutine was not found at compile-time");
#endif // HAVE_BOOST_COROUTINE
}

boost::python::list get_max_cliques_list(GraphInterface& gi)
{
    boost::python::list cliques;
    run_action<>()
        (gi,
         [&](auto& g)
         {
             max_cliques(g,
                         [&](auto& s)
                         {
                             std::vector<size_t> v(s.begin(), s.end());
                             cliques.append(wrap_vector_owned(v));
                         });
         })();
    return cliques;
}


void export_max_cliques()
{
    boost::python::def("max_cliques", &get_max_cliques);
    boost::python::def("max_cliques_list", &get_max_cliques_list);
};
