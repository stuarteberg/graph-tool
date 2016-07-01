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

#include "graph_tool.hh"

#include <boost/python.hpp>
#include "numpy_bind.hh"
#include "hash_map_wrap.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

void collect_vertex_marginals(GraphInterface& gi, boost::any ob,
                              boost::any op, double update)
{
    typedef vprop_map_t<int32_t>::type vmap_t;
    auto b = any_cast<vmap_t>(ob).get_unchecked();

    run_action<>()
        (gi, [&](auto& g, auto p)
         {
             typename property_traits<decltype(p)>::value_type::value_type
                 up = update;
             parallel_vertex_loop
                 (g,
                  [&](auto v)
                  {
                      auto r = b[v];
                      auto& pv = p[v];
                      if (pv.size() <= size_t(r))
                          pv.resize(r + 1);
                      pv[r] += up;
                  });
         },
         vertex_scalar_vector_properties())(op);
}

class BlockPairHist
    : public gt_hash_map<std::pair<int32_t,int32_t>,double>
{
public:

    boost::python::dict get_state()
    {
        boost::python::dict state;
        for (auto& kv : *this)
        {
            auto k = python::make_tuple(kv.first.first,
                                        kv.first.second);
            state[k] = kv.second;
        }
        return state;
    }

    void set_state(boost::python::dict state)
    {
        auto keys = state.keys();
        for (int i = 0; i < python::len(keys); ++i)
        {
            auto k = keys[i];
            int32_t r = python::extract<int32_t>(k[0]);
            int32_t s = python::extract<int32_t>(k[1]);
            double v = python::extract<double>(state[k]);
            (*this)[make_pair(r, s)] = v;
        }
    }

    double get_item(boost::python::object k)
    {
        int32_t r = python::extract<int32_t>(k[0]);
        int32_t s = python::extract<int32_t>(k[1]);
        auto iter = this->find(make_pair(r, s));
        if (iter == this->end())
            return 0;
        return iter->second;
    }

    void set_item(boost::python::object k, double v)
    {
        int32_t r = python::extract<int32_t>(k[0]);
        int32_t s = python::extract<int32_t>(k[1]);
        (*this)[make_pair(r, s)] = v;
    }

};

void collect_edge_marginals(GraphInterface& gi, boost::any ob,
                            boost::any op, double update)
{
    typedef vprop_map_t<int32_t>::type vmap_t;
    auto b = any_cast<vmap_t>(ob).get_unchecked();

    typedef eprop_map_t<python::object>::type
        emap_t;
    auto pe = any_cast<emap_t>(op).get_unchecked(gi.get_edge_index_range());

    run_action<>()
        (gi,
         [&](auto& g)
         {
            parallel_edge_loop
                 (g,
                  [&](const auto& e)
                  {
                      auto u = min(source(e, g), target(e, g));
                      auto v = max(source(e, g), target(e, g));

                      auto r = b[u];
                      auto s = b[v];

                      BlockPairHist& h =
                          boost::python::extract<BlockPairHist&>(pe[e]);

                      h[make_pair(r, s)] += update;
                  });
         })();
}

double mf_entropy(GraphInterface& gi, boost::any opv)
{
    double H=0;
    run_action<>()
        (gi,
         [&](auto& g, auto pv)
         {
             for (auto v : vertices_range(g))
             {
                 double sum = 0;
                 for (auto p : pv[v])
                     sum += p;
                 for (double p : pv[v])
                 {
                     if (p == 0)
                         continue;
                     p /= sum;
                     H -= p * log(p);
                 }
             }
         },
         vertex_scalar_vector_properties())(opv);

    return H;
}

boost::python::tuple bethe_entropy(GraphInterface& gi, boost::any op,
                                   boost::any opv)
{
    typedef vprop_map_t<vector<double>>::type vmap_t;
    vmap_t pv = any_cast<vmap_t>(opv);

    typedef eprop_map_t<boost::python::object>::type emap_t;
    auto pe = any_cast<emap_t>(op).get_unchecked();

    double H=0, Hmf=0;
    run_action<>()
        (gi,
         [&](auto& g)
         {
             for (auto e : edges_range(g))
             {
                 auto u = min(source(e, g), target(e, g));
                 auto v = max(source(e, g), target(e, g));

                 double sum = 0;

                 BlockPairHist& h =
                     boost::python::extract<BlockPairHist&>(pe[e]);

                 for (auto& prs : h)
                 {
                     sum += prs.second;
                     size_t r = prs.first.first;
                     size_t s = prs.first.second;
                     if (r >= pv[u].size())
                         pv[u].resize(r + 1);
                     pv[u][r] += prs.second;
                     if (s >= pv[v].size())
                         pv[v].resize(s + 1);
                     pv[v][s] += prs.second;
                 }

                 for (auto& prs : h)
                 {
                     if (prs.second == 0)
                         continue;
                     double pi = prs.second / sum;
                     H -= pi * log(pi);
                 }
             }

             for (auto v : vertices_range(g))
             {
                 double sum = std::accumulate(pv[v].begin(), pv[v].end(), 0);
                 for (double pi : pv[v])
                 {
                     if (pi == 0)
                         continue;
                     pi /= sum;
                     int kt = 1 - total_degreeS()(v, g);
                     if (kt != 0)
                         H -= kt * (pi * log(pi));
                     Hmf -= pi * log(pi);
                 }
             }
         })();

    return boost::python::make_tuple(H, Hmf);
}

void export_marginals()
{
    using namespace boost::python;

    class_<BlockPairHist>("BlockPairHist")
        .def("__setitem__", &BlockPairHist::set_item)
        .def("__getitem__", &BlockPairHist::get_item)
        .def("__setstate__", &BlockPairHist::set_state)
        .def("__getstate__", &BlockPairHist::get_state).enable_pickling();

    def("vertex_marginals", &collect_vertex_marginals);
    def("edge_marginals", &collect_edge_marginals);
    def("mf_entropy", &mf_entropy);
    def("bethe_entropy", &bethe_entropy);
}
