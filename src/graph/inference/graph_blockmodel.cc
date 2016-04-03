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
#include "random.hh"

#include <boost/python.hpp>

#include "graph_blockmodel_util.hh"
#include "graph_blockmodel.hh"
#include "graph_state.hh"

using namespace boost;
using namespace graph_tool;

GEN_DISPATCH(block_state, BlockState, BLOCK_STATE_params)

python::object make_block_state(boost::python::object ostate,
                                rng_t& rng)
{
    python::object state;
    block_state::make_dispatch(ostate,
                               [&](auto& s){state = python::object(s);},
                               rng);
    return state;
}

degs_map_t get_block_degs(GraphInterface& gi, boost::any ab)
{
    degs_map_t degs;
    vmap_t b = boost::any_cast<vmap_t>(ab);
    run_action<>()(gi,
                   [&](auto& g)
                   {
                       std::vector<gt_hash_map<std::tuple<size_t, size_t>,
                                               size_t>> hist;
                       for (auto v : vertices_range(g))
                       {
                           size_t r = b[v];
                           if (r >= hist.size())
                               hist.resize(r + 1);
                           size_t kin = in_degreeS()(v, g);
                           size_t kout = out_degreeS()(v, g);
                           hist[r][std::make_tuple(kin, kout)]++;
                       }

                       for (size_t r = 0; r < hist.size(); ++r)
                       {
                           auto& deg = degs[r];
                           for (auto& kn : hist[r])
                               deg.emplace_back(get<0>(kn.first),
                                                get<1>(kn.first),
                                                kn.second);
                       }
                   })();
    return degs;
}

degs_map_t get_weighted_block_degs(GraphInterface& gi, degs_map_t& degs,
                                   boost::any ab)
{
    degs_map_t ndegs;
    vmap_t b = boost::any_cast<vmap_t>(ab);
    run_action<>()(gi,
                   [&](auto& g)
                   {
                       std::vector<gt_hash_map<std::tuple<size_t, size_t>,
                                               size_t>> hist;
                       for (auto v : vertices_range(g))
                       {
                           size_t r = b[v];
                           if (r >= hist.size())
                               hist.resize(r + 1);
                           auto& h = hist[r];
                           auto& ks = degs[v];
                           for (auto& k : ks)
                               h[std::make_tuple(get<0>(k), get<1>(k))] += get<2>(k);
                       }

                       for (size_t r = 0; r < hist.size(); ++r)
                       {
                           auto& deg = ndegs[r];
                           for (auto& kn : hist[r])
                               deg.emplace_back(get<0>(kn.first),
                                                get<1>(kn.first),
                                                kn.second);
                       }
                   })();
    return ndegs;
}


template <class Prop>
boost::any get_any(Prop& p)
{
    return any(p);
}

void print_degs(degs_map_t& degs, size_t B)
{
    for (size_t r = 0; r < B; ++r)
    {
        cout << r << ":: ";
        auto& ks = degs[r];
        for (auto& k : ks)
        {
            cout << "(" << get<0>(k) << ", " << get<1>(k) << "): "
                 << get<2>(k) << "  ";
        }
        cout << endl;
    }
}

degs_map_t copy_degs(degs_map_t& degs)
{
    return degs.copy();
}

simple_degs_t copy_simple_degs(simple_degs_t& degs)
{
    return degs;
}

boost::python::tuple bethe_entropy(GraphInterface& gi, size_t B, boost::any op,
                                   boost::any opv)
{
    typedef vprop_map_t<vector<double>>::type vmap_t;
    typedef eprop_map_t<vector<int32_t>>::type emap_t;
    emap_t p = any_cast<emap_t>(op);
    vmap_t pv = any_cast<vmap_t>(opv);

    double H=0, sH=0, Hmf=0, sHmf=0;
    run_action<graph_tool::all_graph_views, boost::mpl::true_>()
        (gi,
         [&](auto& g)
         {
             for (auto v : vertices_range(g))
             {
                 pv[v].resize(B);
                 for (size_t i = 0; i < B; ++i)
                     pv[v][i] = 0;
             }

             H = Hmf = sH = sHmf =  0;

             for (auto e : edges_range(g))
             {
                 auto u = min(source(e, g), target(e, g));
                 auto v = max(source(e, g), target(e, g));

                 double sum = 0;
                 for (size_t r = 0; r < B; ++r)
                     for (size_t s = 0; s < B; ++s)
                     {
                         size_t i = r + B * s;
                         pv[u][r] += p[e][i];
                         pv[v][s] += p[e][i];
                         sum += p[e][i];
                     }

                 for (size_t i = 0; i < B * B; ++i)
                 {
                     if (p[e][i] == 0)
                         continue;
                     double pi = double(p[e][i]) / sum;
                     H -= pi * log(pi);
                     sH += pow((log(pi) + 1) * sqrt(pi / sum), 2);
                 }
             }

             for (auto v : vertices_range(g))
             {
                 double sum = 0;
                 for (size_t i = 0; i < B; ++i)
                     sum += pv[v][i];
                 for (size_t i = 0; i < B; ++i)
                 {
                     if (pv[v][i] == 0)
                         continue;
                     pv[v][i] /= sum;
                     double pi = pv[v][i];
                     double kt = (1 - double(in_degreeS()(v, g)) - double(out_degree(v, g)));
                     if (kt != 0)
                     {
                         H -= kt * (pi * log(pi));
                         sH += pow(kt * (log(pi) + 1) * sqrt(pi / sum), 2);
                     }

                     Hmf -= pi * log(pi);
                     sHmf += pow((log(pi) + 1) * sqrt(pi / sum), 2);
                 }
             }
         })();

    return boost::python::make_tuple(H, sH, Hmf, sHmf);
}


void export_blockmodel_state()
{
    using namespace boost::python;

    block_state::dispatch
        ([&](auto* s)
         {
             typedef typename std::remove_reference<decltype(*s)>::type state_t;

             double (state_t::*virtual_move)(size_t, size_t, bool, bool, bool,
                                             bool, bool) =
                 &state_t::virtual_move;
             size_t (state_t::*sample_block)(size_t, double, vector<size_t>&,
                                             rng_t&)
                 = &state_t::sample_block;
             double (state_t::*get_move_prob)(size_t, size_t, size_t, double,
                                              bool)
                 = &state_t::get_move_prob;
             void (state_t::*merge_vertices)(size_t, size_t)
                 = &state_t::merge_vertices;

             void (state_t::*set_partition)(boost::any&)
                 = &state_t::set_partition;

             class_<state_t> c(name_demangle(typeid(state_t).name()).c_str(),
                               no_init);
             c.def("remove_vertex", &state_t::remove_vertex)
                 .def("add_vertex", &state_t::add_vertex)
                 .def("move_vertex", &state_t::move_vertex)
                 .def("set_partition", set_partition)
                 .def("virtual_move", virtual_move)
                 .def("merge_vertices", merge_vertices)
                 .def("sample_block", sample_block)
                 .def("entropy", &state_t::entropy)
                 .def("get_partition_dl", &state_t::get_partition_dl)
                 .def("get_deg_dl", &state_t::get_deg_dl)
                 .def("get_move_prob", get_move_prob)
                 .def("enable_partition_stats",
                      &state_t::enable_partition_stats)
                 .def("disable_partition_stats",
                      &state_t::disable_partition_stats)
                 .def("is_partition_stats_enabled",
                      &state_t::is_partition_stats_enabled);
         });

    class_<vcmap_t>("unity_vprop_t").def("_get_any", &get_any<vcmap_t>);
    class_<ecmap_t>("unity_eprop_t").def("_get_any", &get_any<ecmap_t>);

    def("make_block_state", &make_block_state);

    class_<std::true_type>("true_type");
    class_<std::false_type>("false_type");

    def("get_block_degs", &get_block_degs);
    def("get_weighted_block_degs", &get_weighted_block_degs);
    class_<degs_map_t>("degs_map_t")
        .def("print", &print_degs)
        .def("copy", &copy_degs);
    class_<simple_degs_t>("simple_degs_t")
        .def("copy", &copy_simple_degs);

    def("bethe_entropy", &bethe_entropy);
}
