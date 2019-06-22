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

#include <boost/python.hpp>

#include "graph_tool.hh"
#include "graph_vertex_similarity.hh"
#include "numpy_bind.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

typedef UnityPropertyMap<uint8_t, GraphInterface::edge_t> ecmap_t;
typedef boost::mpl::push_back<edge_scalar_properties, ecmap_t>::type
        weight_props_t;


void get_dice_similarity(GraphInterface& gi, boost::any as, boost::any weight)
{
    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto& s, auto& w)
         {
             all_pairs_similarity(g, s,
                                  [&](auto u, auto v, auto& mask, auto& w)
                                  {
                                      return dice(u, v, mask, w, g);
                                  }, w);
         },
         all_graph_views(), vertex_floating_vector_properties(),
         weight_props_t())
        (gi.get_graph_view(), as, weight);
}

void get_dice_similarity_pairs(GraphInterface& gi, python::object opairs,
                               python::object osim, boost::any weight)
{
    multi_array_ref<int64_t,2> pairs = get_array<int64_t,2>(opairs);
    multi_array_ref<double,1> sim = get_array<double,1>(osim);

    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto w)
         {
             some_pairs_similarity(g, pairs, sim,
                                   [&](auto u, auto v, auto& mask, auto& w)
                                   {
                                       return dice(u, v, mask, w, g);
                                   }, w);
         },
         all_graph_views(), weight_props_t())
        (gi.get_graph_view(), weight);
}

void get_jaccard_similarity(GraphInterface& gi, boost::any as, boost::any weight)
{
    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto& s, auto w)
         {
             all_pairs_similarity(g, s,
                                  [&](auto u, auto v, auto& mask, auto w)
                                  {
                                      return jaccard(u, v, mask, w, g);
                                  }, w);
         },
         all_graph_views(), vertex_floating_vector_properties(),
         weight_props_t())
        (gi.get_graph_view(), as, weight);
}

void get_jaccard_similarity_pairs(GraphInterface& gi, python::object opairs,
                                  python::object osim, boost::any weight)
{
    multi_array_ref<int64_t,2> pairs = get_array<int64_t,2>(opairs);
    multi_array_ref<double,1> sim = get_array<double,1>(osim);

    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto w)
         {
             some_pairs_similarity(g, pairs, sim,
                                   [&](auto u, auto v, auto& mask, auto w)
                                   {
                                       return jaccard(u, v, mask, w, g);
                                   }, w);
         },
         all_graph_views(), weight_props_t())
        (gi.get_graph_view(), weight);
}

void get_inv_log_weight_similarity(GraphInterface& gi, boost::any as,
                                   boost::any weight)
{
    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto& s, auto w)
         {
             all_pairs_similarity(g, s,
                                  [&](auto u, auto v, auto& mask, auto w)
                                  {
                                      return inv_log_weighted(u, v, mask, w, g);
                                  }, w);
         },
         all_graph_views(), vertex_floating_vector_properties(),
         weight_props_t())
        (gi.get_graph_view(), as, weight);
}

void get_inv_log_weight_similarity_pairs(GraphInterface& gi,
                                         python::object opairs,
                                         python::object osim,
                                         boost::any weight)
{
    multi_array_ref<int64_t,2> pairs = get_array<int64_t,2>(opairs);
    multi_array_ref<double,1> sim = get_array<double,1>(osim);

    if (weight.empty())
        weight = ecmap_t();

    gt_dispatch<>()
        ([&](auto& g, auto w)
         {
             some_pairs_similarity(g, pairs, sim,
                                   [&](auto u, auto v, auto& mask, auto w)
                                   {
                                       return inv_log_weighted(u, v, mask, w, g);
                                   }, w);
         },
         all_graph_views(), weight_props_t())
        (gi.get_graph_view(), weight);
}


void export_vertex_similarity()
{
    python::def("dice_similarity", &get_dice_similarity);
    python::def("dice_similarity_pairs", &get_dice_similarity_pairs);
    python::def("jaccard_similarity", &get_jaccard_similarity);
    python::def("jaccard_similarity_pairs", &get_jaccard_similarity_pairs);
    python::def("inv_log_weight_similarity", &get_inv_log_weight_similarity);
    python::def("inv_log_weight_similarity_pairs",
                &get_inv_log_weight_similarity_pairs);
};
