#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2019 Tiago de Paula Peixoto <tiago@skewed.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from .. import _prop

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_inference as libinference")

from numpy import sqrt

def latent_multigraph(g, epsilon=1e-8, max_niter=0, verbose=False):
    r"""Infer latent Poisson multigraph model given an "erased" simple graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.
    epsilon : ``float`` (optional, default: ``1e-8``)
        Convergence criterion.
    max_niter : ``int`` (optional, default: ``0``)
        Maximum number of iterations allowed (if ``0``, no maximum is assumed).
    verbose : ``boolean`` (optional, default: ``False``)
        If ``True``, display verbose information.

    Returns
    -------
    u : :class:`~graph_tool.Graph`
        Latent graph.
    w : :class:`~graph_tool.EdgePropertyMap`
        Edge property map with inferred multiplicity parameter.

    Examples
    --------
    >>> g = gt.collection.data["as-22july06"]
    >>> gt.scalar_assortativity(g, "out")
    (-0.198384..., 0.001338...)
    >>> u, w = gt.latent_multigraph(g)
    >>> scalar_assortativity(u, "out", eweight=w)
    (-0.048426..., 0.034526...)
    """

    g = g.copy()
    theta_out = g.degree_property_map("out").copy("double")
    theta_out.fa /= sqrt(theta_out.fa.sum())
    if g.is_directed():
        theta_in = g.degree_property_map("in").copy("double")
        theta_in.fa /= sqrt(theta_in.fa.sum())
    else:
        theta_in = theta_out

    w = g.new_ep("double", 1)

    libinference.latent_multigraph(g._Graph__graph,
                                   _prop("e", g, w),
                                   _prop("v", g, theta_out),
                                   _prop("v", g, theta_in),
                                   epsilon, max_niter, verbose)
    return g, w
