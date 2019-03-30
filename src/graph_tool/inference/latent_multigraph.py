#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2018 Tiago de Paula Peixoto <tiago@skewed.de>
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

def latent_multigraph(g, epsilon=1e-8, max_niter=10000):
    r"""
    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be used.

    Returns
    -------

    Notes
    -----

    Examples
    --------
    >>> g = gt.collection.data["football"]
    >>> gt.modularity(g, g.vp.value_tsevans)
    0.5744393497...

    References
    ----------
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
                                   epsilon,
                                   max_niter)
    return g, w
