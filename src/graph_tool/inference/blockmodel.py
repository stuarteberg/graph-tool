#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
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

from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange

from .. import _degree, _prop, Graph, GraphView, libcore, _get_rng, PropertyMap, \
    conv_pickle_state, Vector_size_t, Vector_double, group_vector_property
from .. generation import condensation_graph
from .. stats import label_self_loops
from .. spectral import adjacency
import random
from numpy import *
import numpy
import copy
import collections

from . util import *

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_inference as libinference")

__test__ = False

def set_test(test):
    global __test__
    __test__ = test

def _bm_test():
    global __test__
    return __test__

def get_block_graph(g, B, b, vcount=None, ecount=None):
    if isinstance(ecount, libinference.unity_eprop_t):
        ecount = None
    if isinstance(vcount, libinference.unity_vprop_t):
        vcount = None
    cg, br, vcount, ecount = condensation_graph(g, b,
                                                vweight=vcount,
                                                eweight=ecount,
                                                self_loops=True)[:4]
    cg.vp["count"] = vcount
    cg.ep["count"] = ecount
    cg = Graph(cg, vorder=br)

    cg.add_vertex(B - cg.num_vertices())
    return cg

class BlockState(object):
    r"""The stochastic block model state of a given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge multiplicities (for multigraphs or block graphs).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex multiplicities (for block graphs).
    b : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Initial block labels on the vertices. If not supplied, it will be
        randomly sampled.
    B : ``int`` (optional, default: ``None``)
        Number of blocks (or vertex groups). If not supplied it will be obtained
        from the parameter ``b``.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    pclabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Partition constraint labels on the vertices. This has the same
        interpretation as ``clabel``, but will be used to compute the partition
        description length.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this value, a sparse matrix is used for
        the block graph. Otherwise a dense matrix will be used.
    """

    def __init__(self, g, eweight=None, vweight=None, b=None, B=None,
                 clabel=None, pclabel=None, deg_corr=True, max_BE=1000,
                 **kwargs):

        # initialize weights to unity, if necessary
        if eweight is None or isinstance(eweight, libinference.unity_eprop_t):
            eweight = libinference.unity_eprop_t()
        else:
            eweight = g.own_property(eweight.copy(value_type="int32_t"))
        if vweight is None or isinstance(vweight, libinference.unity_vprop_t):
            vweight = libinference.unity_vprop_t()
        else:
            vweight = g.own_property(vweight.copy(value_type="int32_t"))
        self.eweight = eweight
        self.vweight = vweight

        self.is_edge_weighted = not isinstance(eweight, libinference.unity_eprop_t)
        self.is_vertex_weighted = not isinstance(vweight, libinference.unity_vprop_t)
        self.is_weighted = self.is_edge_weighted or self.is_vertex_weighted

        # configure the main graph and block model parameters
        self.g = g

        self.E = int(self.eweight.fa.sum()) if self.is_edge_weighted else g.num_edges()
        self.N = int(self.vweight.fa.sum()) if self.is_vertex_weighted else g.num_vertices()

        self.deg_corr = deg_corr
        self.overlap = False
        self.degs = kwargs.get("degs", libinference.simple_degs_t())
        if self.degs is None:
            self.degs = libinference.simple_degs_t()

        # ensure we have at most as many blocks as nodes
        if B is not None and b is None:
            B = min(B, self.g.num_vertices())

        if b is None:
            # create a random partition into B blocks.
            if B is None:
                B = get_max_B(self.N, self.E, directed=g.is_directed())
            B = min(B, self.g.num_vertices())
            ba = random.randint(0, B, self.g.num_vertices())
            ba[:B] = arange(B)        # avoid empty blocks
            if B < self.g.num_vertices():
                random.shuffle(ba)
            b = g.new_vp("int")
            b.fa = ba
            self.b = b
        else:
            # if a partition is available, we will incorporate it.
            if isinstance(b, numpy.ndarray):
                self.b = g.new_vp("int")
                self.b.fa = b
            else:
                self.b = b = g.own_property(b.copy(value_type="int32_t"))
            if B is None:
                B = int(self.b.fa.max()) + 1

        if self.b.fa.max() >= B:
            raise ValueError("Maximum value of b is larger or equal to B! (%d vs %d)" %
                             (self.b.fa.max(), B))

        # Construct block-graph
        self.bg = get_block_graph(g, B, self.b, self.vweight, self.eweight)
        self.bg.set_fast_edge_removal()

        self.mrs = self.bg.ep["count"]
        self.wr = self.bg.vp["count"]

        del self.bg.ep["count"]
        del self.bg.vp["count"]

        self.mrp = self.bg.degree_property_map("out", weight=self.mrs)

        if g.is_directed():
            self.mrm = self.bg.degree_property_map("in", weight=self.mrs)
        else:
            self.mrm = self.mrp

        self.B = B

        if pclabel is not None:
            if isinstance(pclabel, PropertyMap):
                self.pclabel = self.g.own_property(pclabel).copy("int")
            else:
                self.pclabel = self.g.new_vp("int")
                self.pclabel.fa = pclabel
        else:
            self.pclabel = self.g.new_vp("int")

        if clabel is not None:
            if isinstance(clabel, PropertyMap):
                self.clabel = self.g.own_property(clabel).copy("int")
            else:
                self.clabel = self.g.new_vp("int")
                self.clabel.fa = clabel
        elif self.pclabel.fa.max() > 0:
            self.clabel = self.pclabel
        else:
            self.clabel = self.g.new_vp("int")

        if not self._check_clabel():
            raise ValueError("provided clabel is inconsistent with node partition")

        self.bclabel = self.get_bclabel()

        self.max_BE = max_BE

        if self.B > self.max_BE:
            self.use_hash = libinference.true_type()
        else:
            self.use_hash = libinference.false_type()

        self.ignore_degrees = kwargs.get("ignore_degrees", None)
        if self.ignore_degrees is None:
            self.ignore_degrees = g.new_vp("bool", False)

        self.merge_map = kwargs.get("merge_map", self.g.vertex_index.copy("int"))

        self.block_list = Vector_size_t()
        self.block_list.extend(arange(self.B, dtype="int"))

        self._abg = self.bg._get_any()
        self._state = libinference.make_block_state(self, _get_rng())


    def __repr__(self):
        return "<BlockState object with %d blocks,%s for graph %s, at 0x%x>" % \
            (self.B, " degree corrected," if self.deg_corr else "", str(self.g),
             id(self))

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        g = copy.deepcopy(self.g, memo)
        return self.copy(g=g)

    def copy(self, g=None, eweight=None, vweight=None, b=None, B=None,
             deg_corr=None, clabel=None, overlap=False, pclabel=None, max_BE=None,
             **kwargs):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor."""

        if not overlap:
            state = BlockState(self.g if g is None else g,
                               eweight=self.eweight if eweight is None else eweight,
                               vweight=self.vweight if vweight is None else vweight,
                               b=self.b.copy() if b is None else b,
                               B=(self.B if b is None else None) if B is None else B,
                               clabel=self.clabel if clabel is None else clabel,
                               pclabel=self.pclabel if pclabel is None else pclabel,
                               deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                               max_BE=self.max_BE if max_BE is None else max_BE,
                               ignore_degrees=kwargs.pop("ignore_degrees", self.ignore_degrees),
                               degs=self.degs.copy(),
                               merge_map=kwargs.get("merge_map", self.merge_map.copy()),
                               **dmask(kwargs, ["ignore_degrees", "merge_map"]))
        else:
            state = OverlapBlockState(self.g if g is None else g,
                                      b=self.b.copy() if b is None else b,
                                      B=(self.B if b is None else None) if B is None else B,
                                      clabel=self.clabel if clabel is None else clabel,
                                      pclabel=self.pclabel if pclabel is None else pclabel,
                                      deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                                      max_BE=self.max_BE if max_BE is None else max_BE,
                                      **kwargs)
        return state


    def __getstate__(self):
        state = dict(g=self.g,
                     eweight=self.eweight if self.is_edge_weighted else None,
                     vweight=self.vweight if self.is_vertex_weighted else None,
                     b=self.b,
                     B=self.B,
                     clabel=self.clabel,
                     pclabel=self.pclabel,
                     deg_corr=self.deg_corr,
                     max_BE=self.max_BE,
                     ignore_degrees=self.ignore_degrees,
                     degs=self.degs if not isinstance(self.degs,
                                                      libinference.simple_degs_t) else None,
                     merge_map=self.merge_map)
        return state

    def __setstate__(self, state):
        conv_pickle_state(state)
        self.__init__(**state)

    def get_block_state(self, b=None, vweight=False, deg_corr=False, **kwargs):
        r"""Returns a :class:`~graph_tool.community.BlockState`` corresponding to the
        block graph (i.e. the blocks of the current state become the nodes). The
        parameters have the same meaning as the in the constructor. If ``vweight
        == True`` the nodes of the block state are weighted with the node
        counts."""

        if deg_corr and vweight:
            if isinstance(self.degs, libinference.simple_degs_t):
                degs = libinference.get_block_degs(self.g._Graph__graph,
                                                   _prop("v", self.g, self.b))
            else:
                degs = libinference.get_weighted_block_degs(self.g._Graph__graph,
                                                            self.degs,
                                                            _prop("v", self.g,
                                                                  self.b))
        else:
            degs = libinference.simple_degs_t()

        state = BlockState(self.bg.copy(), eweight=self.mrs,
                           vweight=self.wr if vweight else self.bg.new_vp("int", 1),
                           b=self.bg.vertex_index.copy("int") if b is None else b,
                           deg_corr=deg_corr,
                           degs=degs,
                           max_BE=self.max_BE, **kwargs)
        return state

    def get_bclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap` corresponding to constraint
        labels for the block graph."""

        bclabel = self.bg.new_vertex_property("int")
        reverse_map(self.b, bclabel)
        pmap(bclabel, self.clabel)
        return bclabel

    def get_bpclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap`` corresponding to partition
        constraint labels for the block graph."""

        bclabel = self.bg.new_vertex_property("int")
        reverse_map(self.b, bclabel)
        pmap(bclabel, self.pclabel)
        return bclabel

    def _check_clabel(self):
        b = self.b.fa + self.clabel.fa * self.B
        b2 = self.b.fa.copy()
        continuous_map(b)
        continuous_map(b2)
        return (b == b2).all()

    def get_blocks(self):
        r"""Returns the property map which contains the block labels for each vertex."""
        return self.b

    def set_blocks(self, b):
        r"""Sets the internal partition of the state."""
        if b.value_type() != "int32_t":
            b = b.copy("int32_t")
        self._state.set_partition(_prop("v", self.g, b))

    def get_bg(self):
        r"""Returns the block graph."""
        return self.bg

    def get_ers(self):
        r"""Returns the edge property map of the block graph which contains the
        :math:`e_{rs}` matrix entries.  For undirected graphs, the diagonal
        values (self-loops) contain :math:`e_{rr}/2`."""
        return self.mrs

    def get_er(self):
        r"""Returns the vertex property map of the block graph which contains the number
        :math:`e_r` of half-edges incident on block :math:`r`. If the graph is
        directed, a pair of property maps is returned, with the number of
        out-edges :math:`e^+_r` and in-edges :math:`e^-_r`, respectively."""
        if self.bg.is_directed():
            return self.mrp, self.mrm
        else:
            return self.mrp

    def get_nr(self):
        r"""Returns the vertex property map of the block graph which contains the block
        sizes :math:`n_r`."""
        return self.wr

    def entropy(self, dl=True, partition_dl=True, degree_dl=True,
                edges_dl=True, dense=False, multigraph=True, dl_ent=False,
                deg_entropy=True, **kwargs):
        r"""Calculate the entropy associated with the current block partition.

        Parameters
        ----------
        dl : ``bool`` (optional, default: ``False``)
            If ``True``, the full description length will be returned.
        partition_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the partition description length
            will be considered.
        degree_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the degree sequence description
            length will be considered.
        edges_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the edge matrix description length
            will be considered.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``False``)
            If ``True``, the multigraph entropy will be used.
        dl_ent : ``bool`` (optional, default: ``False``)
            If ``True``, the description length of the degree sequence will be
            approximated by its entropy.
        deg_entropy : ``bool`` (optional, default: ``True``)
            If ``True``, the degree entropy term that is independent of the
            network partition will be returned.

        Notes
        -----
        For the traditional blockmodel (``deg_corr == False``), the entropy is
        given by

        .. math::

          \mathcal{S}_t &\cong E - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right), \\
          \mathcal{S}^d_t &\cong E - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{n_rn_s}\right),

        for undirected and directed graphs, respectively, where :math:`e_{rs}`
        is the number of edges from block :math:`r` to :math:`s` (or the number
        of half-edges for the undirected case when :math:`r=s`), and :math:`n_r`
        is the number of vertices in block :math:`r` .

        For the degree-corrected variant with "hard" degree constraints the
        equivalent expressions are

        .. math::

            \mathcal{S}_c &\cong -E -\sum_kN_k\ln k! - \frac{1}{2} \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e_re_s}\right), \\
            \mathcal{S}^d_c &\cong -E -\sum_{k^+}N_{k^+}\ln k^+!  -\sum_{k^-}N_{k^-}\ln k^-! - \sum_{rs}e_{rs}\ln\left(\frac{e_{rs}}{e^+_re^-_s}\right),

        where :math:`e_r = \sum_se_{rs}` is the number of half-edges incident on
        block :math:`r`, and :math:`e^+_r = \sum_se_{rs}` and :math:`e^-_r =
        \sum_se_{sr}` are the numbers of out- and in-edges adjacent to block
        :math:`r`, respectively.

        If ``dense == False`` and ``multigraph == True``, the entropy used will
        be of the "Poisson" model, with the additional term:

        .. math::

            {\mathcal{S}_{cm}^{(d)}} = \mathcal{S}_c^{(d)} + \sum_{i>j} \ln A_{ij}! + \sum_i \ln A_{ii}!!


        If ``dl == True``, the description length :math:`\mathcal{L}_t` of the
        model will be returned as well, as described in
        :func:`model_entropy`. Note that for the degree-corrected version the
        description length is

        .. math::

            \mathcal{L}_c = \mathcal{L}_t + \sum_r\min\left(\mathcal{L}^{(1)}_r, \mathcal{L}^{(2)}_r\right),

        with

        .. math::

              \mathcal{L}^{(1)}_r &= \ln{\left(\!\!{n_r \choose e_r}\!\!\right)}, \\
              \mathcal{L}^{(2)}_r &= \ln\Xi_r + \ln n_r! - \sum_k \ln n^r_k!,

        and :math:`\ln\Xi_r \simeq 2\sqrt{\zeta(2)e_r}`, where :math:`\zeta(x)`
        is the `Riemann zeta function
        <https://en.wikipedia.org/wiki/Riemann_zeta_function>`_, and
        :math:`n^r_k` is the number of nodes in block :math:`r` with degree
        :math:`k`. For directed graphs we have instead :math:`k \to (k^-, k^+)`,
        and :math:`\ln\Xi_r \to \ln\Xi^+_r + \ln\Xi^-_r` with :math:`\ln\Xi_r^+
        \simeq 2\sqrt{\zeta(2)e^+_r}` and :math:`\ln\Xi_r^- \simeq
        2\sqrt{\zeta(2)e^-_r}`.

        If ``dl_ent=True`` is passed, this will be approximated instead by

        .. math::

            \mathcal{L}_c \simeq \mathcal{L}_t - \sum_rn_r\sum_kp^r_k\ln p^r_k,

        where :math:`p^r_k = n^r_k / n_r`.

        If the "dense" entropies are requested (``dense=True``), they will be
        computed as

        .. math::

            \mathcal{S}_t  &= \sum_{r>s} \ln{\textstyle {n_rn_s \choose e_{rs}}} + \sum_r \ln{\textstyle {{n_r\choose 2}\choose e_{rr}/2}}\\
            \mathcal{S}^d_t  &= \sum_{rs} \ln{\textstyle {n_rn_s \choose e_{rs}}},

        for simple graphs, and

        .. math::

            \mathcal{S}_m  &= \sum_{r>s} \ln{\textstyle \left(\!\!{n_rn_s \choose e_{rs}}\!\!\right)} + \sum_r \ln{\textstyle \left(\!\!{\left(\!{n_r\choose 2}\!\right)\choose e_{rr}/2}\!\!\right)}\\
            \mathcal{S}^d_m  &= \sum_{rs} \ln{\textstyle \left(\!\!{n_rn_s \choose e_{rs}}\!\!\right)},

        for multigraphs (i.e. ``multigraph == True``). A dense entropy for the
        degree-corrected model is not available, and if requested will raise a
        :exc:`NotImplementedError`.
        """

        if _bm_test() and kwargs.get("test", True):
            args = dict(**locals())
            args.update(**kwargs)
            del args["self"]
            del args["kwargs"]

        xi_fast = kwargs.get("xi_fast", False)
        dl_deg_alt = kwargs.get("dl_deg_alt", True)

        E = self.E
        N = self.N

        S = self._state.entropy(dense, multigraph, deg_entropy)

        if not dense:
            if self.deg_corr:
                S -= E
            else:
                S += E
        if dl:
            if partition_dl:
                S += self._state.get_partition_dl()

            if edges_dl:
                actual_B = (self.wr.a > 0).sum()
                S += model_entropy(actual_B, N, E,
                                   directed=self.g.is_directed(), nr=False)

            if self.deg_corr and degree_dl:
                S_seq = self._state.get_deg_dl(dl_ent, dl_deg_alt, xi_fast)

                S += S_seq

        callback = kwargs.get("callback", None)
        if callback is not None:
            S += callback(self)

        if _bm_test() and kwargs.get("test", True):
            assert not isnan(S) and not isinf(S), \
                "invalid entropy %g (%s) " % (S, str(args))

            args["test"] = False
            Salt = self.copy().entropy(**args)
            assert abs(S - Salt) < 1e-6, \
                "entropy discrepancy after copying (%g %g)" % (S, Salt)

        return S

    def get_matrix(self):
        r"""Returns the block matrix (as a sparse :class:`~scipy.sparse.csr_matrix`),
        which contains the number of edges between each block pair.

        .. warning::

           This corresponds to the adjacency matrix of the block graph, which by
           convention includes twice the amount of edges in the diagonal entries
           if the graph is undirected.

        Examples
        --------

        .. testsetup:: get_matrix

           gt.seed_rng(42)
           np.random.seed(42)
           from pylab import *

        .. doctest:: get_matrix

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=5, deg_corr=True)
           >>> state.mcmc_sweep(niter=1000)
           (...)
           >>> m = state.get_matrix()
           >>> figure()
           <...>
           >>> matshow(m.todense())
           <...>
           >>> savefig("bloc_mat.pdf")

        .. testcleanup:: get_matrix

           savefig("bloc_mat.png")

        .. figure:: bloc_mat.*
           :align: center

           A  5x5 block matrix.

       """

        return adjacency(self.bg, weight=self.mrs)

    def virtual_vertex_move(self, v, s, dense=False, multigraph=True, dl=True,
                            edges_dl=True, partition_dl=True):
        r"""Computes the entropy difference if vertex ``v`` is moved to block ``s``. The
        remaining parameters are the same as in
        :meth:`graph_tool.BlockState.entropy`."""
        return self._state.virtual_move(int(v), s, dense, multigraph, dl,
                                        edges_dl, partition_dl)

    def move_vertex(self, v, s):
        r"""Move vertex ``v`` to block ``s``."""
        self._state.move_vertex(int(v), s)

    def remove_vertex(self, v):
        r"""Remove vertex ``v`` from its current group.

        This optionally accepts a list of vertices to remove.

        .. warning::

           This will leave the state in an inconsistent state before the vertex
           is returned to some other group, or if the same vertex is removed
           twice.
        """
        if isinstance(v, collections.Iterable):
            self._state.remove_vertices(list(v))
        else:
            self._state.remove_vertex(int(v))

    def add_vertex(self, v, r):
        r"""Add vertex ``v`` to block ``r``.

        This optionally accepts a list of vertices and blocks to add.

        .. warning::

           This can leave the state in an inconsistent state if a vertex is
           added twice to the same group.
        """
        if isinstance(v, collections.Iterable):
            self._state.add_vertices(list(v), list(r))
        else:
            self._state.add_vertex(int(v), r)

    def merge_vertices(self, u, v):
        r"""Merge vertex ``u`` into ``v``.

        .. warning::

           This modifies the underlying graph.
        """
        self.move_vertex(u, self.b[v])
        self._state.merge_vertices(int(u), int(v))

    def sample_vertex_move(self, v, c=1., block_list=None):
        r"""Sample block membership proposal of vertex ``v`` according to real-valued
        sampling parameter ``c``: For :math:`c\to 0` the blocks are sampled
        according to the local neighbourhood and their connections; for
        :math:`c\to\infty` the blocks are sampled randomly. If ``block_list`` is
        passed, it should be a list of blocks among which the sampling will be
        constrained if :math:`c = \infty`.
        """
        blist = self.block_list
        if block_list is not None:
            blist = Vector_size_t()
            blist.extend(block_list)
        return self._state.sample_block(int(v), c, blist, _get_rng())

    def get_move_prob(self, v, s, c=1., reverse=False):
        r"""Compute the probability of a move proposal for vertex ``v`` to block ``s``
        according to sampling parameter ``c``, as obtained with
        :meth:`graph_tool.inference.BlockState.sample_vertex_move`. If
        ``reverse == True``, the reverse probability of moving the node back
        from block ``s`` to its current one is obtained.
        """
        return self._state.get_move_prob(int(v), self.b[v], s, c, reverse)

    def get_edges_prob(self, edge_list, missing=True, entropy_args={}):
        """Compute the log-probability of the missing (or spurious if ``missing=False``)
        edges given by ``edge_list`` (a list of ``(source, target)`` tuples, or
        :meth:`~graph_tool.Edge` instances). The values in ``entropy_args`` are
        passed to :meth:`graph_tool.BlockState.entropy()` to calculate the
        log-probability.
        """
        pos = {}
        for u, v in edge_list:
            pos[u] = self.b[u]
            pos[v] = self.b[v]

        Si = self.entropy(**entropy_args)

        self.remove_vertex(pos.keys())

        try:
            if missing:
                new_es = []
                for u, v in edge_list:
                    e = self.g.add_edge(u, v)
                    if self.is_edge_weighted:
                        self.eweight[e] = 1
                    new_es.append(e)
                    self.E += 1
            else:
                old_es = []
                for e in edge_list:
                    u, v = e
                    if isinstance(e, tuple):
                        e = self.g.edge(u, v)
                        if e is None:
                            raise ValueError("edge not found: (%d, %d)" % (int(u),
                                                                           int(v)))

                    if self.is_edge_weighted:
                        self.eweight[e] -= 1
                        if self.eweight[e] == 0:
                            self.g.remove_edge(e)
                    else:
                        self.g.remove_edge(e)
                    old_es.append((u, v))
                    self.E -= 1

            self.add_vertex(pos.keys(), pos.values())

            Sf = self.entropy(**entropy_args)

            self.remove_vertex(pos.keys())

        finally:
            if missing:
                for e in new_es:
                    self.g.remove_edge(e)
                    self.E -= 1
            else:
                for u, v in old_es:
                    if self.is_edge_weighted:
                        e = self.g.edge(u, v)
                        if e is None:
                            e = self.g.add_edge(u, v)
                            self.eweight[e] = 0
                        self.eweight[e] += 1
                    else:
                        self.g.add_edge(u, v)
                    self.E += 1
            self.add_vertex(pos.keys(), pos.values())

        if missing:
            return Si - Sf
        else:
            return Sf - Si

    def _mcmc_sweep_dispatch(self, mcmc_state):
        if (mcmc_state.multigraph and not mcmc_state.dense and
            not isinstance(self.degs, libinference.simple_degs_t)):
            raise ValueError("can only perform mcmc sweep for degree-corrected" +
                             " grouped states with multigraph == True if dense == True")
        return libinference.mcmc_sweep(mcmc_state, self._state,
                                       _get_rng())

    def mcmc_sweep(self, beta=1., c=1., niter=1, entropy_args={},
                   allow_empty=True, sequential=True, parallel=False,
                   vertices=None, block_list=None, verbose=False):
        r"""Perform ``niter`` sweeps of a Metropolis-Hastings rejection sampling MCMC
        to sample network partitions.

        Parameters
        ----------
        beta : ``float`` (optional, default: ``1.``)
            Inverse temperature.
        c : ``float`` (optional, default: ``1.``)
            Sampling parameter ``c`` for move proposals: For :math:`c\to 0` the
            blocks are sampled according to the local neighbourhood of a given
            node and their block connections; for :math:`c\to\infty` the blocks
            are sampled randomly. Note that only for :math:`c > 0` the MCMC is
            guaranteed to be ergodic.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_empty : ``bool`` (optional, default: ``True``)
            Allow movements into empty groups.
        sequential : ``bool`` (optional, default: ``True``)
            If ``sequential == True`` each vertex move attempt is made
            sequentially, where vertices are visited in random order. Otherwise
            the moves are attempted by sampling vertices randomly, so that the
            same vertex can be moved more than once, before other vertices had
            the chance to move.
        parallel : ``bool`` (optional, default: ``False``)
            If ``parallel == True``, vertex movements are attempted in parallel.

            .. warning::

               If ``parallel == True``, the asymptotic exactness of the MCMC
               sampling is not guaranteed.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        block_list : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of blocks which are allowed to be
            sampled if ``c == inf``. Otherwise, all blocks will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.

        References
        ----------
        .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
           greedy heuristic for the inference of stochastic block models", Phys.
           Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
           :arxiv:`1310.4378`
        """

        mcmc_state = DictState(locals())
        entropy_args = overlay(dict(dl=True, partition_dl=True, degree_dl=True,
                                    edges_dl=True, dense=False, multigraph=True),
                               **entropy_args)
        if not entropy_args["dl"]:
            entropy_args = overlay(entropy_args, partition_dl=False,
                                   degree_dl=False, edges_dl=False)
        mcmc_state.update(entropy_args)
        mcmc_state.block_list = self.block_list
        if block_list is not None:
            mcmc_state.block_list = Vector_size_t()
            mcmc_state.block_list.extend(block_list)
        mcmc_state.vlist = Vector_size_t()
        if vertices is None:
            idx = self.g.vertex_index.copy().fa
            if self.is_vertex_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                mcmc_state.vlist.extend(idx[vw > 0])
            else:
                mcmc_state.vlist.extend(idx)
        else:
            mcmc_state.vlist.extend(vertices)
        mcmc_state.E = self.E
        mcmc_state.state = self._state

        if _bm_test():
            assert self._check_clabel(), "invalid clabel before sweep"
            Si = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))

        dS, nmoves = self._mcmc_sweep_dispatch(mcmc_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            Sf = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        return dS, nmoves

    def _gibbs_sweep_dispatch(self, gibbs_state):
        return libinference.gibbs_sweep(gibbs_state, self._state,
                                        _get_rng())

    def gibbs_sweep(self, beta=1., niter=1, entropy_args={}, allow_empty=True,
                    sequential=True, parallel=False, vertices=None,
                    block_list=None, verbose=False):
        r"""Perform ``niter`` sweeps of a rejection-free Gibbs sampling MCMC
        to sample network partitions.

        Parameters
        ----------
        beta : ``float`` (optional, default: ``1.``)
            Inverse temperature.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_empty : ``bool`` (optional, default: ``True``)
            Allow movements into empty groups.
        sequential : ``bool`` (optional, default: ``True``)
            If ``sequential == True`` each vertex move attempt is made
            sequentially, where vertices are visited in random order. Otherwise
            the moves are attempted by sampling vertices randomly, so that the
            same vertex can be moved more than once, before other vertices had
            the chance to move.
        parallel : ``bool`` (optional, default: ``False``)
            If ``parallel == True``, vertex movements are attempted in parallel.

            .. warning::

               If ``parallel == True``, the asymptotic exactness of the MCMC
               sampling is not guaranteed.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        block_list : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of blocks which are allowed to be
            sampled if ``c == inf``. Otherwise, all blocks will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.
        """

        gibbs_state = DictState(locals())
        entropy_args = overlay(dict(dl=True, partition_dl=True, degree_dl=True,
                                    edges_dl=True, dense=False, multigraph=True),
                               **entropy_args)
        if not entropy_args["dl"]:
            entropy_args = overlay(entropy_args, partition_dl=False,
                                   degree_dl=False, edges_dl=False)
        gibbs_state.update(entropy_args)
        gibbs_state.block_list = self.block_list
        if block_list is not None:
            gibbs_state.block_list = Vector_size_t()
            gibbs_state.block_list.extend(block_list)
        gibbs_state.vlist = Vector_size_t()
        if vertices is None:
            idx = self.g.vertex_index.copy().fa
            if self.is_vertex_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                gibbs_state.vlist.extend(idx[vw > 0])
            else:
                gibbs_state.vlist.extend(idx)
        else:
            gibbs_state.vlist.extend(vertices)
        gibbs_state.E = self.E
        gibbs_state.state = self._state

        if _bm_test():
            assert self._check_clabel(), "invalid clabel before sweep"
            Si = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))

        dS, nmoves = self._gibbs_sweep_dispatch(gibbs_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            Sf = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        return dS, nmoves

    def _multicanonical_sweep_dispatch(self, multicanonical_state):
        return libinference.multicanonical_sweep(multicanonical_state,
                                                 self._state, _get_rng())

    def multicanonical_sweep(self, m_state, c=1., niter=1, entropy_args={},
                             allow_empty=True, vertices=None, block_list=None,
                             verbose=False):
        r"""Perform ``niter`` sweeps of a non-Markovian multicanonical sampling using
        the Wang-Landau algorithm.

        Parameters
        ----------
        m_state : :class:`~graph_tool.inference.MulticanonicalState`
            :class:`~graph_tool.inference.MulticanonicalState` instance
            containing the current state of the Wang-Landau run.
        c : ``float`` (optional, default: ``1.``)
            Sampling parameter ``c`` for move proposals: For :math:`c\to 0` the
            blocks are sampled according to the local neighbourhood of a given
            node and their block connections; for :math:`c\to\infty` the blocks
            are sampled randomly. Note that only for :math:`c > 0` the MCMC is
            guaranteed to be ergodic.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_empty : ``bool`` (optional, default: ``True``)
            Allow movements into empty groups.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        block_list : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of blocks which are allowed to be
            sampled if ``c == inf``. Otherwise, all blocks will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.

        References
        ----------
        .. [wang-efficient-2001] Fugao Wang, D. P. Landau, "An efficient, multiple
           range random walk algorithm to calculate the density of states", Phys.
           Rev. Lett. 86, 2050 (2001), :doi:`10.1103/PhysRevLett.86.2050`,
           :arxiv:`cond-mat/0011174`

        .. [belardinelli-wang-2007] R. E. Belardinelli, V. D. Pereyra,
           "Wang-Landau algorithm: A theoretical analysis of the saturation of
           the error", J. Chem. Phys. 127, 184105 (2007),
           :doi:`10.1063/1.2803061`, :arxiv:`cond-mat/0702414`
        """

        niter *= self.g.num_vertices()
        args = dmask(locals(), ["self"])
        multi_state = DictState(args)
        entropy_args = overlay(dict(dl=True, partition_dl=True, degree_dl=True,
                                    edges_dl=True, dense=False, multigraph=True),
                               **entropy_args)
        if not entropy_args["dl"]:
            entropy_args = overlay(entropy_args, partition_dl=False,
                                   degree_dl=False, edges_dl=False)
        multi_state.update(entropy_args)
        multi_state.block_list = self.block_list
        if block_list is not None:
            multi_state.block_list = Vector_size_t()
            multi_state.block_list.extend(block_list)
        multi_state.vlist = Vector_size_t()
        if vertices is None:
            idx = self.g.vertex_index.copy().fa
            if self.is_vertex_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                multi_state.vlist.extend(idx[vw > 0])
            else:
                multi_state.vlist.extend(idx)
        else:
            multi_state.vlist.extend(vertices)
        multi_state.E = self.E
        multi_state.S = self.entropy(xi_fast=True, dl_deg_alt=False,
                                     **dmask(entropy_args,
                                             ["xi_fast", "deg_dl_alt"]))
        multi_state.state = self._state

        multi_state.f = m_state._f
        multi_state.time = m_state._time
        multi_state.refine = m_state._refine
        multi_state.S_min = m_state._S_min
        multi_state.S_max = m_state._S_max
        multi_state.hist = m_state._hist
        multi_state.dens = m_state._density

        S, nmoves, f, time = \
                self._multicanonical_sweep_dispatch(multi_state)

        m_state._f = f
        m_state._time = time

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            Sf = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt"]))
            assert abs(S - Sf) < 1e-6, \
                "inconsistent entropy after sweep %g (%g): %s" % \
                (S, Sf, str(entropy_args))

        return S, nmoves

    def _merge_sweep_dispatch(self, merge_state):
        if not self.is_vertex_weighted or not self.is_edge_weighted:
            raise ValueError("state must be weighted to perform merges")

        if merge_state.multigraph and not merge_state.dense:
            raise ValueError("can only merge multigraphs if dense == True")

        return libinference.merge_sweep(merge_state, self._state,
                                        _get_rng())

    def merge_sweep(self, nmerges=1, niter=10, entropy_args={}, parallel=True,
                    verbose=False):
        r"""Perform ``niter`` merge sweeps, where block nodes are progressively
        merged together in a manner that least increases the entropy.

        Parameters
        ----------
        nmerges : ``int`` (optional, default: ``1``)
            Number block nodes to merge.
        niter : ``int`` (optional, default: ``1``)
            Number of merge attempts to perform for each block node, before the
            best one is selected.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        parallel : ``bool`` (optional, default: ``True``)
            If ``parallel == True``, the merge candidates are obtained in
            parallel.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Notes
        -----

        This function should only be called for block states, obtained from
        :meth:`graph_tool.inference.BlockState.get_block_state`.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices merged.

        References
        ----------
        .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
           greedy heuristic for the inference of stochastic block models", Phys.
           Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
           :arxiv:`1310.4378`
        """
        merge_state = DictState(locals())
        merge_state.E = self.E
        merge_state.state = self._state
        entropy_args = overlay(dict(dl=True, partition_dl=True, degree_dl=True,
                                    edges_dl=True, dense=False, multigraph=True),
                               **entropy_args)
        if not entropy_args["dl"]:
            entropy_args = overlay(entropy_args, partition_dl=False,
                                   degree_dl=False, edges_dl=False)
        merge_state.update(entropy_args)

        if _bm_test():
            self._check_clabel()
            Si = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))

        dS, nmoves = self._merge_sweep_dispatch(merge_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            Sf = self.entropy(xi_fast=True, dl_deg_alt=False,
                              **dmask(entropy_args, ["xi_fast", "deg_dl_alt",
                                                     "callback"]))
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        return dS, nmoves

    def shrink(self, B, **kwargs):
        """Reduces the order of current state by progressively merging groups, until
        only ``B`` are left. All remaining keyword arguments are passed to
        :meth:`graph_tool.inference.BlockState.merge_sweep`.

        This function leaves the current state untouched and returns instead a
        copy with the new partition.
        """
        bstate = self.get_block_state(vweight=True,
                                      clabel=self.get_bclabel(),
                                      deg_corr=self.deg_corr)
        nB = (bstate.wr.a > 0).sum()
        while nB > B:
            bstate.merge_sweep(nB - B, **kwargs)
            nB = (bstate.wr.a > 0).sum()
        b = self.b.copy()
        pmap(b, bstate.merge_map)
        continuous_map(b)
        state = self.copy(b=b)
        if _bm_test():
            nB = (state.wr.a > 0).sum()
            assert nB == B, "wrong number of groups after shrink: %d (should be %d)" % (nB, B)
        return state

    def collect_edge_marginals(self, p=None, update=1.):
        r"""Collect the edge marginal histogram, which counts the number of times
        the endpoints of each node have been assigned to a given block pair.

        This should be called multiple times, e.g. after repeated runs of the
        :meth:`graph_tool.inference.BlockState.mcmc_sweep` function.

        Parameters
        ----------
        p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Edge property map with vector-type values, storing the previous block
            membership counts.  Each vector entry corresponds to ``b[i] + B *
            b[j]``, where ``b`` is the block membership and ``i = min(source(e),
            target(e))`` and ``j = max(source(e), target(e))``. If not provided, an
            empty histogram will be created
        update : float (optional, default: ``1.``)
            Each call increases the current count by the amount given by this
            parameter.

        Returns
        -------
        p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Vertex property map with vector-type values, storing the accumulated
            block membership counts.


        Examples
        --------
        .. testsetup:: collect_edge_marginals

           gt.seed_rng(42)
           np.random.seed(42)

        .. doctest:: collect_edge_marginals

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=4, deg_corr=True)
           >>> pe = None
           >>> state.mcmc_sweep(niter=1000)   # remove part of the transient
           (...)
           >>> for i in range(1000):
           ...     ds, nmoves = state.mcmc_sweep(niter=10)
           ...     pe = state.collect_edge_marginals(pe)
           >>> gt.bethe_entropy(g, state.B, pe)[0]
           4.806605480346...
        """

        if p is None:
            p = self.g.new_ep("vector<double>")

        libinference.edge_marginals(self.g._Graph__graph,
                                    self.B,
                                    _prop("v", self.g, self.b),
                                    _prop("e", self.g, p),
                                    update)
        return p

    def collect_vertex_marginals(self, p=None, update=1.):
        r"""Collect the vertex marginal histogram, which counts the number of times a
        node was assigned to a given block.

        This should be called multiple times, e.g. after repeated runs of the
        :meth:`graph_tool.inference.BlockState.mcmc_sweep` function.

        Parameters
        ----------
        p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Vertex property map with vector-type values, storing the previous block
            membership counts. If not provided, an empty histogram will be created.
        update : float (optional, default: ``1.``)
            Each call increases the current count by the amount given by this
            parameter.

        Returns
        -------
        p : :class:`~graph_tool.PropertyMap`
            Vertex property map with vector-type values, storing the accumulated
            block membership counts.

        Examples
        --------
        .. testsetup:: cvm

           gt.seed_rng(42)
           np.random.seed(42)

        .. doctest:: cvm

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=4, deg_corr=True)
           >>> pv = None
           >>> state.mcmc_sweep(niter=1000)   # remove part of the transient
           (...)
           >>> for i in range(1000):
           ...     ds, nmoves = state.mcmc_sweep(niter=10)
           ...     pv = state.collect_vertex_marginals(pv)
           >>> gt.mf_entropy(g, pv)
           4.989264427809...
           >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie",
           ...               vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.pdf")
           <...>

        .. testcleanup:: cvm

           gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv,
                         output="polbooks_blocks_soft_B4.png")

        .. figure:: polbooks_blocks_soft_B4.*
           :align: center

           "Soft" block partition of a political books network with :math:`B=4`.

        """

        if p is None:
            p = self.g.new_vp("vector<double>")

        libinference.vertex_marginals(self.g._Graph__graph,
                                      _prop("v", self.g, self.b),
                                      _prop("v", self.g, p),
                                      update)
        return p

    def draw(self, **kwargs):
        r"""Convenience wrapper to :func:`~graph_tool.draw.graph_draw` that
        draws the state of the graph as colors on the vertices and edges."""
        gradient = self.g.new_ep("double")
        gradient = group_vector_property([gradient])
        from graph_tool.draw import graph_draw
        return graph_draw(self.g,
                          vertex_fill_color=kwargs.get("vertex_fill_color",
                                                       self.b),
                          vertex_color=kwargs.get("vertex_color", self.b),
                          edge_gradient=kwargs.get("edge_gradient",
                                                   gradient),
                          **dmask(kwargs, ["vertex_fill_color",
                                           "vertex_color",
                                           "edge_gradient"]))

def model_entropy(B, N, E, directed=False, nr=None):
    r"""Computes the amount of information necessary for the parameters of the
    traditional blockmodel ensemble, for ``B`` blocks, ``N`` vertices, ``E``
    edges, and either a directed or undirected graph.

    A traditional blockmodel is defined as a set of :math:`N` vertices which can
    belong to one of :math:`B` blocks, and the matrix :math:`e_{rs}` describes
    the number of edges from block :math:`r` to :math:`s` (or twice that number
    if :math:`r=s` and the graph is undirected).

    For an undirected graph, the number of distinct :math:`e_{rs}` matrices is
    given by,

    .. math::

       \Omega_m = \left(\!\!{\left(\!{B \choose 2}\!\right) \choose E}\!\!\right)

    and for a directed graph,

    .. math::
       \Omega_m = \left(\!\!{B^2 \choose E}\!\!\right)


    where :math:`\left(\!{n \choose k}\!\right) = {n+k-1\choose k}` is the
    number of :math:`k` combinations with repetitions from a set of size
    :math:`n`.

    The total information necessary to describe the model is then,

    .. math::

       \mathcal{L}_t = \ln\Omega_m + \ln\left(\!\!{B \choose N}\!\!\right) +
       \ln N! - \sum_r \ln n_r!,


    where the remaining term is the information necessary to describe the block
    partition, where :math:`n_r` is the number of nodes in block :math:`r`.

    If ``nr`` is ``None``, it is assumed :math:`n_r=N/B`.

    References
    ----------

    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module
       inference in large networks", Phys. Rev. Lett. 110, 148701 (2013),
       :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block
       structures and high-resolution model selection in large networks ",
       Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.

    """

    if directed:
        x = (B * B);
    else:
        x = (B * (B + 1)) / 2;
    if nr is False:
        L = lbinom(x + E - 1, E)
    else:
        L = lbinom(x + E - 1, E) + partition_entropy(B, N, nr)
    return L

def bethe_entropy(g, B, p):
    r"""Compute the Bethe entropy given the edge block membership marginals.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph.
    B : int
        The number of blocks.
    p : :class:`~graph_tool.PropertyMap`
        Edge property map with vector-type values, storing the previous block
        membership counts.  Each vector entry corresponds to ``b[i] + B *
        b[j]``, where ``b`` is the block membership and ``i = min(source(e),
        target(e))`` and ``j = max(source(e), target(e))``.

    Returns
    -------
    H : ``float``
        The Bethe entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_)
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_),
        as would be returned by the :func:`mf_entropy` function.
    pv : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the accumulated
        block membership counts. These are the node marginals, as would be
        returned by the :func:`collect_vertex_marginals` function.

    Notes
    -----

    The Bethe entropy is defined as,

    .. math::

        H = -\sum_{e,(r,s)}\pi_{(r,s)}^e\ln\pi_{(r,s)}^e - \sum_{v,r}(1-k_i)\pi_r^v\ln\pi_r^v,

    where :math:`\pi_{(r,s)}^e` is the marginal probability that the endpoints
    of the edge :math:`e` belong to blocks :math:`(r,s)`, and :math:`\pi_r^v` is
    the marginal probability that vertex :math:`v` belongs to block :math:`r`,
    and :math:`k_i` is the degree of vertex :math:`v` (or total degree for
    directed graphs).

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
       :DOI:`10.1093/acprof:oso/9780198570837.001.0001`
    """
    H = 0
    pv =  g.new_vertex_property("vector<double>")

    H, sH, Hmf, sHmf  = libinference.bethe_entropy(g._Graph__graph, B,
                                                   _prop("e", g, p),
                                                   _prop("v", g, pv))
    return H, Hmf, pv


def mf_entropy(g, p):
    r"""Compute the "mean field" entropy given the vertex block membership marginals.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph.
    p : :class:`~graph_tool.PropertyMap`
        Vertex property map with vector-type values, storing the accumulated block
        membership counts.

    Returns
    -------
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_).

    Notes
    -----

    The "mean field" entropy is defined as,

    .. math::

        H = - \sum_{v,r}\pi_r^v\ln\pi_r^v,

    where :math:`\pi_r^v` is the marginal probability that vertex :math:`v`
    belongs to block :math:`r`.

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
       :DOI:`10.1093/acprof:oso/9780198570837.001.0001`
    """
    H = 0
    for v in g.vertices():
        N = p[v].a.sum()
        if N == 0:
            continue
        pvi = asarray(p[v].a, dtype="float") /  N
        pvi = pvi[pvi > 0]
        H -= (pvi * log(pvi)).sum()
    return H

from . overlap_blockmodel import *