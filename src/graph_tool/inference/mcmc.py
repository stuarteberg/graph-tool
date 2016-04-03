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

from .. import Vector_size_t, Vector_double

import numpy
from . util import *

def mcmc_equilibrate(state, wait=10, nbreaks=2, max_niter=numpy.inf,
                     force_niter=None, epsilon=0, gibbs=False,
                     block_moves=False, mcmc_args={}, entropy_args={},
                     history=False, callback=None, verbose=False):
    r"""Equilibrate a MCMC with a given starting state.

    Parameters
    ----------
    state : Any state class (e.g. :class:`~graph_tool.inference.BlockState`)
        Initial state. This state will be modified during the algorithm.
    wait : ``int`` (optional, default: ``10``)
        Number of iterations to wait for a record-breaking event.
    nbreaks : ``int`` (optional, default: ``2``)
        Number of iteration intervals (of size ``wait``) without record-breaking
        events necessary to stop the algorithm.
    max_niter : ``int`` (optional, default: ``numpy.inf``)
        Maximum number of iterations.
    force_niter : ``int`` (optional, default: ``None``)
        If given, will force the algorithm to run this exact number of
        iterations.
    epsilon : ``float`` (optional, default: ``0``)
        Relative changes in entropy smaller than epsilon will not be considered
        as record-breaking.
    gibbs : ``bool`` (optional, default: ``False``)
        If ``True``, each step will call ``state.gibbs_sweep`` instead of
        ``state.mcmc_sweep``.
    block_moves : ``bool`` (optional, default: ``False``)
        If ``True``, each iteration will be accompanied by a "block move", where
        all vertices of the same group are moved simultaneously.
    mcmc_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to ``state.mcmc_sweep`` (or ``state.gibbs_sweep``).
    history : ``bool`` (optional, default: ``False``)
        If ``True``, a list of tuples of the form ``(iteration, entropy)`` will
        be kept and returned.
    callback : ``function`` (optional, default: ``None``)
        If given, this function will be called after each iteration. The
        function must accept the current state as an argument, and its return
        value must be either `None` or a (possibly empty) list of values that
        will be append to the history, if ``history == True``.
    verbose : ``bool`` or ``tuple`` (optional, default: ``False``)
        If ``True``, progress information will be shown. Optionally, this
        accepts arguments of the type ``tuple`` of the form ``(level, prefix)``
        where ``level`` is a positive integer that specifies the level of
        detail, and ``prefix`` is a string that is prepended to the all output
        messages.

    Notes
    -----

    The MCMC equilibration is attempted by keeping track of the maximum and
    minimum values, and waiting sufficiently long without a record-breaking
    event.

    This function calls ``state.mcmc_sweep`` (or ``state.gibbs_sweep``) at each
    iteration (e.g. :meth:`graph_tool.inference.BlockState.mcmc_sweep` and
    :meth:`graph_tool.inference.BlockState.gibbs_sweep`), and keeps track of
    the value of ``state.entropy(**args)`` with ``args`` corresponding to
    ``mcmc_args["entropy_args"]``.

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
       greedy heuristic for the inference of stochastic block models", Phys.
       Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
       :arxiv:` 1310.4378`
    """

    count = 0
    break_count = 0
    niter = 0
    total_nmoves = 0
    S = state.entropy(**mcmc_args.get("entropy_args", {}))
    min_S = max_S = S
    m_eps = 1e-6
    hist = []
    while count < wait:
        if not gibbs:
            delta, nmoves = state.mcmc_sweep(**mcmc_args)
        else:
            delta, nmoves = state.gibbs_sweep(**mcmc_args)

        if block_moves:
            bstate = state.get_block_state(vweight=True,
                                           clabel=state.get_bclabel())
            if not gibbs:
                ret = bstate.mcmc_sweep(**mcmc_args)
            else:
                ret = bstate.gibbs_sweep(**mcmc_args)

            b = state.b.copy()
            pmap(b, bstate.b)
            state.set_blocks(b)

            delta += ret[0]
            nmoves += ret[1]

        S += delta
        niter += 1
        total_nmoves += nmoves

        if force_niter is not None:
            if niter >= force_niter:
                break
        else:
            if abs(delta) >= (S - delta) * epsilon:
                if S > max_S + m_eps:
                    max_S = S
                    count = 0
                elif S < min_S - m_eps:
                    min_S = S
                    count = 0
                else:
                    count += 1
            else:
                count += 1

            if count >= wait:
                break_count += 1
                if break_count < nbreaks:
                    count = 0
                    min_S = max_S = S

        extra = []
        if callback is not None:
            extra = callback(state)
            if extra is None:
                extra = []

        if check_verbose(verbose):
            print((verbose_pad(verbose) +
                   u"niter: %5d  count: %4d  breaks: %2d  min_S: %#8.8g  " +
                   u"max_S: %#8.8g  S: %#8.8g  ΔS: %#12.6g  moves: %5d %s") %
                   (niter, count, break_count, min_S, max_S, S, delta, nmoves,
                    str(extra) if len(extra) > 0 else ""))

        if history:
            hist.append([S, nmoves] + extra)

        if niter >= max_niter:
            break

    if history:
        return hist
    else:
        return (S, total_nmoves)

def mcmc_anneal(state, beta_range=(1., 10.), niter=100, history=False,
                mcmc_equilibrate_args={}, verbose=False):
    r"""Equilibrate a MCMC at a specified target temperature by performing simulated
    annealing.

    Parameters
    ----------
    state : Any state class (e.g. :class:`~graph_tool.inference.BlockState`)
        Initial state. This state will be modified during the algorithm.
    beta_range : ``tuple`` of two floats (optional, default: ``(1., 10.)``)
        Inverse temperature range.
    niter : ``int`` (optional, default: ``100``)
        Number of steps (in logspace) from the starting temperature to the final
        one.
    history : ``bool`` (optional, default: ``False``)
        If ``True``, a list of tuples of the form ``(iteration, beta, entropy)``
    mcmc_equilibrate_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to :func:`~graph_tool.inference.mcmc_equilibrate`.
    verbose : ``bool`` or ``tuple`` (optional, default: ``False``)
        If ``True``, progress information will be shown. Optionally, this
        accepts arguments of the type ``tuple`` of the form ``(level, prefix)``
        where ``level`` is a positive integer that specifies the level of
        detail, and ``prefix`` is a string that is prepended to the all output
        messages.

    Notes
    -----

    This algorithm employs exponential cooling, where the value of beta is
    multiplied by a constant at each iteration, so that starting from
    `beta_range[0]` the value of `beta_range[1]` is reached after `niter`
    iterations.

    At each iteration, the function
    :func:`~graph_tool.inference.mcmc_equilibrate` is called with the current
    value of `beta` (via the ``mcmc_args`` parameter).

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
       greedy heuristic for the inference of stochastic block models", Phys.
       Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
       :arxiv:` 1310.4378`
    """

    beta = beta_range[0]
    hist = ([], [], [])
    nmoves = 0
    speed = exp((log(beta_range[1]) - log(beta_range[0])) / nsteps)
    mcmc_args = mcmc_equilibrate_args.get("mcmc_args", {})
    while beta < beta_range[1] * speed:
        ret = mcmc_equilibrate(state,
                               **overlay(mcmc_equilibrate_args,
                                         mcmc_args=overlay(mcmc_args,
                                                           beta=beta),
                                         history=history,
                                         verbose=verbose_push(verbose,
                                                              ("β: %#8.6g  " %
                                                               beta))))
        if history:
            ret = list(zip(*ret))
            hist[0].extend([beta] * len(ret[0]))
            hist[1].extend(ret[0])
            hist[2].extend(ret[1])
            S = ret[0][-1]
        else:
            S = ret[0]
            nmoves += ret[1]

        beta *= speed

    if history:
        return list(zip(hist))
    else:
        return S, nmoves


def mcmc_multilevel(state, B, r=2, b_cache=None, anneal=False,
                    mcmc_equilibrate_args={}, anneal_args={}, shrink_args={},
                    verbose=False):
    r"""Equilibrate a MCMC from a starting state with a higher order, by performing
    successive agglomerative initializations and equilibrations until the
    desired order is reached, such that metastable states are avoided.

    Parameters
    ----------
    state : Any state class (e.g. :class:`~graph_tool.inference.BlockState`)
        Initial state. This state will **not** be modified during the algorithm.
    B : ``int``
        Desired model order (i.e. number of groups).
    r : ``int`` (optional, default: ``2``)
        Greediness of agglomeration. At each iteration, the state order will be
        reduced by a factor ``r``.
    b_cache : ``dict`` (optional, default: ``None``)
        If specified, this should be a dictionary with key-value pairs of the
        form ``(B, state)`` that contain pre-computed states of the specified
        order. This dictionary will be modified during the algorithm.
    anneal : ``bool`` (optional, default: ``False``)
        If ``True``, the equilibration steps will use simulated annealing, by
        calling :func:`~graph_tool.inference.mcmc_anneal`, instead of
        :func:`~graph_tool.inference.mcmc_equilibrate`.
    mcmc_equilibrate_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to :func:`~graph_tool.inference.mcmc_equilibrate`.
    mcmc_anneal_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to :func:`~graph_tool.inference.mcmc_anneal`.
    shrink_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to ``state.shrink``
        (e.g. :meth:`graph_tool.inference.BlockState.shrink`).
    verbose : ``bool`` or ``tuple`` (optional, default: ``False``)
        If ``True``, progress information will be shown. Optionally, this
        accepts arguments of the type ``tuple`` of the form ``(level, prefix)``
        where ``level`` is a positive integer that specifies the level of
        detail, and ``prefix`` is a string that is prepended to the all output
        messages.

    Notes
    -----

    This algorithm alternates between equilibrating the MCMC state and reducing
    the state order (via calls to ``state.shrink``,
    e.g. :meth:`graph_tool.inference.BlockState.shrink`).

    This greatly reduces the changes of getting trapped in metastable states if
    the starting point if far away from equilibrium, as discussed in
    [peixoto-efficient-2014]_.

    References
    ----------

    .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
       greedy heuristic for the inference of stochastic block models", Phys.
       Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
       :arxiv:`1310.4378`
    """

    if "mcmc_equilibrate_args" in anneal_args:
        raise ValueError("'mcmc_equilibrate_args' should be passed directly " +
                         "to mcmc_multilevel(), not via 'anneal_args'")
    while state.B > B:
        B_next = max(min(int(round(state.B / r)), state.B - 1), B)

        if b_cache is not None and B_next in b_cache:
            state = b_cache[B_next]
            if check_verbose(verbose):
                print(verbose_pad(verbose) +
                      "shrinking %d -> %d (cached)" % (state.B, B_next))
            continue

        if check_verbose(verbose):
            print(verbose_pad(verbose) +
                  "shrinking %d -> %d" % (state.B, B_next))
        state = state.shrink(B=B_next, **shrink_args)
        if anneal:
            mcmc_anneal(state,
                        **overlay(anneal_args,
                                  mcmc_equilibrate_args=mcmc_equilibrate_args,
                                  verbose=verbose_push(verbose,
                                                       "B=%d  " % state.B)))
        else:
            mcmc_equilibrate(state,
                             **overlay(mcmc_equilibrate_args,
                                       verbose=verbose_push(verbose,
                                                            ("B=%d  " %
                                                             state.B))))
        if b_cache is not None:
            mcmc_args = mcmc_equilibrate_args.get("mcmc_args", {})
            entropy_args = mcmc_args.get("entropy_args", {})
            b_cache[B_next] = (state.entropy(**entropy_args), state)
    return state


class MulticanonicalState(object):
    r"""The density of states of a multicanonical Monte Carlo algorithm. It is used
    by :func:`graph_tool.inference.multicanonical_equilibrate`.

    Parameters
    ----------
    S_min : ``float``
        Minimum energy.
    S_max : ``float``
        Maximum energy.
    nbins : ``int`` (optional, default: ``1000``)
        Number of bins.
    """

    def __init__(self, S_min, S_max, nbins=1000):
        self._S_min = S_min
        self._S_max = S_max
        self._density = Vector_double()
        self._density.resize(nbins)
        self._hist = Vector_size_t()
        self._hist.resize(nbins)

    def __getstate__(self):
        state = [self._S_min, self._S_max,
                 numpy.array(self._density.get_array()),
                 numpy.array(self._hist.get_array())]
        return state

    def __setstate__(self, state):
        S_min, S_max, density, hist = state
        self.__init__(S_min, S_max, len(hist))
        self._density.get_array()[:] = density
        self._hist.get_array()[:] = hist
        return state

    def get_energies(self):
        "Get energy bounds."
        return self._S_min, self._S_max

    def get_allowed_energies(self):
        "Get allowed energy bounds."
        h = self._hist.get_array()
        Ss = self.get_range()
        Ss = Ss[h > 0]
        return Ss[0], Ss[-1]

    def get_range(self):
        "Get energy range."
        return numpy.linspace(self._S_min, self._S_max, len(self._hist))

    def get_density(self):
        """Get density of states, normalized so that the **integral** over the energy
        range is unity."""
        r = numpy.array(self._density.get_array())
        r -= r.max()
        N = len(r)
        dS = (self._S_max - self._S_min) / N
        I = exp(r).sum() * dS
        r -= log(I)
        return r

    def get_prob(self, S):
        r = self.get_density()
        dS = (self._S_max - self._S_min) / N
        j = round((S - self._S_min) / dS)
        return r[j]

    def get_hist(self):
        "Get energy histogram."
        return numpy.array(self._hist.get_array())

    def get_flatness(self):
        "Get energy histogram flatness."
        h = self._hist.get_array()
        if h.sum() == 0:
            return numpy.inf
        h = array(h[h>0], dtype="float")
        h /= h.sum()
        S = -(h * log(h)).sum()
        return len(h) / exp(S) - 1

    def get_mean(self):
        "Get energy mean."
        r = self.get_density()
        N = len(r)
        dS = (self._S_max - self._S_min) / N
        Ss = numpy.linspace(self._S_min, self._S_max, N)
        return (Ss * exp(r) * dS).sum()

    def get_posterior(self):
        "Get posterior mean."
        r = self.get_density()
        N = len(r)
        dS = (self._S_max - self._S_min) / N
        Ss = numpy.linspace(self._S_min, self._S_max, N)
        y = -Ss + r
        y_max = y.max()
        y -= y_max
        return y_max + log((exp(y) * dS).sum())

    def reset_hist(self):
        "Reset energy histogram."
        self._hist.get_array()[:] = 0

def multicanonical_equilibrate(state, m_state, f_range=(1., 1e-6), r=2,
                               flatness=.01, callback=None,
                               multicanonical_args={}, verbose=False):
    r"""Equilibrate a multicanonical Monte Carlo sampling using the Wang-Landau
     algorithm.

    Parameters
    ----------
    state : Any state class (e.g. :class:`~graph_tool.inference.BlockState`)
        Initial state. This state will be modified during the algorithm.
    m_state :  :class:`~graph_tool.inference.MulticanonicalState`
        Initial multicanonical state, where the state density will be stored.
    f_range : ``tuple`` of two floats (optional, default: ``(1., 1e-6)``)
        Range of density updates.
    r : ``float`` (optional, default: ``2.``)
        Greediness of convergence. At each iteration, the density updates will
        be reduced reduced by a factor ``r``.
    flatness : ``float`` (optional, default: ``1e-3``)
        Sufficient histogram flatness threshold used to continue the algorithm.
    callback : ``function`` (optional, default: ``None``)
        If given, this function will be called after each iteration. The
        function must accept the current ``state`` and ``m_state`` as arguments.
    multicanonical_args : ``dict`` (optional, default: ``{}``)
        Arguments to be passed to ``state.multicanonical_sweep`` (e.g.
        :meth:`graph_tool.inference.BlockState.multicanonical_sweep`).
    verbose : ``bool`` or ``tuple`` (optional, default: ``False``)
        If ``True``, progress information will be shown. Optionally, this
        accepts arguments of the type ``tuple`` of the form ``(level, prefix)``
        where ``level`` is a positive integer that specifies the level of
        detail, and ``prefix`` is a string that is prepended to the all output
        messages.

    References
    ----------

    .. [wang-efficient-2001] Fugao Wang, D. P. Landau, "An efficient, multiple
       range random walk algorithm to calculate the density of states", Phys.
       Rev. Lett. 86, 2050 (2001), :doi:`10.1103/PhysRevLett.86.2050`,
       :arxiv:`cond-mat/0011174'
    """

    f = f_range[0]
    while f >= f_range[1]:
        state.multicanonical_sweep(m_state, **overlay(multicanonical_args, f=f))
        hf = m_state.get_flatness()
        if hf < flatness:
            f /= r
            if f >= f_range[1]:
                m_state.reset_hist()

        if callback is not None:
            calback(state, m_state)

        if check_verbose(verbose):
            print(verbose_pad(verbose) + "f: %g  flatness: %g" % (f, hf))

    return m_state