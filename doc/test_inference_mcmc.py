#!/bin/env python

from __future__ import print_function

from pylab import *
from graph_tool.all import *
import numpy.random
from numpy.random import randint
import scipy.stats

numpy.random.seed(43)
seed_rng(43)

verbose = __name__ == "__main__"

g = collection.data["football"]
B = 10


for directed in [True, False]:
    clf()
    g.set_directed(directed)

    state = minimize_blockmodel_dl(g, deg_corr=False, B_min=B, B_max=B)
    state = state.copy(B=B+1)

    c = 0.01
    for v in g.vertices():
        r = state.b[v]
        res = zeros(state.B)
        for s in range(state.B):
            pf = state.get_move_prob(v, s, c=c, reverse=True)
            state.move_vertex(v, s)
            pb = state.get_move_prob(v, r, c=c)
            res[s] = pf - pb
            if abs(res[s]) > 1e-8:
                print("Warning, wrong reverse probability: ", v, r, s, pf, pb,
                      res[s], directed)
            state.move_vertex(v, r)
        plot(res)
    gca().set_ylim(-.1, .1)

    savefig("test_mcmc_reverse_prob_directed%s.pdf" % directed)

for directed in [True, False]:
    clf()
    g.set_directed(directed)

    state = minimize_blockmodel_dl(g, deg_corr=False, B_min=B, B_max=B)
    state = state.copy(B=B+1)

    c = 0.1
    for v in g.vertices():

        # computed probabilities
        mp = zeros(state.B)
        for s in range(state.B):
            mp[s] = state.get_move_prob(v, s, c)

        n_samples = min(int(200 / mp.min()), 1000000)

        # actual samples
        samples = [state.sample_vertex_move(v, c) for i in range(n_samples)]

        # samples from computed distribution
        true_hist = numpy.random.multinomial(n_samples, mp)
        true_samples = []
        for r, count in enumerate(true_hist):
            true_samples.extend([r] * count)

        mp_h = bincount(samples)
        if len(mp_h) < B + 1:
            mp_h = list(mp_h) + [0] * (B + 1 - len(mp_h))
        mp_h = array(mp_h, dtype="float")
        mp_h /= mp_h.sum()
        res = mp - mp_h

        samples = array(samples, dtype="float")
        true_samples = array(true_samples, dtype="float")

        p = scipy.stats.ks_2samp(samples, true_samples)[1]

        if verbose:
            print("directed:", directed, "vertex:", v, "p-value:", p)

        if p < 0.001:
            print(("Warning, move probability for node %d does not " +
                   "match the computed distribution, with p-value: %g") %
                  (v, p))
            clf()
            plot(res)
            savefig("test_mcmc_move_prob_directed%s_v%d.pdf" % (directed, int(v)))

        plot(res)
    gca().set_ylim(-.1, .1)

    savefig("test_mcmc_move_prob_directed%s.pdf" % directed)

g = collection.data["karate"]
B = 3

for directed in [True, False]:
    clf()
    g.set_directed(directed)

    hists = {}

    state = minimize_blockmodel_dl(g, deg_corr=False, B_min=B, B_max=B)
    state = state.copy(B=B+1)

    cs = list(reversed([numpy.inf, 1, 0.1, 0.01, 0.001, "gibbs"]))

    for i, c in enumerate(cs):
        if c != "gibbs":
            mcmc_args=dict(beta=1, c=c, niter=300, allow_empty=True)
        else:
            mcmc_args=dict(beta=1, niter=300, allow_empty=True)
        if i == 0:
            mcmc_equilibrate(state, mcmc_args=mcmc_args, gibbs=c=="gibbs",
                             wait=1000,
                             verbose=(1, "c = %s (t) " % str(c))  if verbose else False)
        hists[c] = mcmc_equilibrate(state, mcmc_args=mcmc_args,
                                    gibbs=c=="gibbs",
                                    force_niter=2000,
                                    verbose=(1, "c = %s " % str(c)) if verbose else False,
                                    history=True)

    for c1 in cs:
        for c2 in cs:
            try:
                if c2 < c1:
                    continue
            except TypeError:
                pass
            Ss1 = array(list(zip(*hists[c1]))[0])
            Ss2 = array(list(zip(*hists[c2]))[0])
            # add very small normal noise, to solve discreetness issue
            Ss1 += numpy.random.normal(0, 1e-6, len(Ss1))
            Ss2 += numpy.random.normal(0, 1e-6, len(Ss2))
            D, p = scipy.stats.ks_2samp(Ss1, Ss2)
            if verbose:
                print("directed:", directed, "c1:", c1, "c2:", c2,
                      "D", D, "p-value:", p)
            if p < .001:
                print(("Warning, distributions for directed=%s (c1, c2) = " +
                       "(%s, %s) are not the same, with a p-value: %g (D=%g)") %
                      (str(directed), str(c1), str(c2), p, D))

    bins = None
    for c in cs:
        hist = hists[c]
        h = histogram(list(zip(*hist))[0], 1000000)
        plot(h[-1][:-1], numpy.cumsum(h[0]), "-", label="c=%s" % str(c))

        if c != numpy.inf:
            hist = hists[numpy.inf]
            h2 = histogram(list(zip(*hist))[0], bins=h[-1])
            res = abs(numpy.cumsum(h2[0]) - numpy.cumsum(h[0]))
            i = res.argmax()
            axvline(h[-1][i], color="grey")

    legend(loc="lower right")
    savefig("test_mcmc_directed%s.pdf" % directed)

print("OK")