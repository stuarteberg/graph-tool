#!/bin/env python

from __future__ import print_function

verbose = __name__ == "__main__"

import os
import sys
if not verbose:
    out = open(os.devnull, 'w')
else:
    out = sys.stdout
from collections import OrderedDict
import itertools
from graph_tool.all import *
import numpy.random
from numpy.random import randint, normal, random

numpy.random.seed(42)
seed_rng(42)

graph_tool.inference.set_test(True)

g = collection.data["football"]
ec = g.new_ep("int", randint(0, 10, g.num_edges()))

rec_p = g.new_ep("double", random(g.num_edges()))
rec_s = g.new_ep("double", normal(0, 10, g.num_edges()))


def gen_state(directed, deg_corr, layers, overlap, rec, allow_empty):
    u = GraphView(g, directed=directed)
    if layers != False:
        state = graph_tool.inference.LayeredBlockState(u, B=u.num_vertices(),
                                                       deg_corr=deg_corr,
                                                       ec=ec.copy(),
                                                       rec=rec,
                                                       overlap=overlap,
                                                       layers=layers == True,
                                                       allow_empty=allow_empty)
    elif overlap:
        state = graph_tool.inference.OverlapBlockState(u, B=2 * u.num_edges(),
                                                       rec=rec,
                                                       deg_corr=deg_corr)
    else:
        state = graph_tool.inference.BlockState(u, B=u.num_vertices(),
                                                rec=rec,
                                                deg_corr=deg_corr,
                                                allow_empty=allow_empty)
    return state


pranges = [("directed", [False, True]),
           ("overlap", [False, True]),
           ("layered", [False, "covariates", True]),
           ("rec", [None, "positive", "signed"]),
           ("deg_corr", [False, True]),
           ("dl", [False, True]),
           ("degree_dl_kind", ["uniform", "distributed", "entropy"]),
           ("allow_empty", [False, True]),
           ("exact", [True, False])]

pranges = OrderedDict(pranges)

def iter_ranges(ranges):
    for vals in itertools.product(*[v for k, v in ranges.items()]):
        yield zip(ranges.keys(), vals)


for pvals in iter_ranges(pranges):
    params = OrderedDict(pvals)

    locals().update(params)

    if not deg_corr and degree_dl_kind != "uniform":
        continue

    if overlap and degree_dl_kind != "distributed":
        continue

    if (rec is not None or layered != False) and not exact:
        continue

    print(params, file=out)

    rec_ = None
    if rec == "positive":
        rec_ = rec_p
    elif rec == "signed":
        rec_ = rec_s

    print("\t mcmc (unweighted)", file=out)
    state = gen_state(directed, deg_corr, layered, overlap, rec_, allow_empty)

    print("\t\t",
          state.mcmc_sweep(beta=0,
                           entropy_args=dict(dl=dl,
                                             degree_dl_kind=degree_dl_kind,
                                             exact=exact)),
          state.get_nonempty_B(), file=out)

    if overlap:
        print("\t\t",
              state.mcmc_sweep(beta=0, bundled=True,
                               entropy_args=dict(dl=dl,
                                                 degree_dl_kind=degree_dl_kind,
                                                 exact=exact)),
              state.get_nonempty_B(), file=out)

    state = gen_state(directed, deg_corr, layered, overlap, rec_, allow_empty)

    if not overlap:
        print("\t mcmc", file=out)
        bstate = state.get_block_state(vweight=True,  deg_corr=deg_corr)

        print("\t\t",
              bstate.mcmc_sweep(beta=0,
                                entropy_args=dict(dl=dl,
                                                  degree_dl_kind=degree_dl_kind,
                                                  exact=exact)),
              bstate.get_nonempty_B(), file=out)

        print("\t\t",
              bstate.mcmc_sweep(beta=0,
                                entropy_args=dict(dl=dl,
                                                  degree_dl_kind=degree_dl_kind,
                                                  exact=exact)),
              bstate.get_nonempty_B(), file=out)

        print("\t\t",
              bstate.gibbs_sweep(beta=0,
                                 entropy_args=dict(dl=dl,
                                                   degree_dl_kind=degree_dl_kind,
                                                   exact=exact)),
              bstate.get_nonempty_B(), file=out)

    print("\t merge", file=out)

    state = gen_state(directed, deg_corr, layered, overlap, rec_, allow_empty)

    if not overlap:
        bstate = state.get_block_state(vweight=True, deg_corr=deg_corr)

        print("\t\t",
              bstate.merge_sweep(50,
                                 entropy_args=dict(dl=dl,
                                                   degree_dl_kind=degree_dl_kind,
                                                   multigraph=False,
                                                   exact=exact)),
              file=out)

        bstate = bstate.copy()

        print("\t\t",
              bstate.mcmc_sweep(beta=0,
                                entropy_args=dict(dl=dl,
                                                  degree_dl_kind=degree_dl_kind,
                                                  exact=exact)),
              file=out)
        print("\t\t",
              bstate.gibbs_sweep(beta=0,
                                 entropy_args=dict(dl=dl,
                                                   degree_dl_kind=degree_dl_kind,
                                                   exact=exact)),
              file=out)
    else:
        print("\t\t",
              state.merge_sweep(50,
                                entropy_args=dict(dl=dl,
                                                  degree_dl_kind=degree_dl_kind,
                                                  multigraph=False,
                                                  exact=exact)),
              file=out)

    print("\t shrink", file=out)

    state = gen_state(directed, deg_corr, layered, overlap, rec_, allow_empty)
    state = state.shrink(B=5, entropy_args=dict(dl=dl,
                                                degree_dl_kind=degree_dl_kind,
                                                multigraph=False,
                                                exact=exact))
    print("\t\t", state.B, "\n", file=out)


pranges = [("directed", [False, True]),
           ("overlap", [False, True]),
           ("layered", [False, "covariates", True]),
           ("rec", [None, "positive", "signed"]),
           ("deg_corr", [False, True]),
           ("degree_dl_kind", ["uniform", "distributed", "entropy"]),
           ("exact", [True])]

pranges = OrderedDict(pranges)

for pvals in iter_ranges(pranges):
    params = OrderedDict(pvals)

    locals().update(params)

    if not deg_corr and degree_dl_kind != "uniform":
        continue

    if overlap and degree_dl_kind != "distributed":
        continue

    print(params, file=out)

    rec_ = None
    if rec == "positive":
        rec_ = rec_p
    elif rec == "signed":
        rec_ = rec_s

    if layered != False:
        state_args = dict(ec=ec, layers=(layered == True), rec=rec_)
    else:
        state_args = dict(rec=rec_)

    entropy_args = dict(exact=exact)

    state = minimize_blockmodel_dl(GraphView(g, directed=directed),
                                   verbose=(1, "\t") if verbose else False,
                                   deg_corr=deg_corr,
                                   overlap=overlap,
                                   layers=layered != False,
                                   state_args=state_args,
                                   mcmc_args=dict(entropy_args=entropy_args))

    print(state.B, state.entropy(), file=out)

    state = minimize_nested_blockmodel_dl(GraphView(g, directed=directed),
                                          verbose=(2, "\t") if verbose else False,
                                          deg_corr=deg_corr,
                                          overlap=overlap,
                                          layers=layered != False,
                                          state_args=state_args,
                                          mcmc_args=dict(entropy_args=entropy_args))
    if verbose:
        state.print_summary()

    print(state.entropy(), "\n", file=out)

graph_tool.inference.set_test(False)

print("OK")