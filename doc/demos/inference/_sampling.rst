.. _sampling:

Sampling from the posterior distribution
----------------------------------------

When analyzing empirical networks, one should be open to the possibility
that there will be more than one fit of the SBM with similar posterior
probabilities. In such situations, one should instead `sample`
partitions from the posterior distribution, instead of simply finding
its maximum. One can then compute quantities that are averaged over the
different model fits, weighted according to their posterior
probabilities.

Full support for model averaging is implemented in ``graph-tool`` via an
efficient `Markov chain Monte Carlo (MCMC)
<https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo>`_ algorithm
[peixoto-efficient-2014]_. It works by attempting to move nodes into
different groups with specific probabilities, and `accepting or
rejecting
<https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm>`_
such moves so that, after a sufficiently long time, the partitions will
be observed with the desired posterior probability. The algorithm is
designed so that its run-time (i.e. each sweep of the MCMC) is linear on
the number of edges in the network, and independent on the number of
groups being used in the model, and hence is suitable for use on very
large networks.

In order to perform such moves, one needs again to operate with
:class:`~graph_tool.inference.blockmodel.BlockState` or
:class:`~graph_tool.inference.nested_blockmodel.NestedBlockState` instances, and calling
their :meth:`~graph_tool.inference.blockmodel.BlockState.mcmc_sweep` methods. For
example, the following will perform 1000 sweeps of the algorithm with
the network of characters in the novel Les Misérables, starting from a
random partition into 20 groups

.. testcode:: model-averaging

   g = gt.collection.data["lesmis"]

   state = gt.BlockState(g, B=20)   # This automatically initializes the state
                                    # with a random partition into B=20
                                    # nonempty groups; The user could
                                    # also pass an arbitrary initial
                                    # partition using the 'b' parameter.

   # Now we run 1,000 sweeps of the MCMC. Note that the number of groups
   # is allowed to change, so it will eventually move from the initial
   # value of B=20 to whatever is most appropriate for the data.

   dS, nattempts, nmoves = state.mcmc_sweep(niter=1000)

   print("Change in description length:", dS)
   print("Number of accepted vertex moves:", nmoves)

.. testoutput:: model-averaging

   Change in description length: -365.317522...
   Number of accepted vertex moves: 38213

.. note::

   Starting from a random partition is rarely the best option, since it
   may take a long time for it to equilibrate. It was done above simply
   as an illustration on how to initialize
   :class:`~graph_tool.inference.blockmodel.BlockState` by hand. Instead, a much
   better option in practice is to start from an approximation to the
   "ground state" obtained with
   :func:`~graph_tool.inference.minimize.minimize_blockmodel_dl`, e.g.

    .. testcode:: model-averaging

       state = gt.minimize_blockmodel_dl(g)
       state = state.copy(B=g.num_vertices())
       dS, nattempts, nmoves = state.mcmc_sweep(niter=1000)

       print("Change in description length:", dS)
       print("Number of accepted vertex moves:", nmoves)

    .. testoutput:: model-averaging

       Change in description length: 1.660677...
       Number of accepted vertex moves: 40461

Although the above is sufficient to implement model averaging, there is a
convenience function called
:func:`~graph_tool.inference.mcmc.mcmc_equilibrate` that is intend to
simplify the detection of equilibration, by keeping track of the maximum
and minimum values of description length encountered and how many sweeps
have been made without a "record breaking" event. For example,

.. testcode:: model-averaging

   # We will accept equilibration if 10 sweeps are completed without a
   # record breaking event, 2 consecutive times.

   gt.mcmc_equilibrate(state, wait=10, nbreaks=2, mcmc_args=dict(niter=10), verbose=True)

will output:

.. testoutput:: model-averaging
    :options: +NORMALIZE_WHITESPACE

    niter:     1  count:    0  breaks:  0  min_S: 706.26857  max_S: 708.14483  S: 708.14483  ΔS:      1.87626  moves:   418 
    niter:     2  count:    0  breaks:  0  min_S: 699.23453  max_S: 708.14483  S: 699.23453  ΔS:     -8.91030  moves:   409 
    niter:     3  count:    0  breaks:  0  min_S: 699.23453  max_S: 715.33531  S: 715.33531  ΔS:      16.1008  moves:   414 
    niter:     4  count:    0  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 723.13301  ΔS:      7.79770  moves:   391 
    niter:     5  count:    1  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 702.93354  ΔS:     -20.1995  moves:   411 
    niter:     6  count:    2  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 706.39029  ΔS:      3.45675  moves:   389 
    niter:     7  count:    3  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 706.80859  ΔS:     0.418293  moves:   404 
    niter:     8  count:    4  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 707.61960  ΔS:     0.811010  moves:   417 
    niter:     9  count:    5  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 706.46577  ΔS:     -1.15383  moves:   392 
    niter:    10  count:    6  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 714.34671  ΔS:      7.88094  moves:   410 
    niter:    11  count:    7  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 706.43194  ΔS:     -7.91477  moves:   383 
    niter:    12  count:    8  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 705.19434  ΔS:     -1.23760  moves:   405 
    niter:    13  count:    9  breaks:  0  min_S: 699.23453  max_S: 723.13301  S: 702.21395  ΔS:     -2.98039  moves:   423 
    niter:    14  count:    0  breaks:  1  min_S: 715.54878  max_S: 715.54878  S: 715.54878  ΔS:      13.3348  moves:   400 
    niter:    15  count:    0  breaks:  1  min_S: 715.54878  max_S: 716.65842  S: 716.65842  ΔS:      1.10964  moves:   413 
    niter:    16  count:    0  breaks:  1  min_S: 701.19994  max_S: 716.65842  S: 701.19994  ΔS:     -15.4585  moves:   382 
    niter:    17  count:    1  breaks:  1  min_S: 701.19994  max_S: 716.65842  S: 715.56997  ΔS:      14.3700  moves:   394 
    niter:    18  count:    0  breaks:  1  min_S: 701.19994  max_S: 719.25577  S: 719.25577  ΔS:      3.68580  moves:   404 
    niter:    19  count:    0  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 723.78811  ΔS:      4.53233  moves:   413 
    niter:    20  count:    1  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 709.77340  ΔS:     -14.0147  moves:   387 
    niter:    21  count:    2  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 714.14891  ΔS:      4.37551  moves:   419 
    niter:    22  count:    3  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 722.05875  ΔS:      7.90984  moves:   399 
    niter:    23  count:    4  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 714.32503  ΔS:     -7.73371  moves:   422 
    niter:    24  count:    5  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 708.53927  ΔS:     -5.78576  moves:   392 
    niter:    25  count:    6  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 714.05889  ΔS:      5.51962  moves:   404 
    niter:    26  count:    7  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 713.93196  ΔS:    -0.126937  moves:   414 
    niter:    27  count:    8  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 709.49863  ΔS:     -4.43333  moves:   410 
    niter:    28  count:    9  breaks:  1  min_S: 701.19994  max_S: 723.78811  S: 707.42167  ΔS:     -2.07696  moves:   397 
    niter:    29  count:    0  breaks:  1  min_S: 699.89982  max_S: 723.78811  S: 699.89982  ΔS:     -7.52185  moves:   388 
    niter:    30  count:    0  breaks:  1  min_S: 698.57305  max_S: 723.78811  S: 698.57305  ΔS:     -1.32677  moves:   391 
    niter:    31  count:    1  breaks:  1  min_S: 698.57305  max_S: 723.78811  S: 706.02629  ΔS:      7.45324  moves:   412 
    niter:    32  count:    2  breaks:  1  min_S: 698.57305  max_S: 723.78811  S: 701.97778  ΔS:     -4.04852  moves:   421 
    niter:    33  count:    3  breaks:  1  min_S: 698.57305  max_S: 723.78811  S: 707.50134  ΔS:      5.52356  moves:   410 
    niter:    34  count:    4  breaks:  1  min_S: 698.57305  max_S: 723.78811  S: 708.56686  ΔS:      1.06552  moves:   424 
    niter:    35  count:    0  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 724.07361  ΔS:      15.5067  moves:   399 
    niter:    36  count:    1  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 723.51969  ΔS:    -0.553915  moves:   384 
    niter:    37  count:    2  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 702.36708  ΔS:     -21.1526  moves:   406 
    niter:    38  count:    3  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 707.60129  ΔS:      5.23420  moves:   405 
    niter:    39  count:    4  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 709.67542  ΔS:      2.07413  moves:   400 
    niter:    40  count:    5  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 714.52753  ΔS:      4.85212  moves:   398 
    niter:    41  count:    6  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 707.86563  ΔS:     -6.66190  moves:   409 
    niter:    42  count:    7  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 718.80926  ΔS:      10.9436  moves:   400 
    niter:    43  count:    8  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 716.37312  ΔS:     -2.43615  moves:   378 
    niter:    44  count:    9  breaks:  1  min_S: 698.57305  max_S: 724.07361  S: 713.76944  ΔS:     -2.60368  moves:   399 
    niter:    45  count:   10  breaks:  2  min_S: 698.57305  max_S: 724.07361  S: 715.29009  ΔS:      1.52066  moves:   421

Note that the value of ``wait`` above was made purposefully low so that
the output would not be overly long. The most appropriate value requires
experimentation, but a typically good value is ``wait=1000``.

The function :func:`~graph_tool.inference.mcmc.mcmc_equilibrate` accepts a
``callback`` argument that takes an optional function to be invoked
after each call to
:meth:`~graph_tool.inference.blockmodel.BlockState.mcmc_sweep`. This function
should accept a single parameter which will contain the actual
:class:`~graph_tool.inference.blockmodel.BlockState` instance. We will use this in
the example below to collect the posterior vertex marginals (via
:class:`~graph_tool.inference.blockmodel.BlockState.collect_vertex_marginals`),
i.e. the posterior probability that a node belongs to a given group:

.. testcode:: model-averaging

   # We will first equilibrate the Markov chain
   gt.mcmc_equilibrate(state, wait=1000, mcmc_args=dict(niter=10))

   pv = None 

   def collect_marginals(s):
      global pv
      pv = s.collect_vertex_marginals(pv)

   # Now we collect the marginals for exactly 100,000 sweeps, at
   # intervals of 10 sweeps:
   gt.mcmc_equilibrate(state, force_niter=10000, mcmc_args=dict(niter=10),
                       callback=collect_marginals)

   # Now the node marginals are stored in property map pv. We can
   # visualize them as pie charts on the nodes:
   state.draw(pos=g.vp.pos, vertex_shape="pie", vertex_pie_fractions=pv,
              edge_gradient=None, output="lesmis-sbm-marginals.svg")

.. figure:: lesmis-sbm-marginals.*
   :align: center
   :width: 450px

   Marginal probabilities of group memberships of the network of
   characters in the novel Les Misérables, according to the
   degree-corrected SBM. The `pie fractions
   <https://en.wikipedia.org/wiki/Pie_chart>`_ on the nodes correspond
   to the probability of being in group associated with the respective
   color.

We can also obtain a marginal probability on the number of groups
itself, as follows.

.. testcode:: model-averaging

   h = np.zeros(g.num_vertices() + 1)

   def collect_num_groups(s):
       B = s.get_nonempty_B()
       h[B] += 1

   # Now we collect the marginals for exactly 100,000 sweeps, at
   # intervals of 10 sweeps:
   gt.mcmc_equilibrate(state, force_niter=10000, mcmc_args=dict(niter=10),
                       callback=collect_num_groups)

.. testcode:: model-averaging
   :hide:

   figure()
   Bs = np.arange(len(h))
   idx = h > 0
   bar(Bs[idx], h[idx] / h.sum(), width=1, color="#ccb974")
   gca().set_xticks([6,7,8,9])
   xlabel("$B$")
   ylabel(r"$P(B|\boldsymbol A)$")
   savefig("lesmis-B-posterior.svg")

.. figure:: lesmis-B-posterior.*
   :align: center

   Marginal posterior probability of the number of nonempty groups for
   the network of characters in the novel Les Misérables, according to
   the degree-corrected SBM.


Hierarchical partitions
+++++++++++++++++++++++

We can also perform model averaging using the nested SBM, which will
give us a distribution over hierarchies. The whole procedure is fairly
analogous, but now we make use of
:class:`~graph_tool.inference.nested_blockmodel.NestedBlockState` instances.

.. note::

    When using :class:`~graph_tool.inference.nested_blockmodel.NestedBlockState` instances
    to perform model averaging, they need to be constructed with the
    option ``sampling=True``.

Here we perform the sampling of hierarchical partitions using the same
network as above.

.. testcode:: nested-model-averaging

   g = gt.collection.data["lesmis"]

   state = gt.minimize_nested_blockmodel_dl(g) # Initialize he Markov
                                               # chain from the "ground
                                               # state"

   # Before doing model averaging, the need to create a NestedBlockState
   # by passing sampling = True.

   # We also want to increase the maximum hierarchy depth to L = 10

   # We can do both of the above by copying.

   bs = state.get_bs()                     # Get hierarchical partition.
   bs += [np.zeros(1)] * (10 - len(bs))    # Augment it to L = 10 with
                                           # single-group levels.

   state = state.copy(bs=bs, sampling=True)

   # Now we run 1000 sweeps of the MCMC

   dS, nattempts, nmoves = state.mcmc_sweep(niter=1000)

   print("Change in description length:", dS)
   print("Number of accepted vertex moves:", nmoves)

.. testoutput:: nested-model-averaging

   Change in description length: 2.371018...
   Number of accepted vertex moves: 56087

Similarly to the the non-nested case, we can use
:func:`~graph_tool.inference.mcmc.mcmc_equilibrate` to do most of the boring
work, and we can now obtain vertex marginals on all hierarchical levels:


.. testcode:: nested-model-averaging

   # We will first equilibrate the Markov chain
   gt.mcmc_equilibrate(state, wait=1000, mcmc_args=dict(niter=10))

   pv = [None] * len(state.get_levels())

   def collect_marginals(s):
      global pv
      pv = [sl.collect_vertex_marginals(pv[l]) for l, sl in enumerate(s.get_levels())]

   # Now we collect the marginals for exactly 100,000 sweeps
   gt.mcmc_equilibrate(state, force_niter=10000, mcmc_args=dict(niter=10),
                       callback=collect_marginals)

   # Now the node marginals for all levels are stored in property map
   # list pv. We can visualize the first level as pie charts on the nodes:
   state_0 = state.get_levels()[0]
   state_0.draw(pos=g.vp.pos, vertex_shape="pie", vertex_pie_fractions=pv[0],
                edge_gradient=None, output="lesmis-nested-sbm-marginals.svg")

.. figure:: lesmis-nested-sbm-marginals.*
   :align: center
   :width: 450px

   Marginal probabilities of group memberships of the network of
   characters in the novel Les Misérables, according to the nested
   degree-corrected SBM. The pie fractions on the nodes correspond to
   the probability of being in group associated with the respective
   color.

We can also obtain a marginal probability of the number of groups
itself, as follows.

.. testcode:: nested-model-averaging

   h = [np.zeros(g.num_vertices() + 1) for s in state.get_levels()]

   def collect_num_groups(s):
       for l, sl in enumerate(s.get_levels()):
          B = sl.get_nonempty_B()
          h[l][B] += 1

   # Now we collect the marginal distribution for exactly 100,000 sweeps
   gt.mcmc_equilibrate(state, force_niter=10000, mcmc_args=dict(niter=10),
                       callback=collect_num_groups)

.. testcode:: nested-model-averaging
   :hide:

   figure()
   f, ax = plt.subplots(1, 5, figsize=(10, 3))
   for i, h_ in enumerate(h[:5]):
       Bs = np.arange(len(h_))
       idx = h_ > 0
       ax[i].bar(Bs[idx], h_[idx] / h_.sum(), width=1, color="#ccb974")
       ax[i].set_xticks(Bs[idx])
       ax[i].set_xlabel("$B_{%d}$" % i)
       ax[i].set_ylabel(r"$P(B_{%d}|\boldsymbol A)$" % i)
       locator = MaxNLocator(prune='both', nbins=5)
       ax[i].yaxis.set_major_locator(locator)
   tight_layout()
   savefig("lesmis-nested-B-posterior.svg")

.. figure:: lesmis-nested-B-posterior.*
   :align: center

   Marginal posterior probability of the number of nonempty groups
   :math:`B_l` at each hierarchy level :math:`l` for the network of
   characters in the novel Les Misérables, according to the nested
   degree-corrected SBM.

Below we obtain some hierarchical partitions sampled from the posterior
distribution.

.. testcode:: nested-model-averaging

   for i in range(10):
       state.mcmc_sweep(niter=1000)
       state.draw(output="lesmis-partition-sample-%i.svg" % i, empty_branches=False)

.. image:: lesmis-partition-sample-0.svg
   :width: 200px
.. image:: lesmis-partition-sample-1.svg
   :width: 200px
.. image:: lesmis-partition-sample-2.svg
   :width: 200px
.. image:: lesmis-partition-sample-3.svg
   :width: 200px
.. image:: lesmis-partition-sample-4.svg
   :width: 200px
.. image:: lesmis-partition-sample-5.svg
   :width: 200px
.. image:: lesmis-partition-sample-6.svg
   :width: 200px
.. image:: lesmis-partition-sample-7.svg
   :width: 200px
.. image:: lesmis-partition-sample-8.svg
   :width: 200px
.. image:: lesmis-partition-sample-9.svg
   :width: 200px
