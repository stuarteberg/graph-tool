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

   Change in description length: -353.848032...
   Number of accepted vertex moves: 37490

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

       Change in description length: 31.622518...
       Number of accepted vertex moves: 43152

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

    niter:     1  count:    0  breaks:  0  min_S: 703.94152  max_S: 730.97213  S: 703.94152  ΔS:     -27.0306  moves:   431 
    niter:     2  count:    1  breaks:  0  min_S: 703.94152  max_S: 730.97213  S: 708.61840  ΔS:      4.67688  moves:   413 
    niter:     3  count:    2  breaks:  0  min_S: 703.94152  max_S: 730.97213  S: 704.60994  ΔS:     -4.00847  moves:   416 
    niter:     4  count:    0  breaks:  0  min_S: 700.85336  max_S: 730.97213  S: 700.85336  ΔS:     -3.75658  moves:   391 
    niter:     5  count:    1  breaks:  0  min_S: 700.85336  max_S: 730.97213  S: 713.22553  ΔS:      12.3722  moves:   387 
    niter:     6  count:    2  breaks:  0  min_S: 700.85336  max_S: 730.97213  S: 703.57357  ΔS:     -9.65196  moves:   434 
    niter:     7  count:    3  breaks:  0  min_S: 700.85336  max_S: 730.97213  S: 715.02440  ΔS:      11.4508  moves:   439 
    niter:     8  count:    0  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 700.68857  ΔS:     -14.3358  moves:   427 
    niter:     9  count:    1  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 717.95725  ΔS:      17.2687  moves:   409 
    niter:    10  count:    2  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 720.02079  ΔS:      2.06354  moves:   435 
    niter:    11  count:    3  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 718.15880  ΔS:     -1.86199  moves:   399 
    niter:    12  count:    4  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 708.06732  ΔS:     -10.0915  moves:   436 
    niter:    13  count:    5  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 712.76007  ΔS:      4.69274  moves:   432 
    niter:    14  count:    6  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 705.60582  ΔS:     -7.15425  moves:   409 
    niter:    15  count:    7  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 704.37333  ΔS:     -1.23249  moves:   434 
    niter:    16  count:    8  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 717.54492  ΔS:      13.1716  moves:   426 
    niter:    17  count:    9  breaks:  0  min_S: 700.68857  max_S: 730.97213  S: 715.05767  ΔS:     -2.48725  moves:   449 
    niter:    18  count:    0  breaks:  1  min_S: 715.77940  max_S: 715.77940  S: 715.77940  ΔS:     0.721731  moves:   448 
    niter:    19  count:    0  breaks:  1  min_S: 708.38072  max_S: 715.77940  S: 708.38072  ΔS:     -7.39868  moves:   447 
    niter:    20  count:    0  breaks:  1  min_S: 705.63447  max_S: 715.77940  S: 705.63447  ΔS:     -2.74625  moves:   441 
    niter:    21  count:    1  breaks:  1  min_S: 705.63447  max_S: 715.77940  S: 707.01766  ΔS:      1.38319  moves:   434 
    niter:    22  count:    2  breaks:  1  min_S: 705.63447  max_S: 715.77940  S: 708.21127  ΔS:      1.19361  moves:   447 
    niter:    23  count:    0  breaks:  1  min_S: 703.12325  max_S: 715.77940  S: 703.12325  ΔS:     -5.08802  moves:   454 
    niter:    24  count:    0  breaks:  1  min_S: 703.05106  max_S: 715.77940  S: 703.05106  ΔS:   -0.0721911  moves:   433 
    niter:    25  count:    1  breaks:  1  min_S: 703.05106  max_S: 715.77940  S: 704.77370  ΔS:      1.72264  moves:   423 
    niter:    26  count:    0  breaks:  1  min_S: 701.61368  max_S: 715.77940  S: 701.61368  ΔS:     -3.16003  moves:   441 
    niter:    27  count:    0  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 721.54373  ΔS:      19.9301  moves:   434 
    niter:    28  count:    1  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 703.33612  ΔS:     -18.2076  moves:   439 
    niter:    29  count:    2  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 710.79425  ΔS:      7.45813  moves:   437 
    niter:    30  count:    3  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 706.35044  ΔS:     -4.44381  moves:   429 
    niter:    31  count:    4  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 713.56014  ΔS:      7.20970  moves:   463 
    niter:    32  count:    5  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 720.16436  ΔS:      6.60422  moves:   445 
    niter:    33  count:    6  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 714.76845  ΔS:     -5.39591  moves:   404 
    niter:    34  count:    7  breaks:  1  min_S: 701.61368  max_S: 721.54373  S: 703.21572  ΔS:     -11.5527  moves:   410 
    niter:    35  count:    0  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 701.53898  ΔS:     -1.67675  moves:   434 
    niter:    36  count:    1  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 708.14043  ΔS:      6.60146  moves:   433 
    niter:    37  count:    2  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 704.07209  ΔS:     -4.06835  moves:   410 
    niter:    38  count:    3  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 704.76811  ΔS:     0.696023  moves:   413 
    niter:    39  count:    4  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 703.54823  ΔS:     -1.21988  moves:   398 
    niter:    40  count:    5  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 713.59891  ΔS:      10.0507  moves:   388 
    niter:    41  count:    6  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 704.40168  ΔS:     -9.19724  moves:   403 
    niter:    42  count:    7  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 707.57723  ΔS:      3.17556  moves:   400 
    niter:    43  count:    8  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 704.09679  ΔS:     -3.48044  moves:   423 
    niter:    44  count:    9  breaks:  1  min_S: 701.53898  max_S: 721.54373  S: 704.64514  ΔS:     0.548354  moves:   419 
    niter:    45  count:   10  breaks:  2  min_S: 701.53898  max_S: 721.54373  S: 715.92329  ΔS:      11.2781  moves:   411 

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

   Change in description length: 15.483135...
   Number of accepted vertex moves: 57684

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
