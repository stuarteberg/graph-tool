.. _reconstruction_dynamics:

Reconstruction from dynamics
++++++++++++++++++++++++++++


.. testsetup:: dynamics

   import os
   try:
      os.chdir("demos/inference")
   except FileNotFoundError:
       pass
   np.random.seed(42)
   gt.seed_rng(44)

.. testcode:: dynamics

   g = gt.collection.data["dolphins"]
  
   ss = []
   for i in range(1000):
       si_state = gt.SIState(g, beta=1)
       s = []
       for j in range(10):
           si_state.iterate_sync()
           s.append(si_state.get_state().copy())
       s = gt.group_vector_property(s)
       ss.append(s)

   u = g.copy()
   u.clear_edges()
   ss = [u.own_property(s) for s in ss]

   rstate = gt.EpidemicsBlockState(u, s=ss, beta=None, r=1e-6, global_beta=.99, 
                                   state_args=dict(B=1), nested=False)

   # Now we collect the marginals for exactly 100,000 sweeps, at
   # intervals of 10 sweeps:

   gm = None
   bm = None

   def collect_marginals(s):
      global gm, bm
      gm = s.collect_marginal(gm)
      b = gt.perfect_prop_hash([s.bstate.b])[0]
      bm = s.bstate.collect_vertex_marginals(bm, b=b)

   gt.mcmc_equilibrate(rstate, force_niter=10000, mcmc_args=dict(niter=10, xstep=0, p=0, h=0),
                       callback=collect_marginals)

   b = bm.new_vp("int", vals=[bm[v].a.argmax() for v in bm.vertices()])

   graph_draw(gm, gm.own_property(g.vp.pos), vertex_shape="square", vertex_color="black",
              vertex_fill_color=b, vertex_pen_width=1,
              edge_pen_width=prop_to_size(gm.ep.eprob, 0, 5), eorder=gm.ep.eprob, output="dolphins")
                       
   eprob = u.ep.eprob
   print("Posterior probability of edge (11, 36):", eprob[u.edge(11, 36)])
   print("Posterior probability of non-edge (15, 73):", eprob[u.edge(15, 73)])
   print("Estimated average local clustering: %g ± %g" % (np.mean(cs), np.std(cs)))
