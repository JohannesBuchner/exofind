ExoFind
==========================================

A radial velocity fitting code based on ExoFit.

This is incomplete work, but the ideas were:

* Make ExoFit callable from Python. Uses fast C code underneath.
* Connect ExoFit with MultiNest instead of MCMC
* Use innovative data-based prior transformations. The FFT of the residuals
  can indicate periods to prefer when searching the parameter space.
  This biasing is cancelled out again in the likelihood (to make Bayesian
  inference consistent), but should make exploring the parameter space simpler.
* Use multi-dimensional exploration to detect planets, expanding and contracting
  the search space as needed. Achieved using differential evolution algorithm
  connected with the Akaike Information Criterion.
* Plotting RV curves.
  
No guarantees.
