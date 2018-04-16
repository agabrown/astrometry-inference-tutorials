# Parallax tutorial 2018

Coryn Bailer-Jones, MPIA Heidelberg (https://mpia.de/homes/calj)

This tutorial concerns inferring distances and velocities from parallaxes and proper motions using the Bayesian approach. It uses R codes in jupyter notebooks, together with simulated data and data from Gaia-DR1 (TGAS). The tutorial is divided into three parts, each with its own notebook:

1. Inference of distance to a single source using its parallax. Includes: simple hierarchical model. 
2. Inference of distance to and size of a cluster using parallax and positions of its members. Includes: naive parallax combination; accommodating correlated measurements; 2D cluster model. See [resources/cluster_inference.pdf](resources/cluster_inference.pdf) for details.
3. Inference of distance to and 2D tangential velocity on the sky of a single source using its parallax and proper motion. Includes: explicit use of MCMC to sample the posterior. See [resources/3D_astrometry_inference.pdf](resources/3D_astrometry_inference.pdf) for details. 

The more generic functions used the tutorials are in the files in the Rcode/ directory.

These tutorials assume that you are familiar with the basic idea of Bayesian inference, and inferring a distance given a parallax and a prior, as described in [Bailer-Jones 2015 (paper 1)](http://adsabs.harvard.edu/abs/2015PASP..127..994B). Additional resources and references are provided below. You should do tutorial 1 before tackling 2 or 3. 

The purpose of these tutorials is just to show how to work with astrometric data in inference problems. Sometimes you will have other relevant information (e.g. colour and apparent magnitude) which can also be used to help infer the distance. And sometimes you won't want to infer a distance at all (e.g. as an intermediate measure; or if you are doing model fitting, which is probably better done in the parallax space where the measurement model is simpler).

Resources and further information:
* [Cluster distance inference](resources/cluster_inference.pdf)
* [Joint inference from parallax and proper motions (CBJ-081)](resources/3D_astrometry_inference.pdf)
* [Bailer-Jones 2015 (paper 1)](http://adsabs.harvard.edu/abs/2015PASP..127..994B)
* [Astraatmadja & Bailer-Jones 2016a (paper 2)](http://adsabs.harvard.edu/abs/2016ApJ...832..137A)
* [Astraatmadja & Bailer-Jones 2016b (paper 3)](http://adsabs.harvard.edu/abs/2016ApJ...833..119A)
* [Gaia Data Release 1](http://adsabs.harvard.edu/abs/2016A%26A...595A...2G)
* [Luri et al. 2018 (tutorial on the use of parallaxes)]()
* [Practical Bayesian Inference (PBI)](http://www.cambridge.org/de/academic/subjects/physics/mathematical-methods/practical-bayesian-inference-primer-physical-scientists?format=PB)

