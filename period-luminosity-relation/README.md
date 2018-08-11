# The Period Luminosity relation

## Description

In this tutorial, we first briefly introduce the Bayesian paradigm for statistical inference and the
concept of hierarchical Bayesian graphical models. We simulate toy problems with R and illustrate
how we can easily represent them and infer their parameters of interest by means of the Stan
programming language and its interface Rstan for R. The second part of the tutorial is devoted to
using the Bayesian hierarchical methodology to derive a set of Period-Luminosity-Metallicity
relations from the Gaia DR1 parallaxes (TGAS) and photometry in several bands. We study the
stability of the solutions and the sensitivity to the prior and hyperprior choices. We reach the
conclusion that there are hidden correlations in the data that prevent the simplified model from
reaching good estimates to the PLZ relations. The full details of the analysis and the solutions to
the problems highlighted in this tutorial can be found in [Delgado et al.
(2018)](http://adsabs.harvard.edu/abs/2018arXiv180301162D). This tutorial accompanies the Gaia DR 2
paper on the use of parallaxes [(Luri et al. 2018, A&A Special Issue for Gaia
DR2)](https://arxiv.org/abs/1804.09376).

## R packages 

The tutorial is written in R and its execution (but not its visualization) requires installation of
the following packages:

* rstan
* ggplo2

For the inline generation of the directed acyclic graphs (DAGs) the following are needed. The
installation of these can be skipped if the [modified version](TutorialPLZ_R.ipynb) of the tutorial
is used.

* DiagrammeR
* DiagrammeRsvg
* magrittr
* svglite
* rsvg
* png

Some of these R packages may require the local installation of libraries outside the R ecosystem.
