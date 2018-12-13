# astrometry-inference-tutorials
This is the repository for the tutorials that accompany the Gaia DR2 paper [Luri et
al. 2018, A&A Special Issue for Gaia DR2](https://doi.org/10.1051/0004-6361/201832964) (also
available on [arxiv.org](https://arxiv.org/abs/1804.09376)) describing recommended practices for the
use of astrometric data (in particular parallaxes) in astronomical data analysis problems.

You can launch the tutorials through [binder](https://mybinder.org). This was installed on 2018.12.12 and is not yet working
properly.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/agabrown/astrometry-inference-tutorials/v1.0)

## Getting things to work

Installation instructions to make sure you have everything needed to run the tutorials are [here](INSTALL.md).

## Tutorials

* [Infer distance to a single source](./single-source) (R)
  * [Bayesian inference of distance from parallax only for a single source](single-source/tutorial/Distance_inference-single_source.ipynb)
  * [Comparison of distance estimators](/single-source/GraphicalUserInterface/Tutorial.ipynb)

* [Infer distance to and size of a cluster](./multiple-source) (R)
  * [Notebook](./multiple-source/Distance_inference-multiple_sources.ipynb)

* [Infer distance and tangential velocity of a source](./3d-distance) (R)
  * [Notebook](./3d-distance/Distance_and_tangential_velocity_inference.ipynb)

* [Luminosity calibration](./luminosity-calibration) (Python)
  * [Explanation of negative parallaxes](./luminosity-calibration/DemoNegativeParallax.ipynb)
  * [Distribution of quantities derived from parallaxes for naive estimators and posteriors based on
    minimal prior information](./luminosity-calibration/Parallax_related_quantities.ipynb)
  * [Basics of handling data truncation](./luminosity-calibration/Handling_Data_Truncation.ipynb)
  * [Simulation of parallax surveys](./luminosity-calibration/Parallax_survey_simulation.ipynb)
  * [Inference of the luminosity of a class of
    stars](./luminosity-calibration/Luminosity_Inference_DistPrior.ipynb)

* [The period luminosity relation](./period-luminosity-relation) (R and Python)
  * [R notebook interfaced to python](./period-luminosity-relation/TutorialPLZ-rp2.ipynb)
  * [R notebook](./period-luminosity-relation/TutorialPLZ_R.ipynb)
