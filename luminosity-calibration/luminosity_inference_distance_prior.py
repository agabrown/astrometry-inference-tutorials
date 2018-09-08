"""
Command line version of the luminosity inference problem discussed in the notebook
Luminosity_Inference_DistPrior.ipynb, part of the tutorials accompanying the paper by Luri et al. (2018)
on the use of Gaia astrometry.

Anthony G.A. Brown Nov 2017 - Sep 2018
<brown@strw.leidenuniv.nl>
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse, json
from scipy.interpolate import interp1d

import pystan
import corner

from stantools import load_stan_code, stan_cache
from parallaxsurveys import UniformDistributionSingleLuminosityHip as udslH
from parallaxsurveys import UniformDistributionSingleLuminosityTGAS as udslT
from parallaxsurveys import showSurveyStatistics

# The following line ensures that this script can also be run remotely, without a display available.
plt.switch_backend('agg')

def run_luminosity_inference(args):
    """
    Generate the simulated parallax survey and then perform the Bayesian inference to estimate the mean
    absolute magnitude of the stars and the variance thereof. 

    Parameters
    ----------

    args - Command line arguments
    """

    distMin = args['distMin']
    distMax = args['distMax']
    plxMin = 1000.0/distMax
    plxMax = 1000.0/distMin
    nstars = args['nstars']
    absMagTrue = args['muM']
    sigmaAbsMagTrue = args['sigmaM']
    magLimit = args['mlim']
    stancontrol = json.loads(args['stancontrol'])

    numChains = 4
    
    print("### Generating parallax survey ... ###")
    if args['cat']=='hip':
        survey = udslH(nstars, distMin, distMax, absMagTrue, sigmaAbsMagTrue, surveyLimit = magLimit)
    else:
        survey = udslT(nstars, distMin, distMax, absMagTrue, sigmaAbsMagTrue, surveyLimit = magLimit)
    survey.setRandomNumberSeed(args['surveyseed'])
    survey.generateObservations()
    print("### ... Done ###")

    showSurveyStatistics(survey, pdfFile="surveyStats.pdf", usekde=False)

    if args['volumecomplete'] or np.isinf(magLimit):
        # Volume limited case
        stan_data = {'minDist':distMin, 'maxDist':distMax, 'N':survey.numberOfStarsInSurvey,
                'obsPlx':survey.observedParallaxes, 'errObsPlx':survey.parallaxErrors,
                'obsMag':survey.observedMagnitudes, 'errObsMag':survey.magnitudeErrors}
        stanmodel = load_stan_code("luminosity_inference_volume_complete_distance_prior.stan")
        sm = stan_cache(stanmodel, model_name="luminosityInferenceDistPriorVC")
        fit = sm.sampling(data = stan_data, pars=['meanAbsMag', 'sigmaAbsMag'], iter=args['staniter'],
                chains=numChains, thin=args['stanthin'], seed=args['stanseed'], control=stancontrol)
    else:
        # Magnitude limited case. Explicit initialization of the true absolute magnitudes is
        # needed. See comments in Stan code.
        maxPossibleAbsMag = survey.apparentMagnitudeLimit-5*np.log10(distMax)+5
        initLow = maxPossibleAbsMag - 4
        initialValuesList = []
        for i in range(numChains):
            initialValuesList.append( dict(absMag=np.random.uniform(initLow, maxPossibleAbsMag,
                size=survey.numberOfStarsInSurvey)) )

        stan_data = {'minDist':distMin, 'maxDist':distMax, 'surveyLimit':survey.apparentMagnitudeLimit,
                'N':survey.numberOfStarsInSurvey, 'obsPlx':survey.observedParallaxes,
                'errObsPlx':survey.parallaxErrors, 'obsMag':survey.observedMagnitudes,
                'errObsMag':survey.magnitudeErrors}
        stanmodel = load_stan_code("luminosity_inference_distance_prior.stan")
        sm = stan_cache(stanmodel, model_name="luminosityInferenceDistPrior")
        fit = sm.sampling(data = stan_data, pars=['meanAbsMag', 'sigmaAbsMag'], iter=args['staniter'],
                chains=numChains, thin=args['stanthin'], seed=args['stanseed'], init=initialValuesList,
                control=stancontrol)
    
    with open('StanSummary.txt','w') as f:
        print(fit.stansummary(), file=f)

    samples = np.vstack([fit.extract()['meanAbsMag'], fit.extract()['sigmaAbsMag']]).transpose()
    
    fig = plt.figure(figsize=(8,8))
    for i in range(1,5):
        fig.add_subplot(2,2,i)
    corner.corner(samples, labels=[r'$\mu_M$', r'$\sigma_M$'], truths=[absMagTrue, sigmaAbsMagTrue],
            truth_color='r', quantiles=[0.16,0.50,0.84], show_titles=True, fig=fig)
    plt.savefig('cornerplot.pdf')

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    defaultDistMin = 1.0
    defaultDistMax = 100.0
    defaultNstars = 50
    defaultMeanAbsMag = 9.0
    defaultSigmaAbsMag = 0.7
    defaultMagLim = np.inf
    defaultCatalogue = 'hip'
    defaultStanIter = 10000
    defaultStanThin = 5
    parser = argparse.ArgumentParser(description="""Run luminosity inference tutorial for distance priors""")
    parser.add_argument("--distMin", help="""Minimum value of distance distribution (default
            {0} pc)""".format(defaultDistMin), default=defaultDistMin, type=float)
    parser.add_argument("--distMax", help="""Maximum value of distance distribution (default
            {0} pc)""".format(defaultDistMax), default=defaultDistMax, type=float)
    parser.add_argument("--nstars", help="""Number of stars in simulated survey (default
            {0})""".format(defaultNstars), default=defaultNstars, type=int)
    parser.add_argument("--muM", help="""Mean true absolute magnitude (default
            {0})""".format(defaultMeanAbsMag), default=defaultMeanAbsMag, type=float)
    parser.add_argument("--sigmaM", help="""Standard deviation true absolute magnitude distribution
            (default {0})""".format(defaultSigmaAbsMag), default=defaultSigmaAbsMag, type=float)
    parser.add_argument("--mlim", help="""Survey limiting magnitude (default
            {0})""".format(defaultMagLim), default=defaultMagLim, type=float)
    parser.add_argument("--cat", help="""Simulated astrometric catalogue (default {0})""".format(defaultCatalogue),
            choices=['hip','tgas'], default=defaultCatalogue, type=str)
    parser.add_argument("--volumecomplete", action="store_true", help="""Use model for volume complete survey""")
    parser.add_argument("--surveyseed", help="""Random number seed for survey simulation (default None)""", type=int, default=None)
    parser.add_argument("--stanseed", help="""Random number seed for Stan MCMC sampler (default None)""", type=int, default=None)
    parser.add_argument("--staniter", help="""Number of iterations per chain for Stan MCMC sampler
            (default {0})""".format(defaultStanIter), type=int, default=defaultStanIter)
    parser.add_argument("--stanthin", help="""Thinning parameter for Stan MCMC sampler (default {0})"""
            .format(defaultStanThin), type=int, default=defaultStanThin)
    parser.add_argument("--stancontrol", help="""Dictionary of parameters to control the sampler’s
            behavior, formulated as JSON string (example: '{"max_treedepth":10}')""", type=str,
            default='{}')
    args=vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    run_luminosity_inference(args)
