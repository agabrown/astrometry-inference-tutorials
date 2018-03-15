"""
Provide classes and methods for simulating simple (magnitude limited) parallax surveys. These consist of
measurements of the parallax and the apparent magnitude of the stars.

Anthony Brown 2011-2018
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc, cm
from sys import stderr
from scipy import isscalar
from scipy.stats import norm, uniform
from scipy.special import erf
from scipy.integrate import quad
from cycler import cycler
from os import path

from plotstyles import useagab, apply_tufte
from distinct_colours import get_distinct
from robuststats import rse
from sklearn.neighbors import KernelDensity

def simDistancesConstantSpaceDensity(numStars, minDistance, maxDistance):
    """
    Simulate distances for stars distributed uniformly in space around the Sun.

    Parameters
    ----------
  
    numStars - number of stars to simulate
    minDistance - lower limit on the distance (pc, closest possible star)
    maxDistance - upper limit on the distance (pc, volume limit of survey)
  
    Returns
    -------

    Vector of distance values
    """
    x=uniform.rvs(loc=0.0, scale=1.0, size=numStars)
    minDCubed=np.power(minDistance,3.0)
    maxDCubed=np.power(maxDistance,3.0)
    return np.power(minDCubed+x*(maxDCubed-minDCubed),1.0/3.0)

def simGaussianAbsoluteMagnitude(numStars, mean, stddev):
    """
    Simulate absolute magnitudes following a Gaussian distribution.

    Parameters
    ----------

    numStars - number of stars to simulate
    mean - mean of the distribution
    stddev - standard deviation of the distribution

    Returns
    -------

    Vector of magnitudes.
    """
    return norm.rvs(loc=mean, scale=stddev, size=numStars)

class ParallaxSurvey:
    """
    Base class for simulating a parallax survey. The survey is assumed to have been conducted for a set
    of stars distributed in space between some minimum aand maximum distance value, according to some
    spatial distribution to be defined by the sub-classes. Likewise the distribution in apparent
    magnitudes is assumed to be specified by the sub-classes.
    """

    def __init__(self, numberOfStars, minDistance, maxDistance, surveyLimit=np.Inf):
        """
        Class constructor/initializer.

        Parameters
        ----------

        numberOfStars - Number of stars to simulate
        minDistance   - Lower limit on the distance (closest possible star, pc)
        maxDistance   - Upper limit on the distance (volume limit of survey, pc)

        Keywords
        --------

        surveyLimit - Apparent magnitude limit of the survey (default: no limit)
        """
        self.numberOfStars = numberOfStars
        self.numberOfStarsInSurvey = numberOfStars
        self.minDistance = minDistance
        self.maxDistance = maxDistance
        self.minParallax = 1000.0/self.maxDistance
        self.maxParallax = 1000.0/self.minDistance
        self.apparentMagnitudeLimit = surveyLimit
        self.seed = None

    def setRandomNumberSeed(self, seed):
        """
        (Re-)Set the random number seed for the simulations. NOTE, also applies to the scipy.stats functions.
  
        Parameters
        ----------
  
        seed - Value of random number seed
        """
        self.seed = seed
        np.random.seed(seed)

    def getRandomNumberSeed(self):
        """
        Get the random number seed that was used for the Monte Carlo simulations.

        Returns
        -------

        seed : int
            The random number seed used. None means the seed was set automatically.
        """
        return self.seed

    def _applyApparentMagnitudeLimit(self):
        """
        Apply the apparent magnitude limit to the simulated survey. The limit is applied to the observed
        apparent magnitude that is actually reported in the survey. This is only an approximation to what
        happens in practice, where a magnitude limit usually arises at the detection stage of the survey.
        The detection errors are normally larger than the catalogued magnitude errors.
        
        Alternatively this simulation can be seen as applying a selection to the catalogued magnitudes.
        """
        indices=(self.observedMagnitudes <= self.apparentMagnitudeLimit)
        self.trueParallaxes=self.trueParallaxes[indices].flatten()
        self.absoluteMagnitudes=self.absoluteMagnitudes[indices].flatten()
        self.apparentMagnitudes=self.apparentMagnitudes[indices].flatten()
        self.parallaxErrors=self.parallaxErrors[indices].flatten()
        self.magnitudeErrors=self.magnitudeErrors[indices].flatten()
        self.observedParallaxes=self.observedParallaxes[indices].flatten()
        self.observedMagnitudes=self.observedMagnitudes[indices].flatten()
        self.numberOfStarsInSurvey=len(self.observedMagnitudes)

class UniformSpaceDistributionSingleLuminosity(ParallaxSurvey):
    """
    Base class for simulated parallax surveys in which the stars are distributed uniformly in space
    between a minimum and maximum distance and with all stars having an absolute magnitude draw from the
    same Normal distribution.
    """
    def __init__(self, numberOfStars, minDistance, maxDistance, meanAbsoluteMagnitude,
            stddevAbsoluteMagnitude, surveyLimit=np.Inf):
        """
        Class constructor/initializer.

        Parameters
        ----------
    
        numberOfStars             - number of stars to simulate
        minDistance               - lower limit on the distance (closest possible star, pc)
        maxDistance               - upper limit on the distance (volume limit of survey, pc)
        meanAbsoluteMagnitude     - Mean of Gaussian absolute magnitude distribution
        stddevAbsoluteMagnitude   - Standard deviation of Gaussian absolute magnitude distribution
  
        Keywords
        --------

        surveyLimit - Apparent magnitude limit of the survey (default: no limit)
        """
        super().__init__(numberOfStars, minDistance, maxDistance, surveyLimit)
        self.meanAbsoluteMagnitude=meanAbsoluteMagnitude
        self.stddevAbsoluteMagnitude=stddevAbsoluteMagnitude

    def generateObservations(self):
        """
        Generate the simulated observations.
        """
        self.trueParallaxes = 1000.0/simDistancesConstantSpaceDensity(self.numberOfStars,
                self.minDistance, self.maxDistance)
        self.absoluteMagnitudes = simGaussianAbsoluteMagnitude(self.numberOfStars,
                self.meanAbsoluteMagnitude, self.stddevAbsoluteMagnitude)
        self.apparentMagnitudes = self.absoluteMagnitudes-5.0*np.log10(self.trueParallaxes)+10.0
        self.parallaxErrors = self._generateParallaxErrors()
        self.magnitudeErrors = self._generateApparentMagnitudeErrors()
        self.observedParallaxes = norm.rvs(loc=self.trueParallaxes, scale=self.parallaxErrors)
        self.observedMagnitudes = norm.rvs(loc=self.apparentMagnitudes, scale=self.magnitudeErrors)
        self._applyApparentMagnitudeLimit()
        self._appMagPdfNorm = 1.0
        #
        # Renormalize if there is a finite apparent magnitude limit.
        #
        if not np.isinf(self.apparentMagnitudeLimit):
            low = self.meanAbsoluteMagnitude - 5*self.stddevAbsoluteMagnitude + 5*np.log10(self.minDistance) - 5.0
            up = np.min([self.meanAbsoluteMagnitude + 5*self.stddevAbsoluteMagnitude +
                5*np.log10(self.maxDistance) - 5.0, self.apparentMagnitudeLimit])
            self._appMagPdfNorm, dummy = quad(self._apparentMagnitude_pdf, low, up)

    def apparentMagnitude_lpdf(self, m):
        """
        Calculate the natural logarithm of the analytical probability density function of the apparent
        magnitudes in the simulated survey. This can be derived from the know parallax and absolute
        magnitude PDFs.

        Parameters
        ----------

        m - The apparent magnitude(s) for which to calculate the PDF.
        """
        c = 0.6*np.log(10.0)
        a = np.power(self.maxDistance, 3.0) - np.power(self.minDistance, 3.0)
        dmH = 5*np.log10(self.maxDistance)-5
        dmL = 5*np.log10(self.minDistance)-5
        zH = (dmH - m + self.meanAbsoluteMagnitude - c*self.stddevAbsoluteMagnitude**2) / \
                (np.sqrt(2.0)*self.stddevAbsoluteMagnitude)
        zL = (dmL - m + self.meanAbsoluteMagnitude - c*self.stddevAbsoluteMagnitude**2) / \
                (np.sqrt(2.0)*self.stddevAbsoluteMagnitude)
        lpdf = np.log(c) - np.log(a) + c*(m - self.meanAbsoluteMagnitude + \
                0.5*c*self.stddevAbsoluteMagnitude**2 + 5.0) + \
                -np.log(2.0) + np.log(erf(zH)-erf(zL))

        return lpdf-np.log(self._appMagPdfNorm)

    def _apparentMagnitude_pdf(self, m):
        """
        Calculate value of the PDF of the apparent magnitudes in the simulated survey.
        """
        return np.exp(self.apparentMagnitude_lpdf(m))

class UniformDistributionSingleLuminosityHip(UniformSpaceDistributionSingleLuminosity):
    """
    Simulate a parallax survey for stars distributed uniformly in space around the sun. The stars all
    have the same luminosity drawn from a Gaussian distribution. The errors on the observed parallaxes
    and apparent magnitudes roughly follow the characteristics of the Hipparcos Catalogue.
    """

    def __init__(self, numberOfStars, minDistance, maxDistance, meanAbsoluteMagnitude,
            stddevAbsoluteMagnitude, surveyLimit=np.Inf):
        """
        Class constructor/initializer.

        Parameters
        ----------
    
        numberOfStars             - number of stars to simulate
        minDistance               - lower limit on the distance (closest possible star, pc)
        maxDistance               - upper limit on the distance (volume limit of survey, pc)
        meanAbsoluteMagnitude     - Mean of Gaussian absolute magnitude distribution
        stddevAbsoluteMagnitude   - Standard deviation of Gaussian absolute magnitude distribution
  
        Keywords
        --------

        surveyLimit - Apparent magnitude limit of the survey (default: no limit)
        """
        super().__init__(numberOfStars, minDistance, maxDistance, meanAbsoluteMagnitude,
                stddevAbsoluteMagnitude, surveyLimit)
        self.parallaxErrorNormalizationMagnitude=4.0
        self.parallaxErrorLogarithmicSlope=0.157
        self.parallaxErrorLogCalibrationFloor=-0.7
        self.magnitudeErrorNormalizationMagnitude=5.0
        self.magnitudeErrorLogarithmicSlope=0.2
        self.magnitudeErrorLogCalibrationFloor=-2.4

    def _generateParallaxErrors(self):
        """
        Generate the parallax errors according to an ad-hoc function of parallax error as a function of
        magnitude.
        """
        errors = self.parallaxErrorLogCalibrationFloor + self.parallaxErrorLogarithmicSlope * \
                (self.apparentMagnitudes - self.parallaxErrorNormalizationMagnitude)
        indices = (errors < self.parallaxErrorLogCalibrationFloor)
        errors[indices] = self.parallaxErrorLogCalibrationFloor
        return np.power(10.0,errors)

    def _generateApparentMagnitudeErrors(self):
        """
        Generate the apparent magnitude errors to an ad-hoc function of magnitude error as a function of
        magnitude.
        """
        errors = self.magnitudeErrorLogCalibrationFloor + self.magnitudeErrorLogarithmicSlope * \
                (self.apparentMagnitudes - self.magnitudeErrorNormalizationMagnitude)
        indices = (errors < self.magnitudeErrorLogCalibrationFloor)
        errors[indices] = self.magnitudeErrorLogCalibrationFloor
        return np.power(10.0,errors)

class UniformDistributionSingleLuminosityTGAS(UniformSpaceDistributionSingleLuminosity):
    """
    Simulate a parallax survey for stars distributed uniformly in space around the sun. The stars all
    have the same luminosity drawn from a Gaussian distribution. The errors on the observed parallaxes
    and apparent magnitudes roughly follow the characteristics of the TGAS Catalogue.
    """

    def __init__(self, numberOfStars, minDistance, maxDistance, meanAbsoluteMagnitude,
            stddevAbsoluteMagnitude, surveyLimit=np.Inf):
        """
        Class constructor/initializer.

        Parameters
        ----------
    
        numberOfStars             - number of stars to simulate
        minDistance               - lower limit on the distance (closest possible star, pc)
        maxDistance               - upper limit on the distance (volume limit of survey, pc)
        meanAbsoluteMagnitude     - Mean of Gaussian absolute magnitude distribution
        stddevAbsoluteMagnitude   - Standard deviation of Gaussian absolute magnitude distribution
  
        Keywords
        --------

        surveyLimit - Apparent magnitude limit of the survey (default: no limit)
        """
        super().__init__(numberOfStars, minDistance, maxDistance, meanAbsoluteMagnitude,
                stddevAbsoluteMagnitude, surveyLimit)
        self.magnitudeErrorNormalizationMagnitude=13.0
        self.magnitudeErrorLogarithmicSlope=0.21
        self.magnitudeErrorLogCalibrationFloor=-3.4
        self.tgasErrorsPdfFile = 'TGAS-parallax-errors-pdf.csv'
        self.tgasErrPdf = None
        if path.isfile(self.tgasErrorsPdfFile):
            self.tgasErrPdf = np.genfromtxt(self.tgasErrorsPdfFile, comments='#', skip_header=1,
                delimiter=',', names=['err','logdens'], dtype=None)
        else:
            print("Cannot find file {0}".format(self.tgasErrorsPdfFile))
            print("exiting")
            exit()

    def _generateParallaxErrors(self):
        """
        Generate the parallax errors according to an ad-hoc function of parallax error as a function of
        magnitude.
        """
        probs = np.exp(self.tgasErrPdf['logdens'])/np.sum(np.exp(self.tgasErrPdf['logdens']))
        errors = np.random.choice(self.tgasErrPdf['err'], size=self.apparentMagnitudes.size, p=probs)
        return errors

    def _generateApparentMagnitudeErrors(self):
        """
        Generate the apparent magnitude errors to an ad-hoc function of magnitude error as a function of
        magnitude.
        """
        errors = self.magnitudeErrorLogCalibrationFloor + self.magnitudeErrorLogarithmicSlope * \
                (self.apparentMagnitudes - self.magnitudeErrorNormalizationMagnitude)
        indices = (errors < self.magnitudeErrorLogCalibrationFloor)
        errors[indices] = self.magnitudeErrorLogCalibrationFloor
        return np.power(10.0,errors)

def showSurveyStatistics(simulatedSurvey, pdfFile=None, pngFile=None, usekde=False):
    """
    Produce a plot with the survey statistics.

    Parameters
    ----------

    simulatedSurvey : Object containing the simulated survey.

    Keywords
    --------

    pdfFile : string
        Name of optional PDF file in which to save the plot.
    pngFile : string
        Name of optional PNG file in which to save the plot.
    usekde  : boolean
        If true use kernel density estimates to show the distribution of survey quantities instead of
        histograms.
    """
    try:
        _ = simulatedSurvey.observedParallaxes.shape
    except AttributeError:
        stderr.write("You have not generated the observations yet!\n")
        return

    parLimitPlot=50.0
    plxSnrLim = 5.0

    positiveParallaxes = (simulatedSurvey.observedParallaxes > 0.0)
    goodParallaxes = (simulatedSurvey.observedParallaxes/simulatedSurvey.parallaxErrors >= plxSnrLim)
    estimatedAbsMags = (simulatedSurvey.observedMagnitudes[positiveParallaxes] +
            5.0*np.log10(simulatedSurvey.observedParallaxes[positiveParallaxes])-10.0)
    relParErr = (simulatedSurvey.parallaxErrors[positiveParallaxes] /
            simulatedSurvey.observedParallaxes[positiveParallaxes])
    deltaAbsMag = estimatedAbsMags - simulatedSurvey.absoluteMagnitudes[positiveParallaxes]

    useagab(usetex=False, fontfam='sans')
    fig = plt.figure(figsize=(18,12))
  
    axA = fig.add_subplot(2,2,1)
    apply_tufte(axA, withgrid=False)
    axA.set_prop_cycle(cycler('color', get_distinct(3)))

    minPMinThird=np.power(simulatedSurvey.minParallax,-3.0)
    maxPMinThird=np.power(parLimitPlot,-3.0)
    x=np.linspace(simulatedSurvey.minParallax,np.min([parLimitPlot,simulatedSurvey.maxParallax]),1001)
    axA.plot(x,3.0*np.power(x,-4.0)/(minPMinThird-maxPMinThird),'--', label='model', lw=3)

    if usekde:
        scatter = rse(simulatedSurvey.trueParallaxes)
        bw = 1.06*scatter*simulatedSurvey.numberOfStarsInSurvey**(-0.2)
        kde = KernelDensity(bandwidth=bw)
        kde.fit(simulatedSurvey.trueParallaxes[:,None])
        samples = np.linspace(simulatedSurvey.trueParallaxes.min(), simulatedSurvey.trueParallaxes.max(), 200)[:,None]
        logdens = kde.score_samples(samples)
        axA.plot(samples, np.exp(logdens), '-', lw=3, label='true')
    else:
        axA.hist(simulatedSurvey.trueParallaxes, bins='auto', density=True, histtype='step', lw=3,
                label='true')

    if usekde:
        scatter = rse(simulatedSurvey.observedParallaxes)
        bw = 1.06*scatter*simulatedSurvey.numberOfStarsInSurvey**(-0.2)
        kde = KernelDensity(bandwidth=bw)
        kde.fit(simulatedSurvey.observedParallaxes[:,None])
        samples = np.linspace(simulatedSurvey.observedParallaxes.min(), simulatedSurvey.observedParallaxes.max(), 200)[:,None]
        logdens = kde.score_samples(samples)
        axA.plot(samples, np.exp(logdens), '-', lw=3, label='observed')
    else:
        axA.hist(simulatedSurvey.observedParallaxes, bins='auto', density=True, histtype='step', lw=3,
                label='observed')

    axA.set_xlabel(r'$\varpi$,  $\varpi_\mathrm{true}$ [mas]')
    axA.set_ylabel(r'$p(\varpi)$, $p(\varpi_\mathrm{true})$')
    leg=axA.legend(loc='best', handlelength=1.0)
    for t in leg.get_texts():
        t.set_fontsize(14)

    axB = fig.add_subplot(2,2,2)
    apply_tufte(axB, withgrid=False)
    axB.set_prop_cycle(cycler('color', get_distinct(3)))

    m = np.linspace(simulatedSurvey.observedMagnitudes.min(), simulatedSurvey.observedMagnitudes.max(), 1000)
    axB.plot(m, np.exp(simulatedSurvey.apparentMagnitude_lpdf(m)), '--', lw=3, label='model')

    if usekde:
        scatter = rse(simulatedSurvey.apparentMagnitudes)
        bw = 1.06*scatter*simulatedSurvey.numberOfStarsInSurvey**(-0.2)
        kde = KernelDensity(bandwidth=bw)
        kde.fit(simulatedSurvey.apparentMagnitudes[:,None])
        samples = np.linspace(simulatedSurvey.apparentMagnitudes.min(), simulatedSurvey.apparentMagnitudes.max(), 200)[:,None]
        logdens = kde.score_samples(samples)
        axB.plot(samples, np.exp(logdens), '-', label='true', lw=3)
    else:
        axB.hist(simulatedSurvey.apparentMagnitudes, bins='auto', density=True, histtype='step', lw=3,
                label='true')

    if usekde:
        scatter = rse(simulatedSurvey.observedMagnitudes)
        bw = 1.06*scatter*simulatedSurvey.numberOfStarsInSurvey**(-0.2)
        kde = KernelDensity(bandwidth=bw)
        kde.fit(simulatedSurvey.observedMagnitudes[:,None])
        samples = np.linspace(simulatedSurvey.observedMagnitudes.min(), simulatedSurvey.observedMagnitudes.max(), 200)[:,None]
        logdens = kde.score_samples(samples)
        axB.plot(samples, np.exp(logdens), '-', label='observed', lw=3)
    else:
        axB.hist(simulatedSurvey.observedMagnitudes, bins='auto', density=True, histtype='step', lw=3,
                label='observed')

    axB.set_xlabel("$m$, $m_\mathrm{true}$")
    axB.set_ylabel("$p(m)$, $p(m_\mathrm{true})$")
    leg=axB.legend(loc='upper left', handlelength=1.0)
    for t in leg.get_texts():
        t.set_fontsize(14)

    axC = fig.add_subplot(2,2,3)
    apply_tufte(axC, withgrid=False)
    axC.set_prop_cycle(cycler('color', get_distinct(3)))

    x=np.linspace(simulatedSurvey.absoluteMagnitudes.min(),
            simulatedSurvey.absoluteMagnitudes.max(), 300)
    axC.plot(x, norm.pdf(x,loc=simulatedSurvey.meanAbsoluteMagnitude,
        scale=simulatedSurvey.stddevAbsoluteMagnitude), '--', lw=3, label='model')

    if usekde:
        scatter = rse(simulatedSurvey.absoluteMagnitudes)
        bw = 1.06*scatter*simulatedSurvey.numberOfStarsInSurvey**(-0.2)
        kde = KernelDensity(bandwidth=bw)
        kde.fit(simulatedSurvey.absoluteMagnitudes[:,None])
        samples = np.linspace(simulatedSurvey.absoluteMagnitudes.min(), simulatedSurvey.absoluteMagnitudes.max(), 200)[:,None]
        logdens = kde.score_samples(samples)
        axC.plot(samples, np.exp(logdens), '-', label='true', lw=3)
    else:
        axC.hist(simulatedSurvey.absoluteMagnitudes, bins='auto', density=True, histtype='step', lw=3,
                label='true')

    if (simulatedSurvey.absoluteMagnitudes[goodParallaxes].size >= 3):
        if usekde:
            scatter = rse(simulatedSurvey.absoluteMagnitudes[goodParallaxes])
            bw = 1.06*scatter*simulatedSurvey.absoluteMagnitudes[goodParallaxes].size**(-0.2)
            kde = KernelDensity(bandwidth=bw)
            kde.fit(simulatedSurvey.absoluteMagnitudes[goodParallaxes][:,None])
            samples = np.linspace(simulatedSurvey.absoluteMagnitudes[goodParallaxes].min(),
                    simulatedSurvey.absoluteMagnitudes[goodParallaxes].max(), 200)[:,None]
            logdens = kde.score_samples(samples)
            axC.plot(samples, np.exp(logdens), '-', label=r'$\varpi/\sigma_\varpi\geq{0:.1f}$'.format(plxSnrLim), lw=3)
        else:
            axC.hist(simulatedSurvey.absoluteMagnitudes[goodParallaxes], bins='auto', density=True, histtype='step', lw=3,
                    label=r'$\varpi/\sigma_\varpi\geq{0:.1f}$'.format(plxSnrLim))
    
    axC.set_xlabel("$M$")
    axC.set_ylabel("$p(M)$")
    leg=axC.legend(loc='best', handlelength=1.0)
    for t in leg.get_texts():
        t.set_fontsize(14)

    axD = fig.add_subplot(2,2,4)
    apply_tufte(axD, withgrid=False)
    axD.set_prop_cycle(cycler('color', get_distinct(3)))
    #if len(relParErr) < 1000:
    axD.plot(simulatedSurvey.trueParallaxes, simulatedSurvey.observedParallaxes-simulatedSurvey.trueParallaxes, '.')
    axD.plot(simulatedSurvey.trueParallaxes[positiveParallaxes],
            simulatedSurvey.observedParallaxes[positiveParallaxes]-simulatedSurvey.trueParallaxes[positiveParallaxes], '.')
    axD.plot(simulatedSurvey.trueParallaxes[goodParallaxes],
            simulatedSurvey.observedParallaxes[goodParallaxes]-simulatedSurvey.trueParallaxes[goodParallaxes], 'o')
    #else:
    #    axD.hexbin(simulatedSurvey.trueParallaxes,
    #            simulatedSurvey.observedParallaxes-simulatedSurvey.trueParallaxes, C=None, cmap=cm.Blues_r, mincnt=1)
    axD.set_xlabel(r"$\varpi_\mathrm{true}$ [mas]")
    axD.set_ylabel("$\\varpi-\\varpi_\\mathrm{true}$ [mas]")

    plt.suptitle("Simulated survey statistics: $N_\\mathrm{{stars}}={0}$, ".format(simulatedSurvey.numberOfStars) +
            "$m_\\mathrm{{lim}}={0}$, ".format(simulatedSurvey.apparentMagnitudeLimit) +
            "$N_\\mathrm{{survey}}={0}$, ".format(simulatedSurvey.numberOfStarsInSurvey) +
            "${0}\\leq\\varpi\\leq{1}$, ".format(simulatedSurvey.minParallax, simulatedSurvey.maxParallax)+
            "$\\mu_M={0}$, ".format(simulatedSurvey.meanAbsoluteMagnitude) + 
            "$\\sigma_M={0:.2f}$".format(simulatedSurvey.stddevAbsoluteMagnitude))
  
    if pdfFile is not None:
        plt.savefig(pdfFile)
    if pngFile is not None:
        plt.savefig(pngFile)
    if (pdfFile is None and pngFile is None):
        plt.show()

def marginal_pdf_distance(r, rmin, rmax, mu, sigma, mlim):
    """
    Calculate the expected marginal distribution of distances given the parallax survey parameters. The
    calculation is only approximate as a the magnitude limit is applied to the error-free true apparent
    magnitude.
    
    Parameters
    ----------
    
    r : float vector
        Values of r for which to calculate p(r).
    rmin : float
        Minimum distance in survey.
    rmax : float
        Maximum distance in survey.
    mu : float
        Mean of the true absolute magnitude distribution.
    sigma : float
        Standard deviation of the true absolute magnitude distribution.
    mlim : float
        Apparent magnitude limit of the survey.
        
    Returns
    -------
    
    p(r) as float vector.
    """
    A = rmax**3-rmin**3
    pdf = lambda x : np.exp(np.log(3) - np.log(A) + 2*np.log(x) + norm.logcdf(mlim-mu-5*np.log10(x)+5, scale=sigma))
    C, dummy = quad(pdf, rmin, rmax)
    return pdf(r)/C

def marginal_pdf_absMag(M, rmin, rmax, mu, sigma, mlim):
    """
    Calculate the expected marginal distribution for the true absolute magnitude given the parallax
    survey parameters. The calculation is only approximate as a the magnitude limit is applied to the
    error-free true apparent magnitude.
    
    Parameters
    ----------
    
    M : float vector
        Values of M_true for which to calculate p(M_true).
    rmin : float
        Minimum distance in survey.
    rmax : float
        Maximum distance in survey.
    mu : float
        Mean of the true absolute magnitude distribution.
    sigma : float
        Standard deviation of the true absolute magnitude distribution.
    mlim : float
        Apparent magnitude limit of the survey.
        
    Returns
    -------
    
    p(M_true) as float vector.
    """
    A = rmax**3-rmin**3
    def _pdf(x):
        rlim = np.power(10,0.2*(mlim-x+5))
        if isscalar(rlim):
            if rlim>rmax:
                rlim=rmax
        else:
            rlim[(rlim>rmax)] = rmax
        return (rlim**3 - rmin**3)/A*norm.pdf(x, loc=mu, scale=sigma)

    C, dummy = quad(_pdf, mu-7*sigma, mu+8*sigma)
    return _pdf(M)/C
