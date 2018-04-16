"""
Quick and dirty fit of TGAS parallax error distribution.
Uses an inverse Gamma distribution.

Anthony G.A. Brown Aug 2016 - Nov 2017
<brown@strw.leidenuniv.nl>
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import argparse
from os import path

from astropy.io import fits
from plotstyles import useagab, apply_tufte

from sklearn.neighbors import KernelDensity
from scipy.stats import invgamma, gamma
from scipy.integrate import simps

def fit_and_plot(args):
    """
    Read in the TGAS parallax error data (or an already calculated KDE of the parallax error
    distribution) and the "fit" the KDE by with a Gamma distribution. Fit quality is simply judged by
    eye!

    Parameters
    ----------

    args - command line arguments.
    """

    offset, alpha, beta = args['distroParams']
    tgasPlxErrPdfFile = 'TGAS-parallax-errors-pdf.csv'

    if path.isfile(tgasPlxErrPdfFile):
        errpdf = np.genfromtxt(tgasPlxErrPdfFile, comments='#', skip_header=1, delimiter=',',
                names=['err','logdens'], dtype=None)
        x = errpdf['err'][:,np.newaxis]
        logdens = errpdf['logdens']
    else:
        tgasData = fits.open('TGAS-allPlxErrorsVsGmag.fits')[1].data
        eplx = tgasData['parallax_error'][:,np.newaxis]
        kde = KernelDensity(kernel='epanechnikov',bandwidth=0.008).fit(eplx)
        x=np.linspace(0,1,1001)[:,np.newaxis]
        logdens=kde.score_samples(x)
        f=open('TGAS-errors-kde.csv', mode='w')
        f.write('#err,logdens\n')
        for xx,yy in zip(x[:,0],logdens):
            f.write("{0},{1}\n".format(xx,yy))
        f.close()

    xx=np.linspace(0.004,1,1000)
    ig=invgamma.pdf(xx,alpha,loc=offset,scale=beta)
    norm = invgamma.cdf(1.0,alpha,loc=offset,scale=beta)
    ig=ig/norm

    useagab(usetex=False, fontfam='sans', sroncolours=False)

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    apply_tufte(ax)
    ax.plot(x[:,0],np.exp(logdens), label='KDE of $\\sigma_\\varpi$ distribution', lw=3)
    ax.plot(xx,ig, label="Fit with InvGamma($x-{0:.3f}|\\alpha={1:.3f}$, $\\beta={2:.3f})$".format(offset,alpha,beta))
    ax.set_xlabel('$\\varpi$ [mas]')
    ax.set_ylabel('pdf')
    ax.legend(loc='upper right', fontsize=12)
    ax.set_title('Parallax error distribution for TGAS', fontsize=14)

    basename = 'FitOfParallaxErrorDistribution'
    if args['pdfOutput']:
        plt.savefig(basename+'.pdf')
    elif args['pngOutput']:
        plt.savefig(basename+'.png')
    else:
        plt.show()

def parseCommandLineArguments():
    """
    Set up command line parsing.
    """
    parser = argparse.ArgumentParser(description="""Basic fit of TGAS parallax error distribution.""")
    parser.add_argument("--distroParams", help="""Parameters of InvGamma distribution, offset in x, alpha,
    beta (default: 0.185, 1.7, 0.2)""", nargs=3, type=float, default=[0.185, 1.7, 0.2])
    parser.add_argument("-p", action="store_true", dest="pdfOutput", help="Make PDF plot")
    parser.add_argument("-b", action="store_true", dest="pngOutput", help="Make PNG plot")
    args=vars(parser.parse_args())
    return args

if __name__ in ('__main__'):
    args=parseCommandLineArguments()
    fit_and_plot(args)
