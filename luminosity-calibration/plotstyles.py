"""
Provides plotting styles.

Anthony Brown Aug 2015 - Jan 2018
"""

import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler

from distinct_colours import get_distinct

def useagab(usetex=True, fontfam='serif', sroncolours=True):
    """
    Configure the plotting style to my liking.

    Parameters
    ----------

    None

    Keywords
    --------

    usetex : boolean
        Whether or not to use LaTeX text (default True).
    fontfam : boolean
        Font family to use (default 'serif')
    sroncolours : boolean
        If true use colour-blind proof distinct colours (https://personal.sron.nl/~pault/).

    Returns
    -------

    Nothing
    """
    line_colours = get_distinct(4)
    if (usetex):
        rc('text', usetex=True)
        rc('text.latex', preamble=r'\usepackage{amsmath}')
    rc('font', family=fontfam, size=18)
    rc('xtick.major', size='6')
    rc('xtick.minor', size='4')
    rc('ytick.major', size='6')
    rc('ytick.minor', size='4')
    rc('lines', linewidth=2.0)
    rc('axes', linewidth=1)
    if (sroncolours):
        rc('axes', prop_cycle=(cycler('color',line_colours)))
    else:
        rc('axes', prop_cycle=(cycler('color',plt.cm.tab10.colors)))
    rc('xtick', direction='out')
    rc('ytick', direction='out')
    rc('grid', color='cbcbcb')
    rc('grid', linestyle='-')
    rc('grid', linewidth=0.5)
    rc('grid', alpha=1.0)
    rc('figure', facecolor='ffffff')
    rc('figure', dpi=80)
    rc('figure.subplot', bottom=0.125)

def apply_tufte(ax, withgrid=False, minorticks=False):
    """
    Apply the "Tufte" style to the plot axes contained in the input axis object. This mimics the sparse
    style advocated by Tufte in his book "The Visual Display of Quantitative Information".

    Parameters
    ----------

    ax - The axis object to configure.

    Keywords
    --------

    withgrid - When true a grid is displayed in the plot background
    minorticks - When true minor tickmarks are drawn.

    Returns
    -------

    Nothing.
    """

    # Move left and bottom spines outward by 5 points
    ax.spines['left'].set_position(('outward', 5))
    ax.spines['bottom'].set_position(('outward', 5))
    # Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params('both', width=1.5, which='both')
    if withgrid:
        ax.grid(True)
    if minorticks:
        ax.minorticks_on()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    ax.tick_params('both', width=1.5, which='major')
    ax.set_facecolor('w')
