# Getting things to work

In principle if a [Python](https://www.python.org) and [R](https://www.r-project.org) environment
are available on your system, including facilities to run [Jupyter Notebooks](https://jupyter.org/),
everything should "just work" (haha, famous last words). The instructions below are to ensure that
your setup contains all the necessary dependencies to run the various tutorials.

## Linux

### Python environment (including Jupyter notebook)

Make your life easy and install the [anaconda](https://www.anaconda.com/download/) distribution,
this can be dropped wherever you like on your hard drive without needing root permissions.
Subsequently install the following:

```
> conda install pystan
> conda install -c astropy corner
> pip install daft (optional)
```

where the installation of `daft` is optional.

### R Environment

There is no equivalent to anaconda for R but installation from source should be fairly painless, or
use the package manager of your Linux distro. Subsequently make sure that your can run R notebooks
by following the [instructions](https://irkernel.github.io/installation/) on
[irkernel.github.io](https://irkernel.github.io) (the "binary" version was used when writing this
INSTALL file). Step 1/2 should be done as "root" for system wide R installations. Step 2/2 should be
done as user. This should result in the `magrittr` R package being installed (needed for some of the
tutorials). Subsequently install the following packages from within R (as root for system wide R
installations):

```
> install.packages("mvtnorm")
> install.packages("PolynomF")
> install.packages("fields")
> install.packages("RColorBrewer")

> install.packages("png")
> install.packages("ggplot2")

> Sys.setenv(MAKEFLAGS = "-j4") 
> install.packages("rstan", type = "source")
```

Note that the [period luminosity relation](./period-luminosity-relation) tutorial requires the
installation of R packages related to graph drawing, as well as the `rpy2` interface between R and
Python. This can be avoided by using the [modified
version](./period-luminosity-relation/TutorialPLZ_R.ipynb) of that tutorial.
