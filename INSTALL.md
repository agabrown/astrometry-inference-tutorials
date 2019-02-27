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
> conda install -c astropy corner
> pip install pygaia
> pip install daft (optional)
```

where the installation of `daft` is optional.

#### PyStan

Follow the instructions [here](https://pystan.readthedocs.io/en/latest/getting_started.html). In
particular, use `pip install pystan`. The necessary compilers can be installed using `conda`.

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
```

Note that the [period luminosity relation](./period-luminosity-relation) tutorial requires the
installation of R packages related to graph drawing, as well as the `rpy2` interface between R and
Python. This can be avoided by using the [modified
version](./period-luminosity-relation/TutorialPLZ_R.ipynb) of that tutorial.

#### RStan

The recommended installation method for RStan (working with RStudio) can be found
[here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

To install from source follow the instructions
[here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Linux) for the installation of
Rstan. Where the last step is to do from within R (run as root for a system wide installation):

```
> Sys.setenv(MAKEFLAGS = "-j4") 
> install.packages("rstan", type = "source")
```

Note that a `Makevars` file in the ~/.R folder with the following contents is useful for
optimization of the compiled Stan models and supressing of irrelevant warnings (see [RStan install
instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Linux)):

```
CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -Wno-builtin-macro-redefined

CXXFLAGS+=-flto -Wno-unused-local-typedefs

CXXFLAGS += -Wno-ignored-attributes -Wno-deprecated-declarations
```

## Windows (only Windows 10 tested)

### Python environment

The [anaconda](https://www.anaconda.com/download/) distribution is again recommended. After
installation open the Anaconda prompt and install the following:

```
> conda install -c astropy corner
> pip install pygaia
> pip install daft (optional)
```

where the installation of `daft` is optional.

__IT IS IMPORTANT__ to now add the location of the anaconda and python executables to your path:
Go to "Settings" (from the windows start menu for example) and then search for "environment" in the
search field at the top of the window. Select "Edit environment variables for your account" and then
edit the "Path" variable, adding the following two paths:

```
C:\Users\....\Anaconda3
C:\Users\....\Anaconda3\Scripts
```

#### PyStan

Installation instructions for pystan can be found
[here](https://pystan.readthedocs.io/en/latest/windows.html#windows). Note in particular that
installation through conda (`conda install pystan`) does not seem to work well. So use pip (after
having installed the necessary build tools as detailed in the instructions).

### R environment

Install the windows version of R and then install the following packages while running R as
"administrator":

```
> install.packages("mvtnorm")
> install.packages("PolynomF")
> install.packages("fields")
> install.packages("RColorBrewer")

> install.packages("png")
> install.packages("ggplot2")
```

Note that the [period luminosity relation](./period-luminosity-relation) tutorial requires the
installation of R packages related to graph drawing, as well as the `rpy2` interface between R and
Python. This can be avoided by using the [modified
version](./period-luminosity-relation/TutorialPLZ_R.ipynb) of that tutorial.

#### RStan

The recommended installation method for RStan (working with RStudio) can be found
[here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

To install from source follow the instructions
[here](https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows) for the
installation of Rstan. 

Subsequently make sure that your can run R notebooks by following the
[instructions](https://irkernel.github.io/installation/) on
[irkernel.github.io](https://irkernel.github.io). Step 1/2 should be done while running as R as
"administrator".

## Mac OsX

### Python environment

Use anaconda for Mac OsX and follow the instruction above for Linux.

#### PyStan

pending...

### R Environment

Use the `.pkg` file from the [CRAN page](https://cran.r-project.org). Install the packages listed
above in the Linux section.

#### RStan

pending...
