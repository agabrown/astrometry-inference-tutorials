# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:41:50 2017

@author: Ariadna
"""

#!/usr/bin/env python
__author__ = "Ariadna Ribes Metidieri"
__copyright__ = ""
__credits__ = ["Xavi Luri, Enrique Utrilla, Alfred Castro, David Orellana, Patricia Ribes"]
__license__ = "GPL v 3.0"
__version__ = "1.0.1"
__maintainer__ = "Ariadna Ribes Metidieri"
__email__ = "aribesmetidieri@gmail.com"
__status__ = "active"

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from scipy.optimize import fsolve
from decimal import *
import matplotlib.pyplot as plt


"""DISTANCE FUNCTIONS"""
#-------------------------------------------POSTERIORS-----------------------------------------------------------------#
def likelihood(x,w,s):
    """Return the likelihood
    
    ===========
    Parameters:
    ===========
        **x**
        : number
            Real distance in kpc units.
        **w**
        : number
            Observed parallax in mas units.
        **s** 
        : number
            Parallaw error in mas units.
    
    ========
    Returns:
    ========
        **number**
            Returns the likelihood.
    """
    
    return 1/(np.sqrt(2*np.pi)*s)*np.exp(-1/(2*s**2)*(w-1/x)**2)

def uniform_distance_posterior(x,w,s,r_lim):    
    """Return the unnormalized probability of a real distance according to the `Uniform Distance Posterior`_.
    
    ===========
    Parameters:
    ===========
        **x**
        : number
            Real distance in kpc units.
        **w**
        : number
            Observed parallax in mas units.
        **s** 
        : number
            Parallaw error in mas units.
        **r_lim**
        : number
            Parameter of the distribution in kpc. For x > r_lim or negative parallaxes returns 0.
    ========
    Returns:
    ========
        **number**
            Unnormalized probability according to the Uniform Distance Posterior.
        
    .. _Uniform Distance Posterior:
        http://iopscience.iop.org/article/10.1086/683116/pdf
    """
    if x > 0 and x <= r_lim:
        return np.float(np.exp(np.float(-(w-1/x)**2/(2*s**2))))/(s*np.sqrt(2*np.pi)*r_lim)
    else:
        return 0
    
def exponentially_decreasing_space_density_posterior(x,w,s,L):
    """Return the unnormalized probability of a real distance according to the  `Exponentially Decreasing Space Density Posterior`_.
    
    ===========
    Parameters:
    ===========
        **x**
        : number
            Real distance in kpc units.
        **w** 
        : number
            Observed parallax in mas units.
        **s**
        : number
            Parallaw error in mas units.
        **r_lim**
        : number
            Parameter of the distribution in kpc. For x > r_lim or negative parallaxes returns 0.
    ========
    Returns:
    ========
        **number**
            Unnormalized probability according to the Exponentially Decreasing Space Density Posterior.
    
    .. _Exponentially Decreasing Space Density Posterior:
        http://iopscience.iop.org/article/10.1086/683116/pdf
    """
    
    if x > 0:
        return x**2*np.exp(-x/L)*np.exp(-(w-1/x)**2/(2*s**2))/(s*np.sqrt(2*np.pi)*2*L**3)
    else:
        return 0


# -------------------------------------------------MODES----------------------------------------------------------------#
def mode_r_uniform(w, r_lim): # Mode according to the uniform distance PDF posterior
    # If the distance computed as the invers of the parallax is not positive or inferior than r_lim the estimated distance is r_lim
    """Returns the mode of the Uniform Distance Probability Distribution Function.
    
    ===========
    Parameters:
    ===========
        **w**
        : number
            Observed parallax in mas units.
        **r_lim**
        : number
            Parameter of the distribution. 
    ========
    Returns:
    ========
        **number**
            If w is negative or `1/w > r_lim` returns r_lim, otherwise returns `1/w`.
    """
    r = 1./w
    if (r <= r_lim) and (r > 0):
        r_mo = r
    else:
        r_mo = r_lim
    return r_mo

def mode_r_exponential(L, w, sig_w): 
    """Returns the mode of the Exponentially Decreasing Space Density Probability Distribution Function.
    
    This method solves the mode equation, obtained through the minimization of the Exponentially Decreasing Space Density Posterior using the numpy.roots function. It also discards the complex roots and selects the real root corresponding to the mode. 
    ===========
    Parameters:
    ===========
        **L**
        : positive number 
            Length scale in kpc. 
        **w**
        : number
            Observed parallax in mas units.
        **sig_w**
        : number
            Parallax error in mas units.
    ========
    Returns:
    ========
        **number**
            If there are 3 solutions and w > 0 it returns the smallest root and if w <0 it  returns the largest one. If there is only one solution, it returns the absolute value.
    """
    a = np.float64(1. / L)
    b = np.float64(-2.)
    c = np.float64(w / sig_w ** 2)
    d = np.float64(- 1. / sig_w ** 2)

    p = np.array([a, b, c, d])

    sol = np.roots(p) #Solving the third degree equation

    l = len(sol)
    solutions = []

    for s in sol:
        if abs(s) == s: #Selecting the real roots
            solutions.append(s)

    # Given the physical limitations s>0 and L>0, only two cases are possible:
    # - All three roots are real and there are two modes and one minimum. If w>0 the solution is the smallest one while if w<0, the greatest one
    # is the solution.
    # - There's only one root.
    if len(solutions) == 3:  # If there are three real roots and the observed parallax is positive, the maximum value will be taken as the solution
        if w >= 0:
            r_mo = np.float(min(solutions))

        else:
            r_mo = np.float(max(solutions))
    else:

        r_mo = np.float64(abs(solutions[0]))
    return r_mo

#---------------------------------------------R*-------------------------------#
def rstar(beta,s,w):
    """Returns the estimate r* of the  `Transformation Method`_ described by Haywood  Smith.
    
    ===========    
    Parameters:
    ===========    
        **beta**
        : number
            Numeric parameter of the transformarion distance estimate. Default is 1.01.
        **s** 
        : number
            Parallax error in arcseconds.
        **w**
        : number
            Observed parallax in arcseconds.
    ========
    Returns:
    ========    
        **number**
            Transformation Distance Estimate in pc.
        
    .. _Transformation Method:
        http://www.astro.ufl.edu/~asthcs/BadH2.pdf
    """
    if w > 0:
        wstar = beta*s/0.8*np.log(1+np.exp(0.8*w/s))
    else:
        wstar = beta*s/0.8*np.log(1+np.exp(0.8*w/s))*np.exp(-0.605*w**2/s**2)
    
    return 1/wstar

#------------------------------------PERCENTILES-------------------------------#
def percentiles(f,r0,r_mode,w,s,par): # Given a distance (i.e. r_mode), it integrates the PDF from r0 to the given distance and returns the unnormalized percentile
    """Returns the percentile corresponding to a given distance r_mode.
    
    This function integrates the function f with parameters w and s from r0 to r_mode using the scipy.integrate.quad function.
    
    ===========
    Parameters:
    ===========    
        **f** 
        : function
            Probability distribution function.
        **r0**
        : number
            Inferior bound in the integration.
        **r_mode**
        : number
            Superio bound in the integration.
        **w** 
        : number
            Parameter of f. 
        **s**
        : number
            Parameter of f.
    ========
    Returns:
    ========
        **number**
            Unnormalized percentile of r_mode.
    """
    p = quad(f,r0,r_mode,args=(w,s,par),epsrel = 1e-12, epsabs = 1e-12)
    return p[0]

def normalized_percentile(p,n):
    """Returns the normalized percentile given the unnormalized percentile and the normalization constant.
    
    If the normalization factor is different from zero, it divides the unnormalized percentile by the normalization factor. 
    
    ===========
    Parameters:
    ===========
        **p**
        : number
            Unnormalized percentile.
        **n**
        : number
            Normalization factor.
    ========        
    Returns:
    ========    
        **float**
            Normalized percentile. Its value ranges from 0 to 1.
    """
    tol = 10**(-25)
    if n == 0:
       
        n = tol
        if p == 0:
            p = tol
   
    return p/n


#-----------------------------------------------NORMALIZATION----------------------------------------------------------#

def normalization(f,par1,w,s,p,r_mode,par2): 
    """Returns the normalization factor of the function f.

    Computes the normalization factor as the summation of the percentile of r_mode plus the integration of the function from r_mode to par2 using the scipy.integrate.quad function. This way numerical errors in the integrration of thin distributions are avoided.
    
    ===========
    Parameters:
    ===========
       **f**
       : function
           Probability distribution function. 
       **par1**
       : number
           Parameter of the function f in kpc.    
       **w**
       : number
           Parameter of the function f in mas.    
       **s**
       : number
           Parameter of the function f in mas.    
       **p**
       : number
           Unnormalized percentile of r_mode.    
       **r_mode**
       : number
           Distance in kpc.      
       **par_2**
       : number
           Superior bound of the integration in kpc.
    ========
    Returns:
    ========
        **number**
            Normalization factor.
   """
    
    N = quad(f, r_mode, par2, args=(w, s, par1),epsrel = 1e-11, epsabs = 1e-11) # We integrate the required posterior from the mode to infinity (par2)
    
    N = N[0] + p # p is the percentile corresponding to the mode, i.e. the integration of the PDF from r0 to r_mode
    
    return N

# -------------------------------------------------MEDIAN---------------------------------------------------------------#

def median(f,w,s,par,par2,r_mode,p,N):
    """Returns the median of the distribution f.
    
    Computes the median of the distribution solving the implicit equation `cdf(x)-0.5 = 0` through the scipy.optimize.brentq function, where `cdf(x)` is the normalized integration from r0 = 0.01 kpc to x using the scipy.integrate.quad function.
    
    ===========
    Parameters:
    ===========    
        **f**
        : function
            Probability distribution function.
        **w**
        :number
            Parameter of the function f in mas.    
        **s**
        : number
            Parameter of the function f in mas.    
        **p** 
        : number
            Unnormalized percentile of r_mode.
        **N**
        : number
            Normalization factor.
    ========
    Returns:
    ========
        **number** 
            Median of the distribution in kpc.
    """
    r0 = 0.01 #kpc
    # The median is the distance which corresponds to the 50% percentile
    # par2 is the 'infinite' of every distribution
    
    if p == 0.5: # If the percentile of the mode = 0.5, the mode and the median coincide
        return r_mode
    
    else: # We try to find the smallest interval that we know contains the median
        if p>0.5: # If the mode is at the right of the median, the median will be between r0 and the mode
            r_inf = r0
            r_sup = r_mode
            
            def function(x):
                I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
                return I[0]/N - 0.5

            r_median = brentq(function, r_inf, r_sup) #Finding the median

            return r_median


        else: # If the mode is at the left of the median, the median will be between the mode and the limit  
            r_inf = r_mode
            r_sup = r_mode
            i = 0
            step = 0.1
            def function(x):
                I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
                return I[0]/N - 0.5
            
            while True:
                i = i + 1
                try:
                    r_median = brentq(function, r_inf, r_sup)
                    return r_median
                
                except ValueError:
                    r_sup = r_sup + step
                    
                if i == 10**3:
                    step = step/10
                    i = 0
                    r_sup = r_mode
            
            

# -------------------------------------------DISTANCE BOUNDS-----------------------------------------------------------#
def distances_from_percentiles(f,N,w,s,stype,p_inf,p_sup,p_ref,r_ref,par,par2):    
    """Returns the quantiles in kpc corresponding to the percentiles p_inf and p_sup.
    
    For stype = 'inf' finds the quantile corresponding to p_inf while otherwise finds  the quantile corresponding to p_sup. It solves the implicit equation `cdf(x)-p = 0` using the scipy.optimize.brentq rootfinding method, where `cdf(x)` is the normalized intregral of f from r0 to x.
    
    ===========
    Parameters:
    ===========
        **f**
        : function
            Probability distribution function.\n
        **N** 
        : number 
            Normalization factor.\n
        **w** 
        : number
            Parameter of the function f in mas.\n    
        **s**
        : number
            Parameter of the function f in mas.\n
        **stype**
        : str
            If stype = 'inf' computes the quantile corresponding to p_inf. Otherwise it computes the quantile corresponding to p_sup.\n
        **p_inf** 
        : float
            Inferior normalized percentile.\n 
        **p_sup**
        : float
            Superior normalized percentile.\n
        **p_ref**
        : float
            Reference percentile corresponding to r_ref.
        **r_ref**
        : float
            Distance reference in kpc corresponding to the percentile p_ref.
        **par**
        : number
            Parameter of f in kpc.\n
        **par2**
        : number
            Superior bound of the integration in kpc.\n
    ========
    Returns:
    ========
        **number**
            Quantile corresponding to the percentile p_inf or p_sup depending on stype.
    """
    r0 = 0.01 #kpc
    if stype == 'inf': # I'm looking for the inferior distance bound
        r_inf = r0
        r_sub = r_ref # The inferior limit will be between 0 and the median
        
        def function(x):
            I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-12, epsabs = 1e-12)
            return I[0] / N - p_inf
        
        r = brentq(function,r_inf,r_sub)
        
        return r

    else: #I'm looking for the superior distance bound
        def function(x):
            I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-12, epsabs = 1e-12)
            return I[0] / N - p_sup
        
        if (quad(f, r0, r_ref, args=(w, s, par), epsrel=1e-12, epsabs=1e-12)[0] / N - p_sup)<0:
            r_inf = r_ref
            
        else:
            r_inf = r0
        
        if (quad(f, r0, r_inf, args=(w, s, par), epsrel=1e-12, epsabs=1e-12)[0] / N - 0.95)*(quad(f, r0, par2, args=(w, s, par), epsrel=1e-11, epsabs=1e-11)[0] / N - 0.95) < 0:
            r_sup = par2
            
            r = brentq(function,r_inf,r_sup)
            return r
            
        else:
            i = 0
            step = 1
            r_sup = r_ref
            while True:
                i = i + 1
                try:
                    r = brentq(function, r_inf, r_sup)
                    return r
                
                except ValueError:
                    r_sup = r_sup + step
                    
                if i == 10**3:
                    step = step/10
                    i = 0
                    r_sup = r_ref
                    
           

def uncertainty_range_tm(w,s,beta,err):
    """Returns the uncertainty range associated to the estimate r* of the `Transformation Method`_.
        
    Returns the interval corresponding to [1/w*(w+err),1/w*(w-err)]. For instance, for err equal to the parallax error, it returns the 68.27 % confidence interval about r* in pc.
    ===========  
    Parameters:
    ===========
        **w**
        : number
            Observed parallax in arcseconds.
        **s**
        : number
            Parallax error in arcseconds.
        **beta**
        : number
            Parameter of the `Transformation Method`_.
        **err**
        : number
            Error interval.
    Returns:
        tuple: 
    .. _Transformation Method:
        http://www.astro.ufl.edu/~asthcs/BadH2.pdf
    """
    
    if w > 0:
        min_wstar = beta*s/0.8*np.log(1+np.exp(0.8*(w+err)/s))
        max_wstar =beta*s/0.8*np.log(1+np.exp(0.8*(w-err)/s))
        
    else:
        min_wstar = beta*s/0.8*np.log(1+np.exp(0.8*(w+err)/s))*np.exp(-0.605*(w+err)**2/s**2)
        max_wstar = beta*s/0.8*np.log(1+np.exp(0.8*(w-err)/s))*np.exp(-0.605*(w-err)**2/s**2)
        
    return (1/min_wstar,1/max_wstar)    




"""DISTANCE MODULUS FUNCTIONS"""
#---------------------------------UNIFORM------------------------------------------------------------------------------#
def likehdmu(x,w,s):
    return 1/(s*np.sqrt(2*np.pi))*np.exp(-(w-10**(-(x+5)/5))**2/(2*s**2))

def dmpdfun(x,w,s,r_lim):
     """ Returns the probability of a given distance modulus using the Uniform Distance Prior.
     
     ===========
     Parameters:
     ===========
         **x**
         : number
             Real distance modulus.
         **w**
         : number
             Observed parallax in arcseconds.
         **s**
         : number
             Parallax error in arcseconds.
         **r_lim**
         : number
             Parameter of the Uniform Distance Prior in pc.
      ========
      Returns:
      ========
          **number**
              Unnormalized probability according to the Uniform Distance Posterior.
     """
     if x >0 and x <=(5*np.log10(r_lim)-5):
         return 10**(x/5)/(5*r_lim)*likehdmu(x,w,s)
     else:
         return 0



def dmod_mode_ud(w,s,r_lim,f):
    """Returns the distance modulus mode of the probability distribution function.
    
    ===========
    Parameters:
    ===========
        **w**
        : number
            Observed parallax in arcseconds.
        **s**
        : number
            Parallax error in arcseconds.
        **r_lim**
        : number
            Parameter of the Uniform Distance prior.
        **f**
        : function
            Probability distribution function.
    ========
    Returns:
    ========
        **number**
        The mode of the distance modulus.
    """
    mu_max = 5*np.log10(r_lim) -5
                   
    if (-4*s**2 + w**2) > 0:
        
        sol = [(w - np.sqrt(-4*s**2 + w**2))/(2*s**2), (w + np.sqrt(-4*s**2 + w**2))/(2*s**2)]
        
        rmin = min(sol)
        rmax = max(sol)
        
        if  (5*np.log10(rmin)-5) >= 0 and (5*np.log10(rmin)-5) <= mu_max and (5*np.log10(rmax)-5) >= 0 and (5*np.log10(rmax)-5)<= mu_max:
            prob = [f(5*np.log10(sol[0])-5,w,s,r_lim),f(5*np.log10(sol[1])-5,w,s,r_lim)]
            
            
            if prob[0] != prob[1]:
                sol1 = [5*np.log10(sol[prob.index(max(prob))])-5,5*np.log10(r_lim) - 5]
                p1 = [f(sol1[0],w,s,r_lim),f(sol1[1],w,s,r_lim)]
                return sol1[p1.index(max(p1))]                
            else:
                return 5*np.log10(r_lim)-5
        
        elif 5*np.log10(rmin)-5>=0 and 5*np.log10(rmin)-5 <= mu_max:
            
            return 5*np.log10(rmin)-5
        
        elif 5*np.log10(rmax)-5>=0 and 5*np.log10(rmax)-5 <= mu_max:
            
            return 5*np.log10(rmax)-5
        
        else:
            
            return 5*np.log10(r_lim)-5     
    else:
        
        return 5*np.log10(r_lim)-5

def dmpdfexp(x,w,s,L):
    """Returns the probability of a given distance modulus using the Expo Prior.
     
     ===========
     Parameters:
     ===========
         **x**
         : number
             Real distance modulus.
         **w**
         : number
             Observed parallax in arcseconds.
         **s**
         : number
             Parallax error in arcseconds.
         **r_lim**
         : number
             Parameter of the Uniform Distance Prior in pc.
      ========
      Returns:
      ========
          **number**
              Unnormalized probability according to the Uniform Distance Posterior.
    """
    #print(x,w,s,L)
    lnp = np.log(1/(L**3*10*np.sqrt(2*np.pi)*s)) + (3*x+10)/5*np.log(10) - (1/(2*s**2)*(w-10**(-(x+5)/5))**2 +10**((x+5)/5)/L)
    #print(lnp)
    ##print(x,w,s,L,lnp)
    return np.exp(lnp)
    
def dmod_mode_exp(w,s,L):
    """Returns the  distance modulus mode of the Exponentially Decreasing Space Density Probabiliyt Distribution Function.
    
    This function solves the third degree equation ´x**3-3*L*x**2+l*w/s**2-L/s**2 = 0´ using numpy.roots and selects the root that represents the mode of the distribution.
    Given the physical limitations s>0 and L>0, only two cases are possible:
    - All three roots are real and there are two modes and one minimum. If w>0 the solution is the smallest one while if w<0, the greatest one is the solution.
    - There's only one real root.
    ===========
    Parameters:
    ===========
         **w**
        : number
            Observed parallax in arcseconds.
        **s**
        : number
            Parallax error in arcseconds.
        **L**
        : possitive number
            Length scale in pc. 
            
    ========
    Returns:
    ========
        **number**
            Mode oof the Exponentially Decreasing Space Density Probability Distribution Function
    """
    a = np.float64(1.)
    b = np.float64(-3*L)
    c = np.float64(L*w/s**2)
    d = np.float64(- L / s** 2)
    
    
    def function(x,a,b,c,d):
        return a*x**3+b*x**2+c*x+d
    
    p = np.array([a, b, c, d])

    sol = np.roots(p) #Solving the third degree equation
    #print('sol',sol)
    solutions = []
    
    for s in sol:
        
        if abs(s) == s:
            
            s2 = abs(s)
        
            solutions.append(s2)
            #print(solutions)
    if len(solutions) != 1 and len(sol) != 0:
        fx = [function(i,a,b,c,d) for i in sol]
        r = float(sol[fx.index(max(fx))])
        #print(5*np.log10(r)-5)
        return 5*np.log10(r)-5
        #return r
    elif len(solutions) == 1:
        r = float(solutions[0])
        #print(5*np.log10(r)-5)
        return 5*np.log10(r)-5
        #return r

   

#---------------------------------------------R*-------------------------------#
def mustar(beta,s,w):
    """Returns the estimate of distance modulus using the `Transformation Method`_ described by Haywood Smith.
    
    ===========
    Parameters:
    ===========
        **beta**
        : number
            Parameter of the estimate.
        **s**
        : number
            Parallax error in arcseconds.
        **w**
        : number
            Observed parallax in arcseconds.
    ========
    Returns:
    ========
        **number**
            Returns the estimate of distance modulus using the Transformatio Method.
    
    .. _Transformation Method:
        http://www.astro.ufl.edu/~asthcs/BadH2.pdf
    """
    if w > 0:
        wstar = beta*s/0.8*np.log(1+np.exp(0.8*w/s))
    else:
        wstar = beta*s/0.8*np.log(1+np.exp(0.8*w/s))*np.exp(-0.605*w**2/s**2)
    
    return 5*np.log10(1/wstar)-5


# -------------------------------------------------MEDIAN---------------------------------------------------------------#

def dmod_median(f,w,s,par,par2,r_mode,p,N,r0):
    """Returns the median distance modulus of the distribution f.
    
    Computes the median of the distribution solving the implicit equation `cdf(x)-0.5 = 0` through the scipy.optimize.brentq function, where `cdf(x)` is the normalized integration from r0 = 0.01 kpc to x using the scipy.integrate.quad function.
    
    ===========
    Parameters:
    ===========    
        **f**
        : function
            Probability distribution function.
        **w**
        :number
            Parameter of the function f in arcseconds.    
        **s**
        : number
            Parameter of the function f in arcseconds.
        **par**
        : number
            Parameter of the function f in pc.
        **par2**
        : number
            Superior distance modulus bound of the integration.
        **p** 
        : number
            Unnormalized percentile of r_mode.
        **N**
        : number
            Normalization factor.
        **r0**
        : number
            Inferior distance modulus bound of the integration.
            
    ========
    Returns:
    ========
        **number** 
            Median of the distribution in kpc.
    """
    if p == 0.5: # If the percentile of the mode = 0.5, the mode and the median coincide
        return r_mode
    
    else: # We try to find the smallest interval that we know contains the median
        
        if p>0.5: # If the mode is at the right of the median, the median will be between r0 and the mode
            r_inf = r0
            r_sup = r_mode
            
            def function(x):
                I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
                return I[0]/N - 0.5
            
            i = 0
            step = 0.1
            while True:
                i = i + 1
               
                try:
                    
                    r_median = brentq(function, r_inf, r_sup)
                   
                    return r_median
                
                except ValueError:
                    r_sup = r_sup + step
                    
                    
                if i == 10**3:
            
                    step = step/10
                    i = 0
                    r_sup = r_mode
            
            
        else: # If the mode is at the left of the median, the median will be between the mode and the limit  
            r_inf = r_mode
            r_sup = r_mode
            
            i = 0
            step = 0.1
            
            def function(x):
                I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
                return I[0]/N - 0.5
            
            while True:
                i = i + 1
                try:
                    
                    r_median = brentq(function, r_inf, r_sup)
                    
                    return r_median
                
                except ValueError:
                    r_sup = r_sup + step
                    
                if i == 10**3:
                  
                    step = step/10
                    i = 0
                    r_sup = r_mode



# -------------------------------------------DISTANCE BOUNDS-----------------------------------------------------------#

def distances_from_percentiles_dmod(f,N,w,s,stype,p_inf,p_sup,p_ref,r_ref,par,par2,r0):
    """Returns the distance modulus quantiles corresponding to the percentiles p_inf and p_sup.
    
    For stype = 'inf' finds the quantile corresponding to p_inf while otherwise finds  the quantile corresponding to p_sup. It solves the implicit equation `cdf(x)-p = 0` using the scipy.optimize.brentq rootfinding method, where `cdf(x)` is the normalized intregral of f from r0 to x.
    
    ===========
    Parameters:
    ===========
        **f**
        : function
            Probability distribution function.\n
        **N** 
        : number 
            Normalization factor.\n
        **w** 
        : number
            Parameter of the function f in arcseconds.\n    
        **s**
        : number
            Parameter of the function f in arcseconds.\n
        **stype**
        : str
            If stype = 'inf' computes the quantile corresponding to p_inf. Otherwise it computes the quantile corresponding to p_sup.\n
        **p_inf** 
        : float
            Inferior normalized percentile.\n 
        **p_sup**
        : float
            Superior normalized percentile.\n
        **p_ref**
        : float
            Reference percentile corresponding to r_ref.
        **r_ref**
        : float
            Distance reference in kpc corresponding to the percentile p_ref.
        **par**
        : number
            Parameter of f in pc.\n
        **par2**
        : number
            Distance modulus superior bound of the integration.\n
        **r0**
        : number
            Distance modulus inferior bound of the integration.
    ========
    Returns:
    ========
        **number**
            Quantile corresponding to the percentile p_inf or p_sup depending on stype.
    """
    if stype == 'inf': # I'm looking for the inferior distance bound
        r_inf = r0
        r_sub = r_ref # The inferior limit will be between 0 and the median
        
        def function(x):
            I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
            return I[0] / N - p_inf
        
        r = brentq(function,r_inf,r_sub)
        
        return r

    else: #I'm looking for the superior distance bound
       
        def function(x):
            I = quad(f, r0, x, args=(w, s, par),epsrel = 1e-11, epsabs = 1e-11)
            return I[0] / N - p_sup
        
        if (quad(f, r0, r_ref, args=(w, s, par), epsrel=1e-11, epsabs=1e-11)[0] / N - p_sup)<0:
            r_inf = r_ref
           
            
        else:
            r_inf = r0
           
        i = 0
        step = 10
        r_sup = r_ref
       
        while True:
            i = i + 1
            try:
                
                r = brentq(function, r_inf, r_sup)
                   
                return r
                
            except ValueError:
                r_sup = r_sup + step
                    
            if i == 10**3:
                   
                step = step/10
                i = 0
                r_sup = r_ref
               


def transformed_parallax(w,s,beta):
    """Returns the estimate of distance modulus using the `Transformation Method`_ described by Haywood Smith.
    
    ===========
    Parameters:
    ===========
        **beta**
        : number
            Parameter of the estimate.
        **s**
        : number
            Parallax error in arcseconds.
        **w**
        : number
            Observed parallax in arcseconds.
    ========
    Returns:
    ========
        **number**
            Returns the estimate of distance modulus using the Transformatio Method.
    
    .. _Transformation Method:
        http://www.astro.ufl.edu/~asthcs/BadH2.pdf
    """
    
    phi = np.log(1+np.exp(5*w/s))/5
    w_trans = beta*s*(1/(np.exp(phi)+np.exp(-5*w/s))+phi)
    return -(5*np.log(w_trans)+5)
    
     
        
        
    
    
def uncertainty_range_mdtm(w,s,beta,err):
    """Returns the uncertainty range associated to the estimate m* of the `Transformation Method`_.
        
    Returns the interval corresponding to [m*(w+err),m*(w-err)]. For instance, for err equal to the parallax error, it returns the 68.27 % confidence interval about m*.
    ===========    
    Parameters:
    ===========    
        **w** 
        : number
            Observed parallax in arcseconds.
        **s**
        : number
            Parallax error in arcseconds.
        **beta** 
        : number
            Parameter of the `Transformation Method`_.
        **err**
        : number
            Error interval.
            
    ========        
    Returns:
    ========    
        **tuple**
            Inferior and superior bound about the estimate m*.
    .. _Transformation Method:
        http://www.astro.ufl.edu/~asthcs/BadH2.pdf
    """
    
    phi_min = np.log(1+np.exp(5*(w-err)/s))/5
    phi_max = np.log(1+np.exp(5*(w+err)/s))/5
    
    w_trans_min = beta*s*(1/(np.exp(phi_min)+np.exp(-5*(w-err)/s))+phi_min)
    w_trans_max = beta*s*(1/(np.exp(phi_max)+np.exp(-5*(w+err)/s))+phi_max)
    
    return (-(5*np.log(w_trans_max)+5),-(5*np.log(w_trans_min)+5))

  
def plot_pdf(r0,r_lim,n_ud,n_exp,r_mode_ud,r_media_ud,r_ud_inf,r_ud_sup,r_mode_exp,r_media_exp,r_exp_inf,r_exp_sup,w,s,par1,par2,f1,f2):
    
    x = np.arange(r0,r_lim,r_lim/500)
    y1 = [f1(i,w,s,par1)/n_ud for i in x]
    y2 = [f2(j,w,s,par2)/n_exp for j in x]
    
    
    x_mode_ud = np.ones(20)*r_mode_ud
    dim_mode_ud = f1(r_mode_ud,w,s,par1)/n_ud
    y_mode_ud = np.arange(0,dim_mode_ud,dim_mode_ud/20)
    
    x_median_ud = np.ones(20)*r_media_ud
    dim_media_ud = f1(r_media_ud,w,s,par1)/n_ud
    y_media_ud = np.arange(0,dim_media_ud,dim_media_ud/20)
    
    x_uncertainty_ud = np.arange(r_ud_inf,r_ud_sup,1)
    y_uncertainty_ud = [f1(k,w,s,par1)/n_ud for k in x_uncertainty_ud]
    
    x_mode_exp = np.ones(20)*r_mode_exp
    dim_mode_exp = f2(r_mode_exp,w,s,par2)/n_exp
    y_mode_exp = np.arange(0,dim_mode_exp,dim_mode_exp/20)
    
    x_median_exp = np.ones(20)*r_media_exp
    dim_media_exp = f2(r_media_exp,w,s,par2)/n_exp
    y_media_exp = np.arange(0,dim_media_exp,dim_media_exp/20)
    
    x_uncertainty_exp = np.arange(r_exp_inf,r_exp_sup,1)
    y_uncertainty_exp = [f2(k,w,s,par2)/n_exp for k in x_uncertainty_exp]
    
    
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=False)
    f.set_size_inches(18.5, 10.5)
    plt.rc('axes', titlesize=25)     
    plt.rc('axes', labelsize=15)
    plt.rc('font',size = 15)
    
    ax1.plot(x, y1)
    ax1.plot(x_mode_ud,y_mode_ud,label = 'Mode',color = "red")
    ax1.plot(x_median_ud,y_media_ud,label = 'Median', color = "black")
    ax1.fill_between(x_uncertainty_ud,0,y_uncertainty_ud,facecolor = 'blue',alpha=0.5,label='90% uncertainty interval')
    ax1.set_title('Uniform Distance PDF')
    ax1.set_xlabel('True distance r (kpc)')
    ax1.set_ylabel('PDF')
    ax1.legend()
    ax2.plot(x, y2)
    ax2.plot(x_mode_exp,y_mode_exp,label = 'Mode',color='red')
    ax2.plot(x_median_exp,y_media_exp,label='Median',color='black')
    ax2.fill_between(x_uncertainty_exp,0,y_uncertainty_exp,color = 'blue',alpha=0.5,label='90% uncertainty interval')
    ax2.set_title('Exponentially Decreasing Space Density PDF')
    ax2.set_xlabel('True distance r (kpc)')
    ax2.set_ylabel('PDF')
    ax2.legend()
    plt.show()
    