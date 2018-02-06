# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 20:03:41 2017

@author: Ariadna
"""

"""Main program using library for the interface application."""

from pyrallaxes import *
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from scipy.optimize import fsolve
from decimal import *
import matplotlib.pyplot as plt

r_lim = np.float64(100) #kpc
r0 = 0.01 #kpc
sigma = 0.341
L = 1.35 #kpc
tol = 10**(-25)
p5 = 0.05
p95 = 0.95
beta = 1.01

def main_un(w,s):
    r_mode_ud = mode_r_uniform(w, r_lim) # Computing the mode of the PDF
                    
    if quad(uniform_distance_posterior, r0, r_lim, args=(w, s, r_lim), epsrel=1e-8, epsabs=1e-8)[0] <= 1e-15 or (s/w >0 and s/w<0.05):
        p_ud = percentiles(uniform_distance_posterior,0.5*r_mode_ud,r_mode_ud,w,s,r_lim)
        n_ud = normalization(uniform_distance_posterior, r_lim, w, s, p_ud, r_mode_ud, 1.5*r_mode_ud)
        p_ud = normalized_percentile(p_ud, n_ud)
        r_media_ud = median(uniform_distance_posterior, w, s, r_lim, 4*r_mode_ud, r_mode_ud, p_ud, n_ud)
        
        r_ud_95 = distances_from_percentiles(uniform_distance_posterior, n_ud, w, s, 'sup', p5, p95, 0.5, r_media_ud, r_lim,1.5*r_mode_ud)
        r_ud_5 = distances_from_percentiles(uniform_distance_posterior, n_ud, w, s, 'inf', p5, p95, 0.5, r_media_ud, r_lim,1.5*r_mode_ud)

    else:
        p_ud = percentiles(uniform_distance_posterior, r0, r_mode_ud, w, s,r_lim)  # Computing the percentile that corresponds to the mode of the PDF
        n_ud = normalization(uniform_distance_posterior, r_lim, w, s, p_ud, r_mode_ud,r_lim)  # Computing the normalization constant of the PDF
        p_ud = normalized_percentile(p_ud, n_ud)  # The percentile of the mode have been normalized p in [0,1]
                
        r_media_ud = median(uniform_distance_posterior,w,s,r_lim,r_lim,r_mode_ud,p_ud,n_ud)
        
        r_ud_95 = distances_from_percentiles(uniform_distance_posterior, n_ud, w, s, 'sup', p5, p95,0.5, r_media_ud, r_lim, r_lim)
        r_ud_5 = distances_from_percentiles(uniform_distance_posterior,n_ud,w,s,'inf',p5,p95,0.5,r_media_ud,r_lim,r_lim)
    
    return r_ud_5*1000,r_mode_ud*1000,r_media_ud*1000,r_ud_95*1000,n_ud

def main_exp(w,s):
   
    r_mode_exp = mode_r_exponential(L, w, s)
                    
    if quad(exponentially_decreasing_space_density_posterior, r0, 1000*L, args=(w, s, L), epsrel=1e-8, epsabs=1e-8)[0] <= 1e-15 or (s/w >0 and s/w<0.05):
       
        p_exp = percentiles(exponentially_decreasing_space_density_posterior,-2*r_mode_exp,r_mode_exp, w, s, L)
        n_exp = normalization(exponentially_decreasing_space_density_posterior, L, w, s, p_exp, r_mode_exp, 1.5*r_mode_exp)
        p_exp = normalized_percentile(p_exp, n_exp)
                    
        r_media_exp = median(exponentially_decreasing_space_density_posterior, w, s, L, 4*r_mode_exp, r_mode_exp, p_exp, n_exp)
            
        r_exp_5 = distances_from_percentiles(exponentially_decreasing_space_density_posterior, n_exp, w, s, 'inf', p5, p95, 0.5, r_media_exp, L,1.5*r_mode_exp)
        r_exp_95 = distances_from_percentiles(exponentially_decreasing_space_density_posterior, n_exp, w, s, 'sup', p5, p95, 0.5, r_media_exp, L,1.5*r_mode_exp)
                    
    else:
        
        p_exp = percentiles(exponentially_decreasing_space_density_posterior,r0,r_mode_exp,w,s,L)
        n_exp = normalization(exponentially_decreasing_space_density_posterior,L,w,s,p_exp,r_mode_exp,r_lim)
        p_exp = normalized_percentile(p_exp,n_exp)
                
        r_media_exp = median(exponentially_decreasing_space_density_posterior,w,s,L,r_lim,r_mode_exp,p_exp,n_exp)
        r_exp_5 = distances_from_percentiles(exponentially_decreasing_space_density_posterior, n_exp, w, s, 'inf', p5, p95,0.5, r_media_exp, L, r_lim)
        r_exp_95 = distances_from_percentiles(exponentially_decreasing_space_density_posterior, n_exp, w, s, 'sup', p5, p95,0.5, r_media_exp, L, r_lim)

    return r_exp_5*1000,r_mode_exp*1000,r_media_exp*1000,r_exp_95*1000,n_exp

def main_trans(w,s):
    stm = s*10**(-3) # Parallax error in arcseconds
    wtm = w*10**(-3) # Parallax in arcseconds
    err = 2*stm
   
    rs = rstar(beta, stm, wtm)
    rsmin,rsmax = uncertainty_range_tm(wtm,stm,beta,err)
    
    return rsmin,rs,rsmax

def main_mun(w,s):
    mu0 = -10
    r_lim = 100*1000 #pc
    mu_lim = 5*np.log(r_lim)-5
    w = w/1000 #arcseconds
    s = s/1000 #arcseconds
    m_mode_un = dmod_mode_ud(w,s,r_lim,dmpdfun)

    if quad(dmpdfun,mu0,m_mode_un,args=(w, s, r_lim),epsrel = 1e-11, epsabs = 1e-11)[0] <= 2.5e-6:
        p_mud = percentiles(dmpdfun,m_mode_un-0.01,m_mode_un,w,s,r_lim) 
        
        n_mud = normalization(dmpdfun, r_lim, w, s, p_mud, m_mode_un, mu_lim)
       
        p_mud = normalized_percentile(p_mud, n_mud)
        
    
        r_media_ud = dmod_median(dmpdfun, w, s, r_lim,m_mode_un+0.01, m_mode_un, p_mud, n_mud,m_mode_un-0.01)
        
        r_ud_95 = distances_from_percentiles_dmod(dmpdfun, n_mud, w, s, 'sup', 0.05, 0.95, 0.5, r_media_ud, r_lim,m_mode_un+0.01,m_mode_un-0.01)
        r_ud_5 = distances_from_percentiles_dmod(dmpdfun, n_mud, w, s, 'inf', 0.05, 0.95, 0.5, r_media_ud, r_lim,m_mode_un+0.01,m_mode_un-0.01)

            
    else:
                    
        p_ud = percentiles(dmpdfun,mu0,m_mode_un,w,s,r_lim)
       
        n_mud = normalization(dmpdfun, r_lim, w, s, p_ud, m_mode_un, mu_lim)
       
        p_ud = normalized_percentile(p_ud, n_mud)
       
        
        r_media_ud = dmod_median(dmpdfun, w, s, r_lim,mu_lim, m_mode_un, p_ud, n_mud,mu0)
       
        r_ud_95 = distances_from_percentiles_dmod(dmpdfun, n_mud, w, s, 'sup', 0.05, 0.95, 0.5, r_media_ud, r_lim,mu_lim,mu0)
        
        r_ud_5 = distances_from_percentiles_dmod(dmpdfun, n_mud, w, s, 'inf', 0.05, 0.95, 0.5, r_media_ud, r_lim,mu_lim,mu0)
        
    return r_ud_5,m_mode_un,r_media_ud,r_ud_95,n_mud
    
def main_mexp(w,s):
    mu0 = -10
    r_lim = 100*1000 #pc
    mu_lim = 5*np.log(r_lim)-5
    L=1.35*1000 #pc
    w = w/1000 #arcseconds
    s = s/1000 #arcseconds
    m = dmod_mode_exp(w,s,L)
    
    p_exp = percentiles(dmpdfexp,m - 10,m, w, s, L)

    n_exp = normalization(dmpdfexp, L, w, s, p_exp, m, m + 10)
    p_exp = normalized_percentile(p_exp, n_exp)
    
    
    r_media_exp = dmod_median(dmpdfexp, w, s, L, m + 10, m, p_exp, n_exp,mu0) 
                 
    r_exp_5 = distances_from_percentiles_dmod(dmpdfexp, n_exp, w, s, 'inf', 0.05, 0.95, 0.5, r_media_exp, L,m + 10,mu0)
   
    r_exp_95 = distances_from_percentiles_dmod(dmpdfexp, n_exp, w, s, 'sup', 0.05, 0.95, 0.5, r_media_exp, L,m + 10,mu0)
    
    
    return r_exp_5,m,r_media_exp,r_exp_95,n_exp

def main_m_trans(w,s):  
    w = w/1000 #arcseconds
    s = s/1000 #arcseconds
    beta = 1.01
    err = 2*s
    mtm = transformed_parallax(w,s,beta)
    mtm_min,mtm_max = uncertainty_range_mdtm(w,s,beta,err)
    return mtm_min,mtm,mtm_max
    

        