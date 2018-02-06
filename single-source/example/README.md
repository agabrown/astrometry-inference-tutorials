#DISTANCE ESTIMATION TOOL
This tool estimates the distance and the distance 
modulus from the trigonometric parallaxes using two
Bayesian Methods with the priors suggested by Coryn 
Bayler Jones (http://iopscience.iop.org/article/
10.1086/683116/pdf) and a frequentist method called 
the Transformation Method described by Haywood Smith 
(http://www.astro.ufl.edu/~asthcs/BadH2.pdf).

This folder contains three executable files __.py,

* Interface: Program to execute in order to obtain
	the graphic interface.

* Library: Set of functions used by the main program.
	All the functions are commented with their 
	inputs, outputs and usage.

* Main: Calculations of the mode, the median and the 
	90% uncertainty interval of confidence using 
	the Bayesian Methods and the r* and m*
	estimates used in the Transformation Method
	with the 90% uncertainty interval.

In order to execute the program it is necessary to 
have installed python 3.6 and to execute the file
interface.py.

Please, take into account that this version is still
being tested and improvements are being developed. 
