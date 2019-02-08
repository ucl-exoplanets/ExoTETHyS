import astropy.io.fits as pyfits
import numpy as np
import os
from scipy.optimize import minimize


def claret4(params, mucut, intscut, weights):
	#make fitting parameters explicit
	c1 = params[0]
	c2 = params[1]
	c3 = params[2]
	c4 = params[3]
	model = 1.0 - c1*(1.0-mucut**0.5) - c2*(1.0-mucut) - c3*(1.0-mucut**1.5) - c4*(1.0-mucut**2.0)
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def power2(params, mucut, intscut, weights):
	#make fitting parameters explicit
	u = params[0]
	a = params[1]
	model = 1.0 - u*(1.0-mucut**a)
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def square_root(params, mucut, intscut, weights):
	#make fitting parameters explicit
	c1 = params[0]
	c2 = params[1]
	model = 1.0 - c1*(1.0-mucut**0.5) - c2*(1.0-mucut)
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def quadratic(params, mucut, intscut, weights):
	#make fitting parameters explicit
	g1 = params[0]
	g2 = params[1]
	model = 1.0 - g1*(1.0-mucut) - g2*(1.0-mucut)**2.0
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def linear(params, mucut, intscut, weights):
	#make fitting parameters explicit
	u = params[0]
	model = 1.0 - u*(1.0-mucut)
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def gen_poly(params, mucut, intscut, weights):
	#make fitting parameters explicit
	g = params.copy()
	n_params = len(g)
	model = 1.0
	for n in range(n_params):
		model -= g[n]*(1.0-mucut**(n+1.0))
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)


def gen_claret(params, mucut, intscut, weights):
	#make fitting parameters explicit
	c = params.copy()
	n_params = len(c)
	model = 1.0
	for n in range(n_params):
		model -= c[n]*(1.0-mucut**((n+1.0)/2))
	return np.sum(weights*((intscut-model)**2)) / np.sum(weights)



limb_darkening_laws_func = {'claret4':claret4, 'square_root':square_root, 'quadratic':quadratic, 'linear':linear, 'gen_claret':gen_claret, 'gen_poly':gen_poly}








