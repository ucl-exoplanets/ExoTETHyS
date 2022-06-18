from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os

#import matplotlib
#if os.environ.get('DISPLAY', '') == '':
#    print('no display found. Using non-interactive Agg backend')
#    matplotlib.use('Agg')
#else:
#    matplotlib.use('TkAgg')

import pkg_resources
import numpy as np
from scipy.optimize import minimize
from collections import OrderedDict
from scipy.interpolate import interp1d
from scipy.integrate import simps
import astropy.units as u

import copy
import pickle

from ._database import Database, databases, sys, urlretrieve, glob, time, shutil

def get_intensities_from_ldcs(mu, coefficients, law, coeffs_type='original', i_rot_angle=None, finite=True, positive=True, monotonically_decreasing=True):
    """
    This function computes the model intensities given the limb-darkening coefficients
    
    :param np.array mu: 1D array with mu values
    :param np.array coefficients: 1D array with the limb-darkening coefficients
    :param str law: name of the limb-darkening law
    :return: the intensities at the given mu values
    :rtype: np.array
    """
    mu = np.atleast_2d( mu )
    coefficients = np.atleast_2d( coefficients )
    #print(mu)
    if np.shape(mu)[0] not in np.shape(coefficients):
        print('ERROR: [In SAIL.get_intensities_from_ldcs] The given mu values and limb-darkening coefficients have incompatible formats.')
        return
    if law=='claret4':
        check_ldcs = check_claret4_ldcs( coefficients, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some claret4 limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        c1, check = make_same_format_2D(coefficients[:,0], mu, first_array='column')
        c2, check = make_same_format_2D(coefficients[:,1], mu, first_array='column')
        c3, check = make_same_format_2D(coefficients[:,2], mu, first_array='column')
        c4, check = make_same_format_2D(coefficients[:,3], mu, first_array='column')
        model = 1.0 - c1*(1.0-mu**0.5) - c2*(1.0-mu) - c3*(1.0-mu**1.5) - c4*(1.0-mu**2.0)
    elif law=='power2':
        check_ldcs = check_power2_ldcs( coefficients, coeffs_type=coeffs_type, finite=finite, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some power2 limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        c_alpha = transform_power2_ldcs( coefficients, coeffs_type, 'original' )
        c, check = make_same_format_2D(c_alpha[:,0], mu, first_array='column')
        alpha, check = make_same_format_2D(c_alpha[:,1], mu, first_array='column')
        #CHECK FORMAT OF C AND ALPHA (ALSO ADJUST transform_power2_ldcs). THE LINE BELOW SHOW HOW IT WORKS WITH (3, 2) FOR c AND (3, 11) FOR mu.
        #prov = 1 - (1-mu**np.transpose(np.atleast_2d(c[:,1])))*np.transpose(np.atleast_2d(c[:,0]))
        model = 1.0 - c*(1.0-mu**alpha)
    elif law=='square_root':
        check_ldcs = check_square_root_ldcs( coefficients, coeffs_type=coeffs_type, i_rot_angle=i_rot_angle, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some square root limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        c1_c2 = transform_square_root_ldcs( coefficients, coeffs_type, 'original', i_rot_angle=i_rot_angle )
        c1, check = make_same_format_2D(c1_c2[:,0], mu, first_array='column')
        c2, check = make_same_format_2D(c1_c2[:,1], mu, first_array='column')
        model = 1.0 - c1*(1.0-mu**0.5) - c2*(1.0-mu)
    elif law=='quadratic':
        check_ldcs = check_quadratic_ldcs( coefficients, coeffs_type=coeffs_type, i_rot_angle=i_rot_angle, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some quadratic limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        u1_u2 = transform_quadratic_ldcs( coefficients, coeffs_type, 'original', i_rot_angle=i_rot_angle )
        u1, check = make_same_format_2D(u1_u2[:,0], mu, first_array='column')
        u2, check = make_same_format_2D(u1_u2[:,1], mu, first_array='column')
        model = 1.0 - u1*(1.0-mu) - u2*(1.0-mu)**2.0
    elif law=='linear':
        check_ldcs = check_linear_ldcs( coefficients, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some linear limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        u, check = make_same_format_2D(coefficients, mu, first_array='column')
        model = 1.0 - u*(1.0-mu)
    elif law=='gen_poly':
        nc = np.shape(coefficients)[1]
        check_ldcs = check_gen_poly_ldcs( coefficients, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some gen_poly', nc, 'limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        model = 1.0
        for n in range(nc):
            cn, check = make_same_format_2D(coefficients[:,n], mu, first_array='column')
            model -= cn*(1.0-mu**(n+1.0))
    elif law=='gen_claret':
        nc = np.shape(coefficients)[1]
        check_ldcs = check_gen_claret_ldcs( coefficients, positive=positive, monotonically_decreasing=monotonically_decreasing)
        if not check_ldcs.all():
            indices = list( np.where( not check_ldcs )[0] )
            print('WARNING: [In SAIL.get_intensities_from_ldcs] Some gen_claret', nc, 'limb-darkening coefficients do not comply with the requirements. The corresponding indices are', indices, '.')
        model = 1.0
        for n in range(nc):
            cn, check = make_same_format_2D(coefficients[:,n], mu, first_array='column')
            model -= cn*(1.0-mu**((n+1.0)/2))
    return  model


def check_claret4_ldcs(params, positive=True, monotonically_decreasing=True):
    """
    This function checks if claret4 limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    c1 = params[:,0]
    c2 = params[:,1]
    c3 = params[:,2]
    c4 = params[:,3]
    check1 = np.isfinite(c1)
    check2 = np.isfinite(c2)
    check3 = np.isfinite(c3)
    check4 = np.isfinite(c4)
    cond = np.column_stack([ check1, check2, check3, check4 ])
    if positive:
        cond1 = (c1 + c2 + c3 + c4 <= 1.0) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = (c1 >= 0.0) #monotonically decreasing
        cond3 = (c1 + 2.0*c2 + 3.0*c3 + 4.0*c4 >= 0.0) #monotonically decreasing
        cond = np.column_stack([cond, cond2, cond3])
    cond = cond.all(axis=1)
    return cond


def get_wrms_claret4(params, mucut, intscut, weights, coeffs_type='original', positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with claret4 limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - c1*(1-mu^(1/2)) - c2*(1-mu) - c3*(1-mu^(3/2)) - c4*(1-mu^2)
    
    :param np.array params: 1D array with four limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_claret4] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_claret4] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:4], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_claret4] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_claret4_ldcs( params, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'claret4')
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)


def transform_power2_ldcs( params, input_type, output_type ):
    """
    This function converts between various formats of power2 limb-darkening coefficients:
    original, (c, alpha) as defined by Hestroffer 1997, A&A, 327, 199
    short_2019, (q1, q2) as defined by Short et al. 2019, RNAAS, 3, 117
    maxted_2018, (h1, h2) as defined by Maxted 2018, A&A, 616, A39
    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :param str output_type: requested format for the limb-darkening coefficients
    :return: the transformed limb-darkening coefficients with the requested output type
    :rtype: float, float
    """
    params = np.asarray( params )
    params = np.atleast_2d( params )
    check = check_2Darray(params, n_col=2)
    if not check:
        print('ERROR: wrong format for the power2 limb-darkening coefficients. It should be an Nx2 array or list of arrays.')
        return
    ldc1 = params[:,0]
    ldc2 = params[:,1]
    #ldc1 = np.atleast_2d( ldc1 )
    #ldc2 = np.atleast_2d( ldc2 )
    #ldc1 = np.transpose( ldc1 )
    #ldc2 = np.transpose( ldc2 )
    #print(ldc1, ldc2)
    if input_type=='original':
        if output_type=='original':
            return np.column_stack([ ldc1, ldc2 ])
        elif output_type=='short_2019':
            q1 = ( 1.0 - ldc1*(2.0**(-ldc2)) )**2
            #q2 = np.random.rand(q1.size)
            q2 = np.zeros_like(q1)
            indices = np.where(q1!=0)[0]
            q2[indices] = (1.0 - ldc1[indices])/np.sqrt(q1[indices])
            return np.column_stack([ q1, q2 ])
        elif output_type=='maxted_2018':
            h1 = 1.0 - ldc1*(1.0-2.0**(-ldc2))
            h2 = ldc1*(2.0**(-ldc2))
            return np.column_stack([ h1, h2 ])
    elif output_type=='original':
        if input_type=='short_2019':
            c = 1.0 - ldc2*np.sqrt(ldc1)
            #alpha = np.random.rand(c.size)
            alpha = np.zeros_like(c)
            indices = np.where(c!=0)
            alpha[indices] = np.log2( c[indices] / (1.0 - np.sqrt(ldc1[indices])) )
            return np.column_stack([ c, alpha ])
        elif input_type=='maxted_2018':
            c = 1.0 - ldc1 + ldc2
            #alpha = np.random.rand(c.size)
            alpha = np.zeros_like(c)
            indices = np.where(c!=0)
            alpha[indices] = np.log2(c[indices]/ldc2[indices])
            return np.column_stack([ c, alpha ])
    else:
        c_alpha = transform_power2_ldcs( params, input_type, 'original' )
        ldcs_tr = transform_power2_ldcs( c_alpha, 'original', output_type )
        return ldcs_tr


#def check_power2_ldcs(params, input_type, warning=True, correct=True):
#    """
#    This function checks if square root limb-darkening coefficients comply with the requirements of
#    (A) positive, (B) monotonically decreasing profile,
#    using the transformation defined by Short et al. 2019, RNAAS, 3, 117.
#    By default, a correction is applied to get valid coefficients.

    
#    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
#    :param str input_type: format of the given limb-darkening coefficients
#    :argument bool correct: if True, apply correction to get valid limb-darkening coefficients
#    :argument bool warning: if True, print warnings
#    :return: if correct, the corrected limb-darkening coefficients, otherwise the same as the input ones
#    :rtype: float, float
#    """
#    q_ldcs = transform_power2_ldcs( params, input_type, 'short_2019' )
#    q_ldcs = np.asarray( q_ldcs )
#    if np.min(q_ldcs)<0.0 or np.max(q_ldcs)>1.0:
#        if correct:
#            q_ldcs[ np.where( q_ldcs<0.0 ) ] = 0.0
#            q_ldcs[ np.where( q_ldcs>1.0 ) ] = 1.0
#            if warning:
#                print('WARNING: There are unphysical power2 limb-darkening coefficients. Fixing them automatically.')
#        elif warning:
#            print('WARNING: There are unphysical power2 limb-darkening coefficients. No corrections applied.')
#    ldc1, ldc2 = transform_power2_ldcs( q_ldcs, 'short_2019', input_type )
#    return ldc1, ldc2

def check_power2_ldcs(params, coeffs_type='original', finite=True, positive=True, monotonically_decreasing=True):
    """
    This function checks if power2 limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    c_alpha = transform_power2_ldcs( params, coeffs_type, 'original' )
    c = c_alpha[:,0]
    alpha = c_alpha[:,1]
    check1 = np.isfinite(c)
    check2 = np.isfinite(alpha)
    cond = np.column_stack([ check1, check2])
    if finite:
        cond1 = (alpha >= 0.0)
        cond = np.column_stack([cond, cond1])
    if positive:
        cond2 = (c <= 1.0) #positive
        cond = np.column_stack([cond, cond2])
    if monotonically_decreasing:
        cond3 = (c >= 0.0) #monotonically decreasing
        cond = np.column_stack([cond, cond2])
    cond = cond.all(axis=1)
    return cond


def get_wrms_power2(params, mucut, intscut, weights, coeffs_type='original', finite=True, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with power2 limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - c*(1-mu^a)
    
    :param np.array params: 1D array with two limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :argument str coeffs_type: format of the given limb-darkening coefficients (default is 'original')
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_power2] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_power2] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:2], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_power2] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_power2_ldcs( params, coeffs_type=coeffs_type, finite=finite, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'power2', coeffs_type=coeffs_type)
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)



def make_same_format_2D(arr1, arr2, first_array='row'):
    arr1 = np.atleast_2d( arr1 )
    arr2 = np.atleast_2d( arr2 )
    check = True
    #check = check_2Darray( arr1 ) and check_2Darray( arr2 )
    #if not check:
    #    return arr1, check
    try:
        r1, c1 = np.shape( arr1 )
        r2, c2 = np.shape( arr2 )
    except:
        check = False
        return arr1, check 
    if r1==r2 and c1==c2:
        return arr1, check
    elif r1==1 and c1==c2 and first_array in ['row', 'whatever']:
        arr1 = np.tile( arr1, (r2,1) )
        return arr1, check
    elif r1==c2 and c1==1 and first_array in ['row', 'whatever']:
        arr1 =  np.transpose( arr1 )
        arr1 = np.tile( arr1, (r2,1) )
        return arr1, check
    elif r1==r2 and c1==1 and first_array in ['column', 'whatever']:
        arr1 = np.tile( arr1, c2 )
        return arr1, check
    elif r1==1 and c1==r2 and first_array in ['column', 'whatever']:
        arr1 =  np.transpose( arr1 )
        arr1 = np.tile( arr1, c2 )
        return arr1, check
    elif r1==1 and c1==1 and first_array=='whatever':
        arr1 = np.tile( arr1, (r2, c2) )
        return arr1, check
    else:
        check = False
        return arr1, check
   


def transform_square_root_ldcs( params, input_type, output_type, i_rot_angle=None, o_rot_angle=None ):
    """
    This function converts between various formats of square root limb-darkening coefficients:
    original, (c1, c2) as defined by Diaz-Cordoves & Gimenez 1992, A&A, 259, 227
    square01_uniform, (q1, q2) as defined by Kipping 2013, MNRAS, 435, 2152
    halfsquare01_uniform, (v1, v2) as defined by Kipping 2013, MNRAS, 435, 2152
    sum_diff, (c1+c2, c1-c2)
    2sum1_1diff2, (2*c1+c2, c1-2*c2) 
    rotated, by applying rotation to the (c1, c2) vector
    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :param str output_type: requested format for the limb-darkening coefficients
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :argument quantity o_rot_angle: angle for the rotated limb-darkening coefficients, requested as output, with respect to the original ones (default is None)
    :return: the transformed limb-darkening coefficients with the requested output type
    :rtype: float, float
    """
    params = np.asarray( params )
    params = np.atleast_2d( params )
    check = check_2Darray(params, n_col=2)
    if not check:
        print('ERROR: wrong format for the square root limb-darkening coefficients. It should be an Nx2 array or list of arrays.')
        return
    ldc1 = params[:,0]
    ldc2 = params[:,1]
    if input_type=='original':
        if output_type=='original':
            return np.column_stack([ ldc1, ldc2 ])
        elif output_type=='square01_uniform':
            q1 = ( ldc1 + ldc2 )**2.0
            q2 = np.zeros_like(q1)
            indices = np.where(q1>0)
            q2[indices] = 0.5*ldc1[indices]/q1[indices]
            return np.column_stack([ q1, q2 ])
        elif output_type=='halfsquare01_uniform':
            v1 = 0.5*ldc1
            v2 = 1.0 - ldc1 - ldc2
            return np.column_stack([ v1, v2 ])
        elif output_type=='sum_diff':
            cplus = ldc1 + ldc2
            cminus = ldc1 - ldc2
            return np.column_stack([ cplus, cminus ])
        elif output_type=='2sum1_1diff2':
            c2sum1 = 2.0*ldc1 + ldc2
            c1diff2 = ldc1 - 2.0*ldc2
            return np.column_stack([ c2sum1, c1diff2 ])
        elif output_type=='rotated':
            cos_t = np.cos(o_rot_angle).value
            sin_t = np.sin(o_rot_angle).value
            w1 = ldc1*cos_t - ldc2*sin_t
            w2 = ldc1*sin_t + ldc2*cos_t
            return np.column_stack([ w1, w2 ])
    elif output_type=='original':
        if input_type=='square01_uniform':
            c1 = 2.0*np.sqrt(ldc1)*ldc2
            c2 = np.sqrt(ldc1)*(1.0-2.0*ldc2)
            return np.column_stack([ c1, c2 ])
        elif input_type=='halfsquare01_uniform':
            c1 = 2.0 * ldc1
            c2 = 1.0 - 2.0*ldc1 - ldc2
            return np.column_stack([ c1, c2 ])
        elif input_type=='sum_diff':
            c1 = 0.5 * ( ldc1 + ldc2 )
            c2 = 0.5 * ( ldc1 - ldc2 )
            return np.column_stack([ c1, c2 ])
        elif input_type=='2sum1_1diff2':
            c1 = 0.2 * ( 2.0*ldc1 + ldc2 )
            c2 = 0.2 * ( ldc1 - 2.0*ldc2 )
            return np.column_stack([ c1, c2 ])
        elif input_type=='rotated':
            cos_t = np.cos(-i_rot_angle).value
            sin_t = np.sin(-i_rot_angle).value
            c1 = ldc1*cos_t - ldc2*sin_t
            c2 = ldc1*sin_t + ldc2*cos_t
            return np.column_stack([ c1, c2 ])
    else:
        c1_c2 = transform_square_root_ldcs( params, input_type, 'original', i_rot_angle=i_rot_angle )
        ldcs_tr = transform_square_root_ldcs( c1_c2, 'original', output_type, o_rot_angle=o_rot_angle )
        return ldcs_tr


def check_square_root_ldcs(params, coeffs_type='original', i_rot_angle=None, positive=True, monotonically_decreasing=True):
    """
    This function checks if square root limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    c1_c2 = transform_square_root_ldcs( params, coeffs_type, 'original', i_rot_angle=i_rot_angle )
    c1 = c1_c2[:,0]
    c2 = c1_c2[:,1]
    check1 = np.isfinite(c1)
    check2 = np.isfinite(c2)
    cond = np.column_stack([ check1, check2])
    if positive:
        cond1 = (c1 + c2 <= 1) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = (c1 >= 0) #monotonically decreasing
        cond3 = (c1 + 2.0*c2 >= 0) #monotonically decreasing
        cond = np.column_stack([cond, cond2, cond3])
    cond = cond.all(axis=1)
    return cond



def get_wrms_square_root(params, mucut, intscut, weights, coeffs_type='original', i_rot_angle=None, o_rot_angle=None, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with square-root limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - c1*(1-mu^(1/2)) - c2*(1-mu)
    
    :param np.array params: 1D array with two limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :argument str coeffs_type: format of the given limb-darkening coefficients (default is 'original')
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_square_root] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_square_root] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:2], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_square_root] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_square_root_ldcs( params, coeffs_type=coeffs_type, i_rot_angle=i_rot_angle, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'square_root', coeffs_type=coeffs_type, i_rot_angle=i_rot_angle)
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)



def transform_quadratic_ldcs( params, input_type, output_type, i_rot_angle=None, o_rot_angle=None ):
    """
    This function converts between various formats of quadratic limb-darkening coefficients:
    original, (u1, u2) as defined by Kopal 1950, HarCi, 454, 1
    square01_uniform, (q1, q2) as defined by Kipping 2013, MNRAS, 435, 2152
    halfsquare01_uniform, (v1, v2) as defined by Kipping 2013, MNRAS, 435, 2152
    sum_diff, (u+, u-) as defined by Brown et al. 2001, 
    2sum1_1diff2, (2*u1+u2, u1-2*u2) as defined by Holman et al. 2006, 
    rotated, as proposed by Pal 2008, 
    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :param str output_type: requested format for the limb-darkening coefficients
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :argument quantity o_rot_angle: angle for the rotated limb-darkening coefficients, requested as output, with respect to the original ones (default is None)
    :return: the transformed limb-darkening coefficients with the requested output type
    :rtype: float, float
    """
    params = np.asarray( params )
    params = np.atleast_2d( params )
    check = check_2Darray(params, n_col=2)
    if not check:
        print('ERROR: wrong format for the quadratic limb-darkening coefficients. It should be an Nx2 array or list of arrays.')
        return
    ldc1 = params[:,0]
    ldc2 = params[:,1]
    if input_type=='original':
        if output_type=='original':
            return np.column_stack([ ldc1, ldc2 ])
        elif output_type=='square01_uniform':
            q1 = ( ldc1 + ldc2 )**2.0
            q2 = np.zeros_like(q1)
            indices = np.where(q1>0)
            q2[indices] = 0.5*ldc1[indices]/q1[indices]
            return np.column_stack([ q1, q2 ])
        elif output_type=='halfsquare01_uniform':
            v1 = 0.5*ldc1
            v2 = 1.0 - ldc1 - ldc2
            return np.column_stack([ v1, v2 ])
        elif output_type=='sum_diff':
            uplus = ldc1 + ldc2
            uminus = ldc1 - ldc2
            return np.column_stack([ uplus, uminus ])
        elif output_type=='2sum1_1diff2':
            u2sum1 = 2.0*ldc1 + ldc2
            u1diff2 = ldc1 - 2.0*ldc2
            return np.column_stack([ u2sum1, u1diff2 ])
        elif output_type=='rotated':
            cos_t = np.cos(o_rot_angle).value
            sin_t = np.sin(o_rot_angle).value
            w1 = ldc1*cos_t - ldc2*sin_t
            w2 = ldc1*sin_t + ldc2*cos_t
            return np.column_stack([ w1, w2 ])
    elif output_type=='original':
        if input_type=='square01_uniform':
            u1 = 2.0*np.sqrt(ldc1)*ldc2
            u2 = np.sqrt(ldc1)*(1.0-2.0*ldc2)
            return np.column_stack([ u1, u2 ])
        elif input_type=='halfsquare01_uniform':
            u1 = 2.0 * ldc1
            u2 = 1.0 - 2.0*ldc1 - ldc2
            return np.column_stack([ u1, u2 ])
        elif input_type=='sum_diff':
            u1 = 0.5 * ( ldc1 + ldc2 )
            u2 = 0.5 * ( ldc1 - ldc2 )
            return np.column_stack([ u1, u2 ])
        elif input_type=='2sum1_1diff2':
            u1 = 0.2 * ( 2.0*ldc1 + ldc2 )
            u2 = 0.2 * ( ldc1 - 2.0*ldc2 )
            return np.column_stack([ u1, u2 ])
        elif input_type=='rotated':
            cos_t = np.cos(-i_rot_angle).value
            sin_t = np.sin(-i_rot_angle).value
            u1 = ldc1*cos_t - ldc2*sin_t
            u2 = ldc1*sin_t + ldc2*cos_t
            return np.column_stack([ u1, u2 ])
    else:
        u1_u2 = transform_quadratic_ldcs( params, input_type, 'original', i_rot_angle=i_rot_angle )
        ldcs_tr = transform_quadratic_ldcs( u1_u2, 'original', output_type, o_rot_angle=o_rot_angle )
        return ldcs_tr


def check_quadratic_ldcs(params, coeffs_type='original', i_rot_angle=None, positive=True, monotonically_decreasing=True):
    """
    This function checks if quadratic limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile,
    using the transformation defined by Kipping 2013, MNRAS, 435, 2152.
    By default, a correction is applied to get valid coefficients.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    u1_u2 = transform_square_root_ldcs( params, coeffs_type, 'original', i_rot_angle=i_rot_angle )
    u1 = u1_u2[:,0]
    u2 = u1_u2[:,1]
    check1 = np.isfinite(u1)
    check2 = np.isfinite(u2)
    cond = np.column_stack([ check1, check2])
    if positive:
        cond1 = (u1 + u2 <= 1) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = (u1 >= 0) #monotonically decreasing
        cond3 = (u1 + 2.0*u2 >= 0) #monotonically decreasing
        cond = np.column_stack([cond, cond2, cond3])
    cond = cond.all(axis=1)
    return cond


def get_wrms_quadratic(params, mucut, intscut, weights, coeffs_type='original', i_rot_angle=None, o_rot_angle=None, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with quadratic limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - u1*(1-mu) - u2*(1-mu)^2
    
    :param np.array params: 1D array with two limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :argument str coeffs_type: format of the given limb-darkening coefficients (default is 'original')
    :argument quantity i_rot_angle: angle for the rotated limb-darkening coefficients, given as input, with respect to the original ones (default is None)
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_quadratic] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_quadratic] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:2], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_quadratic] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_quadratic_ldcs( params, coeffs_type=coeffs_type, i_rot_angle=i_rot_angle, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'quadratic', coeffs_type=coeffs_type, i_rot_angle=i_rot_angle)
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)



def check_linear_ldcs(params, positive=True, monotonically_decreasing=True):
    """
    This function checks if linear limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.
    By default, a correction is applied to get valid coefficients.

    
    :param np.array params: float or 1D array with one limb-darkening coefficients (or N or list)
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    u = copy.deepcopy( params )
    check = np.isfinite(u)
    cond = np.column_stack([ check, check ])[:,:1]
    if positive:
        cond1 = ( u <= 1.0 ) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = ( u >=0.0 ) #monotonically decreasing
        cond = np.column_stack([cond, cond2])
    cond = cond.all(axis=1)
    return cond


def get_wrms_linear(params, mucut, intscut, weights, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with a linear limb-darkening coefficient:
    I(mu)/I(1) = 1.0 - u*(1-mu)
    
    :param np.array params: 1D array with one limb-darkening coefficient
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.column_stack([params, params])[:,:1]
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_linear] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_linear] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:1], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_linear] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_linear_ldcs( params, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'linear')
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)



def check_gen_poly_ldcs(params, positive=True, monotonically_decreasing=True):
    """
    This function checks if gen_poly limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    params = np.atleast_2d( params )
    nc = np.shape(params)[1]
    cond = np.isfinite( params )
    if positive:
        cond1 = (np.sum( params, axis=1 ) <= 1.0) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = (params[:,0] >= 0.0) #monotonically decreasing
        deriv = make_same_format_2D(1+np.arange(nc), params, first_array='row')[0]
        cond3 = (np.sum( deriv*params, axis=1 ) >= 0.0) #monotonically decreasing
        cond = np.column_stack([cond, cond2, cond3])
    cond = cond.all(axis=1)
    return cond


def get_wrms_gen_poly(params, mucut, intscut, weights, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with polynomial limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - g1*(1-mu) - g2*(1-mu^2) - g3*(1-mu^3) - ... - gn*(1-mu^n)
    
    :param np.array params: 1D array with limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    nc = np.shape(params)[1]
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_poly', nc,'] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_poly', nc,'] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:nc], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_poly', nc,'] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_gen_poly_ldcs( params, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'gen_poly')
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)




def check_gen_claret_ldcs(params, positive=True, monotonically_decreasing=True):
    """
    This function checks if gen_claret4 limb-darkening coefficients comply with the requirements of
    (A) positive, (B) monotonically decreasing profile.

    
    :param np.array params: 1D array with two limb-darkening coefficients (or Nx2 array or list)
    :param str input_type: format of the given limb-darkening coefficients
    :return: True, if the limb-darkening coefficients comply with the above requirements; False otherwise
    :rtype: bool (or list of bool)
    """
    params = np.atleast_2d( params )
    nc = np.shape(params)[1]
    cond = np.isfinite( params )
    if positive:
        cond1 = (np.sum( params, axis=1 ) <= 1.0) #positive
        cond = np.column_stack([cond, cond1])
    if monotonically_decreasing:
        cond2 = (params[:,0] >= 0.0) #monotonically decreasing
        deriv = make_same_format_2D(1+np.arange(nc), params, first_array='row')[0]
        cond3 = (np.sum( deriv*params, axis=1 ) >= 0.0) #monotonically decreasing
        cond = np.column_stack([cond, cond2, cond3])
    cond = cond.all(axis=1)
    return cond


def get_wrms_gen_claret(params, mucut, intscut, weights, positive=True, monotonically_decreasing=True, sum_axis=None):
    """
    This function computes the weighted root mean square of residuals
    between the model intensities and the parametrized values
    with claret-n limb-darkening coefficients:
    I(mu)/I(1) = 1.0 - c1*(1-mu^(1/2)) - c2*(1-mu) - c3*(1-mu^(3/2)) - ... - cn*(1-mu^(n/2))
    
    :param np.array params: 1D array with limb-darkening coefficients
    :param np.array mucut: 1D array with mu values
    :param np.array intscut: 1D array with normalized model intensities at the mu values
    :param np.array weights: 1D array of weights for the fitting algorithm
    ..note:: mucut, intscut and weights must have the same size
    :return: the weighted rms of residuals between model and intscut
    :rtype: float
    """
    params = np.atleast_2d(params)
    nc = np.shape(params)[1]
    mucut = np.atleast_2d(mucut)
    intscut = np.atleast_2d(intscut)
    weights = np.atleast_2d(weights)
    mucut, check = make_same_format_2D(mucut, intscut, first_array='row')
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_claret', nc,'] The given mu values and intensities have incompatible formats.')
        return
    weights, check = make_same_format_2D(weights, mucut, first_array='whatever')
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_claret', nc,'] The given weights and mu values or intensities have incompatible formats.')
        return
    params, check = make_same_format_2D( params, intscut[:,:nc], first_array='row' )
    if not check:
        print('ERROR: [In SAIL.get_wrms_gen_claret', nc,'] The given limb-darkening coefficients and mu values or intensities have incompatible formats.')
        return
    check_ldcs = check_gen_claret_ldcs( params, positive=positive, monotonically_decreasing=monotonically_decreasing )
    model = np.zeros_like( intscut ) + np.inf
    if check_ldcs.any():
        model[check_ldcs,:] = get_intensities_from_ldcs(mucut[check_ldcs,:], params[check_ldcs,:], 'gen_claret')
    return np.sum(weights*((intscut-model)**2), axis=sum_axis) / np.sum(weights, axis=sum_axis)

#limb_darkening_laws_func = {'claret4':claret4, 'square_root':square_root, 'quadratic':quadratic, 'linear':linear, 'gen_claret':gen_claret, 'gen_poly':gen_poly}


def str2float(s):
    """
    This function converts a string to a float, if the string is a number.
    
    :param str s:
    :return: the unmodified string or the number corresponding the string
    :rtype: str or float
    """
    try:
        return float(s)
    except ValueError:
        return s


def my_vstack(array1,array2):
    """
    This function is a variant of numpy.vstack that accepts if one of the two arrays is empty.
    
    :param np.array array1:
    :param np.array array2:
    :return: the array formed by stacking the given arrays
    :rtype: np.array
    """
    try:
        return np.vstack((array1,array2))
    except ValueError:
        if len(array1)==0:
            return array2
        elif len(array2)==0:
            return array1
        else:
            return np.vstack((array1,array2))


def check_length(vector, min_length=1, max_length=None):
    """
    This function checks if the length of a vector is within the expected range and returns a boolean value.
    
    :param listOfObjects vector:
    :argument int min_length: minimum length for the vector (default is 1)
    :argument int max_length: maximum length for the vector (default is None)
    :return: True if len(vector)>=min_length and ( len(vector)<=max_length or max_length is None ), False otherwise
    :rtype: bool
    """
    check = True
    if len(vector)<min_length:
        check = False
    elif max_length:
        if len(vector)>max_length:
            check = False
    return check

def check_type_in_list(vector, item_type):
    """
    This function checks that all the elements of a list are of the expected variable type and returns a boolean value.
    
    :param listOfObjects vector:
    :param type item_type: e.g., str, float, int
    :return: True if all the elements in vector are of the required item_type, False otherwise
    :rtype: bool
    """
    check = True
    for item in vector:
        if not isinstance(item, item_type):
            check = False
    return check

def check_integers(vector, min_value=1):
    """
    This function checks that all the elements of a list are integer numbers greater than or equal to min_value and returns a boolean value.
    
    :param listOfObjects vector:
    :argument int min_value: minimum value for the elements (default is 1)
    :return: True if all the elements in vector are integer numbers greater than or equal to min_value, False otherwise
    :rtype: bool
    """
    check = True
    if np.min(vector)<min_value:
        check = False
    vector_float = np.asarray(vector)
    vector_int = np.asarray(vector,int)
    if (vector_float != vector_int).any():
        check = False
    return check

def read_as_numpy_array(path_file):
    """
    This function reads a numpy array from file as numpy.genfromtxt but also checks for the most common errors.
    It returns a numpy array and a boolean value.
    If the boolean value is False, the numpy array will be empty.
    If the error is caused by one or more elements of the array that are not read as numbers, a warning message will inform about the rows where these errors occur.
    
    :param str path_file: absolute or relative path including the file name
    :return: the numpy array contained in the file and True, or an empty array and False
    :rtype: np.array, bool
    """
    check = True
    try:
        file = np.genfromtxt(path_file)
        if not np.isfinite(file).all(): #try using np.argwhere
            for i in range(len(file)):
                if not np.isfinite(file[i]).all():
                    print('WARNING:', path_file, 'invalid number encountered in line', i, '.')
                    check = False
                    return np.array([]), check
        return file, check
    except IOError:
        print('WARNING:', path_file, 'file not found.')
        check = False
        return np.array([]), check
    except ValueError:
        print('WARNING:', path_file, 'file format invalid.')
        check =  False
        return np.array([]), check


def check_2Darray(arr, n_col=None):
    """
    This function checks that an array has 2 dimensions and returns a boolean value.
    Optionally, it may also check the exact number of columns. 
    
    :param np.array arr: 2D array
    :argument int n_col: exact number of columns in arr (default is None)
    :return: True if arr is 2D and it has n_col columns (or n_col is None), False otherwise
    :rtype: bool
    """
    check = True
    if arr.ndim != 2:
        check = False
    elif n_col:
        if np.shape(arr)[1] != n_col:
            check = False
    return check



def read_configuration(filename):
    """
    This function reads the input file line by line and returns a dictionary.
    For each line the first word is a keyword, the following are values (either string or float).
    The lines starting with '#' will be ignored. The values preceded by '!' (without spaces) will be also ignored. 
    
    :param str filename: absolute or relative path including the file name
    :return: the configuration dictionary
    :rtype: dict
    """
    with open(filename, 'r') as file:
        input_dict = {}
        for line in file.readlines():
            content = line.split()
            if line[0] != '#':
                key = content[0]
                value = []
                for item in content[1:]:
                    if item[0] != '!':
                        value += [str2float(item),]
                input_dict[key] = value
    return input_dict


def read_targets_file(filename):
    """
    This function reads the input file line by line and returns a dictionary.
    It is designed specifically to read the target list.
    The lines starting with '#' will be ignored. The values preceded by '!' (without spaces) will be also ignored. 
    
    :param str filename: absolute or relative path including the file name
    :return: the dictionary of targets
    :rtype: dict
    """
    with open(filename, 'r') as file:
        #check = True
        targets_dict = {}
        target_keys = file.readlines()[0].split() #the first line contains the dictionary keys
        n_keys = len(target_keys)
        mandatory_keys = ['star_effective_temperature',] #
        allowed_keys = mandatory_keys + ['target_names', 'star_log_gravity', 'star_metallicity']
        #Checking that all the keywords in the targets input file are valid and not repeated.
        for key in target_keys:
            if key not in allowed_keys:
                print('ERROR:', key, 'is not a valid keyword. In ', filename, '.')
                #check = False
                exit()
            if target_keys.count(key)>1:
                print('ERROR:',key,' entered multiple times in targets_file.')
                exit()
        #Checking that all the mandatory keywords are obtained from the targets input file.
        for key in mandatory_keys:
            if key not in target_keys:
                print('ERROR: mandatory keyword', key, 'is not specified. In ', filename, '.')
                #check = False
                exit()
        #Initialising dictionary with empty lists
        for key in allowed_keys:
            targets_dict[key] = []

    #Now reading file from the second line
    with open(filename, 'r') as file:
        for line in file.readlines()[1:]:
            content = line.split()
            #Checking that each line has the same length
            if line[0] != '#':
                if not check_length(content, min_length=n_keys, max_length=n_keys):
                    print('ERROR: invalid line ', line, 'in ', filename, '(optional). All lines must have the same length.')
                    #check = False
                    exit()
                #Adding values to the target name and/or stellar parameter lists;
                #stellar parameters need to be converted from string to float
                for i in range(len(content)):
                    targets_dict[target_keys[i]] += [str2float(content[i]),]

        #Checking that list of stellar effective temperatures is float 
        if not check_type_in_list(targets_dict['star_effective_temperature'], float):
            print('ERROR: invalid entry for star_effective_temperature. It must be float. In ', filename, '.')
            #check = False
            exit()

        n_targets = len(targets_dict['star_effective_temperature'])

        #If stellar log gravity has not been read from file or omitted for some targets (it must be denoted by 'X' in the file), assign default log gravity.
        if not check_length(targets_dict['star_log_gravity']):
            targets_dict['star_log_gravity'] = [4.5]*n_targets
        else:
            prov = targets_dict['star_log_gravity']
            provcorr = [4.5 if x=='X' else x for x in prov]
            targets_dict['star_log_gravity'] = copy.deepcopy(provcorr)

        #If stellar metallicity has not been read from file or omitted for some targets (it must be denoted by 'X' in the file), assign default metallicity.
        if not check_length(targets_dict['star_metallicity']):
            targets_dict['star_metallicity'] = [0.0]*n_targets
        else:
            prov = targets_dict['star_metallicity']
            provcorr = [0.0 if x=='X' else x for x in prov]
            targets_dict['star_metallicity'] = copy.deepcopy(provcorr)

        #If target names have not been read from file, build name based on the stellar parameters.
        if 'target_names' not in target_keys:
            star_effective_temperature = copy.deepcopy(targets_dict['star_effective_temperature'])
            star_log_gravity = copy.deepcopy(targets_dict['star_log_gravity'])
            star_metallicity = copy.deepcopy(targets_dict['star_metallicity'])
            target_names = []
            for i in range(n_targets):
                target_names += ['teff' + str(star_effective_temperature[i]) + '_logg' + str(star_log_gravity[i]) + '_MH' + str(star_metallicity[i]),]
            targets_dict['target_names'] = copy.deepcopy(target_names)

#    if not check:
#        exit()

    return targets_dict


#def copy_dict(dict):
#    copied_dict = {}
#    for i in dict:
#        copied_dict[i] = copy.deepcopy(dict[i])
#    return copied_dict


def check_configuration(input_dict):
    """
    This function checks and modifies the dictionary obtained from the configuration file.
    It returns a boolean value and the updated dictionary.
    
    :param dict input_dict: 
    :return: a boolean value and the updated configuration dictionary
    :rtype: bool, dict
    """
    check = True
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    mandatory_keys = ['calculation_type', 'stellar_models_grid', 'limb_darkening_laws', 'passbands']
    allowed_keys = mandatory_keys + ['gen_claret_orders', 'gen_poly_orders', 'targets_file', 'target_names', 'star_effective_temperature', 'star_log_gravity', 'star_metallicity', 'star_minimum_effective_temperature', 'star_maximum_effective_temperature', 'star_minimum_log_gravity', 'star_maximum_log_gravity', 'star_minimum_metallicity', 'star_maximum_metallicity', 'wavelength_bins_files', 'user_output', 'passbands_path', 'wavelength_bins_path', 'targets_path', 'output_path']


    #Checking that all the keywords in the input file are valid.
    for key in input_keys:
        if key not in allowed_keys:
            print('ERROR:', key, 'is not a valid keyword.')
            check = False

    #Checking that all the mandatory keywords are obtained from the input file.
    for key in mandatory_keys:
        if key not in input_keys:
            print('ERROR: mandatory keyword', key, 'is not specified.')
            check = False


    #Checking the requested I/O paths (optional).
    passbands_path = ''
    if 'passbands_path' in input_keys:
        passbands_path = input_dict_local['passbands_path']
        if not check_length(passbands_path, max_length=1):
            print('ERROR: invalid length=', len(passbands_path), 'for passbands_path (optional). It must have length=1.')
            check = False
        if not check_type_in_list(passbands_path, str):
            print('ERROR: invalid type for passbands_path. It must be string.')
            check = False
        passbands_path = input_dict_local['passbands_path'][0]

    wavelength_bins_path = ''
    if 'wavelength_bins_path' in input_keys:
        wavelength_bins_path = input_dict_local['wavelength_bins_path']
        if not check_length(wavelength_bins_path, max_length=1):
            print('ERROR: invalid length=', len(wavelength_bins_path), 'for wavelength_bins_path (optional). It must have length=1.')
            check = False
        if not check_type_in_list(wavelength_bins_path, str):
            print('ERROR: invalid type for wavelength_bins_path. It must be string.')
            check = False
        wavelength_bins_path = input_dict_local['wavelength_bins_path'][0]
        if 'wavelength_bins_files' not in input_keys:
            print('WARNING: ignoring keyword wavelength_bins_path as wavelength_bins_files is not provided')
            wavelength_bins_path = ''
    input_dict_local['wavelength_bins_path'] = [wavelength_bins_path]

    targets_path = ''
    if 'targets_path' in input_keys:
        targets_path = input_dict_local['targets_path']
        if not check_length(targets_path, max_length=1):
            print('ERROR: invalid length=', len(targets_path), 'for targets_path (optional). It must have length=1.')
            check = False
        if not check_type_in_list(targets_path, str):
            print('ERROR: invalid type for targets_path. It must be string.')
            check = False
        targets_path = input_dict_local['targets_path'][0]
        if 'targets_file' not in input_keys:
            print('WARNING: ignoring keyword targets_path as targets_file is not provided')
            targets_path = ''
    input_dict_local['targets_path'] = [targets_path]

    output_path = ''
    if 'output_path' in input_keys:
        output_path = input_dict_local['output_path']
        if not check_length(output_path, max_length=1):
            print('ERROR: invalid length=', len(output_path), 'for output_path (optional). It must have length=1.')
            check = False
        if not check_type_in_list(output_path, str):
            print('ERROR: invalid type for output_path. It must be string.')
            check = False
        output_path = input_dict_local['output_path'][0]
    input_dict_local['output_path'] = [output_path]


    #Checking the (mandatory) calculation_type: 'individual' or 'grid'.
    calculation_type = input_dict_local['calculation_type']
    if not check_length(calculation_type, max_length=1):
        print('ERROR: invalid length=', len(calculation_type), 'for calculation_type. It must have length=1.')
        check = False
    else:
        calculation_type = calculation_type[0]
        allowed_calculation_types = ['grid', 'individual']
        if calculation_type not in allowed_calculation_types:
            print('ERROR: invalid calculation_type. It must be either grid or individual.')
            check = False

    #Checking the choice of stellar_models_grid (only one of those available).
    stellar_models_grid = input_dict_local['stellar_models_grid']
    if not check_length(stellar_models_grid, max_length=1):
        print('ERROR: invalid length=', len(stellar_models_grid), 'for stellar_models_grid. It must have length=1.')
        check = False
    else:
        allowed_stellar_models_grid = ['Phoenix_2018', 'Phoenix_2012_13', 'Phoenix_drift_2012', 'Atlas_2000', 'Stagger_2015', 'Stagger_2018']
        stellar_models_grid = stellar_models_grid[0]
        if stellar_models_grid not in allowed_stellar_models_grid:
            print('ERROR:',stellar_models_grid,'is not a valid stellar_models_grid. The allowed names are Phoenix_2018, Phoenix_2012_13, Phoenix_drift_2012, Atlas_2000, Stagger_2018 and Stagger_2015.')
            check = False

    #Checking the choice of limb_darkening_laws; at least one law must be specified in the input file (no default).
    limb_darkening_laws = input_dict_local['limb_darkening_laws']
    if not check_length(limb_darkening_laws):
        print('ERROR: invalid length=', len(limb_darkening_laws), 'for limb_darkening_laws. It must have length>=1.')
        check = False
    allowed_limb_darkening_laws = ['linear', 'quadratic', 'square_root', 'power2', 'claret4', 'gen_claret', 'gen_poly']
    for item in limb_darkening_laws:
        if item not in allowed_limb_darkening_laws:
            print('ERROR:',item,'is not a valid limb_darkening_laws.')
            check = False
        if limb_darkening_laws.count(item)>1:
            print('WARNING:',item,'limb_darkening_laws entered multiple times. Repetitions are ignored.')
        input_dict_local['limb_darkening_laws'] = list(OrderedDict.fromkeys(limb_darkening_laws))
    #Special case gen_claret: generalization of the claret-4 law defined as 1 minus a linear combination of (1-mu)^(n/2) with n integer number.
    #gen_claret_orders defines the list of maximum n to include in the linear combinations.
    if 'gen_claret' in limb_darkening_laws:
        if 'gen_claret_orders' not in input_keys:
            print('ERROR: gen_claret_orders must be given if gen_claret law is required.')
            check = False
        gen_claret_orders = input_dict_local['gen_claret_orders']
        if not check_length(gen_claret_orders):
            print('ERROR: invalid length=', len(gen_claret_orders), 'for gen_claret_orders. It must have length>=1.')
            check = False
        if not check_type_in_list(gen_claret_orders, float):
            print('ERROR: invalid type for item in gen_claret_orders. Items must be positive integer numbers.')
            check = False
        if not check_integers(gen_claret_orders):
            print('ERROR: invalid type for item in gen_claret_orders. Items must be positive integer numbers.')
            check = False
        gen_claret_orders = np.asarray(gen_claret_orders, int)
    if ('gen_claret' not in limb_darkening_laws) and ('gen_claret_orders' in input_keys):
        print('ERROR: invalid keyword gen_claret_orders if gen_claret limb-darkening law is not required.')
        check = False
    #Special case gen_poly: generic polynomial law defined as 1 minus a linear combination of (1-mu)^n with n integer number.
    #gen_poly_orders defines the list of maximum n to include in the linear combinations.
    if 'gen_poly' in limb_darkening_laws:
        if 'gen_poly_orders' not in input_keys:
            print('ERROR: gen_poly_orders must be given if gen_poly law is required.')
            check = False
        gen_poly_orders = input_dict_local['gen_poly_orders']
        if not check_length(gen_poly_orders):
            print('ERROR: invalid length=', len(gen_poly_orders), 'for gen_poly_orders. It must have length>=1.')
            check = False
        if not check_type_in_list(gen_poly_orders, float):
            print('ERROR: invalid type for item in gen_poly_orders. Items must be positive integer numbers.')
            check = False
        if not check_integers(gen_poly_orders):
            print('ERROR: invalid type for item in gen_poly_orders. Items must be positive integer numbers.')
            check = False
        gen_poly_orders = np.asarray(gen_poly_orders, int)
        gen_poly_orders = np.unique(gen_poly_orders)
    if ('gen_poly' not in limb_darkening_laws) and ('gen_poly_orders' in input_keys):
        print('ERROR: invalid keyword gen_poly_orders if gen_poly limb-darkening law is not required.')
        check = False

    #Checking that the passbands list is not empty
    passbands = input_dict_local['passbands']
    if not check_length(passbands):
        print('ERROR: invalid length=', len(passbands), 'for passbands. It must have length>=1.')
        check = False

    #Checking built-in passbands
    passbands_path_package = pkg_resources.resource_filename('exotethys','Passbands/')
    list_passbands = [f for f in os.listdir(passbands_path_package) if f.endswith('.pass')]
    allowed_passbands = [f[:-5] for f in list_passbands] 
    input_dict_local['passbands_ext'] = ['']
    prov = ''.join(passbands)
    #if no input passband_path and input passbands have no extensions, assume that user requests built-in passbands
    if ('passbands_path' not in input_keys) and ('.' not in ''.join(passbands)):
        input_dict_local['passbands_path'] = [passbands_path_package]
        input_dict_local['passbands_ext'] = ['.pass']
        for item in passbands:
            if item not in allowed_passbands:
                print('ERROR:',item,'is not a valid built-in passbands.')
                check = False
    elif ('passbands_path' not in input_keys) and ('.' in ''.join(passbands)):
        input_dict_local['passbands_path'] = ['']
    #input_dict_local['passbands_path'] = [passbands_path]


    #Eliminating repeated passbands
    for item in passbands:
        if passbands.count(item)>1:
            print('WARNING:', item, 'entered multiple times in passbands. Repetitions are ignored.')
        input_dict_local['passbands'] = list(OrderedDict.fromkeys(passbands))


    #Checking the OPTIONAL wavelength_bins_files to split the selected passbands; by DEFAULT the passbands are not split.
    #If wavelength_bins_files is activated, the corresponding line in the input file must contain the wavelength bins file names or the string 'no_bins' (one string per passband).
    passbands = input_dict_local['passbands']
    n_pass = len(passbands)
    if 'wavelength_bins_files' in input_keys:
        wavelength_bins_files = input_dict_local['wavelength_bins_files']
        if not check_length(wavelength_bins_files, min_length=n_pass, max_length=n_pass):
            print('ERROR: invalid length=', len(wavelength_bins_files), 'for wavelength_bins_files (optional). If defined, one file or keyword no_bins per passband is required.')
            check = False
        if not check_type_in_list(wavelength_bins_files, str):
            print('ERROR: invalid type for item in wavelength_bins_files. It must be a file to read or keyword no_bins.')
            check = False


    #Checking the requested user_output (optional): 'basic' or 'complete'. Default is 'basic'.
    user_output = 'basic'
    if 'user_output' in input_keys:
        user_output = input_dict_local['user_output']
        if not check_length(user_output, max_length=1):
            print('ERROR: invalid length=', len(user_output), 'for user_output (optional). It must have length=1.')
            check = False
        user_output = user_output[0]
        allowed_user_output = ['basic', 'complete']
        if user_output not in allowed_user_output:
            print('ERROR: invalid user_output. It must be basic or complete (default is basic).')
            check = False
    input_dict_local['user_output'] = [user_output]

    #Checking keywords for individual/grid calculation_type
    mandatory_keys_individual = ['star_effective_temperature'] 
    allowed_keys_individual = ['target_names', 'star_effective_temperature', 'star_log_gravity', 'star_metallicity']
    mandatory_keys_grid = []
    allowed_keys_grid = ['star_minimum_effective_temperature', 'star_maximum_effective_temperature', 'star_minimum_log_gravity', 'star_maximum_log_gravity', 'star_minimum_metallicity', 'star_maximum_metallicity']

    #Checking that all the mandatory keywords are obtained from the configuration file, if individual calculation_type if without targets_file.
    #Also removing those keywords that will not be used (WARNING messages will appear). 
    if calculation_type == 'individual' and 'targets_file' not in input_keys:
        keys_to_delete = []
        for key in mandatory_keys_individual:
            if key not in input_keys:
                print('ERROR: mandatory keyword', key, 'for individual calculation_type is not specified.')
                check = False
        for key in input_keys:
            if key in allowed_keys_grid:
                print('WARNING:', key, 'is not a valid keyword for individual calculation_type. It will be ignored.')
                keys_to_delete.append(key)
        for i in keys_to_delete:
            del input_dict_local[i]
    elif calculation_type == 'individual' and 'targets_file' in input_keys:
        keys_to_delete = []
        for key in input_keys:
            if key in allowed_keys_individual:
                print('WARNING:', key, 'is not considered for individual calculation_type with targets_file. It will be ignored.')
                keys_to_delete.append(key)
        for i in keys_to_delete:
            del input_dict_local[i]


    #Checking that all the mandatory keywords are obtained from the configuration file, if grid calculation_type.
    #Also removing those keywords that will not be used (WARNING messages will appear). 
    if calculation_type == 'grid':
        keys_to_delete = []
        for key in mandatory_keys_grid:
            if key not in input_keys:
                print('ERROR: mandatory keyword', key, 'for grid calculation_type is not specified.')
                check = False
        for key in input_keys:
            if key in allowed_keys_individual:
                print('WARNING:', key, 'is not a valid keyword for grid calculation_type. It will be ignored.')
                keys_to_delete.append(key)

        for i in keys_to_delete:
            del input_dict_local[i]

    input_keys = list(input_dict_local.keys())

    #Checking the keywords which are specific to one of the alternative calculation_type.
    if 'targets_file' in input_keys:
        targets_file = input_dict_local['targets_file']
        if not check_length(targets_file, max_length=1):
            print('ERROR: invalid length=', len(targets_file), 'for targets_file (optional). It must have length=1.')
            check = False
        if not check_type_in_list(targets_file, str):
            print('ERROR: invalid type for targets_file. It must be string.')
            check = False

    #Checking star_effective_temperature
    if 'star_effective_temperature' in input_keys:
        star_effective_temperature = input_dict_local['star_effective_temperature']
        if not check_length(star_effective_temperature):
            print('ERROR: invalid length=', len(star_effective_temperature), 'for star_effective_temperature. It must have length>=1.')
            check = False
        if not check_type_in_list(star_effective_temperature, float):
            print('ERROR: invalid type for item in star_effective_temperature. Items must be float.')
            check = False

    #Checking star_log_gravity
    if 'star_log_gravity' in input_keys:
        star_log_gravity = input_dict_local['star_log_gravity']
        if not check_length(star_log_gravity):
            print('ERROR: invalid length=', len(star_log_gravity), 'for star_log_gravity (optional). It must have length>=1.')
            check = False
        if not check_type_in_list(star_log_gravity, float):
            print('ERROR: invalid type for item in star_log_gravity (optional). Items must be float.')
            check = False


    #Checking star_metallicity
    if 'star_metallicity' in input_keys:
        star_metallicity = input_dict_local['star_metallicity']
        if not check_length(star_metallicity):
            print('ERROR: invalid length=', len(star_metallicity), 'for star_metallicity (optional). It must have length>=1.')
            check = False
        if not check_type_in_list(star_metallicity, float):
            print('ERROR: invalid type for item in star_metallicity (optional). Items must be float.')
            check = False


    #Checking target_names
    if 'target_names' in input_keys:
        target_names = input_dict_local['target_names']
        if not check_length(target_names):
            print('ERROR: invalid length=', len(target_names), 'for target_names (optional). It must have length>=1.')
            check = False

    #Checking star_minimum_effective_temperature
    if 'star_minimum_effective_temperature' in input_keys:
        star_minimum_effective_temperature = input_dict_local['star_minimum_effective_temperature']
        if not check_length(star_minimum_effective_temperature, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_effective_temperature), 'for star_minimum_effective_temperature. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_effective_temperature, float):
                print('ERROR: invalid type for star_minimum_effective_temperature. It must be float.')
                check = False

    #Checking star_maximum_effective_temperature
    if 'star_maximum_effective_temperature' in input_keys:
        star_maximum_effective_temperature = input_dict_local['star_maximum_effective_temperature']
        if not check_length(star_maximum_effective_temperature, max_length=1):
            print('ERROR: invalid length=', len(star_maximum_effective_temperature), 'for star_maximum_effective_temperature. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_maximum_effective_temperature, float):
                print('ERROR: invalid type for star_maximum_effective_temperature. It must be float.')
                check = False

    #Checking star_minimum_log_gravity
    if 'star_minimum_log_gravity' in input_keys:
        star_minimum_log_gravity = input_dict_local['star_minimum_log_gravity']
        if not check_length(star_minimum_log_gravity, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_log_gravity), 'for star_minimum_log_gravity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_log_gravity, float):
                print('ERROR: invalid type for star_minimum_log_gravity. It must be float.')
                check = False

    #Checking star_maximum_log_gravity
    if 'star_maximum_log_gravity' in input_keys:
        star_maximum_log_gravity = input_dict_local['star_maximum_log_gravity']
        if not check_length(star_maximum_log_gravity, max_length=1):
            print('ERROR: invalid length=', len(star_maximum_log_gravity), 'for star_maximum_log_gravity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_maximum_log_gravity, float):
                print('ERROR: invalid type for star_maximum_log_gravity. It must be float.')
                check = False

    #Checking star_minimum_metallicity
    if 'star_minimum_metallicity' in input_keys:
        star_minimum_metallicity = input_dict_local['star_minimum_metallicity']
        if not check_length(star_minimum_metallicity, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_metallicity), 'for star_minimum_metallicity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_metallicity, float):
                print('ERROR: invalid type for star_minimum_metallicity. It must be float.')
                check = False

    #Checking star_maximum_metallicity
    if 'star_maximum_metallicity' in input_keys:
        star_maximum_metallicity = input_dict_local['star_maximum_metallicity']
        if not check_length(star_maximum_metallicity, max_length=1):
            print('ERROR: invalid length=', len(star_maximum_metallicity), 'for star_maximum_metallicity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_maximum_metallicity, float):
                print('ERROR: invalid type for star_maximum_metallicity. It must be float.')
                check = False

    if check==False:
        exit()

    #input_dict.update(input_dict_local)

    return check, input_dict_local



def get_individual_parameters(input_dict):
    """
    This function extracts and returns the stellar parameters and target names from the input dictionary.
    
    :param dict input_dict: 
    :return: the lists of stellar parameters (effective temperature, log gravity, metallicity) and names to be used for the output files.
    :rtype: list of float, list of float, list of float, list of str
    """
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    #Case: read parameters from targets_file (automatically completed for missing values)
    if 'targets_file' in input_keys:
        targets_file = input_dict_local['targets_file'][0]
        targets_path = ''
        if 'targets_path' in input_keys:
            targets_path = input_dict_local['targets_path'][0]
        targets_dict = read_targets_file(os.path.join(targets_path,targets_file))
        star_effective_temperature = copy.deepcopy(targets_dict['star_effective_temperature'])
        star_log_gravity = copy.deepcopy(targets_dict['star_log_gravity'])
        star_metallicity = copy.deepcopy(targets_dict['star_metallicity'])
        target_names = copy.deepcopy(targets_dict['target_names'])
    #Case: read parameters from entered through the configuration file and eventually complete with default parameters
    else:
        star_effective_temperature = copy.deepcopy(input_dict_local['star_effective_temperature'])
        n_targets = len(star_effective_temperature)
        #
        if 'star_log_gravity' in input_keys:
            star_log_gravity = copy.deepcopy(input_dict_local['star_log_gravity'])
            if len(star_log_gravity)==1 and n_targets>1:
                print('WARNING: assuming star_log_gravity=', star_log_gravity[0], 'for all targets. A single input value was provided.')
                star_log_gravity = star_log_gravity*n_targets
            elif len(star_log_gravity)>1 and len(star_log_gravity)!=n_targets:
                print('ERROR: invalid length for star_log_gravity (optional). It must be equal to the length of star_effective_temperature or a single number.')
                exit()
        elif 'star_log_gravity' not in input_keys:
            print('WARNING: assuming star_log_gravity=4.5 for all targets. No input values were provided.')
            star_log_gravity = [4.5]*n_targets
        #
        if 'star_metallicity' in input_keys:
            star_metallicity = copy.deepcopy(input_dict_local['star_metallicity'])
            if len(star_metallicity)==1 and n_targets>1:
                print('WARNING: assuming star_metallicity=', star_metallicity[0], 'for all targets. A single input value was provided.')
                star_metallicity = star_metallicity*n_targets
            elif len(star_metallicity)>1 and len(star_metallicity)!=n_targets:
                print('ERROR: invalid length for star_metallicity (optional). It must be equal to the length of star_effective_temperature or a single number.')
                exit()
        elif 'star_metallicity' not in input_keys:
            print('WARNING: assuming star_metallicity=0.0 for all targets. No input values were provided.')
            star_metallicity = [0.0]*n_targets
        #
        if 'target_names' in input_keys:
            target_names = input_dict_local['target_names']
            if not check_type_in_list(target_names, str):
                print('ERROR: invalid type for item in target_names. They must be strings.')
                exit()
            if len(target_names)!=n_targets:
                print('ERROR: invalid length for target_names (optional). It must be equal to the length of star_effective_temperature.')
                exit()
        elif 'target_names' not in input_keys:
            target_names = []
            for i in range(n_targets):
                target_names += ['teff' + str(star_effective_temperature[i]) + '_logg' + str(star_log_gravity[i]) + '_MH' + str(star_metallicity[i]),]

    return star_effective_temperature, star_log_gravity, star_metallicity, target_names



def stellar_params_from_file_name(file_name):
    """
    This function extracts and returns the stellar parameters from the file name in the database.
    
    :param str file_name: built-in file name
    :return: a numpy array with the effective temperature, log gravity and metallicity for the input file name.
    :rtype: np.array
    """
    params = os.path.basename(file_name).replace('.pickle', '').split('_')
    teff = float(params[0].replace('teff', ''))
    logg = float(params[1].replace('logg', ''))
    mh = float(params[2].replace('MH', ''))
    return np.array([teff, logg, mh])


def get_grid_parameters(stellar_models_grid):
    """
    This function gets all the file names in the database and the corresponding stellar parameters.
    
    :param str stellar_models_grid:
    :return: the file names and a numpy array with the corresponding stellar parameters (number of files X 3)
    :rtype: list of str, np.array
    """
    files = databases[stellar_models_grid].get_filename_list()
    star_params_grid = np.zeros((len(files),3))
    for i in range(len(files)):
        star_params_grid[i,:] = stellar_params_from_file_name(files[i])
    return files, star_params_grid


def get_subgrid(input_dict):
    """
    This function extracts the file names, indices and stellar parameters from the selected database of stellar model-atmospheres within the selected (optional) parameter ranges.
    
    :param dict input_dict:
    :return: the lists of database file names, the associated stellar parameters and file indices with the (optional) user parameter limits. It exits if the list of files is empty.
    :rtype: list of str, np.array, np.array
    """
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    stellar_models_grid = input_dict_local['stellar_models_grid'][0]
    [files, star_params_grid] = get_grid_parameters(stellar_models_grid)
    subgrid_indices = np.arange(np.shape(star_params_grid)[0])
    #teff_min_indices = np.arange(np.shape(star_params_grid)[0])
    #teff_max_indices = np.arange(np.shape(star_params_grid)[0])
    #logg_min_indices = np.arange(np.shape(star_params_grid)[0])
    #logg_max_indices = np.arange(np.shape(star_params_grid)[0])
    #mh_min_indices = np.arange(np.shape(star_params_grid)[0])
    #mh_max_indices = np.arange(np.shape(star_params_grid)[0])
    check = True
    if 'star_minimum_effective_temperature' in input_keys:
        teff_min = input_dict_local['star_minimum_effective_temperature']
        if teff_min>np.max(star_params_grid[:,0]):
            print('ERROR: star_minimum_effective_temperature is too high.')
            check = False
        elif teff_min<=np.min(star_params_grid[:,0]):
            print('WARNING: star_minimum_effective_temperature is lower than the minimum temperature available in the grid.' )
        teff_min_indices = np.where(star_params_grid[:,0]>=teff_min)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, teff_min_indices)
    if 'star_maximum_effective_temperature' in input_keys:
        teff_max = input_dict_local['star_maximum_effective_temperature']
        if teff_max<np.min(star_params_grid[:,0]):
            print('ERROR: star_maximum_effective_temperature is too low.')
            check = False
        elif teff_max>=np.max(star_params_grid[:,0]):
            print('WARNING: star_maximum_effective_temperature is higher than the maximum temperature available in the grid.' )
        teff_max_indices = np.where(star_params_grid[:,0]<=teff_max)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, teff_max_indices)
    if 'star_minimum_log_gravity' in input_keys:
        logg_min = input_dict_local['star_minimum_log_gravity']
        if logg_min>np.max(star_params_grid[:,1]):
            print('ERROR: star_minimum_log_gravity is too high.')
            check = False
        elif logg_min<=np.min(star_params_grid[:,1]):
            print('WARNING: star_minimum_log_gravity is lower than the minimum log(g) available in the grid.' )
        logg_min_indices = np.where(star_params_grid[:,1]>=logg_min)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, logg_min_indices)
    if 'star_maximum_log_gravity' in input_keys:
        logg_max = input_dict_local['star_maximum_log_gravity']
        if logg_max<np.min(star_params_grid[:,1]):
            print('ERROR: star_maximum_log_gravity is too low.')
            check = False
        elif logg_max>=np.max(star_params_grid[:,1]):
            print('WARNING: star_maximum_log_gravity is higher than the maximum log(g) available in the grid.' )
        logg_max_indices = np.where(star_params_grid[:,1]<=logg_max)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, logg_max_indices)
    if 'star_minimum_metallicity' in input_keys:
        mh_min = input_dict_local['star_minimum_metallicity']
        if mh_min>np.max(star_params_grid[:,2]):
            print('ERROR: star_minimum_metallicity is too high.')
            check = False
        elif mh_min<=np.min(star_params_grid[:,2]):
            print('WARNING: star_minimum_metallicity is lower than the minimum metallicity available in the grid.' )
        mh_min_indices = np.where(star_params_grid[:,2]>=mh_min)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, mh_min_indices)
    if 'star_maximum_metallicity' in input_keys:
        mh_max = input_dict_local['star_maximum_metallicity']
        if mh_max<np.min(star_params_grid[:,2]):
            print('ERROR: star_maximum_metallicity is too low.')
            check = False
        elif mh_max>=np.max(star_params_grid[:,2]):
            print('WARNING: star_maximum_metallicity is higher than the maximum metallicity available in the grid.' )
        mh_max_indices = np.where(star_params_grid[:,2]<=mh_max)[0]
        subgrid_indices = np.intersect1d(subgrid_indices, mh_max_indices)
    if not check:
        exit()
    elif len(subgrid_indices)<1:
        print('ERROR: There are no stellar models in the', stellar_models_grid, 'stellar_models_grid within the selected parameter ranges.')
        exit()
    else:
        star_params_subgrid = star_params_grid[subgrid_indices,:]
        files_subgrid = [files[i] for i in subgrid_indices]
    return files_subgrid, star_params_subgrid, subgrid_indices



def get_neighbour_files_indices(target_name,star_effective_temperature,star_log_gravity,star_metallicity,star_params_grid,stellar_models_grid):
    """
    This function searchs for up to 8 neighbours for a given set of stellar parameters and returns the corresponding indices in the files list (repeated indices are possible).
    
    :param str target_name: name associated with the target to compute
    :param float star_effective_temperature: stellar effective temperature of the target
    :param float star_log_gravity: stellar log gravity of the target
    :param float star_metallicity: stellar metallicity of the target
    :param np.array star_params_grid: 3-column array with stellar effective temperatures, log gravities and metallicities of the database models
    :params str stellar_models_grid: the name of the chosen stellar database
    :return: the indices of the (up to) eight nearest neighbour models (the array will contain 8 indices, but repeated indices are possible).
    :rtype: np.array
    """
    teff_min = np.min(star_params_grid[:,0])
    teff_max = np.max(star_params_grid[:,0])
    logg_min = np.min(star_params_grid[:,1])
    logg_max = np.max(star_params_grid[:,1])
    mh_min = np.min(star_params_grid[:,2])
    mh_max = np.max(star_params_grid[:,2])
    cond1 = (star_effective_temperature<teff_min)
    cond2 = (star_effective_temperature>teff_max)
    cond3 = (star_log_gravity<logg_min)
    cond4 = (star_log_gravity>logg_max)
    cond5 = (star_metallicity<mh_min)
    cond6 = (star_metallicity>mh_max)
    if cond1 or cond2:
        print('WARNING:', target_name, 'cannot be calculated. The star_effective_temperature is out of the range (',teff_min,'-',teff_max,' K) for the stellar_models_grid', stellar_models_grid, '.')
        neigh_indices = np.array([])
    if cond3 or cond4:
        print('WARNING:', target_name, 'cannot be calculated. The star_log_gravity is out of the range (',logg_min,'-',logg_max,') for the stellar_models_grid', stellar_models_grid, '.')
        neigh_indices = np.array([])
    if cond5 or cond6:
        print('WARNING:', target_name, 'cannot be calculated. The star_metallicity is out of the range (',mh_min,'-',mh_max,') for the stellar_models_grid', stellar_models_grid, '.')
        neigh_indices = np.array([])
    if cond1 or cond2 or cond3 or cond4 or cond5 or cond6:
        return neigh_indices
    #
    indices_teff_sup = np.where(star_params_grid[:,0]>=star_effective_temperature)[0]
    indices_teff_inf = np.where(star_params_grid[:,0]<=star_effective_temperature)[0]
    indices_logg_sup = np.where(star_params_grid[:,1]>=star_log_gravity)[0]
    indices_logg_inf = np.where(star_params_grid[:,1]<=star_log_gravity)[0]
    indices_mh_sup = np.where(star_params_grid[:,2]>=star_metallicity)[0]
    indices_mh_inf = np.where(star_params_grid[:,2]<=star_metallicity)[0]
    #
    indices_teff_sup_logg_sup = np.intersect1d( indices_teff_sup, indices_logg_sup )
    indices_teff_sup_logg_sup_mh_sup = np.intersect1d( indices_teff_sup_logg_sup, indices_mh_sup )
    indices_teff_sup_logg_sup_mh_inf = np.intersect1d( indices_teff_sup_logg_sup, indices_mh_inf )
    indices_teff_sup_logg_inf = np.intersect1d( indices_teff_sup, indices_logg_inf )
    indices_teff_sup_logg_inf_mh_sup = np.intersect1d( indices_teff_sup_logg_inf, indices_mh_sup )
    indices_teff_sup_logg_inf_mh_inf = np.intersect1d( indices_teff_sup_logg_inf, indices_mh_inf )
    indices_teff_inf_logg_sup = np.intersect1d( indices_teff_inf, indices_logg_sup )
    indices_teff_inf_logg_sup_mh_sup = np.intersect1d( indices_teff_inf_logg_sup, indices_mh_sup )
    indices_teff_inf_logg_sup_mh_inf = np.intersect1d( indices_teff_inf_logg_sup, indices_mh_inf )
    indices_teff_inf_logg_inf = np.intersect1d( indices_teff_inf, indices_logg_inf )
    indices_teff_inf_logg_inf_mh_sup = np.intersect1d( indices_teff_inf_logg_inf, indices_mh_sup )
    indices_teff_inf_logg_inf_mh_inf = np.intersect1d( indices_teff_inf_logg_inf, indices_mh_inf )
    condd1 = (len(indices_teff_sup_logg_sup_mh_sup)==0)
    condd2 = (len(indices_teff_sup_logg_sup_mh_inf)==0)
    condd3 = (len(indices_teff_sup_logg_inf_mh_sup)==0)
    condd4 = (len(indices_teff_sup_logg_inf_mh_inf)==0)
    condd5 = (len(indices_teff_inf_logg_sup_mh_sup)==0)
    condd6 = (len(indices_teff_inf_logg_sup_mh_inf)==0)
    condd7 = (len(indices_teff_inf_logg_inf_mh_sup)==0)
    condd8 = (len(indices_teff_inf_logg_inf_mh_inf)==0)
    if condd1:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 1 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd2:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 2 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd3:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 3 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd4:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 4 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd5:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 5 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd6:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 6 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd7:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 7 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd8:
        print('WARNING:', target_name, 'cannot be calculated. Neighbour 8 not found for the stellar_models_grid', stellar_models_grid, '.')
    if condd1 or condd2 or condd3 or condd4 or condd5 or condd6 or condd7 or condd8:
        neigh_indices = np.array([])
        return neigh_indices
    #
    indices1t = indices_teff_sup_logg_sup_mh_sup[np.where( star_params_grid[indices_teff_sup_logg_sup_mh_sup,0]-star_effective_temperature == np.min(star_params_grid[indices_teff_sup_logg_sup_mh_sup,0]-star_effective_temperature) )[0] ]
    indices1tg = indices1t[np.where( star_params_grid[indices1t,1]-star_log_gravity == np.min(star_params_grid[indices1t,1]-star_log_gravity) )[0] ]
    index1tgm = indices1tg[np.argmin( star_params_grid[indices1tg,2]-star_metallicity )]
    #
    indices2t = indices_teff_sup_logg_sup_mh_inf[np.where( star_params_grid[indices_teff_sup_logg_sup_mh_inf,0]-star_effective_temperature == np.min(star_params_grid[indices_teff_sup_logg_sup_mh_inf,0]-star_effective_temperature) )[0] ]
    indices2tg = indices2t[np.where( star_params_grid[indices2t,1]-star_log_gravity == np.min(star_params_grid[indices2t,1]-star_log_gravity) )[0] ]
    index2tgm = indices2tg[np.argmin( star_metallicity-star_params_grid[indices2tg,2] )]
    #
    indices3t = indices_teff_sup_logg_inf_mh_sup[np.where( star_params_grid[indices_teff_sup_logg_inf_mh_sup,0]-star_effective_temperature == np.min(star_params_grid[indices_teff_sup_logg_inf_mh_sup,0]-star_effective_temperature) )[0] ]
    indices3tg = indices3t[np.where( star_log_gravity-star_params_grid[indices3t,1] == np.min(star_log_gravity-star_params_grid[indices3t,1]) )[0] ]
    index3tgm = indices3tg[np.argmin( star_params_grid[indices3tg,2]-star_metallicity )]
    #
    indices4t = indices_teff_sup_logg_inf_mh_inf[np.where( star_params_grid[indices_teff_sup_logg_inf_mh_inf,0]-star_effective_temperature == np.min(star_params_grid[indices_teff_sup_logg_inf_mh_inf,0]-star_effective_temperature) )[0] ]
    indices4tg = indices4t[np.where( star_log_gravity-star_params_grid[indices4t,1] == np.min(star_log_gravity-star_params_grid[indices4t,1]) )[0] ]
    index4tgm = indices4tg[np.argmin( star_metallicity-star_params_grid[indices4tg,2] )]
    #
    indices5t = indices_teff_inf_logg_sup_mh_sup[np.where( star_effective_temperature-star_params_grid[indices_teff_inf_logg_sup_mh_sup,0] == np.min(star_effective_temperature-star_params_grid[indices_teff_inf_logg_sup_mh_sup,0]) )[0] ]
    indices5tg = indices5t[np.where( star_params_grid[indices5t,1]-star_log_gravity == np.min(star_params_grid[indices5t,1]-star_log_gravity) )[0] ]
    index5tgm = indices5tg[np.argmin( star_params_grid[indices5tg,2]-star_metallicity )]
    #
    indices6t = indices_teff_inf_logg_sup_mh_inf[np.where( star_effective_temperature-star_params_grid[indices_teff_inf_logg_sup_mh_inf,0] == np.min(star_effective_temperature-star_params_grid[indices_teff_inf_logg_sup_mh_inf,0]) )[0] ]
    indices6tg = indices6t[np.where( star_params_grid[indices6t,1]-star_log_gravity == np.min(star_params_grid[indices6t,1]-star_log_gravity) )[0] ]
    index6tgm = indices6tg[np.argmin( star_metallicity-star_params_grid[indices6tg,2] )]
    #
    indices7t = indices_teff_inf_logg_inf_mh_sup[np.where( star_effective_temperature-star_params_grid[indices_teff_inf_logg_inf_mh_sup,0] == np.min(star_effective_temperature-star_params_grid[indices_teff_inf_logg_inf_mh_sup,0]) )[0] ]
    indices7tg = indices7t[np.where( star_log_gravity-star_params_grid[indices7t,1] == np.min(star_log_gravity-star_params_grid[indices7t,1]) )[0] ]
    index7tgm = indices7tg[np.argmin( star_params_grid[indices7tg,2]-star_metallicity )]
    #
    indices8t = indices_teff_inf_logg_inf_mh_inf[np.where( star_effective_temperature-star_params_grid[indices_teff_inf_logg_inf_mh_inf,0] == np.min(star_effective_temperature-star_params_grid[indices_teff_inf_logg_inf_mh_inf,0]) )[0] ]
    indices8tg = indices8t[np.where( star_log_gravity-star_params_grid[indices8t,1] == np.min(star_log_gravity-star_params_grid[indices8t,1]) )[0] ]
    index8tgm = indices8tg[np.argmin( star_metallicity-star_params_grid[indices8tg,2] )]
    #
    neigh_indices = np.array([index1tgm, index2tgm, index3tgm, index4tgm, index5tgm, index6tgm, index7tgm, index8tgm])
    return neigh_indices


def check_passband_limits(pb_waves, stellar_models_grid):
    """
    This function checks that the wavelengths read from a passband file are within the limits for the chosen stellar_models_grid and returns a boolean value.
    
    :param quantity array pb_waves: 1D array of wavelengths read from the passband file
    :params str stellar_models_grid: the name of the chosen stellar database
    :return: True if the wavelengths are within the limits, False otherwise
    :rtype: bool
    """
    check = True
    if stellar_models_grid == 'Phoenix_2018':
        minimum_wavelength = 500.0 * u.Angstrom
        maximum_wavelength = 25999.0 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Phoenix_2012_13':
        minimum_wavelength = 2500.0 * u.Angstrom
        maximum_wavelength = 99995.0 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Phoenix_drift_2012':
        minimum_wavelength = 10.0 * u.Angstrom
        maximum_wavelength = 9000000.0 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Atlas_2000':
        minimum_wavelength = 90.9 * u.Angstrom
        maximum_wavelength = 1600000.0 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Stagger_2018':
        minimum_wavelength = 1010.0 * u.Angstrom
        maximum_wavelength = 199960.16 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Stagger_2015':
        minimum_wavelength = 2000.172119140625 * u.Angstrom
        maximum_wavelength = 10000.0791015625 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    return check


def read_passband(passbands_path, passbands_ext, passband, stellar_models_grid):
    """
    This function calls another function to read the passband from file and checks its content. It returns the extracted columns and a boolean value.
    
    :param str passbands_path:
    :param str passbands_ext:
    :param str passband:
    ..note:: the above params define the file with the passband to read
    :params str stellar_models_grid: the name of the chosen stellar database
    :return: the wavelengths (in Angstrom), the corresponding pce values (in electron/photon), and a boolean value.
    :rtype: quantity array, quantity array, bool
    """
    [pb, check] = read_as_numpy_array(os.path.join(passbands_path,passband)+passbands_ext)
    if not check:
        print('WARNING: Skipping passband', passband, '.')
        return np.array([]) * u.Angstrom, np.array([]) * u.electron /u.photon, check
    pb = np.atleast_2d(pb)
    if not check_2Darray(pb, n_col=2):
        print('WARNING: invalid format for', passband, 'passband file. It must have 2 columns. Skipping', passband, 'passband.')
        check = False
        return np.array([]) * u.Angstrom, np.array([]) * u.electron /u.photon, check
    if np.min(pb)<0:
        print('WARNING: negative value found in file', passband, 'passband file. Skipping', passband, 'passband.')
        check = False
        return np.array([]) * u.Angstrom, np.array([]) * u.electron /u.photon, check
    pb_waves = pb[:,0] * u.Angstrom
    pb_pce = pb[:,1] * u.electron /u.photon 
    check = check_passband_limits(pb_waves, stellar_models_grid)
    if not check:
        print('WARNING:', passband, 'passband exceeds wavelength range for the', stellar_models_grid, 'stellar model grid. Skipping', passband, 'passband.')
        return np.array([]) * u.Angstrom, np.array([]) * u.electron /u.photon, check
    return pb_waves, pb_pce, check

def read_wavelength_bins(wavelength_bins_path, wavelength_bins_file, pb_waves, passband):
    """
    This function calls another function to read the wavelength bins file and checks its content, unless it is set to no_bins. It returns the extracted array and a boolean value.
    
    :param str wavelength_bins_path: the path to wavelength bins files (without the file name)
    :param str wavelength_bins_file: the name of the wavelength bins file or the string no_bins
    :param quantity array pb_waves: 1D array of wavelengths (in Angstrom)
    ..note:: the above params define the file with the passband to read
    :return: the wavelength bins 2D or empty array (in Angstrom) and a boolean value.
    :rtype: quantity array, bool
    """
    check =  True
    if wavelength_bins_file == 'no_bins':
        check = True
        return np.array([]) * u.Angstrom, check
    wavelength_bins_pathfile = os.path.join(wavelength_bins_path, wavelength_bins_file)
    [wb, check] = read_as_numpy_array(wavelength_bins_pathfile)
    if not check:
        print('WARNING: Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        return np.array([]) * u.Angstrom, check
    wb = np.atleast_2d(wb)
    if not check_2Darray(wb, n_col=2):
        print('WARNING: invalid format for wavelength_bins_file', wavelength_bins_pathfile, '. It must have 2 columns. Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        check = False
        return np.array([]) * u.Angstrom, check
    wb *= u.Angstrom
    if (wb[:,0]>=wb[:,1]).any():
        print('WARNING: invalid line in wavelength_bins_file', wavelength_bins_pathfile, '. The lower limit cannot be greater or equal to the upper limit.')
        check = False
        return np.array([]) * u.Angstrom, check
    if np.min(wb)<np.min(pb_waves) or np.max(wb)>np.max(pb_waves):
        print('WARNING:wavelength_bins_file', wavelength_bins_pathfile, 'exceeds wavelength range for the ', passband, 'passband. Ignoring this file.')
        check = False
        return np.array([]) * u.Angstrom, check
    return wb, check


def get_waves_fromR(lambda1, lambda2, R=1e6):
    """
    This function computes an array of wavelengths with constant spectral resolution, including the minimum and maximum wavelengths.
    The requested resolution can be automatically increased to guarantee the minimum length of 10 for the array of wavelengths.
    
    :param float lambda1: minimum wavelength
    :param float lambda2: maximum wavelength
    :argument float R: spectral resolution
    ..note:: the input params are already checked to avoid errors, if using the boats_calculate_transit or boats_calculate_eclipse functions
    :return: the 1D array of wavelengths
    :rtype: np.array
    """
    log_lambda1 = np.log(lambda1)
    log_lambda2 = np.log(lambda2)
    log_step = np.log( (2*R+1)/(2*R-1) )
    n_steps = np.ceil( (log_lambda2-log_lambda1)/log_step )
    n_steps = np.max([10, n_steps])
    log_waves = np.linspace(log_lambda1, log_lambda2, num=int(n_steps))
    waves = np.exp(log_waves)
    return waves


def get_passband_intensities(model_dict, passbands_dict):

    """
    This function calls the calculation of integrated intensities for various passbands.
    
    :param dict model_dict: dictionary containing the model wavelengths, mu values and specific intensities
    ..note:: dictionary keys are 'mu' (1D numpy array), 'wavelengths' and 'intensities' (quantity array)
    :param dict passbands_dict: dictionary with the spectral responses to use
    ..note:: each passband contains a list with 2 quantity arrays reporting wavelengths and the corresponding pce response
    :return: dictionary with all the passband-integrated normalised intensities and the corresponding mu values
    :rtype: dict
    """
    passbands = list(passbands_dict.keys())
    model_wavelengths = model_dict['wavelengths']
    model_intensities = model_dict['intensities']
    model_mu = model_dict['mu']
    norm_index = np.argmax(model_mu)
    integ_dict = {}
    integ_dict['mu'] = copy.deepcopy(model_mu)
    for passband in passbands:
        my_waves = passbands_dict[passband][0].to(model_wavelengths.unit)
        my_pce = passbands_dict[passband][1].value
        f_interp = interp1d(model_wavelengths.value, model_intensities.value, axis=0, fill_value='extrapolate')
        my_ints = f_interp(my_waves.value)
        my_ints_integ = simps(my_ints*my_pce[:,None]*my_waves.value[:,None], my_waves.value, axis=0)
        integ_dict[passband] = my_ints_integ/my_ints_integ[norm_index]
    return integ_dict



def rescale_and_weights(mu, intensities, stellar_models_grid, user_geometry='pp', user_r_cut=0.99623):
    """
    This function finds the drop-off in the spherical intensity models, removes the values after the drop-off and rescales the other mu/radi values,
    then applies the quasi-spherical cut-off and computes weights for the intensity model-fit.
    
    :param np.array mu: 1D array with mu values calculated for the model-atmosphere
    :param np.array intensities: 1D array with passband-integrated normalised intensities
    :params str stellar_models_grid: the name of the selected stellar database
    :return: mu values, intensities and weights for the limb-darkening fit
    :rtype: np.array, np.array, np.array
    """
    radi = np.sqrt(1.0-mu**2.0)
    if stellar_models_grid in ['Phoenix_2018', 'Phoenix_2012_13', 'Phoenix_drift_2012']:
        #mu increasing, r decreasing
        dint_dr = np.abs( (intensities[1:]-intensities[:-1])/(radi[1:]-radi[:-1]) )
        rmax_dr = 0.5*( radi[np.argmax(dint_dr)+1] + radi[np.argmax(dint_dr)] )
        radicut_dr = radi[np.where(radi<=rmax_dr)]/rmax_dr
        intensitiescut_dr = intensities[np.where(radi<=rmax_dr)]
        radicut_dr_qs = radicut_dr[np.where(radicut_dr<=0.99623)]
        intensitiescut_dr_qs = intensitiescut_dr[np.where(radicut_dr<=0.99623)]
        mucut_dr_qs = np.sqrt(1.0-radicut_dr_qs**2.0)
        weights_dr_qs = np.zeros_like(radicut_dr_qs)
        weights_dr_qs[1:-1] = -0.5*(radicut_dr_qs[2:]-radicut_dr_qs[:-2])
        weights_dr_qs[0] = (1.0-radicut_dr_qs[0]) - 0.5*(radicut_dr_qs[1]-radicut_dr_qs[0])
        weights_dr_qs[-1] = 0.5*radicut_dr_qs[-2]
        return mucut_dr_qs, intensitiescut_dr_qs, weights_dr_qs
    elif stellar_models_grid in ['Atlas_2000']:
        #mu decreasing, r increasing
        weights_dr = np.zeros_like(radi)
        weights_dr[1:-1] = 0.5*(radi[2:]-radi[:-2])
        weights_dr[0] = 0.5*radi[1]
        weights_dr[-1] = (1.0-radi[-1]) + 0.5*(radi[-1]-radi[-2])
        return mu, intensities, weights_dr
    elif stellar_models_grid in ['Stagger_2018', 'Stagger_2015']:
        #mu increasing, r decreasing
        weights_dr = np.zeros_like(radi)
        weights_dr[1:-1] = -0.5*(radi[2:]-radi[:-2])
        weights_dr[0] = (1.0-radi[0]) - 0.5*(radi[1]-radi[0])
        weights_dr[-1] = 0.5*radi[-2]
        return mu, intensities, weights_dr
    elif stellar_models_grid in ['Userfile']:
        return mu, intensities, weights_dr





def get_limb_darkening_coefficients(integ_dict, limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders):
    """
    This function computes the limb-darkening coefficients for the required passbands and limb-darkening laws.
    
    :param dict integ_dict: dictionary with the mu values and the passband-integrated normalised intensities
    :param list of str limb-darkening laws: list of limb-darkening laws chosen by the user
    :param str stellar_models_grid: the name of the selected stellar database
    :param list of int gen_poly_orders: maximum orders of the polynomial limb-darkening laws chosen by the user
    :param list of int gen_claret_orders: maximum orders of the claret-n limb-darkening laws chosen by the user
    :return: dictionary containing the limb-darkening coefficients and the weighted root mean square of the fitting residuals
    :rtype: dict
    """
    mu = integ_dict['mu']
    passbands = [key for key in list(integ_dict.keys()) if key!='mu']
    ldc_dict = {}
    ldc_dict['passbands'] = {}
    for passband in passbands:
        integ_ints = integ_dict[passband]
        ldc_dict['passbands'][passband] = {}
        [res_mu, res_integ_ints, weights] = rescale_and_weights(mu, integ_ints, stellar_models_grid)
        ldc_dict['passbands'][passband]['rescaled_mu'] = res_mu
        ldc_dict['passbands'][passband]['rescaled_intensities'] = res_integ_ints
        ldc_dict['passbands'][passband]['weights'] = weights
        ldc_dict['passbands'][passband]['laws'] = {}
        if 'claret4' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['claret4'] = {}
            conv = False
            params = np.array((0.9, -0.5, 0.9, -0.5))
            while conv==False:
                ldc_claret4 = minimize(get_wrms_claret4, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_claret4.success
                params = ldc_claret4.x
            ldc_dict['passbands'][passband]['laws']['claret4']['coefficients'] = ldc_claret4.x
            ldc_dict['passbands'][passband]['laws']['claret4']['weighted_rms_res'] = np.sqrt(ldc_claret4.fun)
        if 'power2' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['power2'] = {}
            conv = False
            params = np.array((0.5, 1.0))
            while conv==False:
                ldc_power2 = minimize(get_wrms_power2, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_power2.success
                params = ldc_power2.x
            ldc_dict['passbands'][passband]['laws']['power2']['coefficients'] = ldc_power2.x
            ldc_dict['passbands'][passband]['laws']['power2']['weighted_rms_res'] = np.sqrt(ldc_power2.fun)
        if 'square_root' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['square_root'] = {}
            conv = False
            params = np.array((0.9, -0.2))
            while conv==False:
                ldc_square_root = minimize(get_wrms_square_root, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_square_root.success
                params = ldc_square_root.x
            ldc_dict['passbands'][passband]['laws']['square_root']['coefficients'] = ldc_square_root.x
            ldc_dict['passbands'][passband]['laws']['square_root']['weighted_rms_res'] = np.sqrt(ldc_square_root.fun)
        if 'quadratic' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['quadratic'] = {}
            conv = False
            params = np.array((0.9, -0.2))
            while conv==False:
                ldc_quadratic = minimize(get_wrms_quadratic, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_quadratic.success
                params = ldc_quadratic.x
            ldc_dict['passbands'][passband]['laws']['quadratic']['coefficients'] = ldc_quadratic.x
            ldc_dict['passbands'][passband]['laws']['quadratic']['weighted_rms_res'] = np.sqrt(ldc_quadratic.fun)
        if 'linear' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['linear'] = {}
            conv = False
            params = np.array((0.5))
            while conv==False:
                ldc_linear = minimize(get_wrms_linear, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_linear.success
                params = ldc_linear.x
            ldc_dict['passbands'][passband]['laws']['linear']['coefficients'] = ldc_linear.x
            ldc_dict['passbands'][passband]['laws']['linear']['weighted_rms_res'] = np.sqrt(ldc_linear.fun)
        if 'gen_poly' in limb_darkening_laws:
            for n in gen_poly_orders:
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)] = {}
                conv = False
                params = 0.5*np.ones(n)/n
                while conv==False:
                    ldc_gen_poly = minimize(get_wrms_gen_poly, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                    conv = ldc_gen_poly.success
                    params = ldc_gen_poly.x
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)]['coefficients'] = ldc_gen_poly.x
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)]['weighted_rms_res'] = np.sqrt(ldc_gen_poly.fun)
        if 'gen_claret' in limb_darkening_laws:
            for n in gen_claret_orders:
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)] = {}
                conv = False
                params = 0.5*np.ones(n)/n
                while conv==False:
                    ldc_gen_claret = minimize(get_wrms_gen_claret, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xatol':1e-8,'maxfev':10000,'disp':False})
                    conv = ldc_gen_claret.success
                    params = ldc_gen_claret.x
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)]['coefficients'] = ldc_gen_claret.x
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)]['weighted_rms_res'] = np.sqrt(ldc_gen_claret.fun)
    return ldc_dict
        

def interp_ldc(teff, logg, mh, neigh_indices, passband, law, neighbour_limb_darkening_coefficients, grid_files):
    """
    This function interpolates the limb-darkening coefficients by sequential linear interpolation.
    
    :param float teff: stellar effective temperature of the requested model
    :param float logg: stellar surface log gravity of the requested model
    :param float mh: stellar metallicity of the requested model
    :param np.array neigh_indices: indices of the neighbour models in the online database
    :param str passband: name of the passband (eventually, it includes the wavelength bin limits)
    :param str law: name of the limb-darkening law
    :param dict neighbour_limb_darkening_coefficients: dictionary containing the limb-darkening coefficients for the neighbour models
    
    :return: the interpolated limb-darkening coefficients and weighted rms for the requested model
    :rtype: np.array, float
    """
    coeffs_tgm = np.array([])
    wres_tgm = np.atleast_1d(np.array([]))
    star_params_tgm = np.array([])
    for i in neigh_indices:
        label = os.path.splitext(grid_files[i])[0]
        coeffs_tgm = my_vstack(coeffs_tgm, neighbour_limb_darkening_coefficients[label]['passbands'][passband]['laws'][law]['coefficients'])
        wres_tgm = my_vstack(wres_tgm, np.atleast_1d(neighbour_limb_darkening_coefficients[label]['passbands'][passband]['laws'][law]['weighted_rms_res']))
        star_params_tgm = my_vstack(star_params_tgm, neighbour_limb_darkening_coefficients[label]['star_params'])
    #star metallicity interpolation
    coeffs_tg = np.array([])
    wres_tg = np.atleast_1d(np.array([]))
    star_params_tg = np.array([])
    for i in [0, 2, 4, 6]:
        w2 = star_params_tgm[i,2]-mh
        w1 = mh-star_params_tgm[i+1,2]
        if w1+w2==0:
            coeffs_tg = my_vstack(coeffs_tg, coeffs_tgm[i,:])
            wres_tg = my_vstack(wres_tg, np.atleast_1d(wres_tgm[i]))
            star_params_tg = my_vstack(star_params_tg, star_params_tgm[i,:])
        else:
            w_coeffs = (w1*coeffs_tgm[i,:]+w2*coeffs_tgm[i+1,:])/(w1+w2)
            w_wres = (w1*wres_tgm[i]+w2*wres_tgm[i+1])/(w1+w2)
            w_star_params = (w1*star_params_tgm[i,:]+w2*star_params_tgm[i+1,:])/(w1+w2)
            coeffs_tg = my_vstack(coeffs_tg, w_coeffs)
            wres_tg = my_vstack(wres_tg, np.atleast_1d(w_wres))
            star_params_tg = my_vstack(star_params_tg, w_star_params)
    #star log gravity interpolation
    coeffs_t = np.array([])
    wres_t = np.atleast_1d(np.array([]))
    star_params_t = np.array([])
    for i in [0, 2]:
        w2 = star_params_tg[i,1]-logg
        w1 = logg-star_params_tg[i+1,1]
        if w1+w2==0:
            coeffs_t = my_vstack(coeffs_t, coeffs_tg[i,:])
            wres_t = my_vstack(wres_t, np.atleast_1d(wres_tg[i]))
            star_params_t = my_vstack(star_params_t, star_params_tg[i,:])
        else:
            w_coeffs = (w1*coeffs_tg[i,:]+w2*coeffs_tg[i+1,:])/(w1+w2)
            w_wres = (w1*wres_tg[i]+w2*wres_tg[i+1])/(w1+w2)
            w_star_params = (w1*star_params_tg[i,:]+w2*star_params_tg[i+1,:])/(w1+w2)
            coeffs_t = my_vstack(coeffs_t, w_coeffs)
            wres_t = my_vstack(wres_t, np.atleast_1d(w_wres))
            star_params_t = my_vstack(star_params_t, w_star_params)
    #star temperature interpolation
    coeffs = np.array([])
    wres = np.atleast_1d(np.array([]))
    star_params = np.array([])
    w2 = star_params_t[0,0]-teff
    w1 = teff-star_params_t[1,0]
    if w1+w2==0:
        coeffs = my_vstack(coeffs, coeffs_t[0,:])
        wres = my_vstack(wres, np.atleast_1d(wres_t[0]))
    else:
        coeffs = (w1*coeffs_t[0,:]+w2*coeffs_t[1,:])/(w1+w2)
        wres = (w1*wres_t[0]+w2*wres_t[1])/(w1+w2)
    return coeffs, wres



def process_configuration(input_dict):
    """
    This function executes all the operations from the checked input dictionary to create the requested output files.
    
    :param dict: input_dict
    ..note:: all the input dictionary branches end with a list, even if they contain a single element.
    ..note:: the input dictionary may contain more keywords than those specified in the configuration file, as they can be set to default values by the check_configuration function.
    :return: two dictionaries with limb-darkening coefficients (basic output file) and model passband-integrated intensities (saved only if user_output = 'complete') 
    :rtype: dict, dict
    """

    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())

    #Getting the output keywords and calculation type
    output_path = input_dict_local['output_path'][0]
    user_output = input_dict_local['user_output'][0]
    calculation_type = input_dict_local['calculation_type'][0]

    #Getting the limb-darkening laws
    limb_darkening_laws = input_dict_local['limb_darkening_laws']
    gen_poly_orders = None
    gen_claret_orders = None
    if 'gen_poly' in limb_darkening_laws:
        gen_poly_orders = input_dict_local['gen_poly_orders']
        gen_poly_orders = np.asarray(gen_poly_orders, dtype=int)
    if 'gen_claret' in limb_darkening_laws:
        gen_claret_orders = input_dict_local['gen_claret_orders']
        gen_claret_orders = np.asarray(gen_claret_orders, dtype=int)

    #Getting the coupled passband+wavelength_bins keywords
    stellar_models_grid = input_dict_local['stellar_models_grid'][0]
    passbands_path = input_dict_local['passbands_path'][0]
    passbands_ext = input_dict_local['passbands_ext'][0]
    passbands = input_dict_local['passbands']
    n_pass = len(passbands)
    if 'wavelength_bins_files' in input_keys:
        wavelength_bins_files = input_dict_local['wavelength_bins_files']
        wavelength_bins_path = input_dict_local['wavelength_bins_path'][0] 
    else:
        wavelength_bins_files = ['no_bins',]*n_pass
        wavelength_bins_path = ''
    #Computing the photon conversion efficiency for each requested passband + wavelength bins and storing into a dictionary
    passbands_dict = {}
    for i in range(n_pass):
        check = True
        [pb_wavelengths, pb_pce, check] = read_passband(passbands_path, passbands_ext, passbands[i], stellar_models_grid) #user or built-in file
        if check:
            f_interp = interp1d(pb_wavelengths.value, pb_pce.value, fill_value='extrapolate')
            my_waves = get_waves_fromR(np.min(pb_wavelengths.value), np.max(pb_wavelengths.value), R=1e6) * pb_wavelengths.unit #choice of wavelengths
            my_pce = f_interp(my_waves.value) * pb_pce.unit
            passbands_dict[passbands[i]] = [my_waves, my_pce]
            [wb, check] = read_wavelength_bins(wavelength_bins_path, wavelength_bins_files[i], pb_wavelengths, passbands[i]) #includes no_bins case
            if check:
                for wbin in wb:
                    my_waves = get_waves_fromR(wbin[0].value, wbin[1].value, R=1e6) * pb_wavelengths.unit #choice of wavelengths
                    my_pce = f_interp(my_waves.value) * pb_pce.unit
                    passbands_dict[passbands[i]+'_'+str(wbin[0].value)+'_'+str(wbin[1].value)] = [my_waves, my_pce]

    [star_files_grid, star_params_grid] = get_grid_parameters(stellar_models_grid) #Reading file names and stellar parameters from the selected database
    #Individual calculation_type
    if calculation_type=='individual':
        #Obtaining the requested stellar parameters and target names
        [star_effective_temperature, star_log_gravity, star_metallicity, target_names] = get_individual_parameters(input_dict_local)
        n_targets = len(star_effective_temperature)
        #Finding neighbours and eliminating out-of-range requests.
        neighbour_files_indices = np.array([],dtype=int)
        indices_to_delete = []
        for i in range(n_targets):
            neigh_indices = get_neighbour_files_indices(target_names[i],star_effective_temperature[i],star_log_gravity[i],star_metallicity[i],star_params_grid,stellar_models_grid)
            neighbour_files_indices = my_vstack(neighbour_files_indices, neigh_indices)
            if len(neigh_indices)==0:
                indices_to_delete += [i,]
        neighbour_files_indices = np.atleast_2d(neighbour_files_indices)
        star_effective_temperature = np.delete(star_effective_temperature, indices_to_delete)
        star_log_gravity = np.delete(star_log_gravity, indices_to_delete)
        star_metallicity = np.delete(star_metallicity, indices_to_delete)
        target_names = np.delete(target_names, indices_to_delete)
        n_targets = len(star_effective_temperature)
        #Reading neighbours files from the stellar database, and computing the corresponding limb-darkening coefficients
        neighbour_files_indices_list = np.unique(neighbour_files_indices)
        if len(neighbour_files_indices_list)==0:
            print('ERROR: No legal targets to calculate.')
            exit()
        neighbour_model_intensities = {}
        neighbour_passband_intensities = {}
        neighbour_limb_darkening_coefficients = {}
        for i in neighbour_files_indices_list:
            label = os.path.splitext(star_files_grid[i])[0]
            neighbour_model_intensities[label] = databases[stellar_models_grid].get_file_content(dbx_file=star_files_grid[i])
            neighbour_passband_intensities[label] = get_passband_intensities(neighbour_model_intensities[label], passbands_dict)
            neighbour_limb_darkening_coefficients[label] = get_limb_darkening_coefficients(neighbour_passband_intensities[label], limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders)
            neighbour_limb_darkening_coefficients[label]['star_params'] = star_params_grid[i]
        #Checking that the list of calculated passbands is not empty
        final_passbands = list(neighbour_limb_darkening_coefficients[label]['passbands'].keys())
        if len(final_passbands)==0:
            print('ERROR: No legal passbands to calculate.')
            exit()
        #Creating output dictionary with interpolated limb-darkening coefficients
        target_limb_darkening_coefficients = {}
        final_limb_darkening_laws = list(neighbour_limb_darkening_coefficients[label]['passbands'][final_passbands[0]]['laws'].keys())
        for i in range(n_targets):
            target_limb_darkening_coefficients[target_names[i]] = {}
            target_limb_darkening_coefficients[target_names[i]]['star_params'] = np.array([star_effective_temperature[i], star_log_gravity[i], star_metallicity[i]])
            target_limb_darkening_coefficients[target_names[i]]['passbands'] = {}
            for passband in final_passbands:
                target_limb_darkening_coefficients[target_names[i]]['passbands'][passband] = {}
                for law in final_limb_darkening_laws:
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law] = {}
                    [coeffs, w_res] = interp_ldc(star_effective_temperature[i], star_log_gravity[i], star_metallicity[i], neighbour_files_indices[i,:], passband, law, neighbour_limb_darkening_coefficients, star_files_grid)
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law]['coefficients'] = coeffs
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law]['weighted_rms_res'] = w_res
        limb_darkening_coefficients = {}
        limb_darkening_coefficients['neighbour'] = neighbour_limb_darkening_coefficients
        limb_darkening_coefficients['target'] = target_limb_darkening_coefficients
        for i in range(n_targets):
            with open(os.path.join(output_path, target_names[i]+'_ldc.pickle') , 'wb') as out1:
                pickle.dump(target_limb_darkening_coefficients[target_names[i]], out1, protocol=pickle.HIGHEST_PROTOCOL)
        if user_output == 'complete':
            run_name = ''
            for i in range(n_targets):
                run_name += target_names[i]
            with open(os.path.join(output_path, run_name+'_neighbour_ldc.pickle') , 'wb') as out2:
                pickle.dump(neighbour_limb_darkening_coefficients, out2, protocol=pickle.HIGHEST_PROTOCOL)
            with open(os.path.join(output_path, run_name+'_neighbour_intensities.pickle') , 'wb') as out3:
                pickle.dump(neighbour_passband_intensities, out3, protocol=pickle.HIGHEST_PROTOCOL)
        return limb_darkening_coefficients, neighbour_passband_intensities
    #Grid calculation_type
    elif calculation_type=='grid':
        #Obtaining star file names, parameters and file indices for the selected database
        [star_files_subgrid, star_params_subgrid, subgrid_indices] = get_subgrid(input_dict_local)
        n_files = len(star_files_subgrid)
        model_intensities = {}
        passband_intensities = {}
        limb_darkening_coefficients = {}
        #Calculation limb-darkening and output dictionary
        for i in range(n_files):
            #label = 'teff' + str(star_params_subgrid[i,0]) + '_logg' + str(star_params_subgrid[i,1]) + '_MH' + str(star_params_subgrid[i,2])
            label = os.path.splitext(star_files_subgrid[i])[0]
            model_intensities[label] = databases[stellar_models_grid].get_file_content(dbx_file=star_files_grid[subgrid_indices[i]])
            passband_intensities[label] = {}
            passband_intensities[label] = get_passband_intensities(model_intensities[label], passbands_dict)
            limb_darkening_coefficients[label] = get_limb_darkening_coefficients(passband_intensities[label], limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders)
            limb_darkening_coefficients[label]['star_params'] = star_params_subgrid[i]
        with open(os.path.join(output_path, 'grid_ldc.pickle'), 'wb') as out1:
            pickle.dump(limb_darkening_coefficients, out1, protocol=pickle.HIGHEST_PROTOCOL)
        if user_output == 'complete':
            with open(os.path.join(output_path, 'grid_intensities.pickle'), 'wb') as out2:
                pickle.dump(passband_intensities, out2, protocol=pickle.HIGHEST_PROTOCOL)
        return limb_darkening_coefficients, passband_intensities


def ldc_calculate(configuration_file):
    """
    This is the main function to run for calculating the limb-darkening coefficients from the input file.
    It calls the three functions to read, check and process the input file.
    
    :param str configuration_file: absolute or relative path including the file name
    
    :return:
    :rtype:
    """
    input_dict1 = read_configuration(configuration_file) #Reading configuration file
    [check, input_dict] = check_configuration(input_dict1) #Checking configuration file and partial dictionary update
    if check:
        process_configuration(input_dict) #Computing and saving
    else:
        print('Something went wrong with the input.')
