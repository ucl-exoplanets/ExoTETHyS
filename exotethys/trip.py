from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys

if sys.version_info[0] > 2:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve
    input = raw_input

import matplotlib
if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
else:
    matplotlib.use('TkAgg')

import glob
import time
import shutil
import astropy.io.fits as pyfits
import numpy as np
import os
from scipy.optimize import minimize
from scipy.interpolate import LinearNDInterpolator as interp
from scipy.interpolate import interp1d

import copy
import pickle

from ._1databases import *


def str2float(s): #function in common with SAIL
#This function converts a string to a float, if the string is a number.
    try:
        return float(s)
    except ValueError:
        return s

def check_length(vector, min_length=1, max_length=None): #function in common with SAIL
#This function checks if the length of a vector is within the expected range and returns a boolean value.
#Default is length>=1.
    check = True
    if len(vector)<min_length:
        check = False
    elif max_length:
        if len(vector)>max_length:
            check = False
    return check

def check_type_in_list(vector, item_type): #function in common with SAIL
#This function checks that all the elements of a list are of the expected variable type (e.g., string, float, int) and returns a boolean value.
    check = True
    for item in vector:
        if not isinstance(item, item_type):
            check = False
    return check

def check_integers(vector, min_value=1): #function in common with SAIL
#This function checks if all the elements of a vector correspond to integer numbers greater than or equal to min_value and returns a boolean value.
#Default is min_value=1 (positive integers).
    check = True
    if np.min(vector)<min_value:
        check = False
    vector_float = np.asarray(vector)
    vector_int = np.asarray(vector,int)
    if (vector_float != vector_int).any():
        check = False
    return check

def read_as_numpy_array(path_file): #function in common with SAIL
#This function reads a numpy array from file as numpy.genfromtxt but also checks for the most common errors.
#It returns a numpy array and a boolean value.
#If the boolean value is False, the numpy array will be empty.
#If the error is caused by one or more elements of the array that are not read as numbers, a warning message will inform about the rows where these errors occur.
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

def check_2Darray(arr, n_col=None): #function in common with SAIL
#This function checks that an array has 2 dimensions and returns a boolean value.
#Optionally, it may also require the exact number of columns. 
    check = True
    if arr.ndim != 2:
        check = False
    elif n_col:
        if np.shape(arr)[1] != n_col:
            check = False
    return check

def check_1Darray(arr):
#This function checks that an array has 1 dimensions and returns a boolean value.
#Optionally, it may also require the exact number of columns. 
    check = True
    if arr.ndim != 1:
        check = False
    return check




def get_default_value(param):
#This function returns the default value for the requested parameter.
    if param=='rp_over_rs':
        default_value = None
    if param=='sma_over_rs':
        default_value = None
    if param=='inclination':
        default_value = None
    if param=='eccentricity':
        default_value = 0.0
    if param=='arg_pericenter':
        default_value = 0.0
    if param=='period_orbital':
        default_value = None
    if param=='epoch_of_transit':
        default_value = None
    if param=='time_conversion_factor':
        default_value = 1.0
    if param=='n_annuli':
        default_value = 100000
    if param=='interpolation_type':
        default_value = 'linear'
    if param=='interpolation_variable':
        default_value = 'mu'
    if param=='cutting_limb':
        default_value = 'no_cut'
    if param=='rescaling_limb':
        default_value = 'as_cut'
    if param=='rescaling_input_params':
        default_value = 'no'
    if param=='output_path':
        default_value = ''
    if param=='output_filename':
        default_value = None
    if param=='output_fileext':
        default_value = '.pickle'
    if param=='input_limb_path':
        default_value = ''
    if param=='input_series_path':
        default_value = ''
    return default_value



def get_params_limits(param):
#This function returns the limits for the requested parameter (used to compute the projected star-planet separation, z_sep), together with an explanatory string.
    if param=='sma_over_rs':
        min_value = 1.0
        max_value = None
        inf_value = None
        sup_value = None
        sentence = 'sma_over_rs>=1.0'
    if param=='inclination':
        min_value = None
        max_value = 90.0
        inf_value = 0.0
        sup_value = None
        sentence = '0.0<inclination<=90.0'
    if param=='eccentricity':
        min_value = 0.0
        max_value = None
        inf_value = None
        sup_value = 1.0
        sentence = '0.0<=eccentricity<1.0'
    if param=='arg_pericenter':
        min_value = None
        max_value = None
        inf_value = -180.0
        sup_value = 360.0
        sentence = '-180.0<arg_pericenter<360.0'
    if param=='period_orbital':
        min_value = None
        max_value = None
        inf_value = 0.0
        sup_value = None
        sentence = 'period_orbital>0.0'
    if param=='epoch_of_transit':
        min_value = None
        max_value = None
        inf_value = None
        sup_value = None
        sentence = 'float.'
    if param=='time_conversion_factor':
        min_value = None
        max_value = None
        inf_value = 0.0
        sup_value = None
        sentence = 'time_conversion_factor>0.0'
    return min_value, max_value, inf_value, sup_value, sentence


def check_values(vector, min_value=None, max_value=None, inf_value=None, sup_value=None):
#This function checks that all the elements of a list/array are float within given limits (if any). It returns a boolean value.
    check = True
    for item in vector:
        if not isinstance(item, float):
            check = False
        cond1 = (min_value is not None) and item<min_value
        cond2 = (inf_value is not None) and item<=inf_value
        cond3 = (max_value is not None) and item>max_value
        cond4 = (sup_value is not None) and item>=sup_value
        if cond1 or cond2 or cond3 or cond4:
            check = False
    return check




def read_configuration(filename):
#This function reads the input file line by line and returns a dictionary.
#For each line the first word is a keyword, the following are values (either string or float).
#The lines starting with '#' will be ignored. The values preceded by '!' (without spaces) will be also ignored.
    with open(filename, 'r') as file:
        input_dict = {}
        for line in file.readlines():
            content = line.split()
            if line[0] != '#':
                #print(content)
                key = content[0]
                value = []
                for item in content[1:]:
                    if item[0] != '!':
                        value += [str2float(item),]
                input_dict[key] = value
    return input_dict


def check_configuration(input_dict):
#This function checks and modifies the dictionary obtained from the configuration file and; it returns a boolean value and the updated dictionary.
    check = True
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    mandatory_keys = ['input_limb_type', 'input_limb_path', 'input_limb_file', 'input_series_type', 'input_series_path', 'input_series_file', 'rp_over_rs']
    allowed_keys = mandatory_keys + ['output_path', 'output_filename', 'output_fileext', 'sma_over_rs', 'inclination', 'eccentricity', 'arg_pericenter', 'period_orbital', 'epoch_of_transit', 'time_conversion_factor', 'n_annuli', 'interpolation_type', 'interpolation_variable', 'cutting_limb', 'user_cut_mu', 'user_cut_radi', 'rescaling_limb', 'user_rescale_mu', 'user_rescale_radi', 'rescaling_input_params']

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

    #Checking the (mandatory) input_limb_type: 'mu' or 'radi'.
    input_limb_type = input_dict_local['input_limb_type']
    if not check_length(input_limb_type, max_length=1):
        print('ERROR: invalid length=', len(input_limb_type), 'for input_limb_type. It must have length=1.')
        check = False
    else:
        input_limb_type = input_limb_type[0]
        allowed_input_limb_types = ['mu', 'radi']
        if input_limb_type not in allowed_input_limb_types:
            print('ERROR: invalid input_limb_type. It must be either mu or radi.')
            check = False

    #Checking the (mandatory) input_limb_file.
    input_limb_file = input_dict_local['input_limb_file']
    if not check_length(input_limb_file, max_length=1):
        print('ERROR: invalid length=', len(input_limb_file), 'for input_limb_file. It must have length=1.')
        check = False
    if not check_type_in_list(input_limb_file,str):
        print('ERROR: invalid type for input_limb_file. It must be string.')
        check = False

    #Checking the (mandatory) input_series_type: 'phi, 'time' or 'z_sep'.
    input_series_type = input_dict_local['input_series_type']
    if not check_length(input_series_type, max_length=1):
        print('ERROR: invalid length=', len(input_series_type), 'for input_series_type. It must have length=1.')
        check = False
    else:
        input_series_type = input_series_type[0]
        allowed_input_series_types = ['phi', 'time', 'z_sep']
        if input_series_type not in allowed_input_series_types:
            print('ERROR: invalid input_series_type. It must be either phi, time or z_sep.')
            check = False

    #Checking the (mandatory) input_series_file.
    input_series_file = input_dict_local['input_series_file']
    if not check_length(input_series_file, max_length=1):
        print('ERROR: invalid length=', len(input_series_file), 'for input_series_file. It must have length=1.')
        check = False
    if not check_type_in_list(input_series_file,str):
        print('ERROR: invalid type for input_series_file. It must be string.')
        check = False

    #Checking the (mandatory) rp_over_rs.
    rp_over_rs = input_dict_local['rp_over_rs']
    if not check_length(rp_over_rs, max_length=1):
        print('ERROR: invalid length=', len(rp_over_rs), 'for rp_over_rs. It must have length=1.')
        check = False
    if not check_values(rp_over_rs,min_value=0.0):
        print('ERROR: invalid value for rp_over_rs. It must be rp_over_rs>=0.')
        check = False
    elif not check_values(rp_over_rs,max_value=1.0):
        print('WARNING: rp_over_rs>1.0, planet is larger than star.')

    #Checking keywords for input TIME series
    if input_series_type=='time':
        mandatory_keys_time = ['sma_over_rs', 'inclination', 'period_orbital', 'epoch_of_transit']
        forbidden_keys_time = []
        #other_allowed_keys_time = ['eccentricity', 'arg_pericenter', 'time_conversion_factor']
        allowed_keys_time = mandatory_keys_time + ['eccentricity', 'arg_pericenter', 'time_conversion_factor']
        for key in mandatory_keys_time:
            if key not in input_keys:
                print('ERROR: mandatory keyword (for input_series_type=time)', key, 'is not specified.')
                check = False
        for key in forbidden_keys_time:
            if key in input_keys:
                print('ERROR: keyword', key, 'is not valid for input_series_type=time.')
                check = False
        for key in allowed_keys_time:
            if key in input_keys:
                key_prov = input_dict_local[key]
                if not check_length(key_prov, max_length=1):
                    print('ERROR: invalid length=', len(key_prov), 'for', key, '. It must have length=1.')
                    check = False
                [min_key_prov, max_key_prov, inf_key_prov, sup_key_prov, sentence] = get_params_limits(key)
                if not check_values(key_prov, min_value=min_key_prov, max_value=max_key_prov, inf_value=inf_key_prov, sup_value=sup_key_prov):
                    print('ERROR: invalid value for', key, '. It must be '+sentence+'.')
                    check = False
            else:
                input_dict_local[key] = [get_default_value(key)]


    #Checking keywords for input PHI series
    if input_series_type=='phi':
        mandatory_keys_phi = ['sma_over_rs', 'inclination']
        forbidden_keys_phi = ['period_orbital', 'epoch_of_transit', 'time_conversion_factor']
        #other_allowed_keys_phi = ['eccentricity', 'arg_pericenter']
        allowed_keys_phi = mandatory_keys_phi + ['eccentricity', 'arg_pericenter']
        for key in mandatory_keys_phi:
            if key not in input_keys:
                print('ERROR: mandatory keyword (for input_series_type=phi)', key, 'is not specified.')
                check = False
        for key in forbidden_keys_phi:
            if key in input_keys:
                print('ERROR: keyword', key, 'is not valid for input_series_type=phi.')
                check = False
        for key in allowed_keys_phi:
            if key in input_keys:
                key_prov = input_dict_local[key]
                if not check_length(key_prov, max_length=1):
                    print('ERROR: invalid length=', len(key_prov), 'for', key, '. It must have length=1.')
                    check = False
                [min_key_prov, max_key_prov, inf_key_prov, sup_key_prov, sentence] = get_params_limits(key)
                if not check_values(key_prov, min_value=min_key_prov, max_value=max_key_prov, inf_value=inf_key_prov, sup_value=sup_key_prov):
                    print('ERROR: invalid value for', key, '. It must be '+sentence+'.')
                    check = False
            else:
                input_dict_local[key] = [get_default_value(key)]


    #Checking keywords for input Z_SEP series
    if input_series_type=='z_sep':
        mandatory_keys_z_sep = []
        forbidden_keys_z_sep = ['sma_over_rs', 'inclination', 'eccentricity', 'arg_pericenter', 'period_orbital', 'epoch_of_transit', 'time_conversion_factor']
        #other_allowed_keys_z_sep = []
        allowed_keys_z_sep = mandatory_keys_z_sep + []
        for key in mandatory_keys_z_sep:
            if key not in input_keys:
                print('ERROR: mandatory keyword (for input_series_type=z_sep)', key, 'is not specified.')
                check = False
        for key in forbidden_keys_z_sep:
            if key in input_keys:
                print('ERROR: keyword', key, 'is not valid for input_series_type=z_sep.')
                check = False
        for key in allowed_keys_z_sep:
            if key in input_keys:
                key_prov = input_dict_local[key]
                if not check_length(key_prov, max_length=1):
                    print('ERROR: invalid length=', len(key_prov), 'for', key, '. It must have length=1.')
                    check = False
                [min_key_prov, max_key_prov, inf_key_prov, sup_key_prov, sentence] = get_params_limits(key)
                if not check_values(key_prov, min_value=min_key_prov, max_value=max_key_prov, inf_value=inf_key_prov, sup_value=sup_key_prov):
                    print('ERROR: invalid value for', key, '. It must be '+sentence+'.')
                    check = False
            else:
                input_dict_local[key] = [get_default_value(key)]

    #Checking the (optional) input_limb_path.
    if 'input_limb_path' in input_keys:
        input_limb_path = input_dict_local['input_limb_path']
        if not check_length(input_limb_path, max_length=1):
            print('ERROR: invalid length=', len(input_limb_path), 'for input_limb_path. It must have length=1.')
            check = False
        if not check_type_in_list(input_limb_path,str):
            print('ERROR: invalid type for input_limb_path. It must be string.')
            check = False
    else:
        input_dict_local['input_limb_path'] = [get_default_value('input_limb_path')]

    #Checking the (optional) input_series_path.
    if 'input_series_path' in input_keys:
        input_series_path = input_dict_local['input_series_path']
        if not check_length(input_series_path, max_length=1):
            print('ERROR: invalid length=', len(input_series_path), 'for input_series_path. It must have length=1.')
            check = False
        if not check_type_in_list(input_series_path,str):
            print('ERROR: invalid type for input_series_path. It must be string.')
            check = False
    else:
        input_dict_local['input_series_path'] = [get_default_value('input_series_path')]

    #Checking the (optional) n_annuli
    if 'n_annuli' in input_keys:
        n_annuli = input_dict_local['n_annuli']
        if not check_length(n_annuli, max_length=1):
            print('ERROR: invalid length=', len(n_annuli), 'for n_annuli. It must have length=1.')
            check = False
        if not check_integers(n_annuli):
            print('ERROR: invalid type for n_annuli. It must be positive integer. Default is 100000.')
            check = False
    else:
        input_dict_local['n_annuli'] = [get_default_value('n_annuli')]

    #Checking interpolation_type
    if 'interpolation_type' in input_keys:
        interpolation_type = input_dict_local['interpolation_type']
        if not check_length(interpolation_type, max_length=1):
            print('ERROR: invalid length=', len(interpolation_type), 'for interpolation_type. It must have length=1.')
            check = False
        else:
            interpolation_type = input_dict_local['interpolation_type'][0]
            allowed_interpolation_types = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']
            if interpolation_type not in allowed_interpolation_types:
                print('ERROR: invalid interpolation_type. Possible types are', allowed_interpolation_types, '. Default is', get_default_value('interpolation_type'), '.')
                check = False
    else:
        input_dict_local['interpolation_type'] = [get_default_value('interpolation_type')]

    #Checking interpolation_variable
    if 'interpolation_variable' in input_keys:
        interpolation_variable = input_dict_local['interpolation_variable']
        if not check_length(interpolation_variable, max_length=1):
            print('ERROR: invalid length=', len(interpolation_variable), 'for interpolation_variable. It must have length=1.')
            check = False
        else:
            interpolation_variable = input_dict_local['interpolation_variable'][0]
            allowed_interpolation_variables = ['mu', 'radi']
            if interpolation_variable not in allowed_interpolation_variables:
                print('ERROR: invalid interpolation_variable. It must be mu or radi. Default is', get_default_value('interpolation_variable'), '.')
                check = False
    else:
        input_dict_local['interpolation_variable'] = [get_default_value('interpolation_variable')]

    #Checking cutting_limb
    if 'cutting_limb' in input_keys:
        cutting_limb = input_dict_local['cutting_limb']
        if not check_length(cutting_limb, max_length=1):
            print('ERROR: invalid length=', len(cutting_limb), 'for cutting_limb. It must have length=1.')
            check = False
        else:
            cutting_limb = input_dict_local['cutting_limb'][0]
            allowed_cutting_limb = ['no_cut', 'radi_gradient', 'mu_gradient', 'user_cut']
            if cutting_limb not in allowed_cutting_limb:
                print('ERROR: invalid cutting_limb. Possible options are', allowed_cutting_limb, '. Default is', get_default_value('cutting_limb'), '.')
                check = False
            elif cutting_limb=='user_cut':
                if ('user_cut_mu' not in input_keys) and ('user_cut_radi' not in input_keys):
                    print('ERROR: user_cut_mu or user_cut_radi must be given as input if cutting_limb=user_cut.')
                    check = False
                elif ('user_cut_mu' in input_keys) and ('user_cut_radi' in input_keys):
                    print('ERROR: only one between user_cut_mu and user_cut_radi can be given as input.')
                    check = False
            elif cutting_limb!='user_cut':
                if 'user_cut_mu' in input_keys:
                    print('ERROR: user_cut_mu is not a valid keyword unless cutting_limb is user_cut.')
                    check = False
                if 'user_cut_radi' in input_keys:
                    print('ERROR: user_cut_radi is not a valid keyword unless cutting_limb is user_cut.')
                    check = False
    else:
        input_dict_local['cutting_limb'] = [get_default_value('cutting_limb')]

    #Checking user_cut_mu
    if 'user_cut_mu' in input_keys:
        user_cut_mu = input_dict_local['user_cut_mu']
        if not check_length(user_cut_mu, max_length=1):
            print('ERROR: invalid length=', len(user_cut_mu), 'for user_cut_mu. It must have length=1.')
            check = False
        if not check_values(user_cut_mu,inf_value=0.0,sup_value=1.0):
            print('ERROR: invalid value for user_cut_mu. It must be 0.0<user_cut_mu<1.0.')
            check = False
        if check:
            user_cut_mu = input_dict_local['user_cut_mu'][0]
            user_cut_radi = np.sqrt(1.0-user_cut_mu**2.0)
            input_dict_local['user_cut_radi'] = [user_cut_radi]

    #Checking user_cut_radi
    if 'user_cut_radi' in input_keys:
        user_cut_radi = input_dict_local['user_cut_radi']
        if not check_length(user_cut_radi, max_length=1):
            print('ERROR: invalid length=', len(user_cut_radi), 'for user_cut_radi. It must have length=1.')
            check = False
        if not check_values(user_cut_radi,inf_value=0.0,sup_value=1.0):
            print('ERROR: invalid value for user_cut_radi. It must be 0.0<user_cut_radi<1.0.')
            check = False
        if check:
            user_cut_radi = input_dict_local['user_cut_radi'][0]
            user_cut_mu = np.sqrt(1.0-user_cut_radi**2.0)
            input_dict_local['user_cut_mu'] = [user_cut_mu]

    #Checking rescaling_limb
    if 'rescaling_limb' in input_keys:
        rescaling_limb = input_dict_local['rescaling_limb']
        if not check_length(rescaling_limb, max_length=1):
            print('ERROR: invalid length=', len(rescaling_limb), 'for rescaling_limb. It must have length=1.')
            check = False
        else:
            rescaling_limb = input_dict_local['rescaling_limb'][0]
            allowed_rescaling_limb = ['as_cut', 'no_rescale', 'user_rescale']
            if rescaling_limb not in allowed_rescaling_limb:
                print('ERROR: invalid rescaling_limb. Possible options are', allowed_rescaling_limb, '. Default is', get_default_value('rescaling_limb'), '.')
                check = False
            elif rescaling_limb=='user_rescale':
                if ('user_rescale_mu' not in input_keys) and ('user_rescale_radi' not in input_keys):
                    print('ERROR: user_rescale_mu or user_rescale_radi must be given as input if rescaling_limb=user_rescale.')
                    check = False
                elif ('user_rescale_mu' in input_keys) and ('user_rescale_radi' in input_keys):
                    print('ERROR: only one between user_rescale_mu and user_rescale_radi can be given as input.')
                    check = False
            elif rescaling_limb!='user_rescale':
                if 'user_rescale_mu' in input_keys:
                    print('ERROR: user_rescale_mu is not a valid keyword unless rescaling_limb is user_rescale.')
                    check = False
                if 'user_rescale_radi' in input_keys:
                    print('ERROR: user_rescale_radi is not a valid keyword unless rescaling_limb is user_rescale.')
                    check = False
    else:
        input_dict_local['rescaling_limb'] = [get_default_value('rescaling_limb')]

    #Checking user_rescale_mu
    if 'user_rescale_mu' in input_keys:
        user_rescale_mu = input_dict_local['user_rescale_mu']
        if not check_length(user_rescale_mu, max_length=1):
            print('ERROR: invalid length=', len(user_rescale_mu), 'for user_rescale_mu. It must have length=1.')
            check = False
        if not check_values(user_rescale_mu,inf_value=0.0,sup_value=1.0):
            print('ERROR: invalid value for user_rescale_mu. It must be 0.0<user_rescale_mu<1.0.')
            check = False
        if check:
            user_rescale_mu = input_dict_local['user_rescale_mu'][0]
            user_rescale_radi = np.sqrt(1.0-user_rescale_mu**2.0)
            input_dict_local['user_rescale_radi'] = [user_rescale_radi]

    #Checking user_rescale_radi
    if 'user_rescale_radi' in input_keys:
        user_rescale_radi = input_dict_local['user_rescale_radi']
        if not check_length(user_rescale_radi, max_length=1):
            print('ERROR: invalid length=', len(user_rescale_radi), 'for user_rescale_radi. It must have length=1.')
            check = False
        if not check_values(user_rescale_radi,inf_value=0.0,sup_value=1.0):
            print('ERROR: invalid value for user_rescale_radi. It must be 0.0<user_rescale_radi<1.0.')
            check = False
        if check:
            user_rescale_radi = input_dict_local['user_rescale_radi'][0]
            user_rescale_mu = np.sqrt(1.0-user_rescale_radi**2.0)
            input_dict_local['user_rescale_mu'] = [user_rescale_mu]


    #Checking the (optional) rescaling_input_params.
    if 'rescaling_input_params' in input_keys:
        rescaling_input_params = input_dict_local['rescaling_input_params']
        if not check_length(rescaling_input_params, max_length=1):
            print('ERROR: invalid length=', len(rescaling_input_params), 'for rescaling_input_params. It must have length=1.')
            check = False
        else:
            rescaling_input_params = rescaling_input_params[0]
            allowed_rescaling_input_params = ['no', 'yes']
            if rescaling_input_params not in allowed_rescaling_input_params:
                print('ERROR: invalid rescaling_input_params. It must be either no or yes. Default is', get_default_value('rescaling_input_params'), '.')
                check = False
    else:
        input_dict_local['rescaling_input_params'] = [get_default_value('rescaling_input_params')]

    #Checking the (optional) output_path.
    if 'output_path' in input_keys:
        output_path = input_dict_local['output_path']
        if not check_length(output_path, max_length=1):
            print('ERROR: invalid length=', len(output_path), 'for output_path. It must have length=1.')
            check = False
        if not check_type_in_list(output_path,str):
            print('ERROR: invalid type for output_path. It must be string.')
            check = False
    else:
        input_dict_local['output_path'] = [get_default_value('output_path')]

    #Checking the (optional) output_filename.
    if 'output_filename' in input_keys:
        output_filename = input_dict_local['output_filename']
        if not check_length(output_filename, max_length=1):
            print('ERROR: invalid length=', len(output_filename), 'for output_filename. It must have length=1.')
            check = False
        if not check_type_in_list(output_filename,str):
            print('ERROR: invalid type for output_filename. It must be string.')
            check = False
    else:
        input_dict_local['output_filename'] = [get_default_value('output_filename')]


    #Checking the (optional) output_fileext.
    if 'output_fileext' in input_keys:
        output_fileext = input_dict_local['output_fileext']
        if not check_length(output_fileext, max_length=2):
            print('ERROR: invalid length=', len(output_fileext), 'for output_fileext. It must have length<=2.')
            check = False
        else:
            allowed_output_fileext = ['.txt', '.pickle']
            for item in output_fileext:
                if item not in allowed_output_fileext:
                    print('ERROR: invalid output_fileext. It can be .pickle or .txt. Default is', get_default_value('output_fileext'), '.')
                    check = False
    else:
        input_dict_local['output_fileext'] = [get_default_value('output_fileext')]


    if check==False:
        exit()

    return check, input_dict_local



def compute_z_sep(phi, inclination, sma_over_rs, eccentricity, arg_pericenter):
#This function
	theta = 2.0 * np.pi * phi
	inclination = inclination * np.pi / 180.0
	arg_pericenter = arg_pericenter * np.pi / 180.0

	if eccentricity != 0:
		n = len(theta)
		E = np.zeros(n)
		ecc2 = np.sqrt((1.0+eccentricity)/(1.0-eccentricity))
		fref = np.pi/2.0 - omega #setting reference point for the true anomaly
		Eref = 2.0 * np.arctan(1.0/ecc2 * np.tan(fref/2.0))
		if Eref < -np.pi/2.0:
			Eref = Eref + 2.0*np.pi
		Mref = Eref - (ecc * np.sin(Eref))
		for i in range(n):
			Mtmp = theta[i] +  Mref
			Etmp = Mtmp
			for j in range(10):
				Etmp = Etmp + ((Mtmp + ecc*np.sin(Etmp) - Etmp) / (1.0-eccentricity*np.cos(Etmp)))
			E[i] = Etmp
		#calculating true anomaly
		f = 2.0*np.arctan(ecc2*np.tan(E/2.0))
		#calculating distance from true anomaly as fraction
		r_frac = (1.0-eccentricity**2)/(1.0 + eccentricity*np.cos(f))
		#computing z
		z_sep = 1.0 - ((np.sin(inclination)**2.0) * (np.sin(f+arg_pericenter)**2))
		z_sep = sma_over_rs*r_frac*np.sqrt(z_sep)
	
	if eccentricity == 0:
		z_sep = sma_over_rs * np.sqrt( np.ones(len(theta)) - (np.cos(theta) * np.sin(inclination))**2 )

	return z_sep



def get_x_series(input_series_type, input_series_path, input_series_file, sma_over_rs=None, inclination=None, eccentricity=None, arg_pericenter=None, period_orbital=None, epoch_of_transit=None, time_conversion_factor=None):
#This function
    if input_series_type=='z_sep':
        [z_sep, check] = read_as_numpy_array(input_series_path+input_series_file)
        if not check_1Darray(z_sep):
            print('ERROR: invalid format for input_series_file', input_series_path+input_series_file, '. It must have 1 column.')
            exit()
        if np.min(z_sep)<0.0:
            print('ERROR: found negative z_sep in input_series_file', input_series_path+input_series_file, '. Star-planet separation (z_sep) cannot be negative.')
            exit()
        return z_sep

    elif input_series_type=='phi':
        [phi, check] = read_as_numpy_array(input_series_path+input_series_file)
        if not check_1Darray(phi):
            print('ERROR: invalid format for input_series_file', input_series_path+input_series_file, '. It must have 1 column.')
            exit()
        if len(phi)!=len(set(phi)):
            print('ERROR: Duplicate', input_series_type,'value in input_series_file', input_series_path+input_series_file, '.')
            exit()
        z_sep = compute_z_sep(phi, inclination, sma_over_rs, eccentricity, arg_pericenter)
        return phi, z_sep

    elif input_series_type=='time':
        [time, check] = read_as_numpy_array(input_series_path+input_series_file)
        if not check_1Darray(time):
            print('ERROR: invalid format for input_series_file', input_series_path+input_series_file, '. It must have 1 column.')
            exit()
        if len(time)!=len(set(time)):
            print('ERROR: Duplicate', input_series_type,'value in input_series_file', input_series_path+input_series_file, '.')
            exit()
        phi = (time*time_conversion_factor - epoch_of_transit)/period_orbital
        n = round(np.median(phi))
        phi -= n
        z_sep = compute_z_sep(phi, inclination, sma_over_rs, eccentricity, arg_pericenter)
        return time, phi, z_sep





def get_input_limb_model(input_limb_type,input_limb_path,input_limb_file):
    [limb_model, check] = read_as_numpy_array(input_limb_path+input_limb_file)
    if not check_2Darray(limb_model, n_col=2):
        print('ERROR: invalid format for input_limb_file', input_limb_path+input_limb_file, '. It must have 2 columns.')
        exit()
    intensity_model = copy.deepcopy(limb_model[:,1])
    if np.min(intensity_model)<0.0:
        print('ERROR: found negative intensity in input_limb_file', input_limb_path+input_limb_file, '. Intensities cannot be negative.')
        check = False
    if np.min(intensity_model)==0.0 and np.max(intensity_model)==0.0:
        print('ERROR: All intensities are null in input_limb_file', input_limb_path+input_limb_file, '. At least some intensities must be positive.')
        check = False
    muorradi_model = copy.deepcopy(limb_model[:,0])
    if np.min(muorradi_model)<0.0 or np.max(muorradi_model)>1.0:
        print(muorradi_model)
        print('ERROR: Invalid', input_limb_type,'value in input_limb_file', input_limb_path+input_limb_file, '. It must be 0<=', input_limb_type, '<=1.')
        check = False
    if len(muorradi_model)!=len(set(muorradi_model)):
        print('ERROR: Duplicate', input_limb_type,'value in input_limb_file', input_limb_path+input_limb_file, '.')
        check = False
    if check == False:
        exit()
    if input_limb_type=='mu':
        mu_model = copy.deepcopy(limb_model[:,0])
        radi_model = np.sqrt(1.0-mu_model**2.0)
    elif input_limb_type=='radi':
        radi_model = copy.deepcopy(limb_model[:,0])
        mu_model = np.sqrt(1.0-radi_model**2.0)
    sorting_indices = np.argsort(radi_model)
    mu_model = mu_model[sorting_indices]
    radi_model = radi_model[sorting_indices]
    intensity_model = intensity_model[sorting_indices]
    sort_ints = np.argsort(intensity_model)
    if not (sort_ints==np.arange(len(intensity_model),0,-1)-1).all():
        print('WARNING: the model intensities are not radially decreasing. In input_limb_file', input_limb_path+input_limb_file, '.')
    return mu_model, radi_model, intensity_model


def get_max_gradient(muorradi_model,intensity_model):
#This function
    dint_dx = np.abs( (intensity_model[1:]-intensity_model[:-1])/(muorradi_model[1:]-muorradi_model[:-1]) )
    muorradi_maxgrad = 0.5*( muorradi_model[np.argmax(dint_dx)+1] + muorradi_model[np.argmax(dint_dx)] )
    return muorradi_maxgrad


def limbmodel_cut_and_rescale(radi_model,r0cut,r0res,intensity_model):
#This function
    radi_model_cut = radi_model[np.where(radi_model<=r0cut)[0]]
    radi_model_cr = radi_model_cut/r0res
    intensity_model_cr = intensity_model[np.where(radi_model<=r0cut)[0]]
    if np.max(radi_model_cr)>1.0:
        print('WARNING: some radii are larger than 1 after rescaling; these radii are rejected.')
        radi_model_cr = radi_model_cr[np.where(radi_model_cr<=1.0)[0]]
        intensity_model_cr = intensity_model_cr[np.where(radi_model_cr<=1.0)[0]]
    mu_model_cr = np.sqrt(1.0-radi_model_cr**2.0)
    return mu_model_cr, radi_model_cr, intensity_model_cr


def get_limb_grid(mu_model,radi_model,intensity_model,n_annuli,r0res,interpolation_type,interpolation_variable):
#This function
    radi_grid = (0.5+np.arange(n_annuli))/n_annuli
    radi_grid_back = radi_grid*r0res
    if interpolation_variable=='mu':
        mu_grid_back = np.sqrt(1.0-radi_grid_back**2.0)
        get_interp_ints = interp1d(mu_model, intensity_model, kind=interpolation_type, fill_value='extrapolate')
        intensity_grid = get_interp_ints(mu_grid_back)
    elif interpolation_variable=='radi':
        get_interp_ints = interp1d(radi_model, intensity_model, kind=interpolation_type, fill_value='extrapolate')
        intensity_grid = get_interp_ints(radi_grid_back)
    return radi_grid, intensity_grid


def get_star_flux(radi_grid, intensity_grid, n_annuli):
#This function
    dr_grid = 1.0/n_annuli
    flux_grid = 2.0*radi_grid*dr_grid*intensity_grid
    flux_sum = np.sum(flux_grid)
    return flux_grid, flux_sum


def get_occ_star_flux(radi_grid, flux_grid, rp_over_rs, z_sep):
#This function
    F = np.zeros_like(z_sep)
    radi_grid, z_sep = np.meshgrid(radi_grid, z_sep)
    flux_grid, F = np.meshgrid(flux_grid, F)

    cosa = (radi_grid**2.0+z_sep**2.0-rp_over_rs**2.0)/(2.0*radi_grid*z_sep)
    flux = np.where( np.abs( z_sep - radi_grid ) > rp_over_rs, F, np.where( z_sep + radi_grid < rp_over_rs, flux_grid, (flux_grid/np.pi)*np.arccos(cosa) ))
    return np.sum(flux,1)



def process_configuration(input_dict):
#This function
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())

    input_limb_type = input_dict_local['input_limb_type'][0]
    input_limb_path = input_dict_local['input_limb_path'][0]
    input_limb_file = input_dict_local['input_limb_file'][0]
    [mu_model_orig, radi_model_orig, intensity_model_orig] = get_input_limb_model(input_limb_type,input_limb_path,input_limb_file)

    cutting_limb = input_dict_local['cutting_limb'][0]
    if cutting_limb=='radi_gradient':
        r0cut = get_max_gradient(radi_model_orig, intensity_model_orig)
        mu0cut = np.sqrt(1.0-r0cut**2.0)
    elif cutting_limb=='mu_gradient':
        mu0cut = get_max_gradient(mu_model_orig, intensity_model_orig)
        r0cut = np.sqrt(1.0-mu0cut**2.0)
    elif cutting_limb=='user_cut':
        mu0cut = input_dict_local['user_cut_mu'][0]
        r0cut = input_dict_local['user_cut_radi'][0]
    elif cutting_limb=='no_cut':
        mu0cut = 0.0
        r0cut = 1.0

    rescaling_limb = input_dict_local['rescaling_limb'][0]
    if rescaling_limb=='as_cut':
        mu0res = mu0cut
        r0res = r0cut
    elif rescaling_limb=='no_rescale':
        mu0res = 0.0
        r0res = 1.0
    elif rescaling_limb=='user_rescale':
        mu0res = input_dict_local['user_rescale_mu'][0]
        r0res = input_dict_local['user_rescale_radi'][0]

    [mu_model_cr, radi_model_cr, intensity_model_cr] = limbmodel_cut_and_rescale(radi_model_orig,r0cut,r0res,intensity_model_orig)

    n_annuli = input_dict_local['n_annuli'][0]
    n_annuli = int(n_annuli)
    interpolation_type = input_dict_local['interpolation_type'][0]
    interpolation_variable = input_dict_local['interpolation_variable'][0]
    [radi_grid, intensity_grid] = get_limb_grid(mu_model_orig,radi_model_orig,intensity_model_orig,n_annuli,r0res,interpolation_type,interpolation_variable)
    [flux_grid, flux_sum] = get_star_flux(radi_grid, intensity_grid, n_annuli)

    input_series_type = input_dict_local['input_series_type'][0]
    input_series_path = input_dict_local['input_series_path'][0]
    input_series_file = input_dict_local['input_series_file'][0]
    rescaling_input_params = input_dict_local['rescaling_input_params'][0]
    if rescaling_input_params=='no':
        res_par = 1.0
    elif rescaling_input_params=='yes':
        res_par = 1.0/r0res

    if input_series_type=='z_sep':
        z_sep = get_x_series(input_series_type, input_series_path, input_series_file)
        z_sep *= res_par
    elif input_series_type=='phi':
        sma_over_rs = input_dict_local['sma_over_rs'][0]
        sma_over_rs *= res_par
        inclination = input_dict_local['inclination'][0]
        eccentricity = input_dict_local['eccentricity'][0]
        arg_pericenter = input_dict_local['arg_pericenter'][0]
        [phi, z_sep] = get_x_series(input_series_type, input_series_path, input_series_file, sma_over_rs=sma_over_rs, inclination=inclination, eccentricity=eccentricity, arg_pericenter=arg_pericenter)
    elif input_series_type=='time':
        sma_over_rs = input_dict_local['sma_over_rs'][0]
        sma_over_rs *= res_par
        inclination = input_dict_local['inclination'][0]
        eccentricity = input_dict_local['eccentricity'][0]
        arg_pericenter = input_dict_local['arg_pericenter'][0]
        period_orbital = input_dict_local['period_orbital'][0]
        epoch_of_transit = input_dict_local['epoch_of_transit'][0]
        time_conversion_factor = input_dict_local['time_conversion_factor'][0]
        [time, phi, z_sep] = get_x_series(input_series_type, input_series_path, input_series_file, sma_over_rs=sma_over_rs, inclination=inclination, eccentricity=eccentricity, arg_pericenter=arg_pericenter, period_orbital=period_orbital, epoch_of_transit=epoch_of_transit, time_conversion_factor=time_conversion_factor)

    rp_over_rs = input_dict_local['rp_over_rs'][0]
    rp_over_rs *= res_par
    flux_occ = get_occ_star_flux(radi_grid, flux_grid, rp_over_rs, z_sep)
    flux_norm = 1.0 - flux_occ/flux_sum

    output_path = input_dict_local['output_path'][0]
    output_filename = input_dict_local['output_filename'][0]
    output_fileext = input_dict_local['output_fileext']

    if not output_filename:
        output_filename = 'trip_'+os.path.splitext(input_limb_file)[0]#+'_'+os.path.splitext(input_series_file)[0]

    output_dict = {}
    output_dict['configuration_file'] = copy.deepcopy(input_dict_local)
    output_dict['limb_model_original'] = {}
    output_dict['limb_model_original']['mu'] = mu_model_orig
    output_dict['limb_model_original']['radi'] = radi_model_orig
    output_dict['limb_model_original']['intensity'] = intensity_model_orig
    output_dict['limb_model_final'] = {}
    output_dict['limb_model_final']['mu'] = mu_model_cr
    output_dict['limb_model_final']['radi'] = radi_model_cr
    output_dict['limb_model_final']['intensity'] = intensity_model_cr
    output_dict['cut_and_rescale'] = {}
    output_dict['cut_and_rescale']['mu0cut'] = mu0cut
    output_dict['cut_and_rescale']['r0cut'] = r0cut
    output_dict['cut_and_rescale']['mu0res'] = mu0res
    output_dict['cut_and_rescale']['r0res'] = r0res
    output_dict['rescaled_params'] = {}
    output_dict['rescaled_params']['rp_over_rs'] = rp_over_rs
    try:
        output_dict['rescaled_params']['sma_over_rs'] = sma_over_rs
        output_dict['rescaled_params']['inclination'] = inclination
        output_dict['rescaled_params']['eccentricity'] = eccentricity
        output_dict['rescaled_params']['arg_pericenter'] = arg_pericenter
        output_dict['rescaled_params']['period_orbital'] = period_orbital
        output_dict['rescaled_params']['epoch_of_transit'] = epoch_of_transit
    except:
        None
    output_dict['time_series'] = {}
    output_dict['time_series']['z_sep'] = z_sep
    output_dict['time_series']['flux_norm'] = flux_norm
    try:
        output_dict['time_series']['time'] = time
        output_dict['time_series']['phi'] = phi
    except:
        None

    if '.pickle' in output_fileext:
        with open(os.path.join(output_path, output_filename+'.pickle'), 'wb') as outfile:
            pickle.dump(output_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    if '.txt' in output_fileext:
        header_txt = ['z_sep', 'flux_norm']
        output_txt = np.column_stack((z_sep,flux_norm))
        try:
            output_txt = np.column_stack((phi,output_txt))
            header_txt = ['phi']+header_txt
        except:
            None
        try:
            output_txt = np.column_stack((time,output_txt))
            header_txt = ['time']+header_txt
        except:
            None
        np.savetxt(os.path.join(output_path, output_filename+'.txt'), output_txt, fmt='%12.8f', header='    '.join(header_txt))

    return output_dict







def trip_calculate(configuration_file):
#Main function for 
    input_dict1 = read_configuration(configuration_file) #Reading configuration file
    [check, input_dict] = check_configuration(input_dict1) #Checking configuration file and dictionary update
    if check:
        process_configuration(input_dict) #Computing and saving
    else:
        print('Something went wrong with the input.')





















