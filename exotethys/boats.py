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
from astropy.modeling.blackbody import blackbody_lambda
import astropy.units as u
from scipy.interpolate import interp1d
from scipy.integrate import simps

import copy
import pickle

from ._database import Database, databases, sys, urlretrieve, glob, time, shutil


def str2float(s): #function in common with SAIL
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

def check_length(vector, min_length=1, max_length=None): #function in common with SAIL
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

def check_type_in_list(vector, item_type): #function in common with SAIL
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

def check_unit_type(item, ref_unit, others=[]):
    """
    This function checks that a given string corresponds to a physical unit of the expected kind and returns a boolean value.
    
    :param str item: 
    :param quantity ref_unit: reference unit of the expected physical dimension
    :argument listOfObjects others: list of additional valid strings for item, default is empty
    :return: True if item has the same physical type of ref_unit (or item is included in others), False otherwise
    :rtype: bool
    """

    check = True
    if item in others:
        return check
    try:
        u.Unit(item).to(ref_unit)
    except ValueError or u.core.UnitConversionError:
        check = False
    return check


def check_values(vector, min_value=None, max_value=None, inf_value=None, sup_value=None):
    """
    This function checks that all the elements of a list/array are within given limits (if any). It returns a boolean value.
    
    :param list of float vector:
    :argument float or quantity min_value: minimum value for any element of the list (default is None)
    :argument float or quantity max_value: maximum value for any element of the list (default is None)
    :argument float or quantity inf_value: lower limit (not included) for any element of the list (default is None)
    :argument float or quantity sup_value: upper limit (not included) for any element of the list (default is None)
    :return: True if all the vector values are within the limits, False otherwise.
    :rtype: bool
    """
    check = True
    for item in vector:
        cond0 = not np.isfinite(item)
        cond1 = (min_value is not None) and item<min_value
        cond2 = (inf_value is not None) and item<=inf_value
        cond3 = (max_value is not None) and item>max_value
        cond4 = (sup_value is not None) and item>=sup_value
        if cond0 or cond1 or cond2 or cond3 or cond4:
            check = False
    return check



def read_as_numpy_array(path_file): #function in common with SAIL
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

def check_2Darray(arr, n_col=None): #function in common with SAIL
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


def read_configuration(filename): #function in common with SAIL
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
                #print(content)
                key = content[0]
                value = []
                for item in content[1:]:
                    if item[0] != '!':
                        value += [str2float(item),]
                input_dict[key] = value
    return input_dict


def get_transit_duration_T14(rp_over_rs, sma_over_rs, inclination, period):
    """
    This function computes the transit duration between the external contact points.
    
    :param quantity rp_over_rs: ratio of planet and star radii (dimensionless)
    :param quantity sma_over_rs: ratio of orbital semimajor axis and star radius (dimensionless)
    :param quantity inclination: orbital inclination angle
    :param quantity period: orbital period
    :return: the total transit duration, i.e., between the external contact points
    :rtype: quantity
    """
    b_impact = sma_over_rs*np.cos(inclination)
    T14 = np.arcsin( np.sqrt( (1 + rp_over_rs)**2 - b_impact**2 ) / (sma_over_rs*np.sin(inclination)) )
    T14 *= period/(np.pi * u.rad)
    return T14


def get_planet_temperatures(star_effective_temperature, sma_over_rs, albedo, efficiency):
    """
    This function computes the exoplanet day and nightside temperatures, based on Cowan & Agol 2011, ApJ, 729, 54, Equations 4 and 5.
    
    :param quantity star_effective_temperature: the effective temperature of the star (in Kelvin)
    :param quantity sma_over_rs: ratio of orbital semimajor axis and star radius (dimensionless)
    :param quantity albedo: the bond albedo of the exoplanet atmosphere (dimensionless, 0<=albedo<=1)
    :param quantity efficiency: the circulation efficiency of the exoplanet atmosphere (dimensionless, 0<=efficiency<=1)
    :return: the exoplanet day and nightside temperatures and an associated string
    :rtype: quantity
    """
    irradiation_temperature = star_effective_temperature/np.sqrt(sma_over_rs)
    Tday = irradiation_temperature * (1-albedo)**0.25 * ((2.0/3.0) - (5.0/12.0)*efficiency)**0.25
    Tnight = irradiation_temperature * (1-albedo)**0.25 * (efficiency/4.0)**0.25
    label = 'albedo' + str(albedo) + '_efficiency' + str(efficiency)
    return Tday, Tnight, label


def get_planet_albedo_and_efficiency(planet_day_temperature, planet_night_temperature, star_effective_temperature, orbital_semimajor_axis, star_radius):
    """
    This function computes the exoplanet albedo and circulation efficiency from their day and nightside temperatures, based on Cowan & Agol 2011, ApJ, 729, 54, Equations 4 and 5.
    
    :param quantity planet_day_temperature: (in Kelvin)
    :param quantity planet_night_temperature: (in Kelvin)
    :param quantity star_effective_temperature: the effective temperature of the star (in Kelvin)
    :param quantity orbital_semimajor_axis:
    :param quantity star_radius:
    :return: the exoplanet albedo and circulation efficiency
    :rtype: quantity
    """
    sma_over_rs = (orbital_semimajor_axis/star_radius).decompose()
    irradiation_temperature = star_effective_temperature/np.sqrt(sma_over_rs)
    efficiency = (8*planet_night_temperature**4.0) / ( (3*planet_day_temperature**4.0) + (5*planet_night_temperature**4.0) )
    albedo = 1 - ( (3*planet_day_temperature**4.0) + (5*planet_night_temperature**4.0) )/(2*irradiation_temperature**4.0)
    return albedo, efficiency


def check_configuration(input_dict):
    """
    This function checks and modifies the dictionary obtained from the configuration file.
    It returns a boolean value and the updated dictionary.
    
    :param dict input_dict: 
    :return: a bool value and the updated configuration dictionary
    :rtype: bool, dict
    """
    check = True
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    mandatory_keys = ['stellar_models_grid', 'planet_models_grid', 'passbands', 'observing_duration', 'observing_duration_unit', 'telescope_area', 'telescope_area_unit', 'star_radius', 'star_radius_unit', 'orbital_semimajor_axis', 'orbital_semimajor_axis_unit', 'orbital_inclination', 'orbital_inclination_unit', 'orbital_period', 'orbital_period_unit', 'planet_radius', 'planet_radius_unit']
    allowed_keys = mandatory_keys + ['output_path', 'output_filename', 'output_fileext', 'passbands_path', 'wavelength_bins_path', 'wavelength_bins_files', 'star_effective_temperature', 'star_log_gravity', 'star_metallicity', 'planet_bond_albedo', 'planet_circulation_efficiency', 'planet_day_temperature', 'planet_night_temperature', 'system_distance', 'system_distance_unit', 'star_model_path', 'star_model_file', 'rescale_star_flux', 'planet_day_model_path', 'planet_day_model_file', 'planet_night_model_path', 'planet_night_model_file', 'rescale_planet_flux']

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

   #Checking the stellar_models_grid (MANDATORY, 1 ONLY).
    if 'stellar_models_grid' in input_keys:
        stellar_models_grid = input_dict_local['stellar_models_grid']
        if not check_length(stellar_models_grid, max_length=1):
            print('ERROR: invalid length=', len(stellar_models_grid), 'for stellar_models_grid. It must have length=1.')
            check = False
        else:
            allowed_stellar_models_grid = ['Phoenix_2018', 'Phoenix_2012_13', 'Phoenix_drift_2012', 'Atlas_2000', 'Stagger_2015', 'Blackbody', 'Userfile']
            stellar_models_grid = stellar_models_grid[0]
            if stellar_models_grid not in allowed_stellar_models_grid:
                print('ERROR:',stellar_models_grid,'is not a valid stellar_models_grid. The allowed names are Phoenix_2018, Phoenix_2012_13, Phoenix_drift_2012, Atlas_2000, Stagger_2015, Blackbody and Userfile.')
                check = False
        if check:
            #stellar_models_grid = stellar_models_grid[0]
            if stellar_models_grid=='Userfile':
                mandatory_keys_userstar = ['star_model_file', 'rescale_star_flux']
                redundant_keys_userstar = ['star_log_gravity', 'star_metallicity']
                for key in mandatory_keys_userstar:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'is not specified for', stellar_models_grid, 'stellar_models_grid.')
                        check = False
                for key in redundant_keys_userstar:
                    if key in input_keys:
                        print('ERROR:', key, 'is not a valid keyword for', stellar_models_grid, 'stellar_models_grid.')
                        check = False
            elif stellar_models_grid=='Blackbody':
                mandatory_keys_bbstar = ['star_effective_temperature', 'system_distance', 'system_distance_unit']
                redundant_keys_bbstar = ['star_model_path', 'star_model_file', 'rescale_star_flux', 'star_log_gravity', 'star_metallicity']
                for key in mandatory_keys_bbstar:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'is not specified for', stellar_models_grid, 'stellar_models_grid.')
                        check = False
                for key in redundant_keys_bbstar:
                    if key in input_keys:
                        print('ERROR:', key, 'is not a valid keyword for', stellar_models_grid, 'stellar_models_grid.')
                        check = False
            else: #All the others stellar_models_grid
                mandatory_keys_star = ['star_effective_temperature', 'system_distance', 'system_distance_unit']
                redundant_keys_star = ['star_model_path', 'star_model_file', 'rescale_star_flux']
                for key in mandatory_keys_star:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'is not specified for', stellar_models_grid, 'stellar_models_grid.')
                        check = False
                for key in redundant_keys_star:
                    if key in input_keys:
                        print('ERROR:', key, 'is not a valid keyword for', stellar_models_grid, 'stellar_models_grid.')
                        check = False

    #Checking the planet_models_grid (MANDATORY, 1 ONLY).
    if 'planet_models_grid' in input_keys:
        planet_models_grid = input_dict_local['planet_models_grid']
        if not check_length(planet_models_grid, max_length=1):
            print('ERROR: invalid length=', len(planet_models_grid), 'for planet_models_grid. It must have length=1.')
            check = False
        else:
            allowed_planet_models_grid = ['Blackbody', 'Userfile']
            planet_models_grid = planet_models_grid[0]
            if planet_models_grid not in allowed_planet_models_grid:
                print('ERROR:',planet_models_grid,'is not a valid planet_models_grid. The allowed names are Blackbody and Userfile.')
                check = False
        if check:
            #planet_models_grid = planet_models_grid[0]
            if planet_models_grid=='Userfile':
                mandatory_keys_userplanet = ['planet_day_model_file', 'planet_night_model_file', 'rescale_planet_flux']
                redundant_keys_userplanet = ['planet_day_temperature', 'planet_night_temperature', 'planet_circulation_efficiency']
                for key in mandatory_keys_userplanet:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'is not specified for', planet_models_grid, 'planet_models_grid.')
                        check = False
                for key in redundant_keys_userplanet:
                    if key in input_keys:
                        print('ERROR:', key, 'is not a valid keyword for', planet_models_grid, 'planet_models_grid.')
                        check = False
            elif planet_models_grid=='Blackbody':
                mandatory_keys_bbplanet = ['system_distance', 'system_distance_unit']
                redundant_keys_bbplanet = ['planet_day_model_path', 'planet_day_model_file', 'planet_night_model_path', 'planet_night_model_file', 'rescale_planet_flux']
                for key in mandatory_keys_bbplanet:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'is not specified for', planet_models_grid, 'planet_models_grid.')
                        check = False
                for key in redundant_keys_bbplanet:
                    if key in input_keys:
                        print('ERROR:', key, 'is not a valid keyword for', planet_models_grid, 'planet_models_grid.')
                        check = False
                if 'planet_circulation_efficiency' in input_keys:
                    if 'planet_day_temperature' in input_keys:
                        print('ERROR: The planet_circulation_efficiency cannot be used if planet_day_temperature is given.')
                        check = False
                    if 'planet_night_temperature' in input_keys:
                        print('ERROR: The planet_circulation_efficiency cannot be used if planet_night_temperature is given.')
                        check = False
                    if 'star_effective_temperature' not in input_keys:
                        print('ERROR: The star_effective_temperature must be provided if planet_circulation_efficiency is given.') 
                        check = False
                else:
                    if 'planet_day_temperature' and 'planet_night_temperature' not in input_keys:
                        print('ERROR: Either the planet_circulation_efficiency or the planet_day_temperature and planet_night_temperature must be given for', planet_models_grid, 'planet_models_grid.')
                        check = False

    if check and 'rescale_star_flux' in input_keys:
        rescale_star_flux = input_dict_local['rescale_star_flux']
        if not check_length(rescale_star_flux, max_length=1):
            print('ERROR: invalid length=', len(rescale_star_flux), 'for rescale_star_flux. It must have length=1.')
            check = False
        else:
            allowed_rescale_star_flux = ['Yes', 'No']
            rescale_star_flux = rescale_star_flux[0]
            if rescale_star_flux not in allowed_rescale_star_flux:
                print('ERROR:', rescale_star_flux, 'is not a valid entry for rescale_star_flux. The allowed entries are Yes and No.')
                check = False
            elif rescale_star_flux=='Yes':
                mandatory_keys_rescale_star_flux_yes = ['system_distance', 'system_distance_unit']
                for key in mandatory_keys_rescale_star_flux_yes:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'must specified if rescale_star_flux is yes.')
                        check = False

    if check and 'rescale_planet_flux' in input_keys:
        rescale_planet_flux = input_dict_local['rescale_planet_flux']
        if not check_length(rescale_planet_flux, max_length=1):
            print('ERROR: invalid length=', len(rescale_planet_flux), 'for rescale_planet_flux. It must have length=1.')
            check = False
        else:
            allowed_rescale_planet_flux = ['Yes', 'No', 'Star']
            rescale_planet_flux = rescale_planet_flux[0]
            if rescale_planet_flux not in allowed_rescale_planet_flux:
                print('ERROR:', rescale_planet_flux, 'is not a valid entry for rescale_planet_flux. The allowed entries are Yes, No and Star.')
                check = False
            elif rescale_planet_flux=='Yes':
                mandatory_keys_rescale_planet_flux_yes = ['system_distance', 'system_distance_unit']
                for key in mandatory_keys_rescale_planet_flux_yes:
                    if key not in input_keys:
                        print('ERROR: mandatory keyword', key, 'must specified if rescale_planet_flux is yes.')
                        check = False


    #Checking the passbands list (MANDATORY)
    passbands = input_dict_local['passbands']
    if not check_length(passbands):
        print('ERROR: invalid length=', len(passbands), 'for passbands. It must have length>=1.')
        check = False

    #Checking the case of built-in passbands
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

    #Eliminating repeated passbands
    #for item in passbands:
    #    if passbands.count(item)>1:
    #        print('WARNING:', item, 'entered multiple times in passbands. Repetitions are ignored.')
    #    input_dict_local['passbands'] = list(OrderedDict.fromkeys(passbands))

    #Checking the wavelength_bins_files to split the selected passbands (OPTIONAL, by DEFAULT the passbands are not split).
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

    #Checking telescope_area (MANDATORY)
    if 'telescope_area' in input_keys:
        telescope_area = input_dict_local['telescope_area']
        if not check_length(telescope_area, max_length=1):
            print('ERROR: invalid length=', len(telescope_area), 'for telescope_area. It must have length=1.')
            check = False
        if not check_type_in_list(telescope_area, float):
            print('ERROR: invalid type for item in telescope_area. Items must be float.')
            check = False

    #Checking telescope_area_unit (MANDATORY, 1 ONLY)
    if 'telescope_area_unit' in input_keys:
        u_check = True
        telescope_area_unit = input_dict_local['telescope_area_unit']
        if not check_length(telescope_area_unit, max_length=1):
            print('ERROR: invalid length=', len(telescope_area_unit), 'for telescope_area_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(telescope_area_unit, str):
            print('ERROR: invalid type for item in telescope_area_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(telescope_area_unit[0], u.m**2):
                print('ERROR: unrecognised input for telescope_area_unit. It must be a string associated with a unit of area.')
                check = False

    #Checking observing_duration (MANDATORY)
    if 'observing_duration' in input_keys:
        observing_duration = input_dict_local['observing_duration']
        if not check_length(observing_duration):
            print('ERROR: invalid length=', len(observing_duration), 'for observing_duration. It must have length>=1.')
            check = False
        if not check_type_in_list(observing_duration, float):
            print('ERROR: invalid type for item in observing_duration. Items must be float.')
            check = False

    #Checking observing_duration_unit (MANDATORY, 1 ONLY)
    if 'observing_duration_unit' in input_keys:
        u_check = True
        observing_duration_unit = input_dict_local['observing_duration_unit']
        if not check_length(observing_duration_unit, max_length=1):
            print('ERROR: invalid length=', len(observing_duration_unit), 'for observing_duration_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(observing_duration_unit, str):
            print('ERROR: invalid type for item in observing_duration_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(observing_duration_unit[0], u.second, others=['T_14']):
                print('ERROR: unrecognised input for observing_duration_unit. It must be a string associated with a unit of time or T_14.')
                check = False

    #Checking star_effective_temperature (OPTIONAL, 1 ONLY)
    if 'star_effective_temperature' in input_keys:
        star_effective_temperature = input_dict_local['star_effective_temperature']
        if not check_length(star_effective_temperature, max_length=1):
            print('ERROR: invalid length=', len(star_effective_temperature), 'for star_effective_temperature. It must have length=1.')
            check = False
        if not check_type_in_list(star_effective_temperature, float):
            print('ERROR: invalid type for item in star_effective_temperature. Items must be float.')
            check = False

    #Checking star_log_gravity (OPTIONAL, 1 ONLY)
    if 'star_log_gravity' in input_keys:
        star_log_gravity = input_dict_local['star_log_gravity']
        if not check_length(star_log_gravity, max_length=1):
            print('ERROR: invalid length=', len(star_log_gravity), 'for star_log_gravity (optional). It must have length=1.')
            check = False
        if not check_type_in_list(star_log_gravity, float):
            print('ERROR: invalid type for item in star_log_gravity (optional). Items must be float.')
            check = False

    #Checking star_metallicity (OPTIONAL, 1 ONLY)
    if 'star_metallicity' in input_keys:
        star_metallicity = input_dict_local['star_metallicity']
        if not check_length(star_metallicity, max_length=1):
            print('ERROR: invalid length=', len(star_metallicity), 'for star_metallicity (optional). It must have length=1.')
            check = False
        if not check_type_in_list(star_metallicity, float):
            print('ERROR: invalid type for item in star_metallicity (optional). Items must be float.')
            check = False

    #Checking star_radius (MANDATORY, 1 ONLY)
    if 'star_radius' in input_keys:
        star_radius = input_dict_local['star_radius']
        if not check_length(star_radius, max_length=1):
            print('ERROR: invalid length=', len(star_radius), 'for star_radius. It must have length=1.')
            check = False
        if not check_type_in_list(star_radius, float):
            print('ERROR: invalid type for item in star_radius. Items must be float.')
            check = False

    #Checking star_radius_unit (MANDATORY, 1 ONLY)
    if 'star_radius_unit' in input_keys:
        u_check = True
        star_radius_unit = input_dict_local['star_radius_unit']
        if not check_length(star_radius_unit, max_length=1):
            print('ERROR: invalid length=', len(star_radius_unit), 'for star_radius_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(star_radius_unit, str):
            print('ERROR: invalid type for item in star_radius_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(star_radius_unit[0], u.m):
                print('ERROR: unrecognised input for star_radius_unit. It must be a string associated with a unit of length.')
                check = False

    #Checking orbital_semimajor_axis (MANDATORY, 1 ONLY)
    if 'orbital_semimajor_axis' in input_keys:
        orbital_semimajor_axis = input_dict_local['orbital_semimajor_axis']
        if not check_length(orbital_semimajor_axis, max_length=1):
            print('ERROR: invalid length=', len(orbital_semimajor_axis), 'for orbital_semimajor_axis. It must have length=1.')
            check = False
        if not check_type_in_list(orbital_semimajor_axis, float):
            print('ERROR: invalid type for item in orbital_semimajor_axis. Items must be float.')
            check = False

    #Checking orbital_semimajor_axis_unit (MANDATORY, 1 ONLY)
    if 'orbital_semimajor_axis_unit' in input_keys:
        u_check = True
        orbital_semimajor_axis_unit = input_dict_local['orbital_semimajor_axis_unit']
        if not check_length(orbital_semimajor_axis_unit, max_length=1):
            print('ERROR: invalid length=', len(orbital_semimajor_axis_unit), 'for orbital_semimajor_axis_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(orbital_semimajor_axis_unit, str):
            print('ERROR: invalid type for item in orbital_semimajor_axis_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(orbital_semimajor_axis_unit[0], u.m, others=['star_radius']):
                print('ERROR: unrecognised input for orbital_semimajor_axis_unit. It must be a string associated with a unit of length or star_radius.')
                check = False

    #Checking orbital_inclination (MANDATORY, 1 ONLY)
    if 'orbital_inclination' in input_keys:
        orbital_inclination = input_dict_local['orbital_inclination']
        if not check_length(orbital_inclination, max_length=1):
            print('ERROR: invalid length=', len(orbital_inclination), 'for orbital_inclination. It must have length=1.')
            check = False
        if not check_type_in_list(orbital_inclination, float):
            print('ERROR: invalid type for item in orbital_inclination. Items must be float.')
            check = False

    #Checking orbital_inclination_unit (MANDATORY, 1 ONLY)
    if 'orbital_inclination_unit' in input_keys:
        u_check = True
        orbital_inclination_unit = input_dict_local['orbital_inclination_unit']
        if not check_length(orbital_inclination_unit, max_length=1):
            print('ERROR: invalid length=', len(orbital_inclination_unit), 'for orbital_inclination_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(orbital_inclination_unit, str):
            print('ERROR: invalid type for item in orbital_inclination_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(orbital_inclination_unit[0], u.deg):
                print('ERROR: unrecognised input for orbital_inclination_unit. It must be a string associated with a unit of angle.')
                check = False

    #Checking orbital_period (MANDATORY, 1 ONLY)
    if 'orbital_period' in input_keys:
        orbital_period = input_dict_local['orbital_period']
        if not check_length(orbital_period, max_length=1):
            print('ERROR: invalid length=', len(orbital_period), 'for orbital_period. It must have length=1.')
            check = False
        if not check_type_in_list(orbital_period, float):
            print('ERROR: invalid type for item in orbital_period. Items must be float.')
            check = False

    #Checking orbital_period_unit (MANDATORY, 1 ONLY)
    if 'orbital_period_unit' in input_keys:
        u_check = True
        orbital_period_unit = input_dict_local['orbital_period_unit']
        if not check_length(orbital_period_unit, max_length=1):
            print('ERROR: invalid length=', len(orbital_period_unit), 'for orbital_period_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(orbital_period_unit, str):
            print('ERROR: invalid type for item in orbital_period_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(orbital_period_unit[0], u.second):
                print('ERROR: unrecognised input for orbital_period_unit. It must be a string associated with a unit of time.')
                check = False

    #Checking planet_bond_albedo (OPTIONAL)
    if 'planet_bond_albedo' in input_keys:
        planet_bond_albedo = input_dict_local['planet_bond_albedo']
        if not check_length(planet_bond_albedo):
            print('ERROR: invalid length=', len(planet_bond_albedo), 'for planet_bond_albedo (optional). It must have length>=1.')
            check = False
        if not check_type_in_list(planet_bond_albedo, float):
            print('ERROR: invalid type for item in planet_bond_albedo (optional). Items must be float.')
            check = False

    #Checking planet_circulation_efficiency (OPTIONAL)
    if 'planet_circulation_efficiency' in input_keys:
        planet_circulation_efficiency = input_dict_local['planet_circulation_efficiency']
        if not check_length(planet_circulation_efficiency):
            print('ERROR: invalid length=', len(planet_circulation_efficiency), 'for planet_circulation_efficiency (optional). It must have length>=1.')
            check = False
        if not check_type_in_list(planet_circulation_efficiency, float):
            print('ERROR: invalid type for item in planet_circulation_efficiency (optional). Items must be float.')
            check = False

    #Checking planet_day_temperature (OPTIONAL)
    if 'planet_day_temperature' in input_keys:
        planet_day_temperature = input_dict_local['planet_day_temperature']
        if not check_length(planet_day_temperature):
            print('ERROR: invalid length=', len(planet_day_temperature), 'for planet_day_temperature (optional). It must have length>=1.')
            check = False
        if not check_type_in_list(planet_day_temperature, float):
            print('ERROR: invalid type for item in planet_day_temperature (optional). Items must be float.')
            check = False

    #Checking planet_night_temperature (OPTIONAL)
    if 'planet_night_temperature' in input_keys:
        planet_night_temperature = input_dict_local['planet_night_temperature']
        if not check_length(planet_day_temperature, min_length=len(planet_day_temperature), max_length=len(planet_day_temperature)):
            print('ERROR: invalid length=', len(planet_night_temperature), 'for planet_night_temperature (optional). It must have the same length of planet_day_temperature.')
            check = False
        if not check_type_in_list(planet_night_temperature, float):
            print('ERROR: invalid type for item in planet_night_temperature (optional). Items must be float.')
            check = False

    #Checking planet_radius (MANDATORY, 1 ONLY)
    if 'planet_radius' in input_keys:
        planet_radius = input_dict_local['planet_radius']
        if not check_length(planet_radius, max_length=1):
            print('ERROR: invalid length=', len(planet_radius), 'for planet_radius. It must have length=1.')
            check = False
        if not check_type_in_list(planet_radius, float):
            print('ERROR: invalid type for item in planet_radius. Items must be float.')
            check = False

    #Checking planet_radius_unit (MANDATORY, 1 ONLY)
    if 'planet_radius_unit' in input_keys:
        u_check = True
        planet_radius_unit = input_dict_local['planet_radius_unit']
        if not check_length(planet_radius_unit, max_length=1):
            print('ERROR: invalid length=', len(planet_radius_unit), 'for planet_radius_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(planet_radius_unit, str):
            print('ERROR: invalid type for item in planet_radius_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(planet_radius_unit[0], u.m, others=['star_radius']):
                print('ERROR: unrecognised input for planet_radius_unit. It must be a string associated with a unit of length or star_radius.')
                check = False

    #Checking system_distance (MANDATORY, 1 ONLY)
    if 'system_distance' in input_keys:
        system_distance = input_dict_local['system_distance']
        if not check_length(system_distance, max_length=1):
            print('ERROR: invalid length=', len(system_distance), 'for system_distance. It must have length=1.')
            check = False
        if not check_type_in_list(system_distance, float):
            print('ERROR: invalid type for item in system_distance. Items must be float.')
            check = False

    #Checking system_distance_unit (MANDATORY, 1 ONLY)
    if 'system_distance_unit' in input_keys:
        u_check = True
        system_distance_unit = input_dict_local['system_distance_unit']
        if not check_length(system_distance_unit, max_length=1):
            print('ERROR: invalid length=', len(system_distance_unit), 'for system_distance_unit. It must have length=1.')
            check = False
            u_check = False
        if not check_type_in_list(system_distance_unit, str):
            print('ERROR: invalid type for item in system_distance_unit. Items must be string.')
            check = False
            u_check = False
        if u_check:
            if not check_unit_type(system_distance_unit[0], u.m):
                print('ERROR: unrecognised input for system_distance_unit. It must be a string associated with a unit of length.')
                check = False

    #Checking the output_path (OPTIONAL, 1 ONLY).
    if 'output_path' in input_keys:
        output_path = input_dict_local['output_path']
        if not check_length(output_path, max_length=1):
            print('ERROR: invalid length=', len(output_path), 'for output_path. It must have length=1.')
            check = False
        if not check_type_in_list(output_path,str):
            print('ERROR: invalid type for output_path. It must be string.')
            check = False
        if check and not os.path.exists(output_path[0]):
            print('ERROR: the chosen output_path does not exists.')
            check = False
    else:
        input_dict_local['output_path'] = ['']

    #Checking the output_filename (OPTIONAL, 1 ONLY).
    if 'output_filename' in input_keys:
        output_filename = input_dict_local['output_filename']
        if not check_length(output_filename, max_length=1):
            print('ERROR: invalid length=', len(output_filename), 'for output_filename. It must have length=1.')
            check = False
        if not check_type_in_list(output_filename,str):
            print('ERROR: invalid type for output_filename. It must be string.')
            check = False
    else:
        input_dict_local['output_filename'] = ['output_boats']


    #Checking the output_fileext (OPTIONAL).
    if 'output_fileext' in input_keys:
        output_fileext = input_dict_local['output_fileext']
        if not check_length(output_fileext, max_length=2):
            print('ERROR: invalid length=', len(output_fileext), 'for output_fileext. It must have length<=2.')
            check = False
        else:
            allowed_output_fileext = ['.txt', '.pickle']
            for item in output_fileext:
                if item not in allowed_output_fileext:
                    print('ERROR: invalid output_fileext. It can be .pickle or .txt. Default is .pickle.')
                    check = False
    else:
        input_dict_local['output_fileext'] = ['.pickle']


    #Checking the star_model_path (OPTIONAL, 1 ONLY).
    if 'star_model_path' in input_keys:
        star_model_path = input_dict_local['star_model_path']
        if not check_length(star_model_path, max_length=1):
            print('ERROR: invalid length=', len(star_model_path), 'for star_model_path. It must have length=1.')
            check = False
        if not check_type_in_list(star_model_path,str):
            print('ERROR: invalid type for star_model_path. It must be string.')
            check = False
    else:
        input_dict_local['star_model_path'] = ['']


    #Checking the star_model_path (OPTIONAL, 1 ONLY).
    if 'star_model_path' in input_keys:
        star_model_path = input_dict_local['star_model_path']
        if not check_length(star_model_path, max_length=1):
            print('ERROR: invalid length=', len(star_model_path), 'for star_model_path. It must have length=1.')
            check = False
        if not check_type_in_list(star_model_path,str):
            print('ERROR: invalid type for star_model_path. It must be string.')
            check = False
    else:
        input_dict_local['star_model_path'] = ['']

    #Checking the star_model_file (OPTIONAL, 1 ONLY).
    if 'star_model_file' in input_keys:
        star_model_file = input_dict_local['star_model_file']
        if not check_length(star_model_file, max_length=1):
            print('ERROR: invalid length=', len(star_model_file), 'for star_model_file. It must have length=1.')
            check = False
        if not check_type_in_list(star_model_file,str):
            print('ERROR: invalid type for star_model_file. It must be string.')
            check = False

    #Checking the planet_day_model_path (OPTIONAL, 1 ONLY).
    if 'planet_day_model_path' in input_keys:
        planet_day_model_path = input_dict_local['planet_day_model_path']
        if not check_length(planet_day_model_path, max_length=1):
            print('ERROR: invalid length=', len(planet_day_model_path), 'for planet_day_model_path. It must have length=1.')
            check = False
        if not check_type_in_list(planet_day_model_path,str):
            print('ERROR: invalid type for planet_day_model_path. It must be string.')
            check = False
    else:
        input_dict_local['planet_day_model_path'] = ['']

    #Checking the planet_day_model_file (MANDATORY, 1 ONLY).
    if 'planet_day_model_file' in input_keys:
        planet_day_model_file = input_dict_local['planet_day_model_file']
        if not check_length(planet_day_model_file, max_length=1):
            print('ERROR: invalid length=', len(planet_day_model_file), 'for planet_day_model_file. It must have length=1.')
            check = False
        if not check_type_in_list(planet_day_model_file,str):
            print('ERROR: invalid type for planet_day_model_file. It must be string.')
            check = False

    #Checking the planet_night_model_path (OPTIONAL, 1 ONLY).
    if 'planet_night_model_path' in input_keys:
        planet_night_model_path = input_dict_local['planet_night_model_path']
        if not check_length(planet_night_model_path, max_length=1):
            print('ERROR: invalid length=', len(planet_night_model_path), 'for planet_night_model_path. It must have length=1.')
            check = False
        if not check_type_in_list(planet_night_model_path,str):
            print('ERROR: invalid type for planet_night_model_path. It must be string.')
            check = False
    else:
        input_dict_local['planet_night_model_path'] = ['']

    #Checking the planet_night_model_file (MANDATORY, 1 ONLY).
    if 'planet_night_model_file' in input_keys:
        planet_night_model_file = input_dict_local['planet_night_model_file']
        if not check_length(planet_night_model_file, max_length=1):
            print('ERROR: invalid length=', len(planet_night_model_file), 'for planet_night_model_file. It must have length=1.')
            check = False
        if not check_type_in_list(planet_night_model_file,str):
            print('ERROR: invalid type for planet_night_model_file. It must be string.')
            check = False

    return check, input_dict_local




def create_new_dict(input_dict):
    """
    This function checks and modifies the input dictionary returned by the function check_configuration.
    It returns a boolean value and the updated dictionary.
    
    :param dict input_dict: 
    :return: a bool value and the updated configuration dictionary
    :rtype: bool, dict
    """

    check = True
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())

    stellar_models_grid = input_dict_local['stellar_models_grid'][0]
    if not stellar_models_grid in ['Blackbody', 'Userfile']:
        if not 'star_log_gravity' in input_keys:
            input_dict_local['star_log_gravity'] = [4.5]
            print('WARNING: Assuming star_log_gravity = 4.5.')
        if not 'star_metallicity' in input_keys:
            input_dict_local['star_metallicity'] = [0.0]
            print('WARNING: Assuming star_metallicity = 0.0.')

    input_dict_local['telescope_area'] *= u.Unit( input_dict_local['telescope_area_unit'][0] )
    del input_dict_local['telescope_area_unit']
    if not check_values(input_dict_local['telescope_area'], inf_value=0*u.m**2):
        print('ERROR: The telescope_area must be positive.')
        check =  False

    input_dict_local['star_radius'] *= u.Unit( input_dict_local['star_radius_unit'][0] )
    del input_dict_local['star_radius_unit']
    if not check_values(input_dict_local['star_radius'], inf_value=0*u.m):
        print('ERROR: The star_radius must be positive.')
        check =  False

    try:
        input_dict_local['orbital_semimajor_axis'] *= u.Unit( input_dict_local['orbital_semimajor_axis_unit'][0] )
    except ValueError:
        input_dict_local['orbital_semimajor_axis'] *= input_dict_local['star_radius'][0]
    del input_dict_local['orbital_semimajor_axis_unit']
    if not check_values(input_dict_local['orbital_semimajor_axis'], inf_value=input_dict_local['star_radius'][0]):
        print('ERROR: The orbital_semimajor_axis must be longer than the star_radius.')
        check =  False

    input_dict_local['orbital_inclination'] *= u.Unit( input_dict_local['orbital_inclination_unit'][0] )
    del input_dict_local['orbital_inclination_unit']
    if not check_values(input_dict_local['orbital_inclination'], min_value=0*u.deg, max_value=90*u.deg):
        print('ERROR: The orbital_inclination must be in the range 0-90 deg.')
        check =  False

    input_dict_local['orbital_period'] *= u.Unit( input_dict_local['orbital_period_unit'][0] )
    del input_dict_local['orbital_period_unit']
    if not check_values(input_dict_local['orbital_period'], inf_value=0*u.d):
        print('ERROR: The orbital_period must be positive.')
        check =  False

    try:
        input_dict_local['planet_radius'] *= u.Unit( input_dict_local['planet_radius_unit'][0] )
    except ValueError:
        input_dict_local['planet_radius'] *= input_dict_local['star_radius'][0]
    del input_dict_local['planet_radius_unit']
    if not check_values(input_dict_local['planet_radius'], inf_value=0*u.m):
        print('ERROR: The planet_radius must be positive.')
        check =  False

    if 'star_effective_temperature' in input_keys:
        input_dict_local['star_effective_temperature'] *= u.K
        if not check_values(input_dict_local['star_effective_temperature'], inf_value=0*u.K):
            print('ERROR: The star_effective_temperature must be positive.')
            check =  False

    if 'system_distance' in input_keys:
        input_dict_local['system_distance'] *= u.Unit( input_dict_local['system_distance_unit'][0] )
        del input_dict_local['system_distance_unit']
        if not check_values(input_dict_local['system_distance'], min_value=1*u.AU):
            print('ERROR: The system_distance must cannot be shorter than the distance to our Sun.')
            check =  False

    rp_over_rs = ( input_dict_local['planet_radius'][0]/input_dict_local['star_radius'][0] ).decompose()
    sma_over_rs = ( input_dict_local['orbital_semimajor_axis'][0]/input_dict_local['star_radius'][0] ).decompose()
    inclination = input_dict_local['orbital_inclination'][0]
    period = input_dict_local['orbital_period'][0]
    T14 = get_transit_duration_T14(rp_over_rs, sma_over_rs, inclination, period)
    input_dict_local['transit_duration_T14'] = [T14]
    if not check_values(input_dict_local['transit_duration_T14'], inf_value=0*u.h):
        print('ERROR: The planet is not transiting.')
        check =  False

    input_dict_local['observing_duration_labels'] = []
    for value in input_dict_local['observing_duration']:
        label = 'obs_duration_' + str(value) + input_dict_local['observing_duration_unit'][0]
        input_dict_local['observing_duration_labels'] += [label,]

    try:
        input_dict_local['observing_duration'] *= u.Unit( input_dict_local['observing_duration_unit'][0] )
    except ValueError:
        input_dict_local['observing_duration'] *= T14
    del input_dict_local['observing_duration_unit']
    if np.isfinite(T14) and not check_values(input_dict_local['observing_duration'], inf_value=T14, max_value=input_dict_local['orbital_period']-T14):
        print('ERROR: The observing_duration must be longer than the transit and shorter than the orbital period minus the transit duration.')
        check =  False

    if 'planet_bond_albedo' not in input_keys:
        input_dict_local['planet_bond_albedo'] = [0.0]
    if not check_values(input_dict_local['planet_bond_albedo'], min_value=0.0, max_value=1.0):
        print('ERROR: The planet_bond_albedo must be between 0 and 1 (included).')
        check =  False

    if 'planet_circulation_efficiency' in input_keys:
        if not check_values(input_dict_local['planet_circulation_efficiency'], min_value=0.0, max_value=1.0):
            print('ERROR: The planet_circulation_efficiency must be between 0 and 1 (included).')
            check = False
        if check:
            star_effective_temperature = input_dict_local['star_effective_temperature'][0]
            planet_day_temperature = []
            planet_night_temperature = []
            planet_configuration_labels = []
            planet_albedo_adjusted = []
            planet_efficiency_adjusted = []
            #The following for loop consider all possible combinations of albedos and efficiencies
            for albedo in input_dict_local['planet_bond_albedo']:
                for efficiency in input_dict_local['planet_circulation_efficiency']:
                    [Tday, Tnight, label] = get_planet_temperatures(star_effective_temperature, sma_over_rs, albedo, efficiency)
                    planet_day_temperature += [Tday,]
                    planet_night_temperature += [Tnight,]
                    planet_configuration_labels += [label,]
                    planet_albedo_adjusted += [albedo,]
                    planet_efficiency_adjusted += [efficiency,]
            input_dict_local['planet_day_temperature'] = planet_day_temperature
            input_dict_local['planet_night_temperature'] = planet_night_temperature
            input_dict_local['planet_configuration_labels'] = planet_configuration_labels
            input_dict_local['planet_bond_albedo'] = planet_albedo_adjusted
            input_dict_local['planet_efficiency'] = planet_efficiency_adjusted

    if 'planet_day_temperature' in input_keys: #This can happen only if 'planet_circulation_efficiency' was not in the input_keys
        n_temp = len(input_dict_local['planet_day_temperature'])
        input_dict_local['planet_day_temperature'] *= u.K
        input_dict_local['planet_night_temperature'] *= u.K
        if not check_values(input_dict_local['planet_day_temperature'], min_value=0.0 * u.K):
            print('ERROR: The planet_day_temperature cannot be negative.')
            check = False
        if not check_values(input_dict_local['planet_night_temperature'], min_value=0.0 * u.K):
            print('ERROR: The planet_night_temperature cannot be negative.')
            check = False
        if check:
            planet_configuration_labels = []
            planet_albedo_adjusted = []
            planet_Tday_adjusted = []
            planet_Tnight_adjusted = []
            #The following for loop consider all possible combinations of albedos and ordered couples of day+night temperatures
            for albedo in input_dict_local['planet_bond_albedo']:
                for i in range(n_temp):
                    Tday = input_dict_local['planet_day_temperature'][i]
                    Tnight = input_dict_local['planet_night_temperature'][i]
                    label = 'albedo' + str(albedo) + '_Tday' + str(Tday.value) + 'K_Tnight' + str(Tnight.value) + 'K'
                    planet_configuration_labels += [label,]
                    planet_albedo_adjusted += [albedo,]
                    planet_Tday_adjusted += [Tday,]
                    planet_Tnight_adjusted += [Tnight,]
            input_dict_local['planet_configuration_labels'] = planet_configuration_labels
            input_dict_local['planet_bond_albedo'] = planet_albedo_adjusted
            input_dict_local['planet_day_temperature'] = planet_Tday_adjusted
            input_dict_local['planet_night_temperature'] = planet_Tnight_adjusted

    planet_models_grid = input_dict_local['planet_models_grid'][0]
    if planet_models_grid=='Userfile':
        planet_configuration_labels = []
        for albedo in input_dict_local['planet_bond_albedo']:
            label = 'albedo' + str(albedo)
            planet_configuration_labels += [label,]
        input_dict_local['planet_configuration_labels'] = planet_configuration_labels

    return check, input_dict_local


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
    elif stellar_models_grid == 'Stagger_2015':
        minimum_wavelength = 2000.172119140625 * u.Angstrom
        maximum_wavelength = 10000.0791015625 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength or np.max(pb_waves)>maximum_wavelength:
            check = False
    elif stellar_models_grid == 'Blackbody': #This check should be redundant
        minimum_wavelength = 0.0 * u.Angstrom
        if np.min(pb_waves)<minimum_wavelength:
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
    [pb, check] = read_as_numpy_array(os.path.join(passbands_path, passband)+passbands_ext)
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
    [wb, check] = read_as_numpy_array(os.path.join(wavelength_bins_path, wavelength_bins_file))
    if not check:
        print('WARNING: Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        return np.array([]) * u.Angstrom, check
    wb = np.atleast_2d(wb)
    if not check_2Darray(wb, n_col=2):
        print('WARNING: invalid format for wavelength_bins_file', os.path.join(wavelength_bins_path, wavelength_bins_file), '. It must have 2 columns. Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        check = False
        return np.array([]) * u.Angstrom, check
    wb *= u.Angstrom
    if (wb[:,0]>=wb[:,1]).any():
        print('WARNING: invalid line in wavelength_bins_file', os.path.join(wavelength_bins_path, wavelength_bins_file), '. The lower limit cannot be greater or equal to the upper limit.')
        check = False
        return np.array([]) * u.Angstrom, check
    if np.min(wb)<np.min(pb_waves) or np.max(wb)>np.max(pb_waves):
        print('WARNING:wavelength_bins_file', os.path.join(wavelength_bins_path, wavelength_bins_file), 'exceeds wavelength range for the ', passband, 'passband. Ignoring this file.')
        check = False
        return np.array([]) * u.Angstrom, check
    return wb, check


def get_waves_fromR(lambda1, lambda2, R):
    """
    This function computes an array of wavelengths with constant spectral resolution, including the minimum and maximum wavelengths.
    The requested resolution can be automatically increased to guarantee the minimum length of 10 for the array of wavelengths.
    
    :param float lambda1: minimum wavelength
    :param float lambda2: maximum wavelength
    :param float R: spectral resolution
    ..note:: the input params are already checked to avoid errors, if using the boats_calculate_transit or boats_calculate_eclipse functions
    :return: the 1D array of wavelengths
    :rtype: np.array
    """
    log_lambda1 = np.log(lambda1)
    log_lambda2 = np.log(lambda2)
    log_step = np.log( (2*R+1)/(2*R-1) )
    n_steps = np.ceil( (log_lambda2-log_lambda1)/log_step )
    n_steps = np.max([10, n_steps])
    log_waves = np.linspace(log_lambda1, log_lambda2, num=np.int(n_steps))
    waves = np.exp(log_waves)
    return waves


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
    return np.array((teff, logg, mh))


def get_stellar_grid_parameters(stellar_models_grid):
    """
    This function gets all the file names in the database and the corresponding stellar parameters.
    
    :param str file_name: built-in file name
    :return: the file names and a numpy array with the corresponding stellar parameters (number of files X 3)
    :rtype: list of str, np.array
    """
    files = databases[stellar_models_grid].get_filename_list()
    star_params_grid = np.zeros((len(files),3))
    for i in range(len(files)):
        star_params_grid[i,:] = stellar_params_from_file_name(files[i])
    return files, star_params_grid





def get_nearest_file_index(star_effective_temperature, star_log_gravity, star_metallicity, star_params_grid, stellar_models_grid):
    """
    This function returns the index of the nearest neighbour file to the requested stellar parameters in the chosen stellar models grid.
    First it selects the files with the closest temperature; then those with the closest log gravity among the selected ones; finally the one with the closest metallicity.
    
    :param float star_effective_temperature:
    :param float star_log_gravity:
    :param float star_metallicity:
    :param np.array star_params_grid: 3-column array with stellar effective temperatures, log gravities and metallicities of the database models
    :params str stellar_models_grid: the name of the chosen stellar database
    :return: the index of the nearest neighbour
    :rtype: int
    """
    min_teff_diff = np.min(np.abs(star_params_grid[:,0]*u.Kelvin-star_effective_temperature))
    indices_teff = np.where( np.abs(star_params_grid[:,0]*u.Kelvin-star_effective_temperature) == min_teff_diff )[0]
    min_logg_diff = np.min(np.abs(star_params_grid[indices_teff,1]-star_log_gravity))
    indices_logg = np.where( np.abs(star_params_grid[:,1]-star_log_gravity) == min_logg_diff )[0]
    indices_teff_logg = np.intersect1d(indices_teff, indices_logg)
    min_mh_diff = np.min(np.abs(star_params_grid[indices_teff_logg,2]-star_metallicity))
    indices_mh = np.where( np.abs(star_params_grid[:,2]-star_metallicity) == min_mh_diff )[0]
    indices_teff_logg_mh = np.intersect1d(indices_teff_logg, indices_mh)
    index_teff_logg_mh = indices_teff_logg_mh[0]
    return index_teff_logg_mh



def get_model_spectrum(models_grid, params=None, file_to_read=None):
    """
    This function returns the model spectrum from user file, blackbody calculation or built-in dataset
    
    :param str models_grid: the choice of stellar_models_grid or planet_models_grid
    :argument quantity array params: default is None
    :argument str file_to_read: user file to read, default is None
    :return: the model wavelengths (in Angstrom) and the corresponding flux (in erg/(cm^2 s A))
    :rtype: quantity array, quantity array
    """
    if models_grid=='Userfile':
        try:
            model_file = pickle.load(open(file_to_read, 'rb'))
            model_wavelengths = model_file['wavelengths'].to(u.Angstrom)
            if model_file['fluxes'].unit == u.dimensionless_unscaled:
                model_fluxes = model_file['fluxes']
            else:
                model_fluxes = model_file['fluxes'].to(u.erg / (u.cm**2 * u.second * u.Angstrom), equivalencies=u.spectral_density(model_wavelengths))
        except:
            print('ERROR: Something went wrong when reading the file', file_to_read)
            exit()            
    elif models_grid=='Blackbody':
        blackbody_temperature = params[0]
        model_wavelengths = get_waves_fromR(1.0, 2000000.0, 10000.0) * u.Angstrom
        model_fluxes = blackbody_lambda(model_wavelengths, blackbody_temperature)
        model_fluxes *= np.pi * u.sr
    else:
        [star_files_grid, star_params_grid] = get_stellar_grid_parameters(models_grid) #Reading file names and stellar parameters from the selected database
        star_effective_temperature =  params[0]
        star_log_gravity =  params[1]
        star_metallicity =  params[2]
        neighbour_index = get_nearest_file_index(star_effective_temperature, star_log_gravity, star_metallicity, star_params_grid, models_grid)
        print('WARNING: Adopting nearest model in the', models_grid, 'grid: Teff=', star_params_grid[neighbour_index,0]*u.Kelvin, ', logg=', star_params_grid[neighbour_index,1], ', [M/H]=', star_params_grid[neighbour_index,2])
        neighbour_model = databases[models_grid].get_file_content(dbx_file=star_files_grid[neighbour_index])
        model_wavelengths = neighbour_model['wavelengths']
        model_fluxes = neighbour_model['fluxes']
    return model_wavelengths, model_fluxes


def get_photon_spectrum(model_wavelengths, model_fluxes, obj_radius, obj_distance, telescope_area):
    """
    This function computes the model photon fluxes at the telescope primary from the energy fluxes at the given model wavelengths
    
    :param quantity array models_wavelengths:
    :param quantity array models_fluxes: the flux should be expressed in erg/(cm^2 s A)
    :param quantity obj_radius: the radius of the star or planet
    :param quantity obj_distance: the distance of the system from Earth
    :param quantity telescope_area: the collecting area of the telescope
    :return: the photon fluxes in photon/(s A)
    :rtype: quantity array
    """
    scaling_factor = (obj_radius/obj_distance)**2.0 * telescope_area
    scaling_factor = scaling_factor.decompose()
    scaled_fluxes = model_fluxes * scaling_factor
    photon_fluxes = scaled_fluxes.to(u.photon / (u.second * u.Angstrom), equivalencies=u.spectral_density(model_wavelengths))
    return photon_fluxes

def get_relative_spectrum(relative_model_wavelengths, relative_model_fluxes, base_model_wavelengths, base_model_fluxes):
    """
    This function computes the absolute spectral flux given the spectrum relative to a base model, e.g., the planet emission spectrum given the star spectrum and eclipse depths.
    
    :param quantity array relative_model_wavelengths: it should be in Angstrom
    :param quantity array relative_model_fluxes: it should be dimensionless
    :param quantity base_model_wavelengths: it should be in Angstrom
    :param quantity base_model_fluxes:
    :return: the scaled fluxes in absolute units and the corresponding wavelengths
    ..note:: the scaled fluxes are calculated at the base_model_wavelengths
    :rtype: quantity array
    """
    lim1 = np.max([np.min(relative_model_wavelengths.value), np.min(base_model_wavelengths.value)])
    lim2 = np.min([np.max(relative_model_wavelengths.value), np.max(base_model_wavelengths.value)])
    my_waves = get_waves_fromR(lim1, lim2, 10000.0) * u.Angstrom
    f_interp = interp1d(relative_model_wavelengths.value, relative_model_fluxes.value, fill_value='extrapolate')
    relative_interp_fluxes = f_interp(my_waves.value)
    g_interp = interp1d(base_model_wavelengths.value, base_model_fluxes.value, fill_value='extrapolate')
    base_interp_fluxes = g_interp(my_waves.value)
    scaled_fluxes = relative_interp_fluxes * base_interp_fluxes * base_model_fluxes.unit
    return my_waves, scaled_fluxes

def get_passband_fluxes(model_wavelengths, photon_fluxes, passbands_dict):
    """
    This function computes the electrons rate for all the requested instruments/passbands.
    
    :param quantity array model_wavelengths: 1D array with the wavelengths of the precalculated model
    :param quantity array photon_fluxes: 1D array with the corresponding photon fluxes in photon/(s A)
    :param dict passbands_dict: dictionary with the spectral responses to use
    :return: dictionary with all the electron rates for the requested passbands
    :rtype: dict
    """
    passbands = list(passbands_dict.keys())
    photons_dict = {}
    for passband in passbands:
        my_waves = passbands_dict[passband][0]
        my_pce = passbands_dict[passband][1]
        if np.min(my_waves)<np.min(model_wavelengths)-1*u.Angstrom or np.max(my_waves)>np.max(model_wavelengths)+1*u.Angstrom:
            print('WARNING: The passband', passband, 'is out of the model limits. It will be ignored.')
        else:
            f_interp = interp1d(model_wavelengths.value, photon_fluxes.value, fill_value='extrapolate')
            my_photons = f_interp(my_waves.value)
            photons_dict[passband] = (simps(my_photons*my_pce.value, my_waves.value) * photon_fluxes.unit * my_pce.unit * my_waves.unit).decompose()
    return photons_dict


def compute_phase_average(phi1, phi2):
    """
    This function computes the mean integral of a sinusoid between two phases, expressed in units of the orbital period.
    
    :param float phi1: lower bound of the integral
    :param float phi2: upper bound of the integral
    :return: the mean integral of a sinusoid between the two phases
    :rtype: float
    """
    phi1 *= 2*np.pi
    phi2 *= 2*np.pi
    phase_average = 1.0 - ( np.sin(phi2 * u.rad) - np.sin(phi1 * u.rad) )/( phi2 - phi1 )
    phase_average *= 0.5
    return phase_average
    





def process_configuration_transit(input_dict):
    """
    This function executes all the operations from the formatted input dictionary to create the requested output files (for the case of a transit).
    
    :param dict input_dict:
    ..note:: all the input dictionary branches end with a list, even if they contain a single element.
    :return: dictionary containing the information about the input configuration and the calculated output including the transit depth bias
    :rtype: dict
    """

    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())

    #Getting the output keywords
    output_path = input_dict_local['output_path'][0]
    output_filename = input_dict_local['output_filename'][0]
    output_fileext = input_dict_local['output_fileext']

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
            my_waves = get_waves_fromR(np.min(pb_wavelengths.value), np.max(pb_wavelengths.value), 10000.0) * pb_wavelengths.unit #choice of wavelengths
            my_pce = f_interp(my_waves.value) * pb_pce.unit
            passbands_dict[passbands[i]] = [my_waves, my_pce]
            [wb, check] = read_wavelength_bins(wavelength_bins_path, wavelength_bins_files[i], pb_wavelengths, passbands[i]) #includes no_bins case
            if check:
                for wbin in wb:
                    my_waves = get_waves_fromR(wbin[0].value, wbin[1].value, 10000.0) * pb_wavelengths.unit
                    my_pce = f_interp(my_waves.value) * pb_pce.unit
                    passbands_dict[passbands[i]+'_'+str(wbin[0].value)+'_'+str(wbin[1].value)] = [my_waves, my_pce]

    #Computing the electrons rate due to the stellar flux in each passband/wavelength bin and storing into a dictionary:
    #model spectrum (emergent energy flux) -> photon spectrum (at Earth) -> electrons rate (detectors)
    star_radius = input_dict_local['star_radius'][0]
    telescope_area = input_dict_local['telescope_area'][0]
    if stellar_models_grid=='Userfile':
        file_to_read = os.path.join(input_dict_local['star_model_path'][0], input_dict_local['star_model_file'][0])
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, file_to_read=file_to_read)
        if input_dict_local['rescale_star_flux'][0] == 'Yes': #case of emergent flux to be rescaled
                system_distance = input_dict_local['system_distance'][0]
                star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
        else: #case of flux given at Earth
                star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, star_radius, telescope_area)
    elif stellar_models_grid=='Blackbody':
        star_effective_temperature = input_dict_local['star_effective_temperature'][0]
        params = [star_effective_temperature]
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, params=params)
        system_distance = input_dict_local['system_distance'][0]
        star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
    else: #database models
        star_effective_temperature = input_dict_local['star_effective_temperature'][0]
        star_log_gravity = input_dict_local['star_log_gravity'][0]
        star_metallicity = input_dict_local['star_metallicity'][0]
        params = [star_effective_temperature, star_log_gravity, star_metallicity]
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, params=params)
        system_distance = input_dict_local['system_distance'][0]
        star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
    star_electrons_rate_dict = {}
    star_electrons_rate_dict = get_passband_fluxes(star_model_wavelengths, star_photon_spectrum, passbands_dict)

    #Computing the electrons rate due to the planetary day and nightside flux in each passband/wavelength bin and storing into two dictionaries:
    #model spectrum (emergent energy flux) -> photon spectrum (at Earth) -> electrons rate (detectors)
    planet_models_grid = input_dict_local['planet_models_grid'][0]
    planet_radius = input_dict_local['planet_radius'][0]
    planet_day_electrons_rate_dict = {}
    planet_night_electrons_rate_dict = {}
    labels = input_dict_local['planet_configuration_labels']
    n_conf = len(labels)
    if planet_models_grid=='Userfile':
        day_file_to_read = os.path.join(input_dict_local['planet_day_model_path'][0], input_dict_local['planet_day_model_file'][0])
        night_file_to_read = os.path.join(input_dict_local['planet_night_model_path'][0], input_dict_local['planet_night_model_file'][0])
        for i in range(n_conf):
            label = labels[i]
            [planet_day_model_wavelengths, planet_day_model_fluxes] = get_model_spectrum(planet_models_grid, file_to_read=day_file_to_read)
            [planet_night_model_wavelengths, planet_night_model_fluxes] = get_model_spectrum(planet_models_grid, file_to_read=night_file_to_read)
            if input_dict_local['rescale_planet_flux'][0] == 'Star':
                [planet_day_model_wavelengths, planet_day_photon_spectrum] = get_relative_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, star_model_wavelengths, star_photon_spectrum)
                [planet_night_model_wavelengths, planet_night_photon_spectrum] = get_relative_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, star_model_wavelengths, star_photon_spectrum)
            elif input_dict_local['rescale_planet_flux'][0] == 'Yes':
                system_distance = input_dict_local['system_distance'][0]
                planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, system_distance, telescope_area)
                planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, system_distance, telescope_area)
            else:
                planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, planet_radius, telescope_area)
                planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, planet_radius, telescope_area)
            planet_day_electrons_rate_dict[label] = {}
            planet_night_electrons_rate_dict[label] = {}
            planet_day_electrons_rate_dict[label] = get_passband_fluxes(planet_day_model_wavelengths, planet_day_photon_spectrum, passbands_dict)
            planet_night_electrons_rate_dict[label] = get_passband_fluxes(planet_night_model_wavelengths, planet_night_photon_spectrum, passbands_dict)
    elif planet_models_grid=='Blackbody':
        planet_day_temperatures = input_dict_local['planet_day_temperature']
        planet_night_temperatures = input_dict_local['planet_night_temperature']
        system_distance = input_dict_local['system_distance'][0]
        for i in range(n_conf):
            label = labels[i]
            [planet_day_model_wavelengths, planet_day_model_fluxes] = get_model_spectrum(planet_models_grid, params=np.atleast_1d(planet_day_temperatures[i]))
            [planet_night_model_wavelengths, planet_night_model_fluxes] = get_model_spectrum(planet_models_grid, params=np.atleast_1d(planet_night_temperatures[i]))
            planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, system_distance, telescope_area)
            planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, system_distance, telescope_area)
            planet_day_electrons_rate_dict[label] = {}
            planet_night_electrons_rate_dict[label] = {}
            planet_day_electrons_rate_dict[label] = get_passband_fluxes(planet_day_model_wavelengths, planet_day_photon_spectrum, passbands_dict)
            planet_night_electrons_rate_dict[label] = get_passband_fluxes(planet_night_model_wavelengths, planet_night_photon_spectrum, passbands_dict)

    #Getting keywords for planet configurations and observation settings to perform calculations
    planet_albedo = input_dict_local['planet_bond_albedo']
    orbital_semimajor_axis = input_dict_local['orbital_semimajor_axis'][0]
    T14 = input_dict_local['transit_duration_T14'][0]
    orbital_period = input_dict_local['orbital_period'][0]
    observing_duration = input_dict_local['observing_duration']
    labels_obsdur = input_dict_local['observing_duration_labels']
    n_dur = len(observing_duration)
    phi14 = (T14 / orbital_period).decompose()
    phimax = (0.5 * observing_duration / orbital_period).decompose()
    transit_depth = ( (planet_radius / star_radius)**2.0 ).decompose()
    ppm = u.def_unit('ppm')
    passbands_star = list(star_electrons_rate_dict.keys()) #Passbands that were successfully calculated for the star model
    #Performing the final calculation and storing into a dictionary of results
    results_dict = {}
    for i in range(n_conf):
        label = labels[i]
        try:
            planet_day_temperature = input_dict_local['planet_day_temperature'][i]
            planet_night_temperature = input_dict_local['planet_night_temperature'][i]
        except:
            planet_day_temperature = 'not_given'
            planet_night_temperature = 'not_given'
        results_dict[label] = {}
        reflection_factor = planet_albedo[i] * (0.5 * planet_radius / orbital_semimajor_axis)**2.0
        reflection_factor = reflection_factor.decompose()
        passbands_day = list(planet_day_electrons_rate_dict[label].keys()) #Passbands that were successfully calculated for the planet dayside model
        passbands_night = list(planet_night_electrons_rate_dict[label].keys()) #Passbands that were successfully calculated for the planet nightside model
        passbands = [p for p in passbands_star if p in np.intersect1d(passbands_day, passbands_night)]
        for passband in passbands:
            results_dict[label][passband] = {}
            star_electrons_rate = copy.deepcopy(star_electrons_rate_dict[passband])
            planet_day_electrons_rate = copy.deepcopy(planet_day_electrons_rate_dict[label][passband])
            planet_day_electrons_rate += reflection_factor*star_electrons_rate
            planet_night_electrons_rate = copy.deepcopy(planet_night_electrons_rate_dict[label][passband])
            day_phase_average_in = compute_phase_average(0.0, 0.5*phi14)
            night_phase_average_in = 1.0 - day_phase_average_in
            planet_electrons_rate_in = planet_day_electrons_rate*day_phase_average_in + planet_night_electrons_rate*night_phase_average_in
            for j in range(n_dur):
                label2 = labels_obsdur[j]
                results_dict[label][passband][label2] = {}
                day_phase_average_out = compute_phase_average(0.5*phi14, phimax[j])
                night_phase_average_out = 1.0 - day_phase_average_out
                planet_electrons_rate_out = planet_day_electrons_rate*day_phase_average_out + planet_night_electrons_rate*night_phase_average_out
                self_blend = star_electrons_rate / (star_electrons_rate + planet_electrons_rate_out)
                phase_blend = (planet_electrons_rate_out - planet_electrons_rate_in) / (star_electrons_rate + planet_electrons_rate_out)
                transit_depth_biased = self_blend*transit_depth + phase_blend
                results_dict[label][passband][label2] = {}
                results_dict[label][passband][label2]['transit_depth_bias'] = (transit_depth_biased - transit_depth)*1e6 * ppm
                results_dict[label][passband][label2]['phase_blend_bias'] = phase_blend*1e6 * ppm
                results_dict[label][passband][label2]['self_blend_bias'] = results_dict[label][passband][label2]['transit_depth_bias'] - results_dict[label][passband][label2]['phase_blend_bias']
                results_dict[label][passband][label2]['transit_depth'] = transit_depth*1e6 * ppm
                results_dict[label][passband][label2]['transit_duration_T14'] = T14
                results_dict[label][passband][label2]['observing_duration'] = observing_duration[j]
                results_dict[label][passband][label2]['planet_day_temperature'] = planet_day_temperature
                results_dict[label][passband][label2]['planet_night_temperature'] = planet_night_temperature
                results_dict[label][passband][label2]['star_flux'] = star_electrons_rate
                results_dict[label][passband][label2]['planet_day_flux'] = planet_day_electrons_rate
                results_dict[label][passband][label2]['planet_night_flux'] = planet_night_electrons_rate
                results_dict[label][passband][label2]['planet_flux_oot'] = planet_electrons_rate_out
                results_dict[label][passband][label2]['planet_flux_in'] = planet_electrons_rate_in
                #S/N calculation to be added
                number_electrons_out = ( (star_electrons_rate + planet_electrons_rate_out)*(observing_duration[j]-T14) ).decompose()
                number_electrons_in = ( (star_electrons_rate*(1.0-transit_depth) + planet_electrons_rate_in)*T14 ).decompose()
                transit_depth_sigma = (1.0-transit_depth) * np.sqrt((1.0/number_electrons_in.value) + (1.0/number_electrons_out.value))*1e6 * ppm
                results_dict[label][passband][label2]['transit_depth_sigma'] = transit_depth_sigma

    #Saving the final dictionary with input and results
    if '.pickle' in output_fileext:
        final_dict = {}
        final_dict['input_info'] = copy.deepcopy(input_dict_local)
        final_dict['results'] = copy.deepcopy(results_dict)
        with open(os.path.join(output_path, output_filename+'.pickle'), 'wb') as outfile:
            pickle.dump(final_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    if '.txt' in output_fileext:
        text_file = open(os.path.join(output_path, output_filename+'.txt'), 'w')
        for label in labels:
            for label2 in labels_obsdur:
                text_file.write(label+'_'+label2+'\n')
                text_file.write('passband' + '\t' + 'transit_depth' + ' (ppm)\t' + 'transit_depth_sigma' + '\t' + 'transit_depth_bias' + '\t' + 'self_blend_bias' + '\t' + 'phase_blend_bias' + '\n')
                for passband in passbands:
                    text_file.write(passband + '\t' + str(results_dict[label][passband][label2]['transit_depth'].value) + '\t' + str(results_dict[label][passband][label2]['transit_depth_sigma'].value) + '\t' + str(results_dict[label][passband][label2]['transit_depth_bias'].value) + '\t' + str(results_dict[label][passband][label2]['self_blend_bias'].value) + '\t' + str(results_dict[label][passband][label2]['phase_blend_bias'].value) + '\n')

    return final_dict




def process_configuration_eclipse(input_dict):
    """
    This function executes all the operations from the formatted input dictionary to create the requested output files.
    
    :param dict: input_dict
    ..note:: all the input dictionary branches end with a list, even if they contain a single element.
    :return: dictionary containing the information about the input configuration and the calculated output including the transit depth bias
    :rtype: dict
    """

    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())

    #Getting the output keywords
    output_path = input_dict_local['output_path'][0]
    output_filename = input_dict_local['output_filename'][0]
    output_fileext = input_dict_local['output_fileext']

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
            my_waves = get_waves_fromR(np.min(pb_wavelengths.value), np.max(pb_wavelengths.value), 10000.0) * pb_wavelengths.unit #choice of wavelengths
            my_pce = f_interp(my_waves.value) * pb_pce.unit
            passbands_dict[passbands[i]] = [my_waves, my_pce]
            [wb, check] = read_wavelength_bins(wavelength_bins_path, wavelength_bins_files[i], pb_wavelengths, passbands[i]) #includes no_bins case
            if check:
                for wbin in wb:
                    my_waves = get_waves_fromR(wbin[0].value, wbin[1].value, 10000.0) * pb_wavelengths.unit
                    my_pce = f_interp(my_waves.value) * pb_pce.unit
                    passbands_dict[passbands[i]+'_'+str(wbin[0].value)+'_'+str(wbin[1].value)] = [my_waves, my_pce]

    #Computing the electrons rate due to the stellar flux in each passband/wavelength bin and storing into a dictionary:
    #model spectrum (emergent energy flux) -> photon spectrum (at Earth) -> electrons rate (detectors)
    star_radius = input_dict_local['star_radius'][0]
    telescope_area = input_dict_local['telescope_area'][0]
    if stellar_models_grid=='Userfile':
        file_to_read = input_dict_local['star_model_path'][0] + input_dict_local['star_model_file'][0]
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, file_to_read=file_to_read)
        if input_dict_local['rescale_star_flux'][0] == 'Yes': #case of emergent flux to be rescaled
                system_distance = input_dict_local['system_distance'][0]
                star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
        else: #case of flux given at Earth
                star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, star_radius, telescope_area)
    elif stellar_models_grid=='Blackbody':
        star_effective_temperature = input_dict_local['star_effective_temperature'][0]
        params = [star_effective_temperature]
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, params=params)
        system_distance = input_dict_local['system_distance'][0]
        star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
    else: #database models
        star_effective_temperature = input_dict_local['star_effective_temperature'][0]
        star_log_gravity = input_dict_local['star_log_gravity'][0]
        star_metallicity = input_dict_local['star_metallicity'][0]
        params = [star_effective_temperature, star_log_gravity, star_metallicity]
        [star_model_wavelengths, star_model_fluxes] = get_model_spectrum(stellar_models_grid, params=params)
        system_distance = input_dict_local['system_distance'][0]
        star_photon_spectrum = get_photon_spectrum(star_model_wavelengths, star_model_fluxes, star_radius, system_distance, telescope_area)
    star_electrons_rate_dict = get_passband_fluxes(star_model_wavelengths, star_photon_spectrum, passbands_dict)

    #Computing the electrons rate due to the planetary day and nightside flux in each passband/wavelength bin and storing into two dictionaries:
    #model spectrum (emergent energy flux) -> photon spectrum (at Earth) -> electrons rate (detectors)
    planet_models_grid = input_dict_local['planet_models_grid'][0]
    planet_radius = input_dict_local['planet_radius'][0]
    planet_day_electrons_rate_dict = {}
    planet_night_electrons_rate_dict = {}
    labels = input_dict_local['planet_configuration_labels']
    n_conf = len(labels)
    if planet_models_grid=='Userfile':
        day_file_to_read = os.path.join(input_dict_local['planet_day_model_path'][0], input_dict_local['planet_day_model_file'][0])
        night_file_to_read = os.path.join(input_dict_local['planet_night_model_path'][0], input_dict_local['planet_night_model_file'][0])
        for i in range(n_conf):
            label = labels[i]
            [planet_day_model_wavelengths, planet_day_model_fluxes] = get_model_spectrum(planet_models_grid, file_to_read=day_file_to_read)
            [planet_night_model_wavelengths, planet_night_model_fluxes] = get_model_spectrum(planet_models_grid, file_to_read= night_file_to_read)
            if input_dict_local['rescale_planet_flux'][0] == 'Star':
                [planet_day_model_wavelengths, planet_day_photon_spectrum] = get_relative_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, star_model_wavelengths, star_photon_spectrum)
                [planet_night_model_wavelengths, planet_night_photon_spectrum] = get_relative_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, star_model_wavelengths, star_photon_spectrum)
            elif input_dict_local['rescale_planet_flux'][0] == 'Yes':
                system_distance = input_dict_local['system_distance'][0]
                planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, system_distance, telescope_area)
                planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, system_distance, telescope_area)
            else:
                planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, planet_radius, telescope_area)
                planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, planet_radius, telescope_area)
            planet_day_electrons_rate_dict[label] = {}
            planet_night_electrons_rate_dict[label] = {}
            planet_day_electrons_rate_dict[label] = get_passband_fluxes(planet_day_model_wavelengths, planet_day_photon_spectrum, passbands_dict)
            planet_night_electrons_rate_dict[label] = get_passband_fluxes(planet_night_model_wavelengths, planet_night_photon_spectrum, passbands_dict)
    elif planet_models_grid=='Blackbody':
        planet_day_temperatures = input_dict_local['planet_day_temperature']
        planet_night_temperatures = input_dict_local['planet_night_temperature']
        system_distance = input_dict_local['system_distance'][0]
        for i in range(n_conf):
            label = labels[i]
            [planet_day_model_wavelengths, planet_day_model_fluxes] = get_model_spectrum(planet_models_grid, params=np.atleast_1d(planet_day_temperatures[i]))
            [planet_night_model_wavelengths, planet_night_model_fluxes] = get_model_spectrum(planet_models_grid, params=np.atleast_1d(planet_night_temperatures[i]))
            planet_day_photon_spectrum = get_photon_spectrum(planet_day_model_wavelengths, planet_day_model_fluxes, planet_radius, system_distance, telescope_area)
            planet_night_photon_spectrum = get_photon_spectrum(planet_night_model_wavelengths, planet_night_model_fluxes, planet_radius, system_distance, telescope_area)
            planet_day_electrons_rate_dict[label] = {}
            planet_night_electrons_rate_dict[label] = {}
            planet_day_electrons_rate_dict[label] = get_passband_fluxes(planet_day_model_wavelengths, planet_day_photon_spectrum, passbands_dict)
            planet_night_electrons_rate_dict[label] = get_passband_fluxes(planet_night_model_wavelengths, planet_night_photon_spectrum, passbands_dict)

    #Getting keywords for planet configurations and observation settings to perform calculations
    planet_albedo = input_dict_local['planet_bond_albedo']
    orbital_semimajor_axis = input_dict_local['orbital_semimajor_axis'][0]
    T14 = input_dict_local['transit_duration_T14'][0]
    orbital_period = input_dict_local['orbital_period'][0]
    observing_duration = input_dict_local['observing_duration']
    labels_obsdur = input_dict_local['observing_duration_labels']
    n_dur = len(observing_duration)
    phi14 = (T14 / orbital_period).decompose()
    phimax = (0.5 * observing_duration / orbital_period).decompose()
    transit_depth = ( (planet_radius / star_radius)**2.0 ).decompose()
    ppm = u.def_unit('ppm')
    passbands_star = list(star_electrons_rate_dict.keys()) #Passbands that were successfully calculated for the star model
    #Performing the final calculation and storing into a dictionary of results
    results_dict = {}
    for i in range(n_conf):
        label = labels[i]
        try:
            planet_day_temperature = input_dict_local['planet_day_temperature'][i]
            planet_night_temperature = input_dict_local['planet_night_temperature'][i]
        except:
            planet_day_temperature = 'not_given'
            planet_night_temperature = 'not_given'
        results_dict[label] = {}
        reflection_factor = planet_albedo[i] * (0.5 * planet_radius / orbital_semimajor_axis)**2.0
        reflection_factor = reflection_factor.decompose()
        passbands_day = list(planet_day_electrons_rate_dict[label].keys()) #Passbands that were successfully calculated for the planet dayside model
        passbands_night = list(planet_night_electrons_rate_dict[label].keys()) #Passbands that were successfully calculated for the planet nightside model
        passbands = [p for p in passbands_star if p in np.intersect1d(passbands_day, passbands_night)]
        for passband in passbands:
            results_dict[label][passband] = {}
            star_electrons_rate = copy.deepcopy(star_electrons_rate_dict[passband])
            planet_day_electrons_rate = copy.deepcopy(planet_day_electrons_rate_dict[label][passband])
            planet_day_electrons_rate += reflection_factor*star_electrons_rate
            planet_night_electrons_rate = copy.deepcopy(planet_night_electrons_rate_dict[label][passband])
            day_phase_average_in = compute_phase_average(0.5, 0.5*(1+phi14))
            night_phase_average_in = 1.0 - day_phase_average_in
            planet_electrons_rate_in = planet_day_electrons_rate*day_phase_average_in + planet_night_electrons_rate*night_phase_average_in
            for j in range(n_dur):
                label2 = labels_obsdur[j]
                results_dict[label][passband][label2] = {}
                day_phase_average_out = compute_phase_average(0.5*(1+phi14), 0.5+phimax[j])
                night_phase_average_out = 1.0 - day_phase_average_out
                planet_electrons_rate_out = planet_day_electrons_rate*day_phase_average_out + planet_night_electrons_rate*night_phase_average_out
                eclipse_depth_measured = planet_electrons_rate_out / star_electrons_rate
                eclipse_depth_average_in = planet_electrons_rate_in / star_electrons_rate
                eclipse_depth_peak_in = planet_day_electrons_rate / star_electrons_rate
                results_dict[label][passband][label2] = {}
                results_dict[label][passband][label2]['eclipse_depth_measured'] = eclipse_depth_measured*1e6 * ppm
                results_dict[label][passband][label2]['eclipse_depth_average_in'] = eclipse_depth_average_in*1e6 * ppm
                results_dict[label][passband][label2]['eclipse_depth_peak_in'] = eclipse_depth_peak_in*1e6 * ppm
                results_dict[label][passband][label2]['eclipse_duration_T14'] = T14
                results_dict[label][passband][label2]['observing_duration'] = observing_duration[j]
                results_dict[label][passband][label2]['planet_day_temperature'] = planet_day_temperature
                results_dict[label][passband][label2]['planet_night_temperature'] = planet_night_temperature
                results_dict[label][passband][label2]['star_flux'] = star_electrons_rate
                results_dict[label][passband][label2]['planet_day_flux'] = planet_day_electrons_rate
                results_dict[label][passband][label2]['planet_night_flux'] = planet_night_electrons_rate
                results_dict[label][passband][label2]['planet_flux_ooe'] = planet_electrons_rate_out
                results_dict[label][passband][label2]['planet_flux_ine'] = planet_electrons_rate_in
                #S/N calculation to be added
                number_electrons_out = ( (star_electrons_rate + planet_electrons_rate_out)*(observing_duration[j]-T14) ).decompose()
                number_electrons_in = ( star_electrons_rate*T14 ).decompose()
                eclipse_depth_sigma = (1+eclipse_depth_measured) * np.sqrt((1.0/number_electrons_out.value) + (1.0/number_electrons_in.value))*1e6 * ppm
                results_dict[label][passband][label2]['eclipse_depth_sigma'] = eclipse_depth_sigma

    #Saving the final dictionary with input and results
    if '.pickle' in output_fileext:
        final_dict = {}
        final_dict['input_info'] = copy.deepcopy(input_dict_local)
        final_dict['results'] = copy.deepcopy(results_dict)
        with open(os.path.join(output_path, output_filename+'.pickle'), 'wb') as outfile:
            pickle.dump(final_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)

    if '.txt' in output_fileext:
        text_file = open(os.path.join(output_path, output_filename+'.txt'), 'w')
        for label in labels:
            for label2 in labels_obsdur:
                text_file.write(label+'_'+label2+'\n')
                text_file.write('passband' + '\t' + 'eclipse_depth_measured' + ' (ppm)\t' + 'eclipse_depth_sigma' + '\t' + 'eclipse_depth_average_in' + '\t' + 'eclipse_depth_peak_in' + '\n')
                for passband in passbands:
                    text_file.write(passband + '\t' + str(results_dict[label][passband][label2]['eclipse_depth_measured'].value) + '\t' + str(results_dict[label][passband][label2]['eclipse_depth_sigma'].value) + '\t' + str(results_dict[label][passband][label2]['eclipse_depth_average_in'].value) + '\t' + str(results_dict[label][passband][label2]['eclipse_depth_peak_in'].value) + '\n')

    return final_dict













def boats_calculate_transit(configuration_file):
    """
    This is the main function to run for calculating the limb-darkening coefficients from the input file.
    It calls the three functions to read, check and process the input file.
    
    :param str configuration_file: absolute or relative path including the file name
    
    :return:
    :rtype:
    """
    input_dict0 = read_configuration(configuration_file) #Reading configuration file
    [check, input_dict1] = check_configuration(input_dict0) #Checking configuration file and slight dictionary update
    if check:
        [check, input_dict] = create_new_dict(input_dict1)
    if check:
        process_configuration_transit(input_dict) #Computing and saving
    else:
        print('Something went wrong with the input.')


def boats_calculate_eclipse(configuration_file):
    """
    This is the main function to run for calculating the limb-darkening coefficients from the input file.
    It calls the three functions to read, check and process the input file.
    
    :param str configuration_file: absolute or relative path including the file name
    
    :return:
    :rtype:
    """
    input_dict0 = read_configuration(configuration_file) #Reading configuration file
    [check, input_dict1] = check_configuration(input_dict0) #Checking configuration file and slight dictionary update
    if check:
        [check, input_dict] = create_new_dict(input_dict1)
    if check:
        process_configuration_eclipse(input_dict) #Computing and saving
    else:
        print('Something went wrong with the input.')





