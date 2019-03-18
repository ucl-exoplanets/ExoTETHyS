from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *

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


def str2float(s):
#This function convert a string to a float, if the string is a number.
    try:
        return float(s)
    except ValueError:
        return s


def my_vstack(array1,array2):
#This function modifies the numpy.vstack to avoid error if one of the two arrays is empty.
#WARNING: It also works if none of the arrays is empty but they have different lengths, which is not expected to happen in this code.
    try:
        return np.vstack((array1,array2))
    except ValueError:
        if len(array1)<len(array2):
            return array2
        elif len(array1)>len(array2):
            return array1

def check_length(vector, min_length=1, max_length=None):
#This function checks if the length of a vector is within the expected range and returns a boolean value.
#Default is length>=1.
    check = True
    if len(vector)<min_length:
        check = False
    elif max_length:
        if len(vector)>max_length:
            check = False
    return check

def check_type_in_list(vector, item_type):
#This function checks that all the elements of a list are of the expected variable type (e.g., string, float, int) and returns a boolean value.
    check = True
    for item in vector:
        if not isinstance(item, item_type):
            check = False
    return check

def check_integers(vector, min_value=1):
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

def read_as_numpy_array(path_file):
#This function read a numpy array from file as numpy.genfromtxt but also checks for the most common errors.
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


def check_2Darray(arr, n_col=None):
#This function checks that an array has 2 dimensions and returns a boolean value.
#Optionally, it may also require the exact number of columns. 
    check = True
    if arr.ndim != 2:
        check = False
    elif n_col:
        if np.shape(arr)[1] != n_col:
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
                key = content[0]
                value = []
                for item in content[1:]:
                    if item[0] != '!':
                        value += [str2float(item),]
                input_dict[key] = value
    return input_dict


def copy_dict(dict):
    copied_dict = {}

    for i in dict:
        copied_dict[i] = dict[i]

    return copied_dict


def check_configuration(input_dict):
#This function checks whether the dictionary created from the input file is fine and returns a boolean value.
    check = True
    input_dict_local = copy_dict(input_dict)
    input_keys = list(input_dict_local.keys())
    mandatory_keys = ['calculation_type', 'stellar_models_grid', 'limb_darkening_laws', 'passbands']
    allowed_keys = mandatory_keys + ['target_names', 'gen_claret_orders', 'gen_poly_orders', 'star_effective_temperature', 'star_log_gravity', 'star_metallicity', 'star_minimum_effective_temperature', 'star_maximum_effective_temperature', 'star_minimum_log_gravity', 'star_maximum_log_gravity', 'star_minimum_metallicity', 'star_maximum_metallicity', 'wavelength_bins_files', 'user_output', 'output_path']

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

    #Checking the choice of stellar_models_grid among those available. Only one grid can be selected for each execution.
    stellar_models_grid = input_dict_local['stellar_models_grid']
    if not check_length(stellar_models_grid, max_length=1):
        print('ERROR: invalid length=', len(stellar_models_grid), 'for stellar_models_grid. It must have length=1.')
        check = False
    else:
        allowed_stellar_models_grid = ['Phoenix_2018', 'Phoenix_2012_13']
        stellar_models_grid = stellar_models_grid[0]
        if stellar_models_grid not in allowed_stellar_models_grid:
            print('ERROR:',stellar_models_grid,'is not a valid stellar_models_grid. The allowed names are Phoenix_2018, Phoenix_2012_13.')
            check = False

    #Checking the choice of limb_darkening_laws; at least one law must be specified in the input file (no default).
    #TO DO: User-defined limb_darkening_laws should be allowed, but this option is not implemented yet.
    limb_darkening_laws = input_dict_local['limb_darkening_laws']
    if not check_length(limb_darkening_laws):
        print('ERROR: invalid length=', len(limb_darkening_laws), 'for limb_darkening_laws. It must have length>=1.')
        check = False
    allowed_limb_darkening_laws = ['linear', 'quadratic', 'square_root', 'power2', 'claret4', 'gen_claret', 'gen_poly'] #user-defined should be added
    for item in limb_darkening_laws:
        if item not in allowed_limb_darkening_laws:
            print('ERROR:',item,'is not a valid limb_darkening_laws.')
            check = False
        if limb_darkening_laws.count(item)>1:
            print('WARNING:',item,'limb_darkening_laws entered multiple times. Repetitions are ignored.')
        input_dict_local['limb_darkening_laws'] = list(set(limb_darkening_laws))
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

    #Checking that the requested passbands exists in the Passbands folder.
    #The user can add any passband file with the right format. There is no check, at this stage, on the content of the passband files.
    #TO DO: This is not actually the way how the code is currently written. It must be changed.
    passbands = input_dict_local['passbands']
    if not check_length(passbands):
        print('ERROR: invalid length=', len(passbands), 'for passbands. It must have length>=1.')
        check = False
    #allowed_passbands = ['WFC3_G141', 'Kepler'] #user-defined should be added, maybe also rectangle
    list_passbands = [f for f in os.listdir('Passbands/') if f.endswith('.pass')]
    allowed_passbands = [f[:-5] for f in list_passbands]
    for item in passbands:
        if item not in allowed_passbands:
            print('ERROR:',item,'is not a valid passbands.')
            check = False
        if passbands.count(item)>1:
            print('WARNING:', item, 'entered multiple times in passbands. Repetitions are ignored.')
        input_dict_local['passbands'] = list(set(passbands))

    #Checking the input wavelength_bins_files to split the selected passbands.
    #wavelength_bins_files is an optional keyword; by default the passbands are not split.
    #If wavelength_bins_files is activated, the corresponding line in the input file must contain exactly one word for each passband (in the same order).
    #These words can be the names of the files containing the wavelength bins or the string 'no_bins'.
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

    #Checking the requested output_path (optional) where to store the output files.
    #TO DO: Is it a valid path? This should also be checked.
    if 'output_path' in input_keys:
        output_path = input_dict_local['output_path']
        if not check_length(output_path, max_length=1):
            print('ERROR: invalid length=', len(output_path), 'for output_path (optional). It must have length=1.')
            check = False
        if not check_type_in_list(output_path, str):
            print('ERROR: invalid type for output_path. It must be string.')
            check = False
            
    ##The following keywords are specific to one of the alternative calculation_type
    mandatory_keys_individual = ['star_effective_temperature']
    allowed_keys_individual = ['target_names', 'star_effective_temperature', 'star_log_gravity', 'star_metallicity']
    mandatory_keys_grid = []
    allowed_keys_grid = ['star_minimum_effective_temperature', 'star_maximum_effective_temperature', 'star_minimum_log_gravity', 'star_maximum_log_gravity', 'star_minimum_metallicity', 'star_maximum_metallicity']

    #Checking that all the mandatory keywords are obtained from the input file, if the calculation_type is 'individual'.
    #Also removing those keywords that will not be used (WARNING messages will appear). 
    if calculation_type == 'individual':
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

    #Checking that all the mandatory keywords are obtained from the input file, if the calculation_type is 'grid'.
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

    ##Checking the keywords which are specific to one of the alternative calculation_type.
    ##NOTE: There are no checks here that the number of values for the various parameters are consistent. These are checked in the functions get_individual_parameters and get_grid_parameters.
    ##NOTE (part2): This can be simplified if I eliminate the options to add these info from file, which is also not very practical in this way.
    
    #check star_effective_temperature
    if 'star_effective_temperature' in input_keys:
        star_effective_temperature = input_dict_local['star_effective_temperature']
        if not check_length(star_effective_temperature):
            print('ERROR: invalid length=', len(star_effective_temperature), 'for star_effective_temperature. It must have length>=1.')
            check = False
        if len(star_effective_temperature)>1:
            if not check_type_in_list(star_effective_temperature, float):
                print('ERROR: invalid type for item in star_effective_temperature. Items must be float or a single file to read.')
                check = False

    #check star_log_gravity
    if 'star_log_gravity' in input_keys:
        star_log_gravity = input_dict_local['star_log_gravity']
        if not check_length(star_log_gravity):
            print('ERROR: invalid length=', len(star_log_gravity), 'for star_log_gravity (optional). It must have length>=1.')
            check = False
        if len(star_log_gravity)>1:
            if not check_type_in_list(star_log_gravity, float):
                print('ERROR: invalid type for item in star_log_gravity (optional). Items must be float or a single file to read.')
                check = False

    #check star_metallicity
    if 'star_metallicity' in input_keys:
        star_metallicity = input_dict_local['star_metallicity']
        if not check_length(star_metallicity):
            print('ERROR: invalid length=', len(star_metallicity), 'for star_metallicity (optional). It must have length>=1.')
            check = False
        if len(star_metallicity)>1:
            if not check_type_in_list(star_metallicity, float):
                print('ERROR: invalid type for item in star_metallicity (optional). Items must be float or a single file to read.')
                check = False

    #check target_names
    if 'target_names' in input_keys:
        target_names = input_dict_local['target_names']
        if not check_length(target_names):
            print('ERROR: invalid length=', len(target_names), 'for target_names (optional). It must have length>=1.')
            check = False

    #check star_minimum_effective_temperature
    if 'star_minimum_effective_temperature' in input_keys:
        star_minimum_effective_temperature = input_dict_local['star_minimum_effective_temperature']
        if not check_length(star_minimum_effective_temperature, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_effective_temperature), 'for star_minimum_effective_temperature. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_effective_temperature, float):
                print('ERROR: invalid type for star_minimum_effective_temperature. It must be float.')
                check = False

    #check star_maximum_effective_temperature
    if 'star_maximum_effective_temperature' in input_keys:
        star_maximum_effective_temperature = input_dict_local['star_maximum_effective_temperature']
        if not check_length(star_maximum_effective_temperature, max_length=1):
            print('ERROR: invalid length=', len(star_maximum_effective_temperature), 'for star_maximum_effective_temperature. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_maximum_effective_temperature, float):
                print('ERROR: invalid type for star_maximum_effective_temperature. It must be float.')
                check = False

    #check star_minimum_log_gravity
    if 'star_minimum_log_gravity' in input_keys:
        star_minimum_log_gravity = input_dict_local['star_minimum_log_gravity']
        if not check_length(star_minimum_log_gravity, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_log_gravity), 'for star_minimum_log_gravity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_log_gravity, float):
                print('ERROR: invalid type for star_minimum_log_gravity. It must be float.')
                check = False

    #check star_maximum_log_gravity
    if 'star_maximum_log_gravity' in input_keys:
        star_maximum_log_gravity = input_dict_local['star_maximum_log_gravity']
        if not check_length(star_maximum_log_gravity, max_length=1):
            print('ERROR: invalid length=', len(star_maximum_log_gravity), 'for star_maximum_log_gravity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_maximum_log_gravity, float):
                print('ERROR: invalid type for star_maximum_log_gravity. It must be float.')
                check = False

    #check star_minimum_metallicity
    if 'star_minimum_metallicity' in input_keys:
        star_minimum_metallicity = input_dict_local['star_minimum_metallicity']
        if not check_length(star_minimum_metallicity, max_length=1):
            print('ERROR: invalid length=', len(star_minimum_metallicity), 'for star_minimum_metallicity. It must be length=1.')
            check = False
        else:
            if not check_type_in_list(star_minimum_metallicity, float):
                print('ERROR: invalid type for star_minimum_metallicity. It must be float.')
                check = False

    #check star_maximum_metallicity
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

    return check




def get_individual_parameters(input_dict):
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    input_star_effective_temperature = input_dict_local['star_effective_temperature']
    if check_type_in_list(input_star_effective_temperature, str):
        [star_effective_temperature, check] = read_as_numpy_array(input_star_effective_temperature[0])
        if not check:
            exit()
        elif star_effective_temperature.ndim>1:
            print('ERROR: invalid shape for star_effective_temperature read from file. It can be up to 1D vector.')
            exit()
    elif check_type_in_list(input_star_effective_temperature, float):
        star_effective_temperature = input_star_effective_temperature[:]
    n_targets = len(star_effective_temperature)
    #
    if 'star_log_gravity' in input_keys:
        input_star_log_gravity = input_dict_local['star_log_gravity']
        if check_type_in_list(input_star_log_gravity, str):
            [star_log_gravity, check] = read_as_numpy_array(input_star_log_gravity[0])
            if not check:
                exit()
            elif star_log_gravity.ndim>1:
                print('ERROR: invalid shape for star_log_gravity read from file. It can be up to 1D vector.')
                exit()
        elif check_type_in_list(input_star_log_gravity, float):
            star_log_gravity = input_star_log_gravity[:]
        if len(input_star_log_gravity)==1 and n_targets>1:
            print('WARNING: assuming star_log_gravity=', input_star_log_gravity[0], 'for all targets. A single input value was provided.')
            star_log_gravity = input_star_log_gravity*n_targets
        elif len(input_star_log_gravity)>1 and len(input_star_log_gravity)!=n_targets:
            print('ERROR: invalid length for star_log_gravity (optional). It must be equal to the length of star_effective_temperature or a single number.')
            exit()
    elif 'star_log_gravity' not in input_keys:
        print('WARNING: assuming star_log_gravity=4.5 for all targets. No input values were provided.')
        star_log_gravity = [4.5]*n_targets
    #
    if 'star_metallicity' in input_keys:
        input_star_metallicity = input_dict_local['star_metallicity']
        if check_type_in_list(input_star_metallicity, str):
            [star_metallicity, check] = read_as_numpy_array(input_star_metallicity[0])
            if not check:
                exit()
            elif star_metallicity.ndim>1:
                print('ERROR: invalid shape for star_metallicity read from file. It can be up to 1D vector.')
                exit()
        elif check_type_in_list(input_star_metallicity, float):
            star_metallicity = input_star_metallicity[:]
        if len(input_star_metallicity)==1 and n_targets>1:
            print('WARNING: assuming star_metallicity=', input_star_metallicity[0], 'for all targets. A single input value was provided.')
            star_metallicity = input_star_metallicity*n_targets
        elif len(input_star_metallicity)>1 and len(input_star_metallicity)!=n_targets:
            print('ERROR: invalid length for star_metallicity (optional). It must be equal to the length of star_effective_temperature or a single number.')
            exit()
    elif 'star_metallicity' not in input_keys:
        print('WARNING: assuming star_metallicity=0.0 for all targets. No input values were provided.')
        star_metallicity = [0.0]*n_targets
    #
    if 'target_names' in input_keys:
        input_target_names = input_dict_local['target_names']
        if len(input_target_names)==1:
            try:
                with open(input_target_names[0]) as f:
                    target_names = f.read().split()
            except IOError:
                target_names = input_target_names[:]
        elif len(input_target_names)>1:
            target_names = input_target_names[:]
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
    params = file_name.split('/')[-1]
    params = file_name.replace('.pickle', '').split('_')
    teff = float(params[0].replace('teff', ''))
    logg = float(params[1].replace('logg', ''))
    mh = float(params[2].replace('MH', ''))
    return np.array((teff, logg, mh))


def get_grid_parameters(stellar_models_grid):
    files = databases[stellar_models_grid].get_filename_list()
    star_params_grid = np.zeros((len(files),3))
    for i in range(len(files)):
        star_params_grid[i,:] = stellar_params_from_file_name(files[i])
    return files, star_params_grid


def get_subgrid(input_dict):
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
    if stellar_models_grid in ['Phoenix_2018', 'Phoenix_2012_13']:
        indices_teff_sup = np.where(star_params_grid[:,0]>=star_effective_temperature)[0]
        indices_teff_inf = np.where(star_params_grid[:,0]<=star_effective_temperature)[0]
        teff_upper = np.min(star_params_grid[indices_teff_sup,0])
        teff_lower = np.max(star_params_grid[indices_teff_inf,0])
        indices_teff_upper = np.where(star_params_grid[:,0]==teff_upper)[0]
        indices_teff_lower = np.where(star_params_grid[:,0]==teff_lower)[0]
        indices_logg_sup = np.where(star_params_grid[:,1]>=star_log_gravity)[0]
        indices_logg_inf = np.where(star_params_grid[:,1]<=star_log_gravity)[0]
        indices_teff_upper_logg_sup = np.intersect1d( indices_teff_upper, indices_logg_sup )
        indices_teff_upper_logg_inf = np.intersect1d( indices_teff_upper, indices_logg_inf )
        indices_teff_lower_logg_sup = np.intersect1d( indices_teff_lower, indices_logg_sup )
        indices_teff_lower_logg_inf = np.intersect1d( indices_teff_lower, indices_logg_inf )
        logg_upper1 = np.min(star_params_grid[indices_teff_upper_logg_sup,1])
        logg_lower1 = np.max(star_params_grid[indices_teff_upper_logg_inf,1])
        logg_upper2 = np.min(star_params_grid[indices_teff_lower_logg_sup,1])
        logg_lower2 = np.max(star_params_grid[indices_teff_lower_logg_inf,1])
        indices_logg_upper1 = np.where(star_params_grid[:,1]==logg_upper1)[0]
        indices_logg_lower1 = np.where(star_params_grid[:,1]==logg_lower1)[0]
        indices_logg_upper2 = np.where(star_params_grid[:,1]==logg_upper2)[0]
        indices_logg_lower2 = np.where(star_params_grid[:,1]==logg_lower2)[0]
        index_teff_upper_logg_upper = np.intersect1d( indices_teff_upper_logg_sup, indices_logg_upper1 )[0]
        index_teff_upper_logg_lower = np.intersect1d( indices_teff_upper_logg_inf, indices_logg_lower1 )[0]
        index_teff_lower_logg_upper = np.intersect1d( indices_teff_lower_logg_sup, indices_logg_upper2 )[0]
        index_teff_lower_logg_lower = np.intersect1d( indices_teff_lower_logg_inf, indices_logg_lower2 )[0]
        neigh_indices = np.array([index_teff_upper_logg_upper, index_teff_upper_logg_lower, index_teff_lower_logg_upper, index_teff_lower_logg_lower])
        return neigh_indices
            

def check_passband_limits(pb, stellar_models_grid):
    if stellar_models_grid == 'Phoenix_2018':
        minimum_wavelength = 500.0
        maximum_wavelength = 25999.0
        if np.min(pb[:,0])<minimum_wavelength or np.max(pb[:,0])>maximum_wavelength:
            pb = False
        return pb
    elif stellar_models_grid == 'Phoenix_2012_13':
        minimum_wavelength = 2500.0
        maximum_wavelength = 99995.0
        if np.min(pb[:,0])<minimum_wavelength or np.max(pb[:,0])>maximum_wavelength:
            pb = False
        return pb


def get_passband(passband, stellar_models_grid):
    ext = '.pass'
    [pb, check] = read_as_numpy_array('Passbands/'+passband+'.pass')
    if not check:
        print('WARNING: Skipping passband', passband, '.')
        return np.array([]), check
    pb = np.atleast_2d(pb)
    if not check_2Darray(pb, n_col=2):
        print('WARNING: invalid format for', passband, 'passband file. It must have 2 columns. Skipping', passband, 'passband.')
        check = False
        return np.array([]), check
    if np.min(pb)<0:
        print('WARNING: negative value found in file', passband, 'passband file. Skipping', passband, 'passband.')
        check = False
        return np.array([]), check
    pb = check_passband_limits(pb, stellar_models_grid)
    if not check:
        print('WARNING:', passband, 'passband exceeds wavelength range for the', stellar_models_grid, 'stellar model grid. Skipping', passband, 'passband.')
        return np.array([]), check
    return pb, check



def get_wavelength_bins(wavelength_bins_file, pb, passband):
    check =  True
    if wavelength_bins_file == 'no_bins':
        check = True
        return np.array([]), check
    path = 'Wavelength_bins_files/'
    #ext = '.txt'
    [wb, check] = read_as_numpy_array(path+wavelength_bins_file)
    if not check:
        print('WARNING: Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        return np.array([]), check
    wb = np.atleast_2d(wb)
    if not check_2Darray(wb, n_col=2):
        print('WARNING: invalid format for wavelength_bins_file', wavelength_bins_file, '. It must have 2 columns. Ignoring wavelength_bins_file', wavelength_bins_file, '.')
        check = False
        return np.array([]), check
    if (wb[:,0]>=wb[:,1]).any():
        print('WARNING: invalid line in wavelength_bins_file', wavelength_bins_file, '. The lower limit cannot be greater or equal to the upper limit.')
        check = False
        return np.array([]), check
    if np.min(wb)<np.min(pb[:,0]) or np.max(wb)>np.max(pb[:,0]):
        print('WARNING:wavelength_bins_file', wavelength_bins_file, 'exceeds wavelength range for the ', passband, 'passband. Ignoring this passband.')
        check = False
        return np.array([]), check
    return wb, check




def get_response(passband, model_wavelengths):
    resp = np.zeros(len(model_wavelengths))
    passband_wl = np.array((np.min(passband[:,0]),np.max(passband[:,0])))
    for i in range(len(model_wavelengths)):
        if model_wavelengths[i] >= passband_wl[0] and model_wavelengths[i] <= passband_wl[1]:
            prov = np.argmin(np.abs(passband[:,0]-model_wavelengths[i]))
            if passband[prov,0] == model_wavelengths[i]:
                resp[i] = passband[prov,1]
            elif passband[prov,0] < model_wavelengths[i]:
                resp[i] = ( passband[prov,1]*(passband[prov+1,0]-model_wavelengths[i]) + passband[prov+1,1]*(model_wavelengths[i]-passband[prov,0]) ) / (passband[prov+1,0]-passband[prov,0])
            elif passband[prov,0] > model_wavelengths[i]:
                resp[i] = ( passband[prov-1,1]*(passband[prov,0]-model_wavelengths[i]) + passband[prov,1]*(model_wavelengths[i]-passband[prov-1,0]) ) / (passband[prov,0]-passband[prov-1,0])
    return resp


def compute_integrated_intensities(model_wavelengths, model_intensities, response):
    integ_ints = np.zeros(np.shape(model_intensities)[1])
    #model_intensities = model_intensities/np.mean(model_intensities[:,-1]) #arbitrary normalization
    for i in range(1,len(model_wavelengths)-1):
        integ_ints += (model_wavelengths[i+1]-model_wavelengths[i-1])*response[i]*model_wavelengths[i]*model_intensities[i,:]
    integ_ints = integ_ints/integ_ints[-1]
    return integ_ints


def get_integrated_intensities(model_dict, passbands_dict, wavelength_bins_dict):
    model_wavelengths = model_dict['wavelengths']
    model_intensities = model_dict['intensities']
    mu = model_dict['mu']
    passbands = list(passbands_dict.keys())
    integ_dict = {}
    integ_dict['mu'] = mu
    for passband in passbands:
        response = get_response(passbands_dict[passband], model_wavelengths)
        wavelength_bins = wavelength_bins_dict[passband]
        if len(wavelength_bins)==0:
            integ_dict[passband] = compute_integrated_intensities(model_wavelengths, model_intensities, response)
        else:
            for wbin in wavelength_bins:
                wbin_out_indices = np.concatenate(( np.where(model_wavelengths<wbin[0])[0], np.where(model_wavelengths>wbin[1])[0] ))
                wbin_response = response.copy()
                wbin_response[wbin_out_indices] = 0.0
                integ_dict[passband+'_'+str(wbin[0])+'_'+str(wbin[1])] = compute_integrated_intensities(model_wavelengths, model_intensities, wbin_response)
    return integ_dict



def rescale_and_weights(mu, intensities):
    radi = np.sqrt(1.0-mu**2.0)
    dint_dr = np.abs( (intensities[1:]-intensities[:-1])/(radi[1:]-radi[:-1]) )
    rmax_dr = 0.5*( radi[np.argmax(dint_dr)+1] + radi[np.argmax(dint_dr)] )
    radicut_dr = radi[np.where(radi<=rmax_dr)]/rmax_dr
    mucut_dr = np.sqrt(1.0-radicut_dr**2.0)
    intensitiescut_dr = intensities[np.where(radi<=rmax_dr)]
    weights_dr = np.zeros_like(radicut_dr)
    weights_dr[1:-1] = -0.5*(radicut_dr[2:]-radicut_dr[:-2])
    weights_dr[0] = (1.0-radicut_dr[0]) - 0.5*(radicut_dr[1]-radicut_dr[0])
    weights_dr[-1] = 0.5*radicut_dr[-2]
    return mucut_dr, intensitiescut_dr, weights_dr



def get_limb_darkening_coefficients(integ_dict, limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders):
    mu = integ_dict['mu']
    passbands = [key for key in list(integ_dict.keys()) if key!='mu']
    ldc_dict = {}
    ldc_dict['passbands'] = {}
    for passband in passbands:
        integ_ints = integ_dict[passband]
        ldc_dict['passbands'][passband] = {}
        if stellar_models_grid in ['Phoenix_2018', 'Phoenix_2012_13']:
            [res_mu, res_integ_ints, weights] = rescale_and_weights(mu, integ_ints)
            ldc_dict['passbands'][passband]['rescaled_mu'] = res_mu
            ldc_dict['passbands'][passband]['rescaled_intensities'] = res_integ_ints
            ldc_dict['passbands'][passband]['weights'] = weights
            ldc_dict['passbands'][passband]['laws'] = {}
        if 'claret4' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['claret4'] = {}
            conv = False
            params = np.array((0.9, -0.5, 0.9, -0.5))
            while conv==False:
                ldc_claret4 = minimize(claret4, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_claret4.success
                params = ldc_claret4.x
            ldc_dict['passbands'][passband]['laws']['claret4']['coefficients'] = ldc_claret4.x
            ldc_dict['passbands'][passband]['laws']['claret4']['weighted_rms_res'] = np.sqrt(ldc_claret4.fun)
        if 'power2' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['power2'] = {}
            conv = False
            params = np.array((0.5, 1.0))
            while conv==False:
                ldc_power2 = minimize(power2, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_power2.success
                params = ldc_power2.x
            ldc_dict['passbands'][passband]['laws']['power2']['coefficients'] = ldc_power2.x
            ldc_dict['passbands'][passband]['laws']['power2']['weighted_rms_res'] = np.sqrt(ldc_power2.fun)
        if 'square_root' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['square_root'] = {}
            conv = False
            params = np.array((0.9, -0.5))
            while conv==False:
                ldc_square_root = minimize(square_root, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_square_root.success
                params = ldc_square_root.x
            ldc_dict['passbands'][passband]['laws']['square_root']['coefficients'] = ldc_square_root.x
            ldc_dict['passbands'][passband]['laws']['square_root']['weighted_rms_res'] = np.sqrt(ldc_square_root.fun)
        if 'quadratic' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['quadratic'] = {}
            conv = False
            params = np.array((0.9, -0.5))
            while conv==False:
                ldc_quadratic = minimize(quadratic, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_quadratic.success
                params = ldc_quadratic.x
            ldc_dict['passbands'][passband]['laws']['quadratic']['coefficients'] = ldc_quadratic.x
            ldc_dict['passbands'][passband]['laws']['quadratic']['weighted_rms_res'] = np.sqrt(ldc_quadratic.fun)
        if 'linear' in limb_darkening_laws:
            ldc_dict['passbands'][passband]['laws']['linear'] = {}
            conv = False
            params = np.array((0.5))
            while conv==False:
                ldc_linear = minimize(linear, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                conv = ldc_linear.success
                params = ldc_linear.x
            ldc_dict['passbands'][passband]['laws']['linear']['coefficients'] = ldc_linear.x
            ldc_dict['passbands'][passband]['laws']['linear']['weighted_rms_res'] = np.sqrt(ldc_linear.fun)
        if 'gen_poly' in limb_darkening_laws:
            for n in gen_poly_orders:
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)] = {}
                conv = False
                params = np.ones(n)/n
                while conv==False:
                    ldc_gen_poly = minimize(gen_poly, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                    conv = ldc_gen_poly.success
                    params = ldc_gen_poly.x
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)]['coefficients'] = ldc_gen_poly.x
                ldc_dict['passbands'][passband]['laws']['gen_poly'+str(n)]['weighted_rms_res'] = np.sqrt(ldc_gen_poly.fun)
        if 'gen_claret' in limb_darkening_laws:
            for n in gen_claret_orders:
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)] = {}
                conv = False
                params = np.ones(n)/n
                while conv==False:
                    ldc_gen_claret = minimize(gen_claret, params, args=(res_mu, res_integ_ints, weights), method='nelder-mead', options={'xtol':1e-8,'maxfev':10000,'disp':False})
                    conv = ldc_gen_claret.success
                    params = ldc_gen_claret.x
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)]['coefficients'] = ldc_gen_claret.x
                ldc_dict['passbands'][passband]['laws']['gen_claret'+str(n)]['weighted_rms_res'] = np.sqrt(ldc_gen_claret.fun)
    return ldc_dict
        

def interp_ldc(teff, logg, mh, neigh_indices, passband, law, neighbour_limb_darkening_coefficients):
#First interpolate in logg between the two models with the upper temperature value
    if neigh_indices[0]==neigh_indices[1]:
        coeffs_01 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_01 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['passbands'][passband]['laws'][law]['weighted_rms_res']
    else:
        coeffs_0 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_0 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['passbands'][passband]['laws'][law]['weighted_rms_res']
        logg_0 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['star_params'][1]
        coeffs_1 = neighbour_limb_darkening_coefficients[neigh_indices[1]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_1 = neighbour_limb_darkening_coefficients[neigh_indices[1]]['passbands'][passband]['laws'][law]['weighted_rms_res']
        logg_1 = neighbour_limb_darkening_coefficients[neigh_indices[1]]['star_params'][1]
        w0 = (logg_0 - logg)/(logg_0 - logg_1)
        w1 = (logg - logg_1)/(logg_0 - logg_1)
        coeffs_01 = w0*coeffs_0 + w1*coeffs_1
        w_res_01 = w0*w_res_0 + w1*w_res_1
#Then interpolate in logg between the two models with the lower temperature value
    if neigh_indices[2]==neigh_indices[3]:
        coeffs_23 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_23 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['passbands'][passband]['laws'][law]['weighted_rms_res']
    else:
        coeffs_2 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_2 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['passbands'][passband]['laws'][law]['weighted_rms_res']
        logg_2 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['star_params'][1]
        coeffs_3 = neighbour_limb_darkening_coefficients[neigh_indices[3]]['passbands'][passband]['laws'][law]['coefficients']
        w_res_3 = neighbour_limb_darkening_coefficients[neigh_indices[3]]['passbands'][passband]['laws'][law]['weighted_rms_res']
        logg_3 = neighbour_limb_darkening_coefficients[neigh_indices[3]]['star_params'][1]
        w2 = (logg_2 - logg)/(logg_2 - logg_3)
        w3 = (logg - logg_3)/(logg_2 - logg_3)
        coeffs_23 = w2*coeffs_2 + w3*coeffs_3
        w_res_23 = w2*w_res_2 + w3*w_res_3
#Finally interpolate in teff between the two interpolated models with the upper and lower temperatures
    if neigh_indices[0]==neigh_indices[2]:
        coeffs = coeffs_01*1.0
        w_res = w_res_01*1.0
    else:
        teff_0 = neighbour_limb_darkening_coefficients[neigh_indices[0]]['star_params'][0]
        teff_2 = neighbour_limb_darkening_coefficients[neigh_indices[2]]['star_params'][0]
        w01 = (teff_0 - teff)/(teff_0 - teff_2)
        w23 = (teff - teff_2)/(teff_0 - teff_2)
        coeffs = w01*coeffs_01 + w23*coeffs_23
        w_res = w01*w_res_01 + w23*w_res_23
    return coeffs, w_res

            



def process_configuration(input_dict):
    input_dict_local = copy.deepcopy(input_dict)
    input_keys = list(input_dict_local.keys())
    calculation_type = input_dict_local['calculation_type'][0]
    if 'user_output' in input_keys:
        user_output = input_dict_local['user_output'][0]
    else:
        user_output = 'basic'
    if 'output_path' in input_keys:
        output_path = input_dict_local['output_path'][0]
    else:
        output_path = 'Results/'
    stellar_models_grid = input_dict_local['stellar_models_grid'][0]
    limb_darkening_laws = input_dict_local['limb_darkening_laws']
    gen_poly_orders = None
    gen_claret_orders = None
    if 'gen_poly' in limb_darkening_laws:
        gen_poly_orders = input_dict_local['gen_poly_orders']
        gen_poly_orders = np.asarray(gen_poly_orders, dtype=int)
    if 'gen_claret' in limb_darkening_laws:
        gen_claret_orders = input_dict_local['gen_claret_orders']
        gen_claret_orders = np.asarray(gen_claret_orders, dtype=int)
    passbands = input_dict_local['passbands']
    n_pass = len(passbands)
    if 'wavelength_bins_files' in input_keys:
        wavelength_bins_files = input_dict_local['wavelength_bins_files']
    else:
        wavelength_bins_files = ['no_bins',]*n_pass
    passbands_dict = {}
    wavelength_bins_dict = {}
    for i in range(n_pass):
	check = True
        [pb, check] = get_passband(passbands[i], stellar_models_grid)
        if check:
            [wb, check] = get_wavelength_bins(wavelength_bins_files[i], pb, passbands[i])
            if check:
                passbands_dict[passbands[i]] = pb
                wavelength_bins_dict[passbands[i]] = wb
    [star_files_grid, star_params_grid] = get_grid_parameters(stellar_models_grid)
    if calculation_type=='individual':
        [star_effective_temperature, star_log_gravity, star_metallicity, target_names] = get_individual_parameters(input_dict_local)
        n_targets = len(star_effective_temperature)
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
        neighbour_files_indices_list = np.unique(neighbour_files_indices)
        neighbour_model_intensities = {}
        neighbour_integrated_intensities = {}
        neighbour_limb_darkening_coefficients = {}
        for i in neighbour_files_indices_list:
            neighbour_model_intensities[i] = databases[stellar_models_grid].get_file_content(dbx_file=star_files_grid[i])
            neighbour_integrated_intensities[i] = get_integrated_intensities(neighbour_model_intensities[i], passbands_dict, wavelength_bins_dict)
            neighbour_limb_darkening_coefficients[i] = get_limb_darkening_coefficients(neighbour_integrated_intensities[i], limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders)
            neighbour_limb_darkening_coefficients[i]['star_params'] = star_params_grid[i]
        #interp_ldc = get_interp_ldc(neighbour_limb_darkening_coefficients, neighbour_files_indices)
        target_limb_darkening_coefficients = {}
        i0 = neighbour_files_indices_list[0]
        final_passbands = list(neighbour_limb_darkening_coefficients[i0]['passbands'].keys())
        final_limb_darkening_laws = list(neighbour_limb_darkening_coefficients[i0]['passbands'][final_passbands[0]]['laws'].keys())
        for i in range(n_targets):
            target_limb_darkening_coefficients[target_names[i]] = {}
            target_limb_darkening_coefficients[target_names[i]]['star_params'] = np.array([star_effective_temperature[i], star_log_gravity[i], star_metallicity[i]])
            target_limb_darkening_coefficients[target_names[i]]['passbands'] = {}
            for passband in final_passbands:
                target_limb_darkening_coefficients[target_names[i]]['passbands'][passband] = {}
                for law in final_limb_darkening_laws:
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law] = {}
                    [coeffs, w_res] = interp_ldc(star_effective_temperature[i], star_log_gravity[i], star_metallicity[i], neighbour_files_indices[i,:], passband, law, neighbour_limb_darkening_coefficients)
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law]['coefficients'] = coeffs
                    target_limb_darkening_coefficients[target_names[i]]['passbands'][passband][law]['weighted_rms_res'] = w_res
        limb_darkening_coefficients = {}
        limb_darkening_coefficients['neighbour'] = neighbour_limb_darkening_coefficients
        limb_darkening_coefficients['target'] = target_limb_darkening_coefficients
	for i in range(n_targets):
            with open(os.path.join(output_path, target_names[i]+'_ldc.pickle') , 'wb') as out1:
                pickle.dump(target_limb_darkening_coefficients[target_names[i]], out1, protocol=pickle.HIGHEST_PROTOCOL)
        if user_output == 'complete':
            with open(os.path.join(output_path, 'neighbour_ldc.pickle') , 'wb') as out2:
                pickle.dump(neighbour_limb_darkening_coefficients, out2, protocol=pickle.HIGHEST_PROTOCOL)
            with open(os.path.join(output_path, 'neighbour_intensities.pickle') , 'wb') as out3:
                pickle.dump(neighbour_integrated_intensities, out3, protocol=pickle.HIGHEST_PROTOCOL)
        return limb_darkening_coefficients, neighbour_integrated_intensities
    elif calculation_type=='grid':
        [star_files_subgrid, star_params_subgrid, subgrid_indices] = get_subgrid(input_dict_local)
        n_files = len(star_files_subgrid)
        model_intensities = {}
        integrated_intensities = {}
        limb_darkening_coefficients = {}
        for i in range(n_files):
            label = 'teff' + str(star_params_subgrid[i,0]) + '_logg' + str(star_params_subgrid[i,1]) + '_MH' + str(star_params_subgrid[i,2])
            model_intensities[label] = databases[stellar_models_grid].get_file_content(dbx_file=star_files_grid[subgrid_indices[i]])
            integrated_intensities[label] = {}
            integrated_intensities[label] = get_integrated_intensities(model_intensities[label], passbands_dict, wavelength_bins_dict)
            limb_darkening_coefficients[label] = get_limb_darkening_coefficients(integrated_intensities[label], limb_darkening_laws, stellar_models_grid, gen_poly_orders, gen_claret_orders)
            limb_darkening_coefficients[label]['star_params'] = star_params_subgrid[i]
        with open(os.path.join(output_path, 'grid_ldc.pickle'), 'wb') as out1:
            pickle.dump(limb_darkening_coefficients, out1, protocol=pickle.HIGHEST_PROTOCOL)
        if user_output == 'complete':
            with open(os.path.join(output_path, 'grid_intensities.pickle'), 'wb') as out2:
                pickle.dump(integrated_intensities, out2, protocol=pickle.HIGHEST_PROTOCOL)
        return limb_darkening_coefficients, integrated_intensities


def limbcalc(parameters_file):
    input_dict = read_configuration(parameters_file)
    if check_configuration(input_dict):
        process_configuration(input_dict)
    else:
        print('Something went wrong with the input.')
