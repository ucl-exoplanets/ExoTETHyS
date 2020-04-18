import pytest
from exotethys import sail
import numpy as np
import astropy.units as u


testdata = [
    ('ciao', 'ciao'),
    ('5*5', '5*5'),
    ('25', 25.0),
    ('-3.14', -3.14),
]

@pytest.mark.parametrize('toread,reference', testdata)
def test_str2float(toread, reference):
    result = sail.str2float(toread)
    np.testing.assert_equal(result, reference)


#############################


testdata = [
    (np.array([]), np.arange(5), np.arange(5)),
    (np.arange(5), np.array([]), np.arange(5)),
    (np.arange(5), np.arange(5)**2.0, np.array([[0, 1, 2, 3, 4], [0, 1, 4, 9, 16]])),
    (np.ones((2,3)), np.arange(3), np.array([[1, 1, 1], [1, 1, 1], [0, 1, 2]])),
]

@pytest.mark.parametrize('array1,array2,reference', testdata)
def test_my_vstack(array1, array2, reference):
    result = sail.my_vstack(array1,array2)
    np.testing.assert_equal(result, reference)

def test_my_vstack_fails():
    array1 = np.ones((2,3))
    array2 = np.arange(2)
    with pytest.raises(ValueError):
        result = sail.my_vstack(array1,array2)


#############################


testdata = [
    (np.array([]), {},  False),
    (np.array([]), {'min_length': 0}, True),
    (np.array([1]), {'max_length': 1}, True),
    (np.array([1, 2]), {'max_length': 1}, False),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_length(vector, args, reference):
    result = sail.check_length(vector, **args)
    np.testing.assert_equal(result, reference)



#############################

testdata = [
    ([1, '2', 3], float, False),
    ([1, '2', 3], str, False),
    (['1', '2', '3'], str, True),
    (['1', '2', '3'], float, False),
    ([1., 2., 3.], float, True),
]

@pytest.mark.parametrize('vector,item_type,reference', testdata)
def test_check_type_in_list(vector, item_type, reference):
    result = sail.check_type_in_list(vector, item_type)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    (np.arange(6), {},  False),
    (np.arange(6), {'min_value': 0}, True),
    (np.arange(6)+1, {}, True),
    (np.array([1, 2, 3.5]), {}, False),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_integers(vector, args, reference):
    result = sail.check_integers(vector, **args)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    (np.arange(6), {},  False),
    (np.atleast_2d(np.arange(6)), {}, True),
    (np.array([[1, 2, 3], [4, 5, 6]]), {}, True),
    (np.array([[1, 2, 3], [4, 5, 6]]), {'n_col': 3}, True),
    (np.array([[1, 2, 3], [4, 5, 6]]), {'n_col': 2}, False),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_2Darray(vector, args, reference):
    result = sail.check_2Darray(vector, **args)
    np.testing.assert_equal(result, reference)


#############################



testdata = [
    ('teff07500_logg1.5_MH0.0.pickle', (7500, 1.5, 0.0)),
    ('teff18000_logg4.5_MH-0.5.pickle', (18000, 4.5, -0.5)),
]

@pytest.mark.parametrize('filename,reference', testdata)
def test_stellar_params_from_file_name(filename, reference):
    result = sail.stellar_params_from_file_name(filename)
    np.testing.assert_equal(result, reference)


def test_stellar_params_from_file_name_fails():
    filename = 'teff07500__logg1.5_MH0.0.pickle'
    with pytest.raises(ValueError):
        result = sail.stellar_params_from_file_name(filename)
    filename = 'test_data/teff07500_logg1.5_MH0.0.pickle'
    result = sail.stellar_params_from_file_name(filename)
    np.testing.assert_equal(result, (7500, 1.5, 0.0))


#############################

testdata = [
    ('Atlas_2000', [3500.0, 50000.0, 0.0, 5.0, -5.0, 1.0]),
    ('Phoenix_2012_13', [3000.0, 10000.0, 0.0, 6.0, 0.0, 0.0]),
    ('Phoenix_2018', [2300.0, 12000.0, 0.0, 6.0, 0.0, 0.0]),
]

@pytest.mark.parametrize('stellar_models_grid,reference', testdata)
def test_get_grid_parameters(stellar_models_grid, reference):
    files, star_params_grid = sail.get_grid_parameters(stellar_models_grid)
    np.testing.assert_equal(len(files), np.shape(star_params_grid)[0])
    np.testing.assert_equal(np.shape(star_params_grid)[1], 3)
    limits = [np.min(star_params_grid[:,0]), np.max(star_params_grid[:,0]), np.min(star_params_grid[:,1]), np.max(star_params_grid[:,1]), np.min(star_params_grid[:,2]), np.max(star_params_grid[:,2])]
    np.testing.assert_equal(limits, reference)


#############################


testdata = [
    (np.arange(10000, 99996) * u.Angstrom, 'Phoenix_2012_13',  True),
    (np.arange(10000, 99996) * u.Angstrom, 'Phoenix_2018',  False),
    (np.arange(10000, 99996) * u.Angstrom, 'Atlas_2000',  True),
    (np.arange(500, 1000) * u.Angstrom, 'Phoenix_2012_13',  False),
    (np.arange(500, 1000) * u.Angstrom, 'Phoenix_2018',  True),
    (np.arange(500, 1000) * u.Angstrom, 'Atlas_2000',  True),
]

@pytest.mark.parametrize('vector,grid,reference', testdata)
def test_check_passband_limits(vector, grid, reference):
    result = sail.check_passband_limits(vector, grid)
    np.testing.assert_equal(result, reference)


#############################


testdata = [
    (1.0, 2000000.0, 10000.0, [1.0, 2000000.0, 145087]),
    (55000.0, 60000.0, 10000.0, [55000.0, 60000.0, 871]),
    (55000.0, 60000.0, 1.0, [55000.0, 60000.0, 10]),
]

@pytest.mark.parametrize('lambda1,lambda2,R,reference', testdata)
def test_get_waves_fromR(lambda1, lambda2, R, reference):
    vector = sail.get_waves_fromR(lambda1, lambda2, R)
    np.testing.assert_almost_equal([np.min(vector), np.max(vector), len(vector)], reference, decimal=8)


#############################



