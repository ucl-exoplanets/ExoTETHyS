import pytest
from exotethys import boats
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
    result = boats.str2float(toread)
    np.testing.assert_equal(result, reference)


#############################


testdata = [
    (np.array([]), {},  False),
    (np.array([]), {'min_length': 0}, True),
    (np.array([1]), {'max_length': 1}, True),
    (np.array([1, 2]), {'max_length': 1}, False),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_length(vector, args, reference):
    result = boats.check_length(vector, **args)
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
    result = boats.check_type_in_list(vector, item_type)
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
    result = boats.check_2Darray(vector, **args)
    np.testing.assert_equal(result, reference)


#############################



testdata = [
    ('teff07500_logg1.5_MH0.0.pickle', (7500, 1.5, 0.0)),
    ('teff18000_logg4.5_MH-0.5.pickle', (18000, 4.5, -0.5)),
]

@pytest.mark.parametrize('filename,reference', testdata)
def test_stellar_params_from_file_name(filename, reference):
    result = boats.stellar_params_from_file_name(filename)
    np.testing.assert_equal(result, reference)


def test_stellar_params_from_file_name_fails():
    filename = 'teff07500__logg1.5_MH0.0.pickle'
    with pytest.raises(ValueError):
        result = boats.stellar_params_from_file_name(filename)
    filename = 'test_data/teff07500_logg1.5_MH0.0.pickle'
    result = boats.stellar_params_from_file_name(filename)
    np.testing.assert_equal(result, (7500, 1.5, 0.0))


#############################

testdata = [
    ('Atlas_2000', [3500.0, 50000.0, 0.0, 5.0, -5.0, 1.0]),
    ('Phoenix_2012_13', [3000.0, 10000.0, 0.0, 6.0, 0.0, 0.0]),
    ('Phoenix_2018', [2300.0, 12000.0, 0.0, 6.0, 0.0, 0.0]),
    ('Phoenix_drift_2012', [1500.0, 3000.0, 2.5, 5.5, 0.0, 0.0]),
    ('Stagger_2015', [4000.0, 7000.0, 1.5, 5.0, -2.0, 0.5]),
]

@pytest.mark.parametrize('stellar_models_grid,reference', testdata)
def test_get_stellar_grid_parameters(stellar_models_grid, reference):
    files, star_params_grid = boats.get_stellar_grid_parameters(stellar_models_grid)
    np.testing.assert_equal(len(files), np.shape(star_params_grid)[0])
    np.testing.assert_equal(np.shape(star_params_grid)[1], 3)
    limits = [np.min(star_params_grid[:,0]), np.max(star_params_grid[:,0]), np.min(star_params_grid[:,1]), np.max(star_params_grid[:,1]), np.min(star_params_grid[:,2]), np.max(star_params_grid[:,2])]
    np.testing.assert_equal(limits, reference)


#############################


testdata = [
    (np.array([0.0, 0.5, 1.0, 1.5]), {'min_value': 0.0},  True),
    (np.array([0.0, 0.5, 1.0, 1.5]), {'min_value': 0.0, 'max_value': 1.0},  False),
    (np.array([0.0, 0.5, 1.0, 1.5]), {'inf_value': 0.0},  False),
    (np.array([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]), {},  True),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_values(vector, args, reference):
    result = boats.check_values(vector, **args)
    np.testing.assert_equal(result, reference)


#############################


testdata = [
    ('m', u.m, {},  True),
    ('m', u.m**2, {},  False),
    ('AU', u.m, {},  True),
    ('d', u.s, {'others': ['T_14']},  True),
    ('T_14', u.s, {'others': ['T_14']},  True),
    ('T_14', u.s, {'others': []},  False),
]

@pytest.mark.parametrize('item,ref_unit,args,reference', testdata)
def test_check_unit_type(item, ref_unit, args, reference):
    result = boats.check_unit_type(item, ref_unit, **args)
    np.testing.assert_equal(result, reference)



#############################

testdata = [
    (0.15 * u.dimensionless_unscaled, 9.0 * u.dimensionless_unscaled, 90 * u.deg, 3.5 * u.d, 3.425892118713202),
    (0.15 * u.dimensionless_unscaled, 9.0 * u.dimensionless_unscaled, 83 * u.deg, 3.5 * u.d, 1.034871978574793),
    (0.10000001 * u.dimensionless_unscaled, 2.2 * u.dimensionless_unscaled, 60 * u.deg, 1.5 * u.d, 0.0008920940699048776),
]

@pytest.mark.parametrize('rp_over_rs, sma_over_rs, inclination, period, reference', testdata)
def test_get_transit_duration_T14(rp_over_rs, sma_over_rs, inclination, period, reference):
    result = boats.get_transit_duration_T14(rp_over_rs, sma_over_rs, inclination, period)
    np.testing.assert_almost_equal(result.to(u.h).value, reference, decimal=8)



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
    result = boats.check_passband_limits(vector, grid)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    (1.0, 2000000.0, 10000.0, [1.0, 2000000.0, 145087]),
    (55000.0, 60000.0, 10000.0, [55000.0, 60000.0, 871]),
    (55000.0, 60000.0, 1.0, [55000.0, 60000.0, 10]),
]

@pytest.mark.parametrize('lambda1,lambda2,R,reference', testdata)
def test_get_waves_fromR(lambda1, lambda2, R, reference):
    vector = boats.get_waves_fromR(lambda1, lambda2, R)
    np.testing.assert_almost_equal([np.min(vector), np.max(vector), len(vector)], reference, decimal=8)


#############################

waves_test = np.array([1000.0, 10000.0, 50000.0, 100000.0, 500000.0]) * u.Angstrom
fluxes_test = np.array([144.1649012667559, 3741496.1001535463, 19456.223970992392, 1380.7778931393088, 2.4372444918063167]) * u.erg / (u.Angstrom * u.cm**2 * u.s)
fluxes_test_j = np.array([4.808823485037e+12, 1.2480287613351322e+19, 1.6224744362141688e+18, 4.605779285946243e+17, 2.032443134215135e+16]) * u.jansky
photons_Rsun_10pc_25m2 = np.array([9.222862600061873, 2393599.5618332108, 62235.00963957306, 8833.44328455026, 77.96062312384335]) * u.photon / (u.Angstrom * u.s)

testdata = [
    (waves_test, fluxes_test, 1.0 * u.Rsun, 10.0 * u.pc, 25.0 * u.m**2, photons_Rsun_10pc_25m2),
    (waves_test, fluxes_test_j, 1.0 * u.Rsun, 10.0 * u.pc, 25.0 * u.m**2, photons_Rsun_10pc_25m2),
]

@pytest.mark.parametrize('waves,fluxes,radius,distance,area,reference', testdata)
def test_get_photon_spectrum(waves, fluxes, radius, distance, area, reference):
    result = boats.get_photon_spectrum(waves, fluxes, radius, distance, area)
    np.testing.assert_almost_equal(result.value, reference.value, decimal=8)



#############################

testdata = [
    (0, 0.5, 0.5),
    (0, 0.25, 0.1816901138162093),
    (-0.25, 0.25, 0.1816901138162093),
    (0.05, 0.1, 0.056326537753093875),
]

@pytest.mark.parametrize('phi1,phi2,reference', testdata)
def test_compute_phase_average(phi1, phi2, reference):
    result = boats.compute_phase_average(phi1, phi2)
    np.testing.assert_almost_equal(result, reference, decimal=8)





