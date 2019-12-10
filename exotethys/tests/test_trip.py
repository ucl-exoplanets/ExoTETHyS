import pytest
from exotethys import trip
import numpy as np


testdata = [
    ('ciao', 'ciao'),
    ('5*5', '5*5'),
    ('25', 25.0),
    ('-3.14', -3.14),
]

@pytest.mark.parametrize('toread,reference', testdata)
def test_str2float(toread, reference):
    result = trip.str2float(toread)
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
    result = trip.check_length(vector, **args)
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
    result = trip.check_type_in_list(vector, item_type)
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
    result = trip.check_integers(vector, **args)
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
    result = trip.check_2Darray(vector, **args)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    (np.arange(6),  True),
    (np.atleast_2d(np.arange(6)), False),
    (np.array([[1, 2, 3], [4, 5, 6]]), False),
    (np.array([]), False),
]

@pytest.mark.parametrize('vector,reference', testdata)
def test_check_1Darray(vector, reference):
    result = trip.check_1Darray(vector)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    ('rp_over_rs',  None),
    ('eccentricity', 0.0),
]

@pytest.mark.parametrize('query,reference', testdata)
def test_get_default_value(query, reference):
    result = trip.get_default_value(query)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    ('sma_over_rs',  (1.0, None, None, None, 'sma_over_rs>=1.0')),
    ('eccentricity', (0.0, None, None, 1.0, '0.0<=eccentricity<1.0')),
    ('period_orbital', (None, None, 0.0, None, 'period_orbital>0.0')),
    ('time_conversion_factor', (None, None, 0.0, None, 'time_conversion_factor>0.0')),
]

@pytest.mark.parametrize('query,reference', testdata)
def test_get_params_limits(query, reference):
    result = trip.get_params_limits(query)
    np.testing.assert_equal(result, reference)


#############################


testdata = [
    (np.array([0.0, 0.5, 1.0, 1.5]), {'min_value': 0.0},  True),
    (np.array([0.0, 0.5, 1.0, 1.5]), {'min_value': 0.0, 'max_value': 1.0},  False),
    (np.array([0.0, 0.5, 1.0, 1.5]), {'inf_value': 0.0},  False),
    (np.array([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5]), {},  True),
]

@pytest.mark.parametrize('vector,args,reference', testdata)
def test_check_values(vector, args, reference):
    result = trip.check_values(vector, **args)
    np.testing.assert_equal(result, reference)


#############################

testdata = [
    (np.array([0.0, 0.05, 0.1, 0.15, 0.20, 0.25]), 90.0, 9.0, 0.0, 0.0, np.array([0.0, 2.78115295, 5.29006727, 7.28115295, 8.55950865, 9.0])),
    (np.array([0.0, 0.05, 0.1, 0.15, 0.20, 0.25]), 90.0, 9.0, 0.0, 45.0, np.array([0.0, 2.78115295, 5.29006727, 7.28115295, 8.55950865, 9.0])),
    (np.array([0.0, 0.05, 0.1, 0.15, 0.20, 0.25]), 88.0, 9.0, 0.95, 45.0, np.array([0.0183187, 5.65862337, 8.12572219, 9.76760023, 10.9270278, 11.74255736])),
]

@pytest.mark.parametrize('phi, inclination, sma_over_rs, eccentricity, arg_pericenter, reference', testdata)
def test_compute_z_sep(phi, inclination, sma_over_rs, eccentricity, arg_pericenter, reference):
    result = trip.compute_z_sep(phi, inclination, sma_over_rs, eccentricity, arg_pericenter)
    np.testing.assert_almost_equal(result, reference, decimal=8)


#############################

mu_test = np.array([0.016775, 0.018307, 0.019740, 0.021091, 0.022371, 0.023592, 0.024761, 0.025884, 0.026967, 0.028015, 0.029032, 0.030020, 0.030984, 0.031925, 0.032845, 0.033747, 0.034632, 0.035501, 0.036355, 0.037195, 0.038023, 0.038839, 0.039644, 0.040438, 0.041224, 0.042002, 0.042772, 0.043534, 0.044288, 0.045033, 0.045764, 0.046473, 0.047149, 0.047773, 0.048326, 0.048784, 0.049176, 0.049434, 0.049640, 0.049846, 0.050065, 0.050300, 0.050557, 0.050839, 0.051149, 0.051490, 0.051866, 0.052279, 0.052732, 0.071677, 0.090622, 0.109570, 0.128510, 0.147460, 0.166400, 0.185350, 0.204290, 0.223240, 0.242190, 0.261130, 0.280080, 0.299020, 0.317970, 0.336910, 0.355860, 0.374800, 0.393750, 0.412690, 0.431640, 0.450580, 0.469530, 0.488480, 0.507420, 0.526370, 0.545310, 0.564260, 0.583200, 0.602150, 0.621090, 0.640040, 0.658980, 0.677930, 0.696870, 0.715820, 0.734760, 0.753710, 0.772660, 0.791600, 0.810550, 0.829490, 0.848440, 0.867380, 0.886330, 0.905270, 0.924220, 0.943160, 0.962110, 0.981050, 1.000000])

ints_test = np.array([3.89376100e-04, 5.52602300e-04, 7.19631400e-04, 9.32646700e-04, 1.21701850e-03, 1.60428810e-03, 2.13953270e-03, 2.88093200e-03, 3.91488780e-03, 5.36872860e-03, 7.40813830e-03, 1.02938394e-02, 1.43877708e-02, 2.02052327e-02, 2.84783043e-02, 4.01805184e-02, 5.65845803e-02, 7.92191135e-02, 1.09642158e-01, 1.48959723e-01, 1.96961847e-01, 2.50804630e-01, 3.04507265e-01, 3.49998079e-01, 3.81146735e-01, 3.98119111e-01, 4.06390960e-01, 4.11325778e-01, 4.15210461e-01, 4.18602985e-01, 4.21629051e-01, 4.24335277e-01, 4.26738768e-01, 4.28827729e-01, 4.30590390e-01, 4.31995561e-01, 4.33158100e-01, 4.33906669e-01, 4.34494763e-01, 4.35076343e-01, 4.35684123e-01, 4.36330132e-01, 4.37023837e-01, 4.37771585e-01, 4.38579980e-01, 4.39454066e-01, 4.40398814e-01, 4.41414485e-01, 4.42504435e-01, 4.75967106e-01, 4.99535830e-01, 5.19261846e-01, 5.36853958e-01, 5.53055524e-01, 5.68261032e-01, 5.82708589e-01, 5.96553651e-01, 6.09905216e-01, 6.22844561e-01, 6.35432544e-01, 6.47716067e-01, 6.59731522e-01, 6.71508719e-01, 6.83071221e-01, 6.94436103e-01, 7.05618665e-01, 7.16628347e-01, 7.27476743e-01, 7.38169314e-01, 7.48711531e-01, 7.59109611e-01, 7.69365426e-01, 7.79483695e-01, 7.89466096e-01, 7.99315152e-01, 8.09032478e-01, 8.18619964e-01, 8.28080297e-01, 8.37414092e-01, 8.46623246e-01, 8.55708593e-01, 8.64672516e-01, 8.73517101e-01, 8.82242088e-01, 8.90850847e-01, 8.99343430e-01, 9.07722899e-01, 9.15990356e-01, 9.24146685e-01, 9.32194839e-01, 9.40135070e-01, 9.47969489e-01, 9.55700367e-01, 9.63330127e-01, 9.70858186e-01, 9.78288037e-01, 9.85620210e-01, 9.92857448e-01, 1.00000000e+00])    

testdata = [
    (mu_test, ints_test, 0.0392415),
]

@pytest.mark.parametrize('array1, array2, reference', testdata)
def test_get_max_gradient(array1, array2, reference):
    result = trip.get_max_gradient(array1, array2)
    np.testing.assert_almost_equal(result, reference, decimal=8)


#############################

radi_test = np.sqrt(1.0-mu_test**2.0)

testdata = [
    (radi_test, 1.0, 1.0, ints_test, [mu_test, radi_test, ints_test]),
    (radi_test, 0.9992297557, 1.0, ints_test, [mu_test[22:], radi_test[22:], ints_test[22:]]),
    (radi_test, 0.9992297557, 0.9992297557, ints_test, [np.sqrt(1.0-(radi_test[22:]/0.9992297557)**2.0), radi_test[22:]/0.9992297557, ints_test[22:]]),
]

@pytest.mark.parametrize('array1, r0cut, r0res, array2, reference', testdata)
def test_limbmodel_cut_and_rescale(array1, r0cut, r0res, array2, reference):
    result = trip.limbmodel_cut_and_rescale(array1, r0cut, r0res, array2)
    np.testing.assert_almost_equal(result, reference, decimal=10)

#############################





