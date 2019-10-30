import pytest
from exotethys import sail
import numpy as np

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
