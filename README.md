# ExoTETHyS
Beta version -- writing subpackage tutorials

This repository stores codes to help modeling exoplanetary transits and eclipsing binaries.

If you use this code for your research, please consider citing
Morello et al. 2019, ... (arXiv)

INSTALLATION
python setup.py install

LIST OF SUBPACKAGES
sail: this subpackage can provide sets of stellar limb-darkening coefficients with
      - continuous ranges of the stellar parameters (effective temperature, surface log gravity, scaled solar metallicity);
      - built-in or user-defined passbands (with or without spectroscopic bins);
      - built-in or user-defined limb-darkening laws;
      - using different databases of stellar model-atmospheres.

how to run:
```
>>> from exotethys import sail  
>>> sail.ldc_calculate('examples/sail_example1.txt')  
```
      
trip: this subpackage can compute transit light-curves by using stellar specific intensities rather than (approximate) limb-darkening coefficients.
how to run:
```
>>> from exotethys import trip  
>>> trip.trip_calculate('examples/trip_example.txt')  
```
