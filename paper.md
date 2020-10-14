---
title: 'ExoTETHyS II: SAIL on a TRIP with BOATS'
tags:
  - astronomy
  - exoplanets
  - exoplanetary transits
  - galactic dynamics
  - stellar limb-darkening
authors:
- name: Giuseppe Morello
  orcid: 0000-0002-4262-5661
  affiliation: "1, 2, 3"
- name: Marine Martin-Lagarde
  affiliation: 1  
- name: Christophe Cossou
  affiliation: 1  
- name: Tiziano Zingales
  affiliation: "4, 3"  
- name: Rene Gastaud
  affiliation: 1  
- name: Pierre-Olivier Lagage
  affiliation: 1  
- name: Andrea Chiavassa
  affiliation: 5  
affiliations:
- name: "AIM, CEA, CNRS, Université Paris-Saclay, Université de Paris, F-91191 Gif-sur-Yvette, France"  
  index: 1   
- name: "Dept. of Physics & Astronomy, University College London, Gower Street, WC1E 6BT London, UK"  
  index: 2  
- name: "INAF-Osservatorio Astronomico di Palermo, Piazza del Parlamento 1, 90134 Palermo, Italy"  
  index: 3  
- name: "Laboratoire d'astrophysique de Bordeaux, Univ. Bordeaux, CNRS, B18N, allée Geoffroy Saint-Hilaire, 33615 Pessac, France"  
  index: 4  
- name: "Université Côte d'Azur, Observatoire de la Côte d'Azur, CNRS, Lagrange, CS 34229, Nice, France"  
  index: 5  
date: 8 October 2020
bibliography: paper.bib
aas-doi: 
aas-journal: The Astronomical Journal
---

# Summary

We present here the second release of ExoTETHyS: Tools for Exoplanetary Transits around Host Stars.
ExoTETHyS is a Python package which aims to provide a stand-alone set of tools for modeling spectro-photometric observations of transiting exoplanets. All the code is written in Python3 and is backward compatible with Python2. Here we provide an overview of the subpackages available within the current version:

1. SAIL (Stellar Atmosphere Intensity Limb), i.e., a calculator of stellar limb-darkening coefficients that outperforms the existing software by one order of magnitude in terms of light-curve model accuracy, i.e., down to <10 parts per million (ppm);  
2. TRIP (Transit Ring-Integrated Profile), which can compute an exact transit light-curve by direct integration of the occulted stellar flux from the model intensities, without using a parameterization (limb-darkening law) to approximate the stellar intensity profile;
3. BOATS (Bias in the Occultation Analysis of Transiting Systems), which offers a way to test common approximations in transit and eclipse light-curve fits;
4. manage_database, to deal with the storage of files that can be read by other subpackages.

We refer to the relevant papers for description of the scientific applications and performances of SAIL and TRIP (@morello2020a), and BOATS (@morello2020b,@martin-lagarde2020).
Technical usage details can be found on GitHub.

# Statement of need

The SAIL subpackage of ExoTETHyS has been used in the analysis of spectroscopic observations of exoplanetary transits (e.g., @houyip2020a, @changeat2020), as well as in future instrument simulators (@sarkar2020), and other technical studies (@houyip2020b).
The new BOATS subpackage includes several functions that are especially useful in the preparation of observing proposals, and to account for some less known effects in the analysis of exoplanetary transits and eclipses. The high significance of such effects has been demonstrated for some targets to be osberved with the upcoming James Webb Space Telescope.

# Differences with the previous version

We summarize the list of changes:

- BOATS and manage_database are brand-new subpackages of this release;
- The online database has been extended and changed in format;
- New functions and technical changes have been implemented in pre-existing subpackages.

## Database updates

Some calculations make use of precalculated stellar model files. These files are stored privately, but they are automatically downloaded if requested by running programs.
Several grids of stellar-atmosphere models were already available, i.e., ATLAS\_2000 (@claret2000), PHOENIX\_2012\_13 (@claret2012,@claret2013), and PHOENIX\_2018 (@claret2018).
We have now added a Stagger\_2015 grid (@magic2015,@chiavassa2018).

Each file contains a dictionary that includes the specific intensities emerging from several points on the stellar disk at various wavelengths. In the previous version of the database files, the wavelength and intensities were arranged on numpy arrays, their units were specified in the accompanying software documentation.
The new database files incorporate information on physical units in quantity arrays by means of the astropy.units library.
Additionally, they also include a quantity vector with the disk-integrated flux as obtained from the intensity profiles at each wavelength, in order to speed up the new algorithms that make use of stellar spectra (without spatial/angular resolution).


## SAIL updates

Some of the SAIL functions were slightly modified to deal with the new format of the database files.
Another modification concerns the calculation of the product between the model spectrum of stellar intensities and the instrumental response.
In Exotethys version 1, the instrumental response was interpolated on the wavelengths of the stellar model, the latter usually having a higher spectral resolution. However, this approach is not always valid. For example, the Atlas\_2000 spectra are very sparsely sampled in the infrared at wavelengths longer than 10 $\mu$m, so that instrumental passbands or spectral bins narrower than the distance between two consecutive wavelengths of the model were zeroed.
In Exotethys version 2, both the instrumental response and model spectra are interpolated on a common wavelength grid with R=10000 or higher covering the requested wavelength range. We checked that this approach overcomes the issue explained above, otherwise leading to results identical to those obtained with the previous method in all non-problematic cases. The new method is also more efficient, as it limits the interpolation interval to that of interest instead of using the whole spectral range of the models.


## TRIP updates

The only change for the TRIP subpackage is the reduction of the number of annuli used by default to compute the occulted stellar flux from 100000 to 20000. We checked that this number is more than sufficient to guarantee the numerical precision well below 1 ppm, while significantly reducing the memory usage and computing time.


# Acknowledgments  

This work was supported by the LabEx P2IO, the French ANR contract 05-BLANNT09-573739, the European Unions Horizon 2020 Research and Innovation Programme (Grant Agreement N 776403). T.Z. has been founded by the European Research Council (ERC) under the European Union's Horizon 2020 Research and Innovation Programme (Grant Agreement N 679030/WHIPLASH).


# References
