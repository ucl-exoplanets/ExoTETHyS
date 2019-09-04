# ExoTETHyS

[//]: # (![MacDown logo](http://macdown.uranusjr.com/static/images/logo-160.png))

Beta version -- README tutorial in preparation

This repository stores codes to help modeling exoplanetary transits and eclipsing binaries.

If you use this code for your research, please consider citing Morello et al. 2019, ... (arXiv:1908.09599)

[//]: # (## Dependencies)

[//]: # (The code is consistent with python2/3. It makes use of)
[//]: # (- numpy, scipy, os, sys, glob, time, shutil, copy, pickle, astropy.io.fits, matplotlib)

## Download and installation
1. Go to <https://github.com/ucl-exoplanets/ExoTETHyS/> and click the green button "Clone or download", then click "Download ZIP" to download the whole repository.
2. Open a terminal window and access the ExoTETHyS-master folder (you may need to unzip first). The exact path and instructions depends on which operating systems you run.
3. After accessing the ExoTETHyS-master folder from terminal, type ```
	python setup.py install``` to install the package.  
	Alternatively you could import the package without installation, if you run python from the ExoTETHyS-master folder.
4. To test the installation, you can open a python shell and try the following commands:

    ```
    >>> import exotethys  
    >>> from exotethys import sail  
    >>> from exotethys import trip  
    ```

## List of subpackages
1. SAIL (Stellar Atmosphere Intensity Limb)  
   This subpackage can provide sets of stellar limb-darkening coefficients with
      - continuous ranges of the stellar parameters (effective temperature, surface log gravity, scaled solar metallicity);
      - built-in or user-defined passbands (with or without spectroscopic bins);
      - built-in or user-defined limb-darkening laws;
      - using different databases of stellar model-atmospheres.

   How to run:

    ```
    >>> from exotethys import sail  
    >>> sail.ldc_calculate('examples/sail_example1.txt')   
    ```

2. TRIP (Transit Ring-Integrated Profile)  
   This subpackage can compute transit light-curves by using stellar specific intensities rather than (approximate) limb-darkening coefficients.
   
   How to run:
   
    ```
    >>> from exotethys import trip  
    >>> trip.trip_calculate('examples/trip_example.txt')  
    ```

## SAIL configuration file
The SAIL configuration file is a text file in which each line begins with a keyword followed by one or more values associated with the keyword. The lines starting with \# are ignored; keyword values preceded by \! are also ignored. Examples of configuration files can be found in the "examples" folder.

Below I describe the available keywords:

**calculation\_type** (MANDATORY)  
values: "individual" (for specific targets) or "grid" (for database models)

**stellar\_models\_grid** (MANDATORY)  
values: "Phoenix_2018", "Phoenix_2012_13" or "Atlas" (choose only one of them)

**limb\_darkening\_laws** (MANDATORY)  
values: "claret4", "power2", "square_root", "quadratic", "linear", "gen\_claret" and/or "gen\_poly"

**gen\_claret\_orders** (only if gen\_claret, mandatory)  
values: positive integer numbers

**gen\_poly\_orders** (only if gen\_poly, mandatory)  
values: positive integer numbers

**passbands\_path** (OPTIONAL)  
values: path to user passbands (it should end with /)

**passbands** (MANDATORY)  
values: built-in passband names if not passbands\_path ("Kepler", "TESS", "WFC3\_G141", "STIS\_G430L", "STIS\_G750L", "irac1subarray", "irac2subarray", "irac3subarray", "irac4subarray", "uniform\_2012\_13") or user file names

**wavelength\_bins\_path** (OPTIONAL)  
values: path to wavelength bins files (it should end with /)

**wavelength\_bins\_files** (OPTIONAL)  
values: user file names or "no\_bins"

**user\_output** (OPTIONAL)  
values: "basic" (default) or "complete"

**output\_path** (OPTIONAL)  
values: path to where to store the results (it should end with /)

**targets\_path** (only if individual calculation\_type, optional)  
values: path to target file (it should end with /)

**targets\_file** (only if individual calculation\_type, optional)  
values: user file name (examples are in the "target\_list" folder)

**target\_names** (only if individual calculation\_type, if not targets\_file, optional)  
values: string type for target names

**star\_effective\_temperature** (only if individual calculation\_type, if not targets\_file, mandatory)  
values: float type for stellar temperatures

**star\_log\_gravity** (only if individual calculation\_type, if not targets\_file, optional)  
values: float type for stellar log(g) (default is 4.5)

**star\_metallicity** (only if individual calculation\_type, if not targets\_file, optional)  
values: float type for stellar \[M/H\] (default is 0.0)

**star\_minimum\_effective\_temperature** (only if grid calculation\_type, optional)  
values: single float value

**star\_maximum\_effective\_temperature** (only if grid calculation\_type, optional)  
values: single float value

**star\_minimum\_log\_gravity** (only if grid calculation\_type, optional)  
values: single float value

**star\_maximum\_log\_gravity** (only if grid calculation\_type, optional)  
values: single float value

**star\_minimum\_metallicity** (only if grid calculation\_type, optional)  
values: single float value

**star\_maximum\_metallicity** (only if grid calculation\_type, optional)  
values: single float value




