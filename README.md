# ExoTETHyS

# ![ExoTETHyS logo](https://github.com/ucl-exoplanets/ExoTETHyS/blob/master/logo.png)

Version 1.0.1

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
values: "Phoenix\_2018", "Phoenix\_2012\_13" or "Atlas\_2000" (choose only one of them)

**limb\_darkening\_laws** (MANDATORY)  
values: "claret4", "power2", "square_root", "quadratic", "linear", "gen\_claret" and/or "gen\_poly"

**gen\_claret\_orders** (only if gen\_claret, mandatory)  
values: positive integer numbers

**gen\_poly\_orders** (only if gen\_poly, mandatory)  
values: positive integer numbers

**passbands\_path** (OPTIONAL)  
values: path to user passbands (it should end with /)

**passbands** (MANDATORY)  
values: built-in passband names if not passbands\_path ("Kepler", "TESS", "WFC3\_G141", "STIS\_G430L", "STIS\_G750L", "irac1subarray", "irac2subarray", "irac3subarray", "irac4subarray", "uniform\_phoenix\_2012\_13", "uniform\_phoenix\_2018", "uniform\_atlas\_2000") or user file names

**wavelength\_bins\_path** (OPTIONAL)  
values: path to wavelength bins files (it should end with /)

**wavelength\_bins\_files** (OPTIONAL)  
values: user file names or "no\_bins"  
This refers to a text file with 2 columns reporting the lower and upper limits of the desired wavelength bins (each row defines a wavelength bin).

**user\_output** (OPTIONAL)  
values: "basic" (default) or "complete"

**output\_path** (OPTIONAL)  
values: path to where to store the results (it should end with /)

**targets\_path** (only if individual calculation\_type, optional)  
values: path to target file (it should end with /)

**targets\_file** (only if individual calculation\_type, optional)  
values: user file name (examples are in the "target\_list" folder)  
This refers to a text file with up to 4 columns. The first row contains the keywords in arbitrary order: target\_names (optional), star\_effective\_temperature (mandatory), star\_log\_gravity (optional) and star\_metallicity (optional). The other rows contain the corresponding key values for the requested targets. If one of the optional parameters is not available, it should be reported as "X", so that all rows have the same number of elements.

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


## TRIP configuration file
The TRIP configuration file is a text file in which each line begins with a keyword followed by one or more values associated with the keyword. The lines starting with \# are ignored; keyword values preceded by \! are also ignored. Examples of configuration files can be found in the "examples" folder.

Below I describe the available keywords:

**input\_limb\_type** (MANDATORY)  
values: "mu" or "radi" (specify first column of the input limb-darkening file)

**input\_limb\_path** (OPTIONAL)  
values: path to input limb-darkening file (it should end with /)

**input\_limb\_file** (MANDATORY)  
values: user file name (examples are in the "aux\_files" folder)
This refers to a text file with 2 columns. The first column reports the mu or radi coordinates that define the positions on the projected stellar disk (0 <= mu, radi <= 1). The second column reports the specific intensities with dimension [energy]/([area][time][wavelength][solid angle]).

**input\_series\_type** (MANDATORY)  
values: "phi" (orbital phase values - integer numbers correspond to inferior conjunction) or "time" (Julian Dates or any other time units) or "z\_sep" (sky-projected star-planet separations in units of the the stellar radius)

**input\_series\_path** (OPTIONAL)  
values: path to input timeseries file (it should end with /)

**input\_series\_file** (MANDATORY)  
values: user file name (examples are in the "aux\_files" folder)  
This refers to a text file with 1 column reporting the points at which to compute the light-curve.

**rp\_over\_rs** (MANDATORY)
values: float type (non-negative)

**sma\_over\_rs** (only if phi or time input\_series\_type, mandatory)
values: float type (sma\_over\_rs >= 1)

**inclination** (only if phi or time input\_series\_type, mandatory)
values: float type (0 <= inclination <= 90)

**eccentricity** (only if phi or time input\_series\_type, optional)
values: float type (0 <= eccentricity < 1.0, default is 0.0)

**arg\_pericenter** (only if phi or time input\_series\_type, optional)
values: float type (-180.0 < arg\_pericenter < 360.0, default is 0.0)

**period\_orbital** (only if time input\_series\_type, mandatory)
values: float type (positive)

**epoch\_of\_transit** (only if time input\_series\_type, mandatory)
values: float type  
It must be in the same units of the input timeseries.

**time\_conversion\_factor** (only if time input\_series\_type, optional)
values: float type (positive, default is 1)  
It is the ratio between timeseries and orbital period units (e.g., time\_conversion\_factor = 1.157407407407407e-05 if timeseries in seconds and orbital period in days).

**n\_annuli** (OPTIONAL)
values: positive integer number (default is 100000)

**interpolation\_type** (OPTIONAL)
values: "linear" (default), "nearest", "zero", "slinear", "quadratic", "cubic", "previous" or "next"  
Check scipy.interpolate.interp1d function for details.

**interpolation\_variable** (OPTIONAL)
values: "mu" (default) or "radi"

**cutting\_limb** (OPTIONAL)
values: "no\_cut" (default), "radi\_gradient", "mu\_gradient" or "user\_cut"  
Different options for truncating the input limb-darkening model.

**user\_cut\_mu** (only if user\_cut cutting\_limb and no user\_cut\_radi, mandatory)
values: float number (0 < user\_cut\_mu < 1)

**user\_cut\_radi** (only if user\_cut cutting\_limb and no user\_cut\_mu, mandatory)
values: float number (0 < user\_cut\_radi < 1)

**rescaling\_limb** (OPTIONAL)
values: "as\_cut" (default), "no\_rescale" or "user\_rescale"  
Different options for rescaling the mu or radi coordinates of the input limb-darkening model.

**user\_rescale\_mu** (only if user\_rescale rescaling\_limb and no user\_rescale\_radi, mandatory)
values: float number (0 < user\_rescale\_mu < 1)

**user\_rescale\_radi** (only if user\_rescale rescaling\_limb and no user\_rescale\_mu, mandatory)
values: float number (0 < user\_rescale\_radi < 1)

**rescaling\_input\_params** (OPTIONAL)
values: "no" (default) or "yes"  
If yes, the rp\_over\_rs and sma\_over\_rs values are divided by the rescaling radius.

**output\_path** (OPTIONAL)  
values: path to where to store the results (it should end with /)

**output\_filename** (OPTIONAL)  
values: string type (without extension)

**output\_fileext** (OPTIONAL)  
values: ".pickle" (default) and/or ".txt"







