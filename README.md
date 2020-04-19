# ExoTETHyS

# ![ExoTETHyS logo](https://github.com/ucl-exoplanets/ExoTETHyS/blob/master/logo.png)

Version 2.x.x [![Build Status](https://travis-ci.org/ucl-exoplanets/ExoTETHyS.svg?branch=master)](https://travis-ci.org/ucl-exoplanets/ExoTETHyS) [![DOI](https://zenodo.org/badge/169268509.svg)](https://zenodo.org/badge/latestdoi/169268509)

TO BE COMPLETED  
ExoTETHyS is an open-source package for modeling exoplanetary transits, eclipsing binaries and related phenomena.

If you use this code for your research, please consider citing Morello et al. 2019, ... (arXiv:1908.09599)

## Dependencies

The code is consistent with python2/3. It makes use of
- numpy, scipy

## Download and installation
### From Pypi  (NOT YET FOR THIS VERSION)
Type on the terminal:
```pip install exotethys```  
This command installs the latest stable version of the package.

### From GitHub  
The most updated version of the package is available on GitHub. As this version is constantly under development, stability is not guaranteed.  
1. Go to <https://github.com/ucl-exoplanets/ExoTETHyS/tree/v2-beta> and click the green button "Clone or download", then click "Download ZIP" to download the whole repository. Alternatively type `git clone --single-branch --branch v2-beta https://github.com/ucl-exoplanets/ExoTETHyS` in a terminal window. /
2. Access the root folder from terminal (you may need to unzip first).
3. After accessing the root folder from terminal, type 
```pip install .```  
to install the package. Otherwise, you could import the package without installation, if you run python from the root folder.
4. To test the installation, you can type:

    ```
    pytest PATH_TO_ROOT/exotethys/tests/test_sail.py  
    pytest PATH_TO_ROOT/exotethys/tests/test_trip.py 
    pytest PATH_TO_ROOT/exotethys/tests/boats_trip.py 
    ```
NOTE: The root folder name depends on the download process. It appears to be "ExoTETHyS-2-beta" if downloaded from the web browser interface, "ExoTETHyS" if git cloned from terminal.

### Before using  
If this is the first time that you are using ExoTETHyS, you can skip this subsection.  
If you had already installed/used an older version of ExoTETHyS, you should delete the old database folder to avoid incompatibility issues. If a file from the old database is mistakenly used with this new version, the run will fail raising an error message. The database folder is named ".exotethys" and should be located in your home. You should be able locate, then remove this folder.  
This operation can also be performed by using the manage_database subpackage of ExoTETHyS, as follows:  
WARNING: The following operation is irreversible. It is highly recommended that you read more about the manage_database subpackage before deciding to perform this operation.

```
>>> from exotethys import manage_database as mdb 
>>> mdb.rm_database() 
Are you sure that you want to delete the directory /Users/pepe/.exotethys? [y/N]: y 
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
    >>> sail.ldc_calculate('PATH_TO_ROOT/examples/sail_example1.txt')   
    ```

2. TRIP (Transit Ring-Integrated Profile)  
   This subpackage can compute transit light-curves by using stellar specific intensities rather than (approximate) limb-darkening coefficients.
   
   How to run:
   
    ```
    >>> from exotethys import trip  
    >>> trip.trip_calculate('PATH_TO_ROOT/examples/trip_example.txt')  
    ```
WARNING: running this trip\_example will consume a lot of memory (>10 GB), because by default TRIP uses 100000 annuli to compute the integral stellar flux. The user can set a different number of annuli by uncommenting the line with "n\_annuli" in the trip\_example.txt and changing the relevant number (5000-10000 should be sufficient to get a nice looking light-curve, but the absolute precision is not guaranteed).

3. BOATS (Bias in the Occultation Analysis of Transiting Systems)  
   This subpackage can compute the potential bias in transit/eclipse depth due to neglecting the exoplanetary flux and/or its variation with the orbital phase (common approximations), along with the photon noise limited error bars. It can also be used as a general Exposure Time Calculator. 
   
   How to run:
   
    ```
    >>> from exotethys import boats  
    >>> boats.boats_calculate_transit('PATH_TO_ROOT/examples/boats_example4.txt')  
    >>> boats.boats_calculate_eclipse('PATH_TO_ROOT/examples/boats_example5.txt')  
    ```
NOTE: The examples are written to be launched from root directory level. Alternatively, the paths in the examples need to be personalized by the user.

4. manage_database  
   This subpackage can be used to manage the folder ".exotethys" that is created in your home the first time that a file is downloaded to perform a calculation. It contains 3 functions to list, copy and remove the content of this folder.
   
   How to run:
   
    ```
    >>> from exotethys import manage_database as mdb  
    
    >>> path, content = mdb.ls_database()  
    #path should be '/your_home/.exotethys'  
    #content is the list of folders that are inside, e.g.,  ['Atlas_2000', 'Phoenix_2012_13', 'Phoenix_2018']  
    >>> path, filenames = mdb.ls_database(grid='Atlas_2000')  
    #path should be '/your_home/.exotethys/Atlas_2000'  
    #filenames is the list of filenames in this folder, e.g.,  
    #['teff04500_logg4.5_MH0.0.pickle', 'teff04750_logg4.5_MH0.0.pickle', 'teff05000_logg4.0_MH1.0.pickle']  
    >>> path, filenames = mdb.ls_database(grid='Atlas_2000', starts='teff04', ends='.pickle', contains=['logg4.', 'MH0.0'])  
    #path should be '/your_home/.exotethys/Atlas_2000'  
    #filenames is the list of filenames in this folder starting, ending and containing the specified strings
    
    >>> mdb.cp_database('all', final_path)
    #This command copies all the folders contained in .exotethys to final_path/ExoTETHyS_database/  
    >>> mdb.cp_database('Atlas_2000', final_path)  
    #This command copies the folder Atlas_2000 to final_path (the entire folder is copied, not just the files in it)  
    >>> mdb.cp_database('Atlas_2000', final_path, files=filenames)
    #This command copies the files in the filenames list from the Atlas_2000 folder to final_path 
    #(no new Atlas_2000 folder is created at the final_path)  
    
    >>> mdb.rm_database()
    #This command will remove the whole database (irreversible). Yes/No confirmation will be asked to proceed.  
    >>> mdb.rm_database(grid='Atlas_2000')  
    #This command will remove the Atlas_2000 folder from the database (irreversible).  
    #Yes/No confirmation will be asked to proceed. 
    >>> mdb.rm_database(grid='Atlas_2000', files=filenames)  
    #This command will remove the files in the filenames list from the Atlas_2000 folder (irreversible).  
    #Yes/No confirmation will be asked to remove each file. 
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

NOTE: The parameters must be within the parameter space of the selected stellar\_models\_grid. Further information is provided in the "Description of examples". 

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

## Description of examples
**sail\_example1**: This example is to compute the limb-darkening coefficients for a single stellar target and photometric passband. It creates a file named "teff6065.0_logg4.36_MH0.0_ldc.pickle".

**sail\_example2**: This example is to test the complete output, including the stellar intensity profile and coefficients for the neighbour models. It creates three files named "teff28300.0_logg4.35_MH-0.23_ldc.pickle", "teff28300.0_logg4.35_MH-0.23_neighbour_ldc.pickle" and "teff28300.0_logg4.35_MH-0.23_neighbour_intensities.pickle".

**sail\_example3**: This example contains stellar parameters just outside the covered parameter space. If running this example, the code will exit with 2 warnings (one for each neighbour not found) and one error:
```
WARNING: teff28300.0_logg2.6_MH-0.23 cannot be calculated. Neighbour 3 not found for the stellar_models_grid Atlas_2000 .
WARNING: teff28300.0_logg2.6_MH-0.23 cannot be calculated. Neighbour 4 not found for the stellar_models_grid Atlas_2000 .
ERROR: No legal targets to calculate.
```
In particular, the code searches for the 8 nearest neighbours that defines a volume containing the requested point (28000, 2.6, -0.23 ) in the following order:  
neighbour 1: (+, +, +)  
neighbour 2: (+, +, -)  
neighbour 3: (+, -, +)  
neighbour 4: (+, -, -)  
neighbour 5: (-, +, +)  
neighbour 6: (-, +, -)  
neighbour 7: (-, -, +)  
neighbour 8: (-, -, -)  
where "+" and "-" denote a parameter value higher or lower than in the requested point. In this case the code fails to find neighbours with higher effective temperature and lower surface gravity, as they do not exist in the Atlas\_2000 database. If more targets were requested, the code would just skip this target and move to the next one. Given there are no other requested targets, the code prints the error message and exits without producing any output.  
The user can extract the information about all the available models in the grids by typing the following commands:  
```
>>> from exotethys import sail  
>>> [files_Atlas_2000, params_Atlas_2000] = sail.get_grid_parameters('Atlas_2000') 
>>> [files_Phoenix_2012_13, params_Phoenix_2012_13] = sail.get_grid_parameters('Phoenix_2012_13') 
>>> [files_Phoenix_2018, params_Phoenix_2018] = sail.get_grid_parameters('Phoenix_2018') 
```
The first variable is the list of file names in the database, the second variable is the 3-column numpy array with the corresponding stellar parameters.

**sail\_example4**: This example computes the limb-darkening coefficients for a grid of models. It creates a file named "grid_ldc.pickle".

**sail\_example5**: This example computes the limb-darkening coefficients for a single stellar target over multiple spectroscopic bins within an instrument passband. It creates a file named "teff6065.0_logg4.36_MH0.0_ldc.pickle".

**sail\_example6**: This example computes the limb-darkening coefficients for a single stellar target over multiple spectroscopic bins within an instrument passband and for another photometric passband. It creates a file named "teff6065.0_logg4.36_MH0.0_ldc.pickle".

**sail\_example7**: This example computes the limb-darkening coefficients for two targets with names over multiple spectroscopic bins within an instrument passband. It creates two files named "HD209458b_ldc.pickle" and "WASP43b_ldc.pickle".

**sail\_example8**: This example computes the limb-darkening coefficients for three targets read from file over multiple spectroscopic bins within an uniform passband. It creates three files named "HD209458b_ldc.pickle", "Sun_ldc.pickle" and "Cool_ldc.pickle".

**trip\_example**: This example computes an "exact" transit light-curve based on some auxiliary input files. It creates two files named "trip_ld_Teff6100.0_logg4.5_MH0.0_TESS.pickle" and "trip_ld_Teff6100.0_logg4.5_MH0.0_TESS.txt".  







