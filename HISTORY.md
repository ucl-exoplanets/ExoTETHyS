## History of versions

**version 1.0.0**  
First online release

**version 1.0.1**  
Version under revision by the [Journal of Open Source Software (JOSS)](https://joss.theoj.org/)
* Added new database Atlas\_2000  
* Phoenix\_2012\_13 and Phoenix\_2018 databases migrated to a different online directory  
* Added check to passband wavelength limits  
* Added error message if all required targets cannot be calculated  

**version 1.0.2**  
Version fully revised by the JOSS reviewers  
* Added documentation with docstrings  
* Added ``pytest`` tests and [Travis](https://docs.travis-ci.com/) continuous integration  
* Fixed small issues and bugs  
* Added history and community guidelines files  
* Made repository pip-installable  

**version 2.0.0**  
Second online release  
* Added BOATS and manage_database subpackages  
* Added Phoenix_drift_2012 and Stagger_2015 stellar models grids in the online database  
* Using quantity arrays with astropy.units instead of numpy arrays in the files of the stellar model grids  
* Including disk-integrated stellar spectra in addition to the specific intensities in the files of the stellar model grids  
* Added new precalculated passbands  
* Slightly revised algorithms for the SAIL subpackage, but results are unaffected  
* Reduced default number of annuli for the TRIP calculation, after testing precision  
* Provided detailed user manuals for each subpackage and a revised README.adoc  
* Added FAQ section in the wiki  
* Added new examples and auxiliary files  
* New organization of the GitHub repository
