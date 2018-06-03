# BayesAss3-SNPs
Modification of BayesAss 3.0.4 to allow handling of large SNP datasets generated via methods such as RADseq protocols.

This is an enhancement of the code produced by the Rannala Research Group: http://www.rannala.org/inference-of-recent-migration/

Most changes will not be obvious to the average user.  However, the most important changes will now allow the user to input large SNP datasets that previously would have caused BayesAss to crash with a "segmentation fault" error.  The upper limit of 420 loci has been removed.

## Installation:

To compile, you must have a recent version of the GNU Science Library (GSL) and the C++ Boost Libraries installed on your computer.  Installation of these libraries under Ubuntu is best handled through your system's package manager.  Both can be installed with the command:

`sudo apt-get install libgsl-dev libboost-dev libboost-program-options-dev`

Download the source code, change directories into the BayesAss3-SNPs folder, and compile by issuing the command:

`Make`

You may receive several warnings when compiling.  I am currently working through these to eliminate them.

### Important note for Mac users
Compilation of this program will fail if GSL has been installed via Homebrew.  Please instead download the GSL source code from https://www.gnu.org/software/gsl/, compile, and install.  

### Important note about input files
Prior to running the program, your input files should be filtered to remove any locus for which all data is missing in all individuals.  The program will crash with an error related to the GSL random number generator if this has not been done.

## Using BayesAss3-SNPs

The program functions the same as the original with the following exceptions:
* The input immanc formatted file must be specified using the -F command line option.  Previously this was not a requirement for declaring the input file on the command line.
* The number of loci must now be specified using the -l command line option.  Previously this was not a requirement.
* I have added a file converter that will convert the two-line per sample Structure file format to immanc format in my file_converters repository (https://github.com/smussmann82/file_converters).  In full disclosure, this converter has not been robustly tested.


## List of changes:
2018-03-30:
* Added a warning message for when a data file contains no data for any individuals at one locus.

2017-07-14:
* The files BA3indiv.txt and BA3trace.txt are now named based upon the input file.  For example, running the program with a file named input.immanc would result in outputs named input.indiv.txt and input.trace.txt.

2017-07-12:
* Fixed a bug that was created when the Boost Program_options library was implemented for command line option handling - the "debug" and "trace" options had accidentally switched with one another.
* Default values of "false" are now assigned to the command line boolean switches.
* Further reduced the number of global variables
* New required command line option added: --loci or -l must be used to specify the number of loci in the input file.
* Data structures for which the amount of required memory is dependent upon the number of loci are now dynamically allocated.
* I now consider the program to be usable, but plan to make additional upgrades and modify code to fit my preferences.

2017-04-26:
* Converted the indiv struct into a custom class

2017-04-21: 
* Parsing of command line options and setting of default values is now handled by the Boost Program_options library.  This has significantly reduced and simplified the code required for this purpose.
* The help menu now displays default values for all program options.
* Boolean switches operated by command line options were changed from int to bool
* The program now has fewer global variables
