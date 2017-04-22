# BayesAss3-SNPs
Modification of BayesAss 3.0.4 to allow handling of large SNP datasets generated via methods such as RADseq protocols.

This is an enhancement of the code produced by the Rannala Research Group: http://www.rannala.org/inference-of-recent-migration/

Most changes will not be obvious to the average user.  However, the most important changes will now allow the user to input large SNP datasets that previously would have caused BayesAss to crash with a "segmentation fault" error.  The upper limit of 420 loci has been removed.

To compile, you must have a recent version of the GNU Science Library (GSL) and the C++ Boost Libraries installed on your computer.  Download the source code, change directories into the BayesAss3-SNPs folder, and compile by issuing the command:

`Make`

## Using BayesAss3-SNPs

The program functions the same as the original with the following exceptions:
* The input immanc formatted file must be specified using the -F command line option.  Previously this was not a requirement for declaring the input file on the command line.

## List of changes:

2017-04-21: 
* Parsing of command line options and setting of default values is now handled by the Boost Program_options library.  This has significantly reduced and simplified the code required for this purpose.
* The help menu now displays default values for all program options.
* Boolean switches operated by command line options were changed from int to bool
* The program now has fewer global variables
