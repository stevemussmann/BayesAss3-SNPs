[![DOI](https://zenodo.org/badge/89026759.svg)](https://zenodo.org/badge/latestdoi/89026759)

# BayesAss3-SNPs
Modification of BayesAss 3.0.4 to allow handling of large SNP datasets generated via methods such as RADseq protocols.

This is an enhancement of the code produced by the Rannala Research Group: http://www.rannala.org/inference-of-recent-migration/

Most changes will not be obvious to the average user.  However, the most important changes will now allow the user to input large SNP datasets that previously would have caused BayesAss to crash with a "segmentation fault" error.  The upper limit of 420 loci has been removed.

## Citing BA3-SNPS

There are two publications you must cite for this program. Firstly, all credit for the BayesAss algorithm goes to the original authors.  Please cite: **G.A. Wilson and B. Rannala 2003. Bayesian inference of recent migration rates using multilocus genotypes. *Genetics* 163: 1177-1191.**

Secondly, our manuscript for BA3-SNPs has been accepted in Methods in Ecology and Evolution.  Please also cite: **S.M. Mussmann, M.R. Douglas, T.K. Chafin, and M.E. Douglas 2019. BA3‐SNPs: Contemporary migration reconfigured in BayesAss for next‐generation sequence data. *Methods in Ecology and Evolution* https://doi.org/10.1111/2041-210X.13252.**

## Manual Installation:

If you would prefer to install via Docker container, skip to the next section. If you would prefer to install the program manually, then keep reading. To compile, you must have a recent version of the GNU Science Library (GSL) and the C++ Boost Libraries installed on your computer.  Installation of these libraries under Ubuntu is best handled through your system's package manager.  Both can be installed with the command:

`sudo apt-get install libgsl-dev libboost-dev libboost-program-options-dev`

Download the source code, change directories into the BayesAss3-SNPs folder, and compile by issuing the command:

`make`

To install the program to your path (/usr/local/bin), use the following command:

`sudo make install`

You may receive several warnings when compiling.  I am currently working through these to eliminate them.

### Important note for Mac users
Compilation of this program will fail if GSL has been installed via Homebrew.  Please instead download the GSL source code from https://www.gnu.org/software/gsl/, compile, and install.  

## Docker Installation:

I now provide a Docker container to streamline the installation process. If you want to avoid potential issues with dependencies or manual setup, then this is now the recommended method for running the program. To get started, pull the Docker image using the following command:

```
docker pull mussmann/ba3snps:1.0
```

Launch the container by placing the "runDocker.sh" script in the folder from which you want to run the container, then executing it as shown below. This script can be found in the "Docker" folder of the github repository provided above.

```
./runDocker.sh
```

This script creates a folder named "data" in the directory on your machine from which you launched the Docker container. You can put any input files for ba3snps into this folder and they will be accessible inside the container (in /app/data/). Any outputs written to this folder and any of its subdirectories will still be accessible after you exit the container. If you write any output to other locations inside the container, they will be lost upon exit. In addition to BA3-SNPs, the countLociImmanc.sh script, [BA3-SNPS-autotune](https://github.com/stevemussmann/BA3-SNPS-autotune), and file conversion scripts ([stacksStr2Immanc.pl and pyradStr2Immanc.pl](https://github.com/stevemussmann/file_converters)) are also installed in your $PATH in this container. 

## Using BayesAss3-SNPs

The program functions similarly to the original code. Users should refer to the original BayesAss user manual, available from: https://github.com/brannala/BA3/blob/master/doc/BA3Manual.pdf. The instructions within the original documentation will help you run the program and understand the MCMC mixing parameters. There are, however, a few notable changes that users must be aware of to make BA3-SNPs function properly:
* The input immanc formatted file must be specified using the -F command line option.  Previously this was not a requirement for declaring the input file on the command line.
* The number of loci must now be specified using the -l command line option.  Previously this was not a requirement.
* Output files BA3indiv.txt and BA3trace.txt are now named based upon your input file name.  For example, running BA3-SNPs on a file named input.immanc would result in outputs named input.indiv.txt (=BA3indiv.txt) and input.trace.txt (=BA3trace.txt).


### Important note about input files
Prior to running the program, your input files should be filtered to remove any locus for which all data is missing in all individuals.  The program will crash with an error related to the GSL random number generator if this has not been done.

### Helper scripts
I have added file converters that will convert the two-line per sample Structure file format to immanc format in my file_converters repository (https://github.com/stevemussmann/file_converters).  There are separate file converters based upon popular RAD data assembly pipelines (Stacks and pyRAD). 

I now also provide a script that will count the number of loci present in .immanc formatted files. Syntax is as follows:
```
./countLociImmanc.sh -f filename
```

## List of changes:
2022-02-19:
* Added countLociImmanc.sh which will report the number of loci present in an immanc-formatted file.
* Created Docker container for executing BA3-SNPS, BA3-SNPS-autotune.py, countLociImmanc.sh, stacksStr2immanc.pl, and pyradStr2immanc.pl

2019-07-11:
* Updated version number to 1.1.
* Fixed immanc file parsing bug in which individuals would sometimes be assigned the incorrect population. This occurred when all samples from the same population did not occur consecutively within the immanc file. This bug carried over from the original BA3 source code, so you should repeat any analyses you conducted in any version of BayesAss if your input immanc file was not sorted according to population.
* Added buffer flush when printing to MCMC trace file. The trace file should update more regularly now.
* Fixed indenting in file parsing source code.
* Updated program banner.

2019-02-16:
* Fixed casting int as char to eliminate compiler warning.

2018-06-03:
* Added install command to makefile.
* Updated installation instructions.
* Updated program banner

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
