/*
 *  BA3.h
 *  BA3
 *
 *  Created by Bruce Rannala on 11/29/10.
 *  Copyright 2010 University of California Davis. All rights reserved.
 *
 */

#ifndef BA3_HEADERS
#define BA3_HEADERS
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <map>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <cmath>

#include "indiv.h"

using namespace std;

//const int MAXLOCI=10000;
const int MAXALLELE=300;
const int MAXPOPLN=100;
const int MAXINDIV=5000;
const int MAXLINELENGTH=1000000;
/*
struct indiv
{
	int genotype[MAXLOCI][2];
	vector<int> missingGenotypes;
	unsigned int samplePopln;
	unsigned int migrantPopln;
	unsigned int migrantAge;
	double logL;
};
*/
struct ancestryProbs
{
	double poplnPostProb[MAXPOPLN][3];
};

void printBanner(void);
void readInputFile(Indiv **sampleIndiv, unsigned int &noIndiv, unsigned int &noLoci, unsigned int &noPopln, unsigned int *noAlleles, string &infileName, int MAXLOCI);
void getEmpiricalAlleleFreqs(double ***alleleFreqs, Indiv **sampleIndiv, unsigned int *noAlleles, unsigned int noPopln, unsigned int noLoci, unsigned int noIndiv, bool debug);
void fillMigrantCounts(Indiv **sampleIndiv, long int ***migrantCounts, unsigned int noIndiv, unsigned int noPopln);
double migCountLogProb(long int ***migrantCounts, double **migrationRates, unsigned int noPopln);
double logLik(Indiv *Indiv, double ***alleleFreqs, double *FStat, unsigned int noLoci);
double oneLocusLogLik(Indiv *Indiv, double ***alleleFreqs, double *FStat, int chosenLocus);
// void proposeMigrantAncDrop(int &migrantPopln, int &migrantAge, int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncDrop(unsigned int &migrantPopln, unsigned int &migrantAge, unsigned int samplePopln, int noPopln, long int ***migrantCounts);
void proposeMigrantAncAdd(unsigned int &migrantPopAdd, unsigned int &migrantAgeAdd,unsigned int migrantPopDrop, unsigned int migrantAgeDrop, 
						  unsigned int samplePopln, int noPopln);
void parseComLine(int argc, char **argv, string &infileName, int &seed, unsigned int &mciter, unsigned int &sampling, unsigned int &burnin, string &outfileName,
         double &deltaM, double &deltaA, double &deltaF, bool &verbose, bool& settings,
         bool &genotypes, bool &trace, bool &debug, bool &nolikelihood, int &MAXLOCI);
std::string remove_extension(const std::string& filename);
#endif
