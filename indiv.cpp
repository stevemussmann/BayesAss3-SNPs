#include "indiv.h"

Indiv::Indiv(int maxloci)
{
	genotype.resize(maxloci);
	for(int i=0; i<maxloci; i++){
		genotype[i].resize(2);
	}
}

Indiv::Indiv(const Indiv& orig)
{

}

Indiv::~Indiv()
{

}

void Indiv::setAllele(int l, int m, int allele)
{
	genotype[l][m] = allele;
}

void Indiv::setSamplePopln(int currPoplnID)
{
	samplePopln = currPoplnID;
}

void Indiv::setMigrantAge(unsigned int a)
{
	migrantAge = a;
}

void Indiv::setMigrantPopln(unsigned int a)
{
	migrantPopln = a;
}

void Indiv::setLogL(double l)
{
	logL = l;
}

void Indiv::pushMissing(unsigned int j)
{
	missingGenotypes.push_back(j);
}

int Indiv::getAllele(int l, int m)
{
	return genotype[l][m];
}

int Indiv::getMissingGenotype(long unsigned int i)
{
	return missingGenotypes.at(i);
}

unsigned int Indiv::getSamplePopln()
{
	return samplePopln;
}

unsigned int Indiv::getMigrantPopln()
{
	return migrantPopln;
}

unsigned int Indiv::getMigrantAge()
{
	return migrantAge;
}

unsigned int Indiv::getMissingGenSize()
{
	return missingGenotypes.size();
}

const std::vector<int>& Indiv::getMissingVector() const
{
	return missingGenotypes;
}

double Indiv::getLogL()
{
	return logL;
}
