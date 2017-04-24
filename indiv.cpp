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


