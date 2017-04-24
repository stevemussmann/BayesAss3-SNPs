#ifndef INDIV_H
#define INDIV_H

#include <vector>

class Indiv {
	public:
		Indiv(int maxloci);
		Indiv(const Indiv& orig);
		virtual ~Indiv();
	protected:
		std::vector<std::vector<int> > genotype;
		std::vector<int> missingGenotypes;
		unsigned int samplePopln;
		unsigned int migrantPopln;
		unsigned int migrantAge;
		double logL;
};


#endif
