#ifndef INDIV_H
#define INDIV_H

#include <vector>

class Indiv {
	public:
		Indiv(int maxloci);
		Indiv(const Indiv& orig);
		virtual ~Indiv();
		void setAllele(int l, int m, int allele);
		void setSamplePopln(int currPoplnID);
		void setMigrantAge(unsigned int a);
		void setMigrantPopln(unsigned int a);
		void setLogL(double l);
		void pushMissing(unsigned int j);
		int getAllele(int l, int m);
		int getMissingGenotype(long unsigned int i);
		unsigned int getSamplePopln();
		unsigned int getMigrantPopln();
		unsigned int getMigrantAge();
		unsigned int getMissingGenSize();
		const std::vector<int>& getMissingVector() const;
		double getLogL();
	private:
		std::vector<std::vector<int> > genotype;
		std::vector<int> missingGenotypes;
		unsigned int samplePopln;
		unsigned int migrantPopln;
		unsigned int migrantAge;
		double logL;
};


#endif
