#ifndef HF_SCF_H
#define HF_SCF_H  

#include <vector>
#include "Molecule.h"
#include "STO.h"



class SCF
{
public:
	void add_orbital (int ngauss, int L, double* alphai, double* ai, std::vector<Atom>* coords);
	void print_all_orbitals (const char* filename);
	double compute_sto_overlap(STO* s1, STO * s2);
	double compute_gto_overlap(GTO* g1, GTO * g2);
	double od_overlap(int i, int j, double alpha, double beta, double xa, double xb);
	

private:
	std::vector <STO> orbitals;

};




#endif