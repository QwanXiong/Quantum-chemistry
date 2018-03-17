#ifndef HF_SCF_H
#define HF_SCF_H  

#include <vector>
#include "Molecule.h"
#include "STO.h"
#include "Algorithm.h"



class SCF : public Algorithm
{
public:
	// SCF(int _N){
	// 	if (N > 12)
	// 	{
	// 		throw std::range_error("N exceeds 12 in class Algorithm!");
	// 	}
	// 	bincoef =  new int [(N+2)*(N+1)/2];
	// 	//printf("Size: %d\n",(N+2)*(N+1)/2 );
	// 	fill_bincoef();
	// };
	SCF(){};
	setN(int _N)
	{
		N = _N;
		if (N > 12)
		{
			throw std::range_error("N exceeds 12 in class Algorithm!");
		}
		bincoef =  new int [(N+2)*(N+1)/2];
		//printf("Size: %d\n",(N+2)*(N+1)/2 );
		fill_bincoef();
	}
	void add_orbital (int ngauss, int L, double* alphai, double* ai, std::vector<Atom>* coords);
	void print_all_orbitals (const char* filename);
	
	


	std::vector <STO> orbitals;

	private:

};




#endif