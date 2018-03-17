#include <stdio.h>
#include <vector>

#include "HF_SCF.h"
#include "Molecule.h"
#include "STO.h"
#include "GTO.h"

using std::vector;

void SCF::add_orbital(int ngauss, int L, double* alphai, double* ai, vector<Atom>* coords)// указатель на вектор (один)
{


	int func_counter = 0;
	printf("size of coords %d\n", coords->size());
	for ( int nat = 0; nat < coords->size(); ++nat) //number of atoms
	{
		for (int l = 0; l <= L; ++l)
		{
			for (int k = 0; k <= (L-l); ++k)
			{
				//primitives[counter].set_data(alphai[ng], l, k, L-l-k);
				//printf("alpha from SCF %f\n", alphai[0]);
				
				//orbitals.push_back(STO(ngauss, l, k, L-l-k, alphai, ai, &((*coords)[nat]) ));

			 STO tempSTO(ngauss, l, k, L-l-k, alphai, ai, &((*coords)[nat]) );
				 orbitals.push_back(tempSTO);
				
				//printf("alpha from SCF %f\n", orbitals[0].alphai[0]);
				++ func_counter;
			}// end k
		}//end l
				//printf("%f\n",ai[ng]);
				//printf("number of functions: %d\n", nfunc);
	}//end nat
	 
	//printf("total number of functions for a given atom: %d\n", func_counter);
	//printf("orbitals size %d\n", orbitals.size());
}

void SCF::print_all_orbitals (const char* filename)
{
	FILE* fil = fopen(filename, "w");
	fprintf(fil, "Total number of orbitals: %d\n", orbitals.size());
	for (int i = 0; i < orbitals.size(); ++i)
	{
		fprintf(fil, "Orbital # %d\n", i+1);
		fprintf(fil, "quantum numbers: %d  %d  %d\n",orbitals[i].l, orbitals[i].m, orbitals[i].n );
		for (int ng = 0; ng < orbitals[i].ngauss; ++ng)
		{
			fprintf(fil, "alpha[%d] = %f   a[%d] = %f   anorm[%d] = %f  norm[%d] = %f\n",
				ng, orbitals[i].alphai[ng], ng, orbitals[i].ai[ng], ng, orbitals[i].primitives[ng].a, 
				ng, sqrt(orbitals[i].primitives[ng].norm) );

		}
		fprintf(fil, "Coordinates:\n");
		Atom temp_at( orbitals[i].primitives[0].get_coords());
		fprintf(fil,"x = %f\n",temp_at.x);
		fprintf(fil,"y = %f\n",temp_at.y);
		fprintf(fil,"z = %f\n",temp_at.z);
		fprintf(fil, "Full norm: %f\n",
			sqrt(compute_sto_overlap_v2(&orbitals[i], &orbitals[i])));
		fprintf(fil, "\n");
	}
	fclose(fil);
}

