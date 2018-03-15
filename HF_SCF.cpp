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
			fprintf(fil, "alpha[%d] = %f   a[%d] = %f   anorm[%d] = %f\n",
				ng, orbitals[i].alphai[ng], ng, orbitals[i].ai[ng], ng, orbitals[i].primitives[ng].a );

		}
		fprintf(fil, "Coordinates:\n");
		Atom temp_at( orbitals[i].primitives[0].get_coords());
		fprintf(fil,"x = %f\n",temp_at.x);
		fprintf(fil,"y = %f\n",temp_at.y);
		fprintf(fil,"z = %f\n",temp_at.z);
		fprintf(fil, "\n");
	}
	fclose(fil);
}

double SCF::compute_sto_overlap(STO* s1, STO * s2)
{
	double norm = 0.0;
	for (int i = 0; i < s1->ngauss; ++i)
	{
		for (int j = i; j < s2->ngauss; ++j)
		{
			double coef = (i == j) ? 1.0 : 2.0;
			norm += coef*(s1->primitives[i].a)*(s2->primitives[j].a)* 
				compute_gto_overlap(&(s1->primitives[i]), &(s2->primitives[j]));
				//printf("%f  %f  %f\n", primitives[i].a, primitives[j].a, compute_overlap(&primitives[i], &primitives[j]));
		}
	}
	return norm;
}

double SCF::compute_gto_overlap(GTO* g1, GTO * g2)
{
	double xa = g1->gto_coords.x;
	double xb = g2->gto_coords.x;
	double ya = g1->gto_coords.y;
	double yb = g2->gto_coords.y;
	double za = g1->gto_coords.z;
	double zb = g2->gto_coords.z;

	double alpha = g1 -> alpha;
	double beta = g2 -> alpha;

	double rsq = (xa-xb)*(xa-xb)+(ya-yb)*(ya-yb) + (za-zb)*(za-zb);


	double xov = od_overlap(g1->l, g2->l, alpha, beta, xa, xb );
	double yov = od_overlap(g1->m, g2->m, alpha, beta, ya, yb );
	double zov = od_overlap(g1->n, g2->n, alpha, beta, za, zb );


	return xov*yov*zov*exp(-alpha*beta*(rsq)/(alpha+beta))*pow(M_PI/(alpha+beta), 3.0/2.0);
}

double SCF::od_overlap(int i, int j, double alpha, double beta, double xa, double xb)
{
	double eta = alpha + beta;
	double xpb = alpha/eta*(xa-xb);
	double xpa = -beta/eta*(xa-xb);

	// int curr_i = 0;
	// int curr_j = 0;

	// double sij = pow(M_PI/eta, 1.0/2.0);
	// double sim1j = 0;
	// double sijm1 = 0;

	if ((i < 0) || (j < 0))
		return 0;
	else if ((i == 0) && (j == 0))
	{
		return 1.0;//pow(M_PI/eta, 1.0/2.0)*exp(-alpha*beta*(xa-xb)*(xa-xb)/eta);
	}
	else if (i == 0)
	{
		return xpb*od_overlap(i, j-1, alpha, beta, xa, xb) + 1/2.0/eta*(j-1)*od_overlap(i, j-2, alpha, beta, xa, xb);
	}
	else
	{
		return xpa*od_overlap(i-1, j, alpha, beta, xa, xb) + 
				1/2.0/eta*((i-1)*od_overlap(i-2, j, alpha, beta, xa, xb) + j*od_overlap(i-1, j-1, alpha, beta, xa, xb));
	}

	double res = od_overlap(i ,j, alpha, beta, xa, xb);





}