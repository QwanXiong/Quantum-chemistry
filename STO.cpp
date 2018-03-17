#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "STO.h"
#include "GTO.h"


#include "Molecule.h"

using std::vector;


STO::STO(const STO& csto): ngauss(csto.ngauss), l(csto.l), m(csto.m), n(csto.n)
{
	alphai = new double [ngauss];
	ai = new double [ngauss];
	primitives = new GTO[ngauss];
	for (int i = 0; i < ngauss; ++i)
	{
			alphai[i] = csto.alphai[i];
			// if (L == -1)
			// {
			// 	ai[2*i+0] = _ai[2*i+0];
			// 	ai[2*i+1] = _ai[2*i+1];
			// }
			//else
			ai[i] = csto.ai[i];
			primitives[i].copy_from(csto.primitives[i]);
		
	}

	//printf("copying!!!\n");
}

STO::STO (int _ngauss, int _l, int _m, int _n, double * _alphai, double* _ai, Atom* coords):
 l(_l), m(_m), n(_n), ngauss(_ngauss)
 {


		alphai = new double [ngauss];

	// //	if (L == -1)
	// 		//ai = new double [2*ngauss];
	// //	else
		ai = new double [ngauss];

		//printf("ngauss: %d\n", ngauss);

		for (int i = 0; i < ngauss; ++i)
		{
			alphai[i] = _alphai[i];
			// if (L == -1)
			// {
			// 	ai[2*i+0] = _ai[2*i+0];
			// 	ai[2*i+1] = _ai[2*i+1];
			// }
			//else
			ai[i] = _ai[i];
		}

		
		// nang = (L + 2)*(L + 1)/2;


		// if (L == -1)
		// 	nfunc = 4*ngauss;
		// else
		// 	nfunc = ngauss*nang;
		
 		primitives = new GTO[ngauss];
		 


		// if ( L == -1)
		// {
		// 	for (int ng = 0; ng < ngauss; ++ng)
		// 	{
		// 		printf("a unnormalized: %f  %f\n",ai[2*ng+0], ai[2*ng+1]);
		// 		primitives[4*ng+0].set_data(alphai[ng], 0, 0, 0);
		// 		primitives[4*ng+0].calc_norm();

		// 		primitives[4*ng+1].set_data(alphai[ng], 1, 0, 0);
		// 		primitives[4*ng+2].set_data(alphai[ng], 0, 1, 0);
		// 		primitives[4*ng+3].set_data(alphai[ng], 0, 0, 1);
				
		// 		primitives[4*ng+1].calc_norm();
						
		// 		//ai[2*ng+0] /= sqrt(primitives[4*ng+0].norm);
		// 		//ai[2*ng+1] /= sqrt(primitives[4*ng+1].norm);

		// 		primitives[4*ng+0].set_a(ai[2*ng+0]/sqrt(primitives[4*ng+0].norm));

		// 		for (int i = 1; i < 4; ++i)
		// 		{
		// 			primitives[4*ng+i].set_a(ai[2*ng+1]/sqrt(primitives[4*ng+1].norm));
		// 		}
						
		// 				//printf("%d  %d  %d  %d  %d\n",counter, ng, l, k, L-l-k);

		// 		//++counter;
				 
		// 		printf("a normalized: %f  %f\n",ai[2*ng+0], ai[2*ng+1]);
		// 	}// end ng


		// }
		//else
		{
			for (int counter = 0; counter < ngauss; ++counter)
			{
			//primitives[counter].set_data(alphai[ng], l, k, L-l-k);
				primitives[counter].set_data(alphai[counter], l, m, n, coords);
				primitives[counter].calc_norm();

						
						 
				primitives[counter].set_a(ai[counter]/sqrt(primitives[counter].norm));
				//primitives[counter].set_a(ai[counter]);
				 
						
						//printf("%d  %d  %d  %d  %d\n",counter, ng, l, k, L-l-k);

				//++counter;
			}

			//printf("calculating\n");

		}// end else
		//printf("alpha[0] from STO: %f\n", alphai[0]);
	}// end constructor


	// double STO::compute_overlap (GTO * g1, GTO* g2)
	// {
	// 	int lsum = (g1 -> l+g2 -> l);
	// 	int msum = (g1 -> m+g2 -> m);
	// 	int nsum = (g1 -> n+g2 -> n);

	// 	if ((lsum % 2 == 1) || (msum % 2 == 1) || (nsum % 2 == 1))
	// 		return 0;

	// 	Atom temp_coords = g1->gto_coords+g2->gto_coords;

	// 	GTO gsum((g1->alpha+g2->alpha)/2.0, lsum/2, msum/2, nsum/2, &temp_coords);
	// 	gsum.calc_norm();
	// 	return gsum.norm;
	// }

	// double STO::compute_STO_norm ()
	// {
	// 	double norm = 0.0;
	// 	for (int i = 0; i < ngauss; ++i)
	// 	{
	// 		for (int j = i; j < ngauss; ++j)
	// 		{
	// 			double coef = (i == j) ? 1.0 : 2.0;
	// 			norm += coef*primitives[i].a*primitives[j].a* compute_overlap(&primitives[i], &primitives[j]);
	// 			//printf("%f  %f  %f\n", primitives[i].a, primitives[j].a, compute_overlap(&primitives[i], &primitives[j]));
	// 		}
	// 	}
	// 	return norm;
	// }