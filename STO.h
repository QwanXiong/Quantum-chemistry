#ifndef STO_H
#define STO_H 

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <vector>

#include "GTO.h"
#include "Molecule.h"


class STO
{
public:
	//STO (int _ngauss, int _L, double * _alphai, double* _ai, double _r0); 
	STO (int _ngauss, int _l, int _m, int _n, double * _alphai, double* _ai, Atom* coords); 
	STO(const STO& csto);
	

	~STO()
	{
 
		delete [] alphai;
		delete [] ai;
		 
		
		delete [] primitives;
	}

	int   l, m ,n, ngauss   ;
//	double r0;
	double * alphai;
	double * ai;
	GTO* primitives;
	//std::vector<GTO> primitives;

	//исходный контракционный коэффициент необходимо умножать на norm, если 
	//"STO" нормирована с изначальными коэффициентами и примитивами as is
	//и делить на norm, если изначально она нормирована с нормированными примитивами!

	// double compute_overlap (GTO * g1, GTO* g2);
	 

	// double compute_STO_norm ();
 


};





#endif