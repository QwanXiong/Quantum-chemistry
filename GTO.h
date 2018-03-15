#ifndef GTO_H
#define GTO_H 

#include <cmath>
#include "Molecule.h"





class GTO
{
public:
	GTO(double _alpha, int _l, int _m, int _n, Atom* coords);
	
	GTO(){}
	~GTO(){}

	

	//for odd number n!!!
	int dfact (int n)
	{
		int res = 1;
		for (int i = 0; i < (n+1)/2; ++i)
			res *= (n-i*2);
		return res;
	}


	void copy_from(const GTO& cgto)  
	{
		l = cgto.l;
		m = cgto.m;
		n = cgto.n;

		alpha = cgto.alpha;
		a = cgto.a;
		norm = cgto.norm;

		 gto_coords = cgto.gto_coords;
		// gto_coords.x = cgto.gto_coords.x;
		// gto_coords.y = cgto.gto_coords.y;
		// gto_coords.z = cgto.gto_coords.z;

	}




	void set_data(double _alpha, int _l, int _m, int _n, Atom* coords);

	void set_a(double _a)
	{
		a = _a;
	}
	Atom gto_coords;
	double  a;
	int l;
	double alpha;
	int m, n;
	

	
	double norm;
	//double x0, y0, z0;
	

	//уже с учетом квадратного корня!
	void calc_norm();
	

	Atom get_coords()
	{
		return gto_coords;
	}

};



#endif