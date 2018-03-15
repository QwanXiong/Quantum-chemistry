

#include <cmath>
#include <stdio.h>
#include "GTO.h"
#include "Molecule.h"

GTO::GTO(double _alpha, int _l, int _m, int _n, Atom* coords): alpha(_alpha), l(_l), m(_m), n(_n) {
	 // 	x0 = coords->x;
		// y0 = coords->y;
		// z0 = coords->z;
	gto_coords.x = coords->x;
	gto_coords.y = coords->y;
	gto_coords.z = coords->z;

}



void GTO::set_data(double _alpha, int _l, int _m, int _n, Atom* coords)
{
	alpha = _alpha;
	l = _l;
	m = _m;
	n = _n;
		// x0 = coords->x;
		// y0 = coords->y;
		// z0 = coords->z;
	gto_coords.x = coords->x;
	gto_coords.y = coords->y;
	gto_coords.z = coords->z;
}

void GTO::calc_norm()
{
	double xf = (double)dfact(2*l-1);
	double yf = (double)dfact(2*m-1);
	double zf = (double)dfact(2*n-1);
		// /norm =  xf*yf*zf / pow(2.0*alpha, (l+m+n+1.5)) / pow(2.0, (l+m+n)) * pow(M_PI , 1.5);
	norm = (xf*yf*zf / pow(2.0*alpha, (l+m+n+1.5)) / pow(2.0, (l+m+n)) * pow(M_PI , 1.5));
	//printf("normalization: %f\n", norm);
}
	