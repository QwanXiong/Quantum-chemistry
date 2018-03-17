#ifndef ALGORITHM_H
#define ALGORITHM_H 

#include "STO.h"
#include "GTO.h"

#include <stdexcept>

class Algorithm
{
public:
	Algorithm(int _N): N(_N)
	{
		if (N > 12)
		{
			throw std::range_error("N exceeds 12 in class Algorithm!");
		}
		bincoef =  new int [(N+2)*(N+1)/2];
		//printf("Size: %d\n",(N+2)*(N+1)/2 );
		fill_bincoef();
		
	}
	Algorithm(){};

	void setN(int _N)
	{
		if (_N > 12)
		{
			throw std::range_error("N exceeds 12 in class Algorithm!");
		}
		bincoef =  new int [(_N+2)*(_N+1)/2];
		//printf("Size: %d\n",(N+2)*(N+1)/2 );
		fill_bincoef();
	}

	~Algorithm()
	{
		delete [] bincoef;
	}

	int dfact (int n)
	{
		int res = 1;
		for (int i = 0; i < (n+1)/2; ++i)
			res *= (n-i*2);
		return res;
	}
	double compute_sto_overlap_v1(STO* s1, STO * s2);
	double compute_sto_overlap_v2(STO* s1, STO * s2);
	double compute_gto_overlap_v1(GTO* g1, GTO * g2);
	double compute_gto_overlap_v2(GTO* g1, GTO * g2);
	//double compute_gto_overlap_v2( );
	double od_overlap(int i, int j, double alpha, double beta, double xa, double xb);
	int get_bincoef(int n, int m);// C_n^m
	double Boysf (int n, double x);
	double compute_gto_coulomb_1e(GTO * g1, GTO * g2, Atom* C);
	//double compute_gto_coulomb_1e();

	double f(int j, int l, int m, double a, double b);
	//void compute_f_arr( );
	void compute_f_arr(int l1, int l2, int m1, int m2, int n1, int n2, Atom * PA,Atom * PB );
	double A_1e(int i, int r, int u, int l1, int l2, double ax, double ay, double az, double gamma);
	double compute_gto_coulomb_2e(GTO *g1, GTO* g2, GTO * g3, GTO *g4);
	//double plan_1e_calculation(GTO* g1, GTO *g2, Atom*);
	//double finish_1e_calculation();
	double fill_CI(int l1, int l2, int l3, int l4, 
		double PAx, double PBx, double QCx, double QDx, double QPx,double gamma1, double gamma2, double* CI );
	double fill_H(int l1, int l2, double a, double b, double gamma, double *H);
protected:
	int N;
	int* bincoef;
	void fill_bincoef();
	
	double * fx_arr;//f up to l1+l2;
	double * fy_arr;//f up to m1+m2;
	double * fz_arr;//f up to n1+n2;

	double * HL;//Lmax = l1+l2;
	double * HM;//Mmax = l3+l4;

	double * CI;//Imax = Lmax+Mmax = l1+l2+l3+l4
	double * CJ;
	double * CK;
	//double  fill_f_for_sto(GTO* g1, GTO *g2 );

	//Atom A,B,P, PA, PB, C, CP;
	//double alpha1, alpha2, gamma;
	//int l1, l2, m1, m2, n1, n2;

	
};

#endif