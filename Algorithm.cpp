

#include <stdio.h>
#include <vector>
#include <gsl/gsl_sf_gamma.h>

#include "Algorithm.h"
 
#include "Molecule.h"
#include "STO.h"
#include "GTO.h"

#include <stdio.h>



void Algorithm::fill_bincoef()
{
	int nfact = 1;
	int mfact = 1;
	int counter = 0;
	int nmfact = 1;
	for (int n = 0; n <= N; ++n)
	{
		for (int i = 1; i <= n; ++i)
		{
			nmfact *= i;
		}
		
		for (int m = 0; m <= n; ++m)
		{
			bincoef[counter] = nfact/mfact/nmfact;
			//printf("%d\n", bincoef[counter]);
			//printf("n: %d, m: %d, C_n^m %d\n",n,m, bincoef[counter] );
			mfact *= (m+1);
			(m == n) ? (nmfact = 1) :(nmfact /= (n-m));
			//;
			++counter;	
			//printf("%d\n", counter);
		}
		mfact = 1;
		nmfact = 1;
		nfact *= (n+1);
	}

	printf("BINCOEF SUCCESSFULLY FILLED\n");
}
//TODO: Как-нибудь вместить заполнение массива прямо в функцию f, чтобы не происходило лишних вызовов
void Algorithm::compute_f_arr(int l1, int l2, int m1, int m2, int n1, int n2, Atom * PA,Atom * PB)
{
	for (int i = 0; i <= l1+l2; ++i)
	{
		 fx_arr[i] = f(i,l1,l2,PA->x,PB->x);
		// printf("computing %f\n",f(i,l1,l2,PA->x,PB->x));
	}

	for (int i = 0; i <= m1+m2; ++i)
	{
		 fy_arr[i] = f(i,m1,m2,PA->y,PB->y);
	}

	for (int i = 0; i <= n1+n2; ++i)
	{
		 fz_arr[i] = f(i,n1,n2,PA->z,PB->z);
	}
}

int Algorithm::get_bincoef(int n, int m)
{
	return bincoef[(n+1)*n/2+m];
}

double Algorithm::Boysf(int n, double x)
{
	return 0.5*pow(x,-0.5-n)*(gsl_sf_gamma(0.5+n) - gsl_sf_gamma_inc(0.5+n,x));
}

double Algorithm::compute_sto_overlap_v1(STO* s1, STO * s2)
{
	double norm = 0.0;
	for (int i = 0; i < s1->ngauss; ++i)
	{
		for (int j = i; j < s2->ngauss; ++j)
		{
			double coef = (i == j) ? 1.0 : 2.0;
			//double coef = 1.0;
			norm += coef*(s1->primitives[i].a)*(s2->primitives[j].a)* 
				compute_gto_overlap_v1(&(s1->primitives[i]), &(s2->primitives[j]));
			// norm += coef* 
			// 	compute_gto_overlap(&(s1->primitives[i]), &(s2->primitives[j]));
				//printf("%s\n", );
				// printf("%f  %f  %e\n", s1->primitives[i].a, s2->primitives[j].a, 
				// 	compute_gto_overlap(&(s1->primitives[i]), &(s2->primitives[j])));
		}
	}
	return norm;
}



double Algorithm::compute_gto_overlap_v1(GTO* g1, GTO * g2)
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

	//printf("xov, yox, zov : %f %f %f\n", xov, yov, zov);
	double res = xov*yov*zov*exp(-alpha*beta*(rsq)/(alpha+beta))*pow(M_PI/(alpha+beta), 3.0/2.0);
	// printf("res: %e\n", exp(-alpha*beta*(rsq)/(alpha+beta)));
	// printf("res: %e\n", exp(-alpha*beta*(rsq)/(alpha+beta)));

	return res ;
}

double Algorithm::compute_sto_overlap_v2(STO* s1, STO * s2)
{
	double norm = 0.0;
	for (int i = 0; i < s1->ngauss; ++i)
	{
		for (int j = i; j < s2->ngauss; ++j)
		{
			//fill_f_for_sto(&(s1->primitives[i]),&(s2->primitives[j]));//new inside!!!
			double coef = (i == j) ? 1.0 : 2.0;
			//double coef = 1.0;
			norm += coef*(s1->primitives[i].a)*(s2->primitives[j].a)* 
				compute_gto_overlap_v2( &(s1->primitives[i]), &(s2->primitives[j]) );
				//finish_1e_calculation();// delete inside!!!
			// norm += coef* 
			// 	compute_gto_overlap(&(s1->primitives[i]), &(s2->primitives[j]));
				//printf("%s\n", );
				// printf("%f  %f  %e\n", s1->primitives[i].a, s2->primitives[j].a, 
				// 	compute_gto_overlap(&(s1->primitives[i]), &(s2->primitives[j])));
		}
	}
	return norm;
}

// double Algorithm::plan_1e_calculation(GTO* g1, GTO *g2, Atom* nuc)
// {

// 	alpha1 = g1->alpha;
// 	alpha2 = g2->alpha;
// 	gamma = alpha1+alpha2;
// 	A = g1->gto_coords;
// 	B = g2->gto_coords;

// 	C = *nuc;

// 	P = (alpha1*A+alpha2*B)/gamma;
// 	PA = P - A;
// 	PB = P - B;
// 	CP = P-C;
	
// 	l1 = g1->l;
// 	l2 = g2->l;

// 	m1 = g1->m;
// 	m2 = g2->m;

// 	n1 = g1->n;
// 	n2 = g2->n;

// 	fx_arr = new double [l1+l2+1];
// 	fy_arr = new double [m1+m2+1];
// 	fz_arr = new double [n1+n2+1];
// 	compute_f_arr( );
// }

// double Algorithm::fill_f_for_sto(GTO* g1, GTO *g2 )
// {

// 	alpha1 = g1->alpha;
// 	alpha2 = g2->alpha;
// 	gamma = alpha1+alpha2;
// 	A = g1->gto_coords;
// 	B = g2->gto_coords;

	 

// 	P = (alpha1*A+alpha2*B)/gamma;
// 	PA = P - A;
// 	PB = P - B;
	 
	
// 	l1 = g1->l;
// 	l2 = g2->l;

// 	m1 = g1->m;
// 	m2 = g2->m;

// 	n1 = g1->n;
// 	n2 = g2->n;

// 	fx_arr = new double [l1+l2+1];
// 	fy_arr = new double [m1+m2+1];
// 	fz_arr = new double [n1+n2+1];
// 	compute_f_arr( );
// }

// double Algorithm::finish_1e_calculation()
// {
// 	delete [] fx_arr;
// 	delete [] fy_arr;
// 	delete [] fz_arr;
// }
//сгенерировать fx, fy, fz прямо в sto_overlap, ведь для каждого примитива расстояние одинаково
//более того, эти 3 массива надо генерировать сразу в одной функции, чтобы потом дальше эти же 3 матрицы
//использовать как для матрицы перекрывания, так и для 1е кулоновских интегралов
double Algorithm::compute_gto_overlap_v2(GTO* g1, GTO * g2)
//double Algorithm::compute_gto_overlap_v2(  )
{
	double alpha1 = g1->alpha;
	double alpha2 = g2->alpha;
	  double gamma = alpha1+alpha2;

	// P = (alpha1*A+alpha2*B)/gamma;
	// PA = P - A;
	// PB = P - B;
	// CP = P-C;
	// 
	Atom P = (alpha1*g1->gto_coords+alpha2*g2->gto_coords)/gamma;
	Atom PA = P - g1->gto_coords;
	Atom PB = P - g2->gto_coords;
	//Atom CP = P -*C;

	fx_arr = new double [g1->l+g2->l+1];
	fy_arr = new double [g1->m+g2->m+1];
	fz_arr = new double [g1->n+g2->n+1];

	compute_f_arr(g1->l,g2->l,g1->m,g2->m,g1->n,g2->n,&PA, &PB);
	 
	double sumx = fx_arr[0], sumy = fy_arr[0], sumz = fz_arr[0];
	double tim1df = 1.0; //double factorial of (2i-1)
	for (int i = 1; i <= (g1->l+g2->l)/2; ++i)
	{
		tim1df *= (2*i-1);
		sumx += fx_arr[2*i]*tim1df/pow(2*gamma,i);

	}

	double tjm1df = 1.0; //double factorial of (2i-1)
	for (int j = 1; j <= (g1->m+g2->m)/2; ++j)
	{
		tjm1df *= (2*j-1);
		sumy += fy_arr[2*j]*tjm1df/pow(2*gamma,j);

	}
 
 	double tkm1df = 1.0; //double factorial of (2i-1)
	for (int k = 1; k <= (g1->n+g2->n)/2; ++k)
	{
		tkm1df *= (2*k-1);
		sumz += fz_arr[2*k]*tkm1df/pow(2*gamma,k);

	}

	
	delete [] fx_arr;
	delete [] fy_arr;
	delete [] fz_arr;

	double rsq =g1->gto_coords.dist_sq(g2->gto_coords);


	// double xov = od_overlap(g1->l, g2->l, alpha, beta, xa, xb );
	// double yov = od_overlap(g1->m, g2->m, alpha, beta, ya, yb );
	// double zov = od_overlap(g1->n, g2->n, alpha, beta, za, zb );

	//printf("xov, yox, zov : %f %f %f\n", xov, yov, zov);
	double res = sumx*sumy*sumz*exp(-alpha1*alpha2*(rsq)/gamma)*pow(M_PI/gamma, 3.0/2.0);
	// printf("res: %e\n", exp(-alpha*beta*(rsq)/(alpha+beta)));
	// printf("res: %e\n", exp(-alpha*beta*(rsq)/(alpha+beta)));

	return res ;
}

double Algorithm::od_overlap(int i, int j, double alpha, double beta, double xa, double xb)
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

double Algorithm::compute_gto_coulomb_1e(GTO * g1, GTO * g2, Atom* C)
//double Algorithm::compute_gto_coulomb_1e(  )
{
	double alpha1 = g1->alpha;
	double alpha2 = g2->alpha;
	double gamma = alpha1+alpha2;
	Atom P = (alpha1*g1->gto_coords+alpha2*g2->gto_coords)/gamma;
	Atom PA = P - g1->gto_coords;
	Atom PB = P - g2->gto_coords;
	Atom CP = P -*C;

	 //double cp_sq = P.dist_sq(*C);
	double cp_sq = CP.x*CP.x+CP.y*CP.y + CP.z*CP.z;
	double ab_sq = g1->gto_coords.dist_sq(g2->gto_coords);
	//double ab_sq=A.dist_sq(B);

	double sumx, sumy, sumz , sum = 0.0  ;
	int nu;
	for (int i = 0; i <= g1->l+g2->l; ++i)
	{
		 for (int r = 0; r <= i/2; ++r)
		 {
		 	for (int u = 0; u <= (i-2*r)/2; ++u)
		 	{
		 		for (int j = 0; j <= g1->m+g2->m; ++j)
				{
					 for (int s = 0; s <= j/2; ++s)
					 {
					 	for (int v = 0; v <= (j-2*s)/2; ++v)
					 	{
					 		
					 		for (int k = 0; k <=  g1->n+ g2->n; ++k)
							{
								 for (int t = 0; t <= k/2; ++t)
								 {
								 	for (int w = 0; w <= (k-2*t)/2; ++w)
								 	{
								 		sumz = A_1e(k, t, w, g1->n,g1->n, PA.z, PB.z, CP.z, gamma);
								 		sumy = A_1e(j, s, v, g1->m, g2->m, PA.y, PB.y, CP.y, gamma);
					 					sumx = A_1e(i, r, u, g1->l, g2->l, PA.x, PB.x, CP.x, gamma);
					 					nu = i + j + k - 2*(r + s + t) - (u + v +w);
					 					sum  += sumx*sumy*sumz*Boysf(nu,cp_sq*gamma);
					 					printf("calculating sum!\n");
								 	}
								 }
							}

					 	}
					 }
				}
		 		 
		 	}
		 }
	}
 	
 // 	printf("Boys: %f\n", Boysf(1,3.14));
	// return  2.0*M_PI/gamma*exp(-alpha1*alpha2*ab_sq/gamma)*
	// (alpha1/gamma*(g1->gto_coords.x-g2->gto_coords.x)*Boysf(0,cp_sq*gamma)+CP.x*Boysf(1,cp_sq*gamma));
	return sum*2.0*M_PI/gamma*exp(-alpha1*alpha2*ab_sq/gamma);


}

double Algorithm::A_1e(int i, int r, int u, int l1, int l2, double pa, double pb, double cp, double gamma)
{
	double sign = pow(-1.0,i+u);
	double num = gsl_sf_fact(i)*pow(cp, i-2*r-2*u)*pow(1/4.0/gamma,r+u);
	double den  = gsl_sf_fact(r)*gsl_sf_fact(u)*gsl_sf_fact(i-2*r-2*u);

	return sign*num/den*f(i,l1,l2,pa,pb);
}

double Algorithm::f(int j, int l, int m, double a, double b)
{
	double res = 0.0;
	int Ninterm = l + m -j;
	//Ninterm = (Ninterm > l) ? l : Ninterm;

	for (int k = 0; k <= Ninterm; ++k)
	{
		 
		if ((Ninterm - k) > m)
			continue;
		if (k > l)
			break;
		 
		res += (double)get_bincoef(l,k)* (double)get_bincoef(m,Ninterm-k)*pow(a, k)*pow(b, Ninterm-k);
		// double bin1 = gsl_sf_fact(l)/gsl_sf_fact(k)/gsl_sf_fact(l-k);
		// double bin2 = gsl_sf_fact(m)/gsl_sf_fact(Ninterm-k)/gsl_sf_fact(m-Ninterm+k);
		// res += bin1*bin2*pow(a, k)*pow(b, Ninterm-k);
		//printf("calc\n");
		//printf("k = %d  Ninterm = %d a^k = %f b^(N-k) = \n", );
	}
	return res;
}

double Algorithm::compute_gto_coulomb_2e(GTO *g1, GTO* g2, GTO * g3, GTO *g4)
{
	int l1, l2, l3, l4;
	int m1, m2, m3, m4;
	int n1, n2, n3, n4;

	l1 = g1->l;
	l2 = g2->l;
	l3 = g3->l;
	l4 = g4->l;

	m1 = g1->m;
	m2 = g2->m;
	m3 = g3->m;
	m4 = g4->m;

	n1 = g1->n;
	n2 = g2->n;
	n3 = g3->n;
	n4 = g4->n;

	double alpha1 = g1->alpha;
	double alpha2 = g2->alpha;
	double alpha3 = g3->alpha;
	double alpha4 = g4->alpha;

	double gamma1 = alpha1+alpha2;
	double gamma2 = alpha3+alpha4;

	Atom P = (alpha1*g1->gto_coords+alpha2*g2->gto_coords)/gamma1; 
	Atom Q = (alpha3*g3->gto_coords+alpha4*g4->gto_coords)/gamma2; 
	Atom QP = Q-P;
	Atom PA = P-g1->gto_coords;
	Atom PB = P-g2->gto_coords;
	Atom QC = Q-g3->gto_coords;
	Atom QD = Q-g4->gto_coords;


	double pq_sq = Q.dist_sq(P);
	double ab_sq = g1->gto_coords.dist_sq(g2->gto_coords);
	double cd_sq = g3->gto_coords.dist_sq(g4->gto_coords);

	double delta = 1/gamma1/4.0 + 1/gamma2/4.0; 

	
	int numax = l1+l2+l3+l4+m1+m2+m3+m4+n1+n2+n3+n4;
	int Imax = l1+l2+l3+l4;
	int Jmax = m1+m2+m3+m4;
	int Kmax = n1+n2+n3+n4;

	CI = new double[Imax+1];
	CJ = new double[Jmax+1];
	CK = new double[Kmax+1];


	 fill_CI(l1, l2, l3, l4, PA.x, PB.x, QC.x, QD.x, QP.x, gamma1, gamma2,CI);
	 fill_CI(m1, m2, m3, m4, PA.y, PB.y, QC.y, QD.y, QP.y, gamma1, gamma2,CJ);
	 fill_CI(n1, n2, n3, n4, PA.z, PB.z, QC.z, QD.z, QP.z, gamma1, gamma2,CK);


	double sumC = 0.0;
	for (int nu = 0; nu <= numax; ++nu)
	{
		for (int I = 0; I <= Imax; ++I)
		{
			for (int J = 0; J <= Jmax; ++J)
			{
				int K = nu - I- J;
				 if ((K < 0) || (K>Kmax))
				 	continue;

				// printf("calculating CI %f\n",CI[I]);
				 sumC += CI[I]*CJ[J]*CK[K]*Boysf(nu, pq_sq/4.0/delta);
			}
		}
	}

	return 2*M_PI*M_PI/gamma2/gamma1*sqrt(M_PI/(gamma1+gamma2))*sumC*
			exp(-alpha1*alpha2*ab_sq/gamma1-alpha3*alpha4*cd_sq/gamma2);

	delete [] CI;
	delete [] CJ;
	delete [] CK;
}

double Algorithm::fill_CI(int l1, int l2, int l3, int l4,
 double PAx, double PBx, double QCx, double QDx, double QPx,double gamma1, double gamma2, double* CI )
{

	int Imax = l1+l2+l3+l4;
	int Lmax = l1+l2;
	int Mmax = l3+l4;
	
	HL = new double[Lmax+1];
	HM = new double[Mmax+1];

	fill_H(l1,l2,PAx,PBx,gamma1,HL);
	fill_H(l3,l4,QCx,QDx,gamma2,HM);

	double delta = 1/gamma1/4.0 + 1/gamma2/4.0; 

	double sum = 0.0;
	for (int I = 0; I <= Imax; ++I)
	{
		for (int L = 0; L <= Lmax; ++L)
		{
			for (int M = 0; M <= Mmax; ++M)
			{
				int Umax = (L + M)/2;
				int U = L + M - I;
				if ((U < 0 )|| (U > Umax))
					continue;


				sum += HL[L]*pow(-1,M)*HM[M]*
				gsl_sf_fact(L+M)*pow(-1,U)*pow(QPx,L+M-2*U)/
				(gsl_sf_fact(U)*gsl_sf_fact(L+M-2*U)*pow(delta, L+M-U));

			}
		}
		CI[I] = sum;
		sum = 0.0;
	}
	delete [] HM;
	delete [] HL;

	//return sum;


}

double Algorithm::fill_H(int l1, int l2, double a, double b, double gamma, double *H)
{
	if ((l1 + l2) == 0)
	{
		 H[0] = 1.0;
	}	
	else if ( (l1 == 1) && (l2 == 1) )
	{
		H[0] = a*b + 1/2.0/gamma;
		H[1] = (a+b)/4.0/gamma;
		H[2] = (1/4.0/gamma)*(1/4.0/gamma);
	}
	else if( ((l1 == 1) && (l2 == 0)) ||  ((l1 == 0) && (l2 == 1)) )
	{
		H[0] = a;
		H[1] = 1/4.0/gamma;
	}
	else
		printf("SHIT!\n");
}