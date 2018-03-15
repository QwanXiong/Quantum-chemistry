#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <map>
#include <string>
#include <vector>

using std::ifstream;
using std::vector;


//for odd number n!!!
int dfact (int n)
{
	int res = 1;
	for (int i = 0; i < (n+1)/2; ++i)
		res *= (n-i*2);
	return res;
}


class GTO
{
public:
	GTO(double _alpha, int _l, int _m, int _n): alpha(_alpha), l(_l), m(_m), n(_n) {

	}
	GTO() 
	{

	}

	void set_data(double _alpha, int _l, int _m, int _n)
	{
		alpha = _alpha;
		l = _l;
		m = _m;
		n = _n;
	}

	void set_a(double _a)
	{
		a = _a;
	}
	int l, m, n;
	double alpha, a;
	double norm;

	//уже с учетом квадратного корня!
	void calc_norm()
	{
		double xf = (double)dfact(2*l-1);
		double yf = (double)dfact(2*m-1);
		double zf = (double)dfact(2*n-1);
		// /norm =  xf*yf*zf / pow(2.0*alpha, (l+m+n+1.5)) / pow(2.0, (l+m+n)) * pow(M_PI , 1.5);
		norm = (xf*yf*zf / pow(2.0*alpha, (l+m+n+1.5)) / pow(2.0, (l+m+n)) * pow(M_PI , 1.5));
	}

};

class STO
{
public:
	STO (int _ngauss, int _L, double * _alphai, double* _ai, double _r0) : L(_L), ngauss(_ngauss), r0(_r0)
	{


		alphai = new double [ngauss];

		if (L == -1)
			ai = new double [2*ngauss];
		else
			ai = new double [ngauss];

		for (int i = 0; i < ngauss; ++i)
		{
			alphai[i] = _alphai[i];
			if (L == -1)
			{
				ai[2*i+0] = _ai[2*i+0];
				ai[2*i+1] = _ai[2*i+1];
			}
			else
				ai[i] = _ai[i];
		}

		
		nang = (L + 2)*(L + 1)/2;


		if (L == -1)
			nfunc = 4*ngauss;
		else
			nfunc = ngauss*nang;
		
		primitives = new GTO[nfunc];
		int counter = 0;


		if ( L == -1)
		{
			for (int ng = 0; ng < ngauss; ++ng)
			{
				printf("a unnormalized: %f  %f\n",ai[2*ng+0], ai[2*ng+1]);
				primitives[4*ng+0].set_data(alphai[ng], 0, 0, 0);
				primitives[4*ng+0].calc_norm();

				primitives[4*ng+1].set_data(alphai[ng], 1, 0, 0);
				primitives[4*ng+2].set_data(alphai[ng], 0, 1, 0);
				primitives[4*ng+3].set_data(alphai[ng], 0, 0, 1);
				
				primitives[4*ng+1].calc_norm();
						
				//ai[2*ng+0] /= sqrt(primitives[4*ng+0].norm);
				//ai[2*ng+1] /= sqrt(primitives[4*ng+1].norm);

				primitives[4*ng+0].set_a(ai[2*ng+0]/sqrt(primitives[4*ng+0].norm));

				for (int i = 1; i < 4; ++i)
				{
					primitives[4*ng+i].set_a(ai[2*ng+1]/sqrt(primitives[4*ng+1].norm));
				}
						
						//printf("%d  %d  %d  %d  %d\n",counter, ng, l, k, L-l-k);

				//++counter;
				 
				printf("a normalized: %f  %f\n",ai[2*ng+0], ai[2*ng+1]);
			}// end ng


		}
		else
		{
			for (int ng = 0; ng < ngauss; ++ng)
			{
				for (int l = 0; l <= L; ++l)
				{
					for (int k = 0; k <= (L-l); ++k)
					{
						primitives[counter].set_data(alphai[ng], l, k, L-l-k);
						primitives[counter].calc_norm();
						
						//ai[ng] /= sqrt(primitives[counter].norm);
						printf("norm: %f\n",primitives[counter].norm);
						primitives[counter].set_a(ai[ng]/sqrt(primitives[counter].norm));
						
						//printf("%d  %d  %d  %d  %d\n",counter, ng, l, k, L-l-k);

						++counter;
					}// end k
				}//end l
				//printf("%f\n",ai[ng]);
				//printf("number of functions: %d\n", nfunc);
			}// end ng

		}// end else
	}// end constructor

	~STO()
	{
		delete [] alphai;
		delete [] ai;
	}

	int L, ngauss, nang, nfunc;
	double r0;
	double * alphai;
	double * ai;
	GTO* primitives;

	//исходный контракционный коэффициент необходимо умножать на norm, если 
	//"STO" нормирована с изначальными коэффициентами и примитивами as is
	//и делить на norm, если изначально она нормирована с нормированными примитивами!

	double compute_overlap (GTO * g1, GTO* g2)
	{
		int lsum = (g1 -> l+g2 -> l);
		int msum = (g1 -> m+g2 -> m);
		int nsum = (g1 -> n+g2 -> n);

		if ((lsum % 2 == 1) || (msum % 2 == 1) || (nsum % 2 == 1))
			return 0;

		GTO gsum((g1->alpha+g2->alpha)/2.0, lsum/2, msum/2, nsum/2);
		gsum.calc_norm();
		return gsum.norm;
	}

	double compute_STO_norm ()
	{
		double norm = 0.0;
		for (int i = 0; i < nfunc; ++i)
		{
			for (int j = i; j < nfunc; ++j)
			{
				double coef = (i == j) ? 1.0 : 2.0;
				norm += coef*primitives[i].a*primitives[j].a* compute_overlap(&primitives[i], &primitives[j]);
				//printf("%f  %f  %f\n", primitives[i].a, primitives[j].a, compute_overlap(&primitives[i], &primitives[j]));
			}
		}
		return norm;
	}



};

// void reader()
// {

// }

void read_basis(const char* filename)
{
	ifstream infile(filename);
	char buffer[256];
	char* buffer2;

	int ngauss;
	int L;
	double * alpha;
	double * a;

	infile.getline(buffer, 256);
	while(!infile.eof())
	{

		
		buffer2 = strtok(buffer, " ");
		/*buffer2 = strtok(NULL, " ");
		buffer2 = strtok(NULL, " ");
		buffer2 = strtok(NULL, " ");*/
		if ((buffer2 != NULL) && !strcmp(buffer2,"S"))
			L = 0;

		if ((buffer2 != NULL) && !strcmp(buffer2,"P"))
			L = 1;

		if ((buffer2 != NULL) && !strcmp(buffer2,"D"))
			L = 2;

		if ((buffer2 != NULL) && !strcmp(buffer2,"F"))
			L = 3;

		if ((buffer2 != NULL) && !strcmp(buffer2,"L"))
			L = -1;

		buffer2 = strtok(NULL, " ");
		if (buffer2 != NULL)
			ngauss = atoi(buffer2);

		alpha = new double[ngauss];
		
		if (L == -1)
			a = new double [ngauss*2];
		else
			a = new double[ngauss];


		for (int i = 0; i < ngauss; ++i)
		{
			infile.getline(buffer, 256);
			buffer2 = strtok(buffer, " ");
			buffer2 = strtok(NULL, " ");
			//printf("%s\n", buffer2);
			if (L == -1)
			{
				if (buffer2 != NULL)
					alpha[i] = atof(buffer2);
				 buffer2 = strtok(NULL, " ");
				 
				if (buffer2 != NULL)
					a[2*i+0] = atof(buffer2);
				buffer2 = strtok(NULL, " ");
				 
				if (buffer2 != NULL)
					a[2*i+1] = atof(buffer2);
			}
			else
			{
				if (buffer2 != NULL)
					alpha[i] = atof(buffer2);
				 buffer2 = strtok(NULL, " ");
				// printf("printing\n");
				// printf("%f\n", atof(buffer2));
				if (buffer2 != NULL)
					a[i] = atof(buffer2);
			}
			

		}

		infile.getline(buffer, 256);

		 
	}

	for (int i = 0; i < ngauss; ++i)
	{
		printf("alpha[%d] = %f\n",i, alpha[i]);
	}

	if (L == -1)
	{
		for (int i = 0; i < 2*ngauss; ++i)
		{
			printf("a[%d] = %f\n",i, a[i]);
		}
	}
	else
	{
		for (int i = 0; i < ngauss; ++i)
		{
			printf("a[%d] = %f\n",i, a[i]);
		}
	}
	
	printf("L = %d\n", L);
	printf("ngauss = %d\n", ngauss);

	STO orb(ngauss, L, alpha, a, 0.0);

	printf("norm = %f\n", orb.compute_STO_norm());

	infile.close();
	delete[] a;
	delete[] alpha;
}




struct Atom
{
	Atom(){}
	Atom(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
	double x, y, z;

};

class Molecule
{
public:
	Molecule():curr_pos(0) {}
	void add_atom(int charge, double x , double y, double z);
	void get_coords(int charge, vector<Atom> * atom_coords);
	void get_positions()
	{
		for (auto it = position_in_array.begin(); it != position_in_array.end(); ++it)
		{
			 printf("position: %d  %d\n",it -> first, position_in_array[6] );
		}
	}
	void get_all_coords()
	{
		for (auto it = set_of_atoms.begin(); it < set_of_atoms.end(); ++it)
		{
			for (auto itt = it->begin(); itt < it->end(); ++itt)
			{
			 	printf("x: %f  y: %f  z: %f\n", itt->x, itt->y, itt-> z);
			 	//printf("\n");
			 	
			}
		}
	}


private:
	
	//Atom ** set_of_atoms;
	vector<vector<Atom> > set_of_atoms;
	bool check_charge(int charge);
	std::map<int, int> position_in_array;
	std::map<int,int> number_of_given_atoms;
	std::map<int, int>::iterator it;
	int curr_pos; //позиция, в которую будет вписываться новый массив атомов одного вида




};

void Molecule::add_atom(int charge, double x , double y, double z)
{

	if ((it = position_in_array.find(charge)) != position_in_array.end())
	{
		number_of_given_atoms[charge]++; //увеличиваем на 1 количество атомов этого типа
		set_of_atoms[it -> second].emplace_back(x,y,z); // добавляем координаты данного типа атомов

	}
	else
	{
		//printf("%d  %f  %f  %f\n", curr_pos, x, y, z);
		number_of_given_atoms[charge] = 1; //создаем новый член в map
		//vector<Atom> curr_at(1,Atom(x,y,z));
		//curr_at.emplace_back(x,y,z);
		//set_of_atoms.push_back(curr_at);
		set_of_atoms.emplace_back(1,Atom(x,y,z));// добавляем в набор атомов вектор, который инициализируем одним элементом Atom
		//set_of_atoms[curr_pos].emplace_back(x,y,z);
		position_in_array[charge] = curr_pos;
		++curr_pos;

	//	printf("sizes: %d  %d\n",number_of_given_atoms.size(), position_in_array.size() );
	//	printf("vec size: %d\n", set_of_atoms[0].size());
	}
}

void Molecule::get_coords(int charge, vector<Atom> * atom_coords)
{

	atom_coords->assign(set_of_atoms[ position_in_array[charge] ].begin(), set_of_atoms[ position_in_array[charge] ].end());
}


void read_coords(const char* filename, Molecule* curr_mol)
{

	ifstream infile(filename);
	char buffer[256];
	char* buffer2;

	 
	 
	int charge;
	double x, y, z;

	
	while(!infile.eof())
	{
		infile.getline(buffer, 256);
		
		buffer2 = strtok(buffer, " ");
		if (buffer2 == NULL)
			continue;
		
		buffer2 = strtok(NULL, " ");
		/*buffer2 = strtok(NULL, " ");
		buffer2 = strtok(NULL, " ");
		buffer2 = strtok(NULL, " ");*/
		 

		if (buffer2 != NULL)
			charge = atoi(buffer2);

		buffer2 = strtok(NULL, " ");

		if (buffer2 != NULL)
			x = atof(buffer2);

		buffer2 = strtok(NULL, " ");

		if (buffer2 != NULL)
			y = atof(buffer2);

		buffer2 = strtok(NULL, " ");

		if (buffer2 != NULL)
			z = atof(buffer2);	


		//printf("%d  %f  %f  %f\n", charge, x, y, z);

		curr_mol->add_atom(charge, x, y, z);

		vector<Atom> atoms ;
		 curr_mol->get_coords(charge, &atoms);
		 for (auto it = atoms.begin(); it < atoms.end(); ++it)
		 {
		 	printf("x: %f  y: %f  z: %f\n", it->x, it->y, it-> z);
		 	//printf("\n");
		 	
		 }
		 printf(" \n");

		//infile.getline(buffer, 256);
	} 
	 

	 curr_mol->get_positions();

	 printf("\n\n");
	 curr_mol->get_all_coords();
	 
	 

 

	infile.close();
	 

}

int main(int argc, char const *argv[])
{
	// printf("%d\n", dfact(13));

	// GTO gto1(2.432879, 0, 0, 0);
	// GTO gto2(0.433051, 0, 0, 0);

	// GTO * objects = new GTO[2];

	// objects[0].set_data(2.432879, 0, 0, 0);
	// objects[1].set_data(0.433051, 0, 0, 0);

	// objects[0].calc_norm();



	// gto1.calc_norm();
	// gto2.calc_norm();
	// printf("normalization: %f\n", gto1.norm );
	// printf("normalization: %f\n", objects[0].norm );
	// printf("%f\n", pow(0.5643,3+1+4));

	// delete [] objects;


	// double * alpha = new double[2];
	// alpha[0] = 2.432879;
	// alpha[1] = 0.433051;

	// double * a = new double[2];
	// a[0] = 0.59717;
	// a[1] = 0.258303;

	// STO orb(2, 0, alpha, a, 0.0);

	// delete[] a;
	// delete[] alpha;


	//read_basis("basis.txt");
	Molecule current_molecule;
	read_coords("coords.txt",&current_molecule);
	

	

	system("pause");
	return 0;
}