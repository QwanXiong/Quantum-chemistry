#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <map>
#include <string>
#include <vector>

#include "Algorithm.h"
#include "Molecule.h"
#include "Elements.h"
#include "GTO.h"
#include "STO.h"
#include "HF_SCF.h"

using std::ifstream;
using std::vector;



void read_basis(const char* filename, SCF* scf_operator, Molecule* curr_mol);
 


void read_coords(const char* filename, Molecule* curr_mol);

int main(int argc, char const *argv[])
{

	
	char basis_filename[256];
	char coord_filename[256];

	printf("Write the name of basis file:\n");
	scanf("%s", basis_filename);
	printf("Write the name of coordinates file:\n");
	scanf("%s", coord_filename);



	 Molecule current_molecule;
	 read_coords(coord_filename,&current_molecule);
	// read_coords("coords.txt",&current_molecule);

	SCF scf_operator ;
	scf_operator.setN(10);

	
	read_basis(basis_filename, &scf_operator, &current_molecule);
	// read_basis("basis.txt", &scf_operator, &current_molecule);

	 scf_operator.print_all_orbitals("orbitals.txt");


	 



	

	system("pause");
	return 0;
}


	
 /*
	Atom at1(0.3, 0.0, -0.0);
	 Atom at2(0.0, 0.4, 0.5);


	  GTO g1(0.19,1, 1, 1, &at1);
	  GTO g2(0.19, 0, 1, 0, &at1);
	   GTO g3(0.19, 1, 1, 0, &at2);
	    GTO g4(0.19, 0, 1, 0, &at2);
	 Atom nuc(0.2, 0.3, 1.2);

 
	
	for (int i = 0; i < 100000; ++i)
	{
		scf_operator.compute_gto_coulomb_2e(&g1, &g2, &g3, &g4);
	}
	printf("integral: %f\n", scf_operator.compute_gto_coulomb_2e(&g1, &g2, &g3, &g4));

	  */


	/*printf("overlap: %f\n",scf_operator.od_overlap(8,2, 0.12, 0.43, -1.2, 0.14));
	scf_operator.plan_1e_calculation(&g1, &g2,&nuc);
	printf("overlap: %f\n",scf_operator.compute_gto_overlap_v1(&g1, &g2));
	printf("overlap: %f\n",scf_operator.compute_gto_overlap_v2( &g1, &g2));
	scf_operator.finish_1e_calculation();*/

	


 

	
	/* printf("overlap between STOs: %f\n",
	 	scf_operator.compute_sto_overlap_v2(&scf_operator.orbitals[0], &scf_operator.orbitals[1]));

	 printf("calculating\n");
	 for (int i = 0; i < 1000000; ++i)
	 {
	 	scf_operator.compute_sto_overlap_v1(&scf_operator.orbitals[0], &scf_operator.orbitals[1]);
	 }
*/


void read_basis(const char* filename, SCF* scf_operator, Molecule* curr_mol)
{
	ChemElements Elements;
	ifstream infile(filename);
	char buffer[256];
	char* buffer2;

	int ngauss;
	int L;
	double * alpha;
	double * a;
	int charge, curr_charge;

	std::vector<Atom> curr_atoms;

	
	while(!infile.eof())
	{
		infile.getline(buffer, 256);
		
		buffer2 = strtok(buffer, " ");
		if (buffer2 == NULL)
			continue;

		if ((charge = Elements.if_get_charge(std::string(buffer2)) ) != -1)//если нашли имя элемента
		{
			infile.getline(buffer, 256);//читаем дальше
			buffer2 = strtok(buffer, " ");
			curr_charge = charge;
		}
		else// если не нашли, не делаем ничего -- может, это что-то нужное
		{
			// infile.getline(buffer, 256);
			// continue;
		}

		

		//printf("buffer: %s\n", buffer);

		//buffer2 = strtok(buffer, " ");
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

		 curr_mol->get_coords(curr_charge, &curr_atoms);
		 //printf("coordintes from main: %f  %f\n", curr_atoms[0].x, curr_atoms[1].x);
		 scf_operator->add_orbital(ngauss, L, alpha, a,&curr_atoms );
		// printf("size of curr_atoms: %d\n", curr_atoms.size());
	
		printf("charge: %d  %d\n", curr_charge, charge);

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

		delete[] a;
		delete[] alpha;

		 
	}



	//STO orb(ngauss, L, alpha, a, 0.0);
	

	//printf("norm = %f\n", orb.compute_STO_norm());

	infile.close();
	
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
		 // for (auto it = atoms.begin(); it < atoms.end(); ++it)
		 // {
		 // 	printf("x: %f  y: %f  z: %f\n", it->x, it->y, it-> z);
		 // 	//printf("\n");
		 	
		 // }
		// printf(" \n");

		//infile.getline(buffer, 256);
	} 
	 

	 //curr_mol->get_positions();

	 //printf("\n\n");
	 curr_mol->get_all_coords();
	 curr_mol->print_charge_coords();
	 curr_mol -> nuclear_repulsion_energy();
	 
	 

 

	infile.close();
	 

}