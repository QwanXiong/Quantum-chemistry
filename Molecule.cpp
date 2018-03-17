
#include <stdio.h>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include "Molecule.h"

using std::vector;

Atom operator * ( const double alpha,const Atom& at)
{
	Atom temp_at(at.x*alpha, at.y*alpha, at.z*alpha );
	return temp_at;

}


void Molecule::get_positions()
{
	for (auto it = position_in_array.begin(); it != position_in_array.end(); ++it)
	{
		 printf("position: %d  %d\n",it -> first, position_in_array[it->first] );
	}
}

void Molecule::get_all_coords()
{
	printf("Listing all read coordinates:\n");
	for (auto it = set_of_atoms.begin(); it < set_of_atoms.end(); ++it)
	{
		for (auto itt = it->begin(); itt < it->end(); ++itt)
		{
				printf("x: %f  y: %f  z: %f\n", itt->x, itt->y, itt-> z);
		 	//printf("\n");
			 	
		}
	}
	// printf("number of atoms: %d  %d\n", number_of_given_atoms[1], number_of_given_atoms[6]);
	// printf("x-coord: %f   %f\n",set_of_atoms[position_in_array[1]][0].x, set_of_atoms[position_in_array[1]][1].x );
	// printf("\n");
}

void Molecule::print_charge_coords()
{
	printf("Listing all read coordinates and corresponding charges:\n");
	for (auto it = position_in_array.begin(); it != position_in_array.end(); ++it)
	{
		for (auto itt = set_of_atoms[it -> second].begin(); itt < set_of_atoms[it -> second].end(); ++itt)
		{
				printf("Z: %d x: %f  y: %f  z: %f\n",it->first ,itt->x, itt->y, itt-> z);
		 	//printf("\n");
			 	
		}
	}
	// printf("number of atoms: %d  %d\n", number_of_given_atoms[1], number_of_given_atoms[6]);
	// printf("x-coord: %f   %f\n",set_of_atoms[position_in_array[1]][0].x, set_of_atoms[position_in_array[1]][1].x );
	// printf("\n");
}

double Molecule::nuclear_repulsion_energy()
{
	double rep =0.0;
	double R =0.0;
	printf("Calculating nuclear repusion energy:\n");
	for (auto it1 = position_in_array.begin(); it1 != position_in_array.end(); ++it1)
	{
		for (auto itt1 = set_of_atoms[it1 -> second].begin(); itt1 < set_of_atoms[it1 -> second].end(); ++itt1)
		{
				for (auto it2 = position_in_array.begin(); it2 != position_in_array.end(); ++it2)
				{
					for (auto itt2 = set_of_atoms[it2 -> second].begin(); itt2 < set_of_atoms[it2-> second].end(); ++itt2)
					{
						if (itt2 != itt1)
						{
							//printf("calculation of sum\n");
							R = sqrt(pow((itt1->x-itt2->x),2.0)+pow((itt1->y-itt2->y),2.0)+pow((itt1->z-itt2->z),2.0));
							//printf("distance: %f\n",R );
							rep += (it1->first)*(it2->first)/2.0/R;
							//printf("charge: %d\n",it1->first );
							//printf("rep: %f\n",(it1->first)*(it2->first)/2.0/R );
						}
					}
				}
		 	//printf("\n");
			 	
		}
	}
	// printf("number of atoms: %d  %d\n", number_of_given_atoms[1], number_of_given_atoms[6]);
	// printf("x-coord: %f   %f\n",set_of_atoms[position_in_array[1]][0].x, set_of_atoms[position_in_array[1]][1].x );
	// printf("\n");
	nuc_rep_en = rep;
	printf("Nuclear repulsion energy is: %f\n", rep);
	return rep;
}

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
		unique_charges.push_back(charge);
		++curr_pos;

	//	printf("sizes: %d  %d\n",number_of_given_atoms.size(), position_in_array.size() );
	//	printf("vec size: %d\n", set_of_atoms[0].size());
	}
}

void Molecule::get_coords(int charge, vector<Atom> * atom_coords)
{

	atom_coords->assign(set_of_atoms[ position_in_array[charge] ].begin(), set_of_atoms[ position_in_array[charge] ].end());
}
