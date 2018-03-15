
#include <stdio.h>
#include <map>
#include <string>
#include <vector>
#include "Molecule.h"

using std::vector;


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
