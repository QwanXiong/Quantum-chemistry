#ifndef MOLECULE_H
#define MOLECULE_H 


 

#include <map>
#include <string>
#include <vector>

using std::vector;


struct Atom
{
	Atom():x(0.0), y(0.0), z(0.0){}
	Atom(double _x, double _y, double _z): x(_x), y(_y), z(_z) {}
	Atom (const Atom& at)
	{
		x = at.x;
		y = at.y;
		z = at.z;
	}

	~Atom(){}
	double x, y, z;
	Atom operator + (const Atom & at)
	{
		Atom temp_at(this->x+at.x, this->y + at.y, this->z + at.z );
		return temp_at;

	}

	Atom& operator=(const Atom& at) {
        //проверка на самоприсваивание
        if (this == &at) {
            return *this;
        }
        x= at.x;
        y = at.y;
        z = at.z;
        return *this;
    }

};


class Molecule
{
public:
	Molecule():curr_pos(0) {}
	void add_atom(int charge, double x , double y, double z);
	void get_coords(int charge, vector<Atom> * atom_coords);
	void get_positions();
 
	void get_all_coords();
	 


private:
	
	//Atom ** set_of_atoms;
	vector<vector<Atom> > set_of_atoms;
	bool check_charge(int charge);
	std::map<int, int> position_in_array;
	std::map<int,int> number_of_given_atoms;
	std::map<int, int>::iterator it;
	int curr_pos; //позиция, в которую будет вписываться новый массив атомов одного вида




};



#endif