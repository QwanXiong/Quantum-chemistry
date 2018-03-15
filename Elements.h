#ifndef ELEMENTS_H
#define ELEMENTS_H 

#include <map>

class ChemElements
{
public:
	ChemElements()
	{
		elements["HYDROGEN"] = 1;
		elements["HELIUM"] = 2;
		elements["CARBON"] = 6;
		elements["NITROGEN"] = 7;
		elements["OXYGEN"] = 8;
		elements["FLUORINE"] = 9;
		elements["NEON"] = 10;

	}

	int if_get_charge(std::string name)
	{
		//return elements[name];
		if ((it = elements.find(name)) != elements.end())
		{
			return it->second;

		}
		else
			return -1;
	}

private:
	std::map<std::string, int> elements;
	std::map<std::string, int>::iterator it;

};



#endif