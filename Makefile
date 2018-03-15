all: main.exe clean

CXXFLAGS := -std=gnu++11

main.exe: main.o  STO.o GTO.o HF_SCF.o Molecule.o
	g++ $(CXXFLAGS)  main.o STO.o GTO.o HF_SCF.o Molecule.o -o main.exe


main.o: 
	g++ $(CXXFLAGS) -c main.cpp

STO.o: 
	g++ $(CXXFLAGS) -c  STO.cpp 

GTO.o: 
	g++ $(CXXFLAGS) -c  GTO.cpp

HF_SCF.o:
	g++ $(CXXFLAGS) -c  HF_SCF.cpp

Molecule.o:
	g++ $(CXXFLAGS) -c  Molecule.cpp

 


clean:
	rm *.o