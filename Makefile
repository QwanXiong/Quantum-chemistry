all: main.exe clean

CXXFLAGS := -std=gnu++11 
LIBS := -lgsl -lgslcblas

main.exe: main.o  STO.o GTO.o HF_SCF.o Molecule.o Algorithm.o
	g++ $(CXXFLAGS)  main.o STO.o GTO.o HF_SCF.o Molecule.o Algorithm.o $(LIBS) -o testfiles/main.exe


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

Algorithm.o:
	g++ $(CXXFLAGS) -c  Algorithm.cpp

 


clean:
	rm *.o